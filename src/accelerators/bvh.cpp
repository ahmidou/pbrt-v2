
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// accelerators/bvh.cpp*
#include "stdafx.h"
#include "accelerators/bvh.h"
#include "probes.h"
#include "paramset.h"
#include "timer.h"

#include <algorithm>
#include <cilk/cilk.h>
#include "radixSort.h"

#define AAC_PARALLEL_GRAN 6
#define AAC_DELTA 20 // 20 for HQ, 4 for fast
#define AAC_EPSILON 0.1f //0.1f for HQ, 0.2 for fast
#define AAC_ALPHA (0.5f-AAC_EPSILON)
//#define AAC_C (0.5f * pow(AAC_DELTA, 0.5f + AAC_EPSILON))
#define MORTON_CODE_END 0
#define MORTON_CODE_START 29

#pragma mark - morton code
/* morton code */

Timer combineTimer;
Timer leafsTimer;
Timer dataMergeTimer;
Timer appendTimer;
Timer baseCaseTimer;
Timer bitchTimer;

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
uint32_t expandBits(uint32_t v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
uint32_t morton3D(float x, float y, float z)
{
//    Warning("morton3D: x == %f", x);
    
    x = min(max(x * 1024.0f, 0.0f), 1023.0f);
    y = min(max(y * 1024.0f, 0.0f), 1023.0f);
    z = min(max(z * 1024.0f, 0.0f), 1023.0f);
    
//    Warning("morton3D: new x == %u", (unsigned int)x);
    
    uint32_t xx = expandBits((uint32_t)x);
//    Warning("morton3D: xx == %u", xx);
    uint32_t yy = expandBits((uint32_t)y);
    uint32_t zz = expandBits((uint32_t)z);
    return ((zz << 2) | (yy << 1) | xx);
}

static inline int getDimForMorton3D(uint32_t idx)
{
    return idx % 3;
}

static inline float AAC_C() {
    return (0.5f * powf(AAC_DELTA, 0.5f + AAC_EPSILON));
}

static inline uint32_t AAC_F(uint32_t x) {
    return (uint32_t) (ceil(AAC_C() * powf(x, AAC_ALPHA)));
}

/* end of morton code stuff*/


// BVHAccel Local Declarations

struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() { }
    BVHPrimitiveInfo(int pn, const BBox &b)
        : primitiveNumber(pn), bounds(b) {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }
    int primitiveNumber;
    Point centroid;
    BBox bounds;
};


struct BVHBuildNode {
    // BVHBuildNode Public Methods
    BVHBuildNode() { children[0] = children[1] = NULL; }
    ~BVHBuildNode() { }
    void InitLeaf(uint32_t first, uint32_t n, const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
    void InitInterior(uint32_t axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
    BVHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};


struct CompareToMid {
    CompareToMid(int d, float m) { dim = d; mid = m; }
    int dim;
    float mid;
    bool operator()(const BVHPrimitiveInfo &a) const {
        return a.centroid[dim] < mid;
    }
};


struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BVHPrimitiveInfo &a,
                    const BVHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};


struct CompareToBucket {
    CompareToBucket(int split, int num, int d, const BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const BVHPrimitiveInfo &p) const;

    int splitBucket, nBuckets, dim;
    const BBox &centroidBounds;
};


bool CompareToBucket::operator()(const BVHPrimitiveInfo &p) const {
    int b = nBuckets * ((p.centroid[dim] - centroidBounds.pMin[dim]) /
            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
    if (b == nBuckets) b = nBuckets-1;
    Assert(b >= 0 && b < nBuckets);
    return b <= splitBucket;
}


struct LinearBVHNode {
    BBox bounds;
    union {
        uint32_t primitivesOffset;    // leaf
        uint32_t secondChildOffset;   // interior
    };

    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};


static inline bool IntersectP(const BBox &bounds, const Ray &ray,
        const Vector &invDir, const uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}



// BVHAccel Method Definitions
BVHAccel::BVHAccel(const vector<Reference<Primitive> > &p,
                   uint32_t mp, const string &sm) {
    maxPrimsInNode = min(255u, mp);
    for (uint32_t i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    if (sm == "sah")         splitMethod = SPLIT_SAH;
    else if (sm == "aac")    splitMethod = SPLIT_AAC;
    else if (sm == "middle") splitMethod = SPLIT_MIDDLE;
    else if (sm == "equal")  splitMethod = SPLIT_EQUAL_COUNTS;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                sm.c_str());
        splitMethod = SPLIT_SAH;
    }
    
    Info("BVH split method \"%s\"", sm.c_str());

    if (primitives.size() == 0) {
        nodes = NULL;
        return;
    }
    // Build BVH from _primitives_
    PBRT_BVH_STARTED_CONSTRUCTION(this, primitives.size());

    // Initialize _buildData_ array for primitives
    vector<BVHPrimitiveInfo> buildData;
    buildData.reserve(primitives.size());
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        BBox bbox = primitives[i]->WorldBound();
        buildData.push_back(BVHPrimitiveInfo(i, bbox));
    }

    // Recursively build BVH tree for primitives
    MemoryArena buildArena;
    uint32_t totalNodes = 0;
    vector<Reference<Primitive> > orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode *root;
    if (splitMethod == SPLIT_AAC) {
        buildAAC(buildData, 0,
                              primitives.size(), &totalNodes,
                              orderedPrims);
    } else {
        root = recursiveBuild(buildArena, buildData, 0,
                              primitives.size(), &totalNodes,
                              orderedPrims);
        
        primitives.swap(orderedPrims);

        // Compute representation of depth-first traversal of BVH tree
        nodes = AllocAligned<LinearBVHNode>(totalNodes);
        for (uint32_t i = 0; i < totalNodes; ++i)
            new (&nodes[i]) LinearBVHNode;
        uint32_t offset = 0;
        flattenBVHTree(root, &offset);
        Assert(offset == totalNodes);
    }
    Info("BVH created with %d nodes for %d primitives (%.2f MB)", totalNodes,
         (int)primitives.size(), float(totalNodes * sizeof(LinearBVHNode))/(1024.f*1024.f));
    PBRT_BVH_FINISHED_CONSTRUCTION(this);
}


BBox BVHAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : BBox();
}
#pragma mark - AAC

uint32_t findBestMatch(std::vector<BVHBuildNode*> &clusters, uint32_t i) {
    float closestDist = INFINITY;
    uint32_t idx = i;
    for (uint32_t j = 0; j < clusters.size(); ++j) {
        if (i == j) continue;
        
        BBox combined = Union(clusters[i]->bounds, clusters[j]->bounds);
        float d = combined.SurfaceArea();
        if (d < closestDist) {
            closestDist = d;
            idx = j;
        }
    }
    
    return idx;
}

std::vector<BVHBuildNode*> combineCluster(std::vector<BVHBuildNode*> &clusters, uint32_t endSize,
                                          uint32_t *totalNodes, int dim) {
    
    std::vector<uint32_t> closest(clusters.size(), 0);
    
    for (uint32_t i = 0; i < clusters.size(); ++i) {
        closest[i] = findBestMatch(clusters, i);
    }
    
    while (clusters.size() > endSize) {
        float bestDist = INFINITY;
        uint32_t leftIdx = 0;
        uint32_t rightIdx = 0;
        
        for (uint32_t i = 0; i < clusters.size(); ++i) {
            BBox combined = Union(clusters[i]->bounds, clusters[closest[i]]->bounds);
            float d = combined.SurfaceArea();
            if (d < bestDist) {
                bestDist = d;
                leftIdx = i;
                rightIdx = closest[i];
            }
        }
        
        (*totalNodes)++;
        BVHBuildNode *node = new BVHBuildNode();
        
        node->InitInterior(dim,
                           clusters[leftIdx],
                           clusters[rightIdx]);
        clusters[leftIdx] = node;
        clusters[rightIdx] = clusters.back();
        closest[rightIdx] = closest.back();
        clusters.pop_back();
        closest.pop_back();
        closest[leftIdx] = findBestMatch(clusters, leftIdx);
        
        for (uint32_t i = 0; i < clusters.size(); ++i) {
            if (closest[i] == leftIdx || closest[i] == rightIdx)
                closest[i] = findBestMatch(clusters, i);
            else if (closest[i] == closest.size()) {
                closest[i] = rightIdx;
            }
        }
    }
    
    
    return clusters;
}

uint32_t makePartition(std::pair<uint32_t, uint32_t> *sortedMC,
                       uint32_t start, uint32_t end, int partitionbit) {
    
    uint32_t curFind = (1 << partitionbit);
    
    if (((sortedMC[start].first & curFind) == (sortedMC[end-1].first & curFind))
        || (partitionbit < MORTON_CODE_END)) {
        if (partitionbit < MORTON_CODE_END)
            Info("makePartition spliting in half. SUPPOSED TO BE RARE");
        return start + (end-start)/2;
    }
    
    uint32_t lower = start;
    uint32_t upper = end;
    
    while (lower < upper) {
        uint32_t mid = lower + (upper-lower)/2;
        if ((sortedMC[mid].first & curFind) == 0) {
            lower = mid+1;
        } else {
            upper = mid;
        }
    }
    
    return lower;
}

void findBestMatch(AAC_Data *aacData, AAC_DataCoord *dc, int n, int cIdx) {
//    float closestDist = INFINITY;
//    uint32_t idx = 0;
    int start = dc->start;
    
    
    aacData->minArea[cIdx] = INFINITY;
    for (int j=start; j<cIdx; j++)
        if (aacData->area[cIdx][j] < aacData->minArea[cIdx]) {
            aacData->minArea[cIdx] = aacData->area[cIdx][j];
            aacData->minPos[cIdx] = j;
        }
    for (int j=cIdx+1; j<start+n; j++)
        if (aacData->area[j][cIdx] < aacData->minArea[cIdx]) {
            aacData->minArea[cIdx] = aacData->area[j][cIdx];
            aacData->minPos[cIdx] = j;
        }
}

std::vector<BVHBuildNode*> combineCluster(std::vector<BVHBuildNode*> &clusters, uint32_t endSize,
                        uint32_t *totalNodes, int dim, AAC_Data *aacData, AAC_DataCoord *dc) {
    
//    std::vector<uint32_t> closest(clusters.size(), 0);
    
//    for (uint32_t i = 0; i < clusters.size(); ++i) {
//        closest[i] = findBestMatch(clusters, i);
//    }
    
//    Warning("combineCluster START");
    
    while (clusters.size() > endSize) {
        float bestDist = INFINITY;
        uint32_t leftIdx = 0;
        uint32_t rightIdx = 0;
        
        for (uint32_t i = 0; i < clusters.size(); ++i) {
            int adjustedI = i + dc->start;
            float d = aacData->minArea[adjustedI];
            if (d < bestDist) {
                bestDist = d;
                leftIdx = i;
                rightIdx = aacData->minPos[adjustedI] - dc->start;
            }
        }
        
        (*totalNodes)++;
        BVHBuildNode *node = new BVHBuildNode();
        
        node->InitInterior(dim,
                           clusters[leftIdx],
                           clusters[rightIdx]);
        clusters[leftIdx] = node;
        clusters[rightIdx] = clusters.back();
//        closest[rightIdx] = closest.back();
        
        uint32_t ai = leftIdx+dc->start;
        for (uint32_t j=dc->start; j<ai; j++) {
            BBox combined = Union(clusters[leftIdx]->bounds, clusters[j-dc->start]->bounds);
            aacData->area[ai][j] = combined.SurfaceArea();
        }
        for (uint32_t j=ai+1; j<dc->start+clusters.size(); j++) {
            BBox combined = Union(clusters[leftIdx]->bounds, clusters[j-dc->start]->bounds);
            aacData->area[j][ai] = combined.SurfaceArea();
        }
        
        clusters.pop_back();
//        closest.pop_back();
        
        uint32_t bi = rightIdx+dc->start;
        for (uint32_t j=dc->start; j<bi; j++)
            aacData->area[bi][j] = aacData->area[dc->start+clusters.size()][j];
        for (uint32_t j=bi+1; j<dc->start+clusters.size(); j++)
            aacData->area[j][bi] = aacData->area[dc->start+clusters.size()][j];
        
        aacData->minArea[bi] = aacData->minArea[dc->start+clusters.size()];
        aacData->minPos[bi] = aacData->minPos[dc->start+clusters.size()];
        
        findBestMatch(aacData, dc, clusters.size(), ai);
        
        for (uint32_t i = 0; i < clusters.size(); ++i) {
            uint32_t newI = i+dc->start;
            if (uint32_t(aacData->minPos[newI]) == ai || uint32_t(aacData->minPos[newI]) == bi) {
                findBestMatch(aacData, dc,clusters.size(), newI);
            } else if (uint32_t(aacData->minPos[newI]) == dc->start+clusters.size()) {
                aacData->minPos[newI] = bi;
            }
        }
    }
    
    dc->end = dc->start + clusters.size();
//    Warning("combineCluster END");
    return clusters;
}

//Set leaves.
void aacSetLeaf(std::vector<BVHBuildNode*> clusters, AAC_Data *aacData, AAC_DataCoord dc) {
//    Warning("aacSetLeaf START");
    int numleaf = clusters.size();
    int start = dc.start;
    for (int i=0; i<numleaf; i++) {
        aacData->minArea[start+i] = INFINITY;
    }
    
    for (int i=start; i<start+numleaf; i++)  {
        for (int j=start; j<i; j++) {
            BBox combined = Union(clusters[i-start]->bounds, clusters[j-start]->bounds);
            float d = combined.SurfaceArea();
            aacData->area[i][j] = d;
            if (aacData->area[i][j] < aacData->minArea[i]) {
                aacData->minArea[i] = aacData->area[i][j];
                aacData->minPos[i] = j;
            }
            if (aacData->area[i][j] < aacData->minArea[j]) {
                aacData->minArea[j] = aacData->area[i][j];
                aacData->minPos[j] = i;
            }
        }
    }
//    Warning("aacSetLeaf END");
}

//Combine two subtaskes.
void aacDataMerge(std::vector<BVHBuildNode*> clusters, AAC_Data *aacData, AAC_DataCoord leftC, AAC_DataCoord rightC) {
//    Warning("aacDataMerge START");
    int start = leftC.start;
    int len1 = leftC.end- leftC.start;
    int len2 = rightC.end-rightC.start;
    for (int i=start; i<start+len1; i++)
        for (int j=start+len1; j<start+len1+len2; j++) {
            // align i and j to zero for clusters idx
//            Warning("aacDataMerge Getting Surface Area");
            BBox combined = Union(clusters[i-start]->bounds, clusters[j-start]->bounds);
            float d = combined.SurfaceArea();
//            Warning("aacDataMerge Done Getting Surface Area");

            aacData->area[j][i] = d;
            if (aacData->area[j][i] < aacData->minArea[i]) {
                aacData->minArea[i] = aacData->area[j][i];
                aacData->minPos[i] = j;
            }
            if (aacData->area[j][i] < aacData->minArea[j]) {
                aacData->minArea[j] = aacData->area[j][i];
                aacData->minPos[j] = i;
            }
        }
    
//    Warning("aacDataMerge END");
}

void recursiveBuildAAC(vector<BVHPrimitiveInfo> &buildData,
                                std::pair<uint32_t, uint32_t> *mortonCodes,
                                uint32_t start,
                                uint32_t end, uint32_t *totalNodes, int partitionBit,
                                std::vector<BVHBuildNode*> *clusterData,
                                AAC_Data *aacData, AAC_DataCoord *coord) {
    
    if (end-start == 0) {
        return;
    }
    
    int dim = getDimForMorton3D(partitionBit);
    
    if (end-start < AAC_DELTA) {
//        baseCaseTimer.Start();
        Assert(coord->start + (end-start) < aacData->size);
        
        std::vector<BVHBuildNode*> clusters;
        (*totalNodes) += (end-start);
        for (uint32_t i = start; i < end; ++i) {
            // Create leaf _BVHBuildNode_
            BVHBuildNode *node = new BVHBuildNode();
            uint32_t primIdx = mortonCodes[i].second;
            node->InitLeaf(primIdx, 1, buildData[primIdx].bounds);
            // deal with firstPrimOffset later with DFS
            clusters.push_back(node);
        }
        
//        baseCaseTimer.Stop();
//        
//        leafsTimer.Start();
        aacSetLeaf(clusters, aacData, *coord);
//        leafsTimer.Stop();
        
//        combineTimer.Start();
        *clusterData = combineCluster(clusters, AAC_F(AAC_DELTA), totalNodes, dim, aacData, coord);
//        combineTimer.Stop();
        return;
    }
    
    
    uint32_t splitIdx = makePartition(mortonCodes, start, end, partitionBit);
    
    int newPartionBit = partitionBit - 1;
    std::vector<BVHBuildNode*> leftC;
    std::vector<BVHBuildNode*> rightC;
    uint32_t rightTotalnodes = 0;
    
//    cilk_spawn recursiveBuildAAC(buildData, mortonCodes, start, splitIdx, totalNodes, newPartionBit, &leftC);
//    recursiveBuildAAC(buildData, mortonCodes, splitIdx, end, &rightTotalnodes, newPartionBit, &rightC);
//    cilk_sync;
    AAC_DataCoord leftDC = AAC_DataCoord(*coord);
    recursiveBuildAAC(buildData, mortonCodes, start, splitIdx, totalNodes, newPartionBit, &leftC, aacData, &leftDC);
    AAC_DataCoord rightDC = AAC_DataCoord(*coord);
    rightDC.start = leftDC.end;
    recursiveBuildAAC(buildData, mortonCodes, splitIdx, end, &rightTotalnodes, newPartionBit, &rightC, aacData, &rightDC);
    
    (*totalNodes) += rightTotalnodes;
    
//    appendTimer.Start();
    leftC.insert( leftC.end(), rightC.begin(), rightC.end() );
//    appendTimer.Stop();
    
    
//    dataMergeTimer.Start();
    aacDataMerge(leftC, aacData, leftDC, rightDC);
//    dataMergeTimer.Stop();
    
//    combineTimer.Start();
    *clusterData = combineCluster(leftC, AAC_F(end-start), totalNodes, dim, aacData, coord);
//    combineTimer.Stop();
}

void spawnRecursiveBuild(vector<BVHPrimitiveInfo> &buildData,
                         std::pair<uint32_t, uint32_t> *mortonCodes,
                         uint32_t start,
                         uint32_t end, uint32_t *totalNodes, int partitionBit,
                         TaskTree* &tt) {
    
    if (end-start == 0) {
        return;
    }
    
    if (MORTON_CODE_START-partitionBit > AAC_PARALLEL_GRAN) {
        AAC_Data *aac_data = new AAC_Data;
        float minSize = float(AAC_DELTA/2);
        int dataSize = 4*int(minSize * pow((end-start), 0.5-AAC_EPSILON/2) + 1e-5);
        aac_data->init(dataSize);
        
        AAC_DataCoord coord = AAC_DataCoord(0, dataSize);
        std::vector<BVHBuildNode*> *clusters = new std::vector<BVHBuildNode*>();
        
        recursiveBuildAAC(buildData, mortonCodes, start, end, totalNodes,
                          MORTON_CODE_START, clusters, aac_data, &coord);
        
        
        
        tt = new TaskTree();
        
        tt->bvhNodes = clusters;
        tt->aacData = aac_data;
        tt->aacDataC = coord;
        tt->numPrims = end-start;
        
        return;
    }
    
    
    uint32_t splitIdx = makePartition(mortonCodes, start, end, partitionBit);
    
    int newPartionBit = partitionBit - 1;
    uint32_t rightTotalnodes = 0;
    
    tt = new TaskTree();
    
    cilk_spawn spawnRecursiveBuild(buildData, mortonCodes, start, splitIdx, totalNodes, newPartionBit, (tt->leftNode));
    spawnRecursiveBuild(buildData, mortonCodes, splitIdx, end, &rightTotalnodes, newPartionBit, (tt->rightNode));
    cilk_sync;
    
    (*totalNodes) += rightTotalnodes;
}

void sequentialBuild(uint32_t *numPrims, uint32_t *totalNodes, int partitionBit,
                     std::vector<BVHBuildNode*> *clusterData,
                     AAC_Data *aacData, AAC_DataCoord *coord,
                     TaskTree* tt) {
    
    int dim = getDimForMorton3D(partitionBit);
    int newPartionBit = partitionBit - 1;
    
    if (!(tt->leftNode) && !(tt->rightNode)) {
        // push back clusterData
        *clusterData = *(tt->bvhNodes);
        
        AAC_Data *srcData = tt->aacData;
        AAC_DataCoord sdc = tt->aacDataC;
        
        Assert(coord->start + (sdc.end - sdc.start) <= coord->end);
        
        // fill in coordinates
        coord->end = coord->start + (sdc.end - sdc.start);
        
        // fill in aacData
        for (int i=sdc.start; i<sdc.end; i++)  {
            int destI = coord->start + (i - sdc.start);
            for (int j=sdc.start; j<i; j++) {
                int destJ = coord->start + (j - sdc.start);
                aacData->area[destI][destJ] = srcData->area[i][j];
            }
            
            aacData->minArea[destI] = srcData->minArea[i];
            aacData->minPos[destI] = srcData->minPos[i] + coord->start;
        }
        
        (*numPrims) += tt->numPrims;
        
        delete tt;
        
        return;
    } else if (!(tt->rightNode)) {
        sequentialBuild(numPrims, totalNodes, newPartionBit, clusterData, aacData, coord, tt->leftNode);
        delete tt;
        return;
    } else if (!(tt->leftNode)) {
        sequentialBuild(numPrims, totalNodes, newPartionBit, clusterData, aacData, coord, tt->rightNode);
        delete tt;
        return;
    }
    
    std::vector<BVHBuildNode*> leftC;
    std::vector<BVHBuildNode*> rightC;
    
    // do left
    AAC_DataCoord leftDC = AAC_DataCoord(*coord);
    
    sequentialBuild(numPrims, totalNodes, newPartionBit, &leftC, aacData, &leftDC, tt->leftNode);
    
    // do right
    AAC_DataCoord rightDC = AAC_DataCoord(*coord);
    rightDC.start = leftDC.end;
    sequentialBuild(numPrims, totalNodes, newPartionBit, &rightC, aacData, &rightDC, tt->rightNode);
    
    leftC.insert( leftC.end(), rightC.begin(), rightC.end() );
    
    aacDataMerge(leftC, aacData, leftDC, rightDC);
    
    *clusterData = combineCluster(leftC, AAC_F(*numPrims), totalNodes, dim, aacData, coord);
    delete tt;
}

void hybridRecursiveBuild(vector<BVHPrimitiveInfo> &buildData,
                         std::pair<uint32_t, uint32_t> *mortonCodes,
                         uint32_t start,
                         uint32_t end, uint32_t *totalNodes, int partitionBit,
                         TaskTree* &tt) {
    
    if (end-start == 0) {
//        Warning("hybridRecursiveBuild: NOTHING HERE");
        return;
    }
    
    if (MORTON_CODE_START-partitionBit > AAC_PARALLEL_GRAN) {
//        Warning("hybridRecursiveBuild: BASE CASE START");
        AAC_Data *aac_data = new AAC_Data;
        float minSize = float(AAC_DELTA/2);
        int dataSize = 4*int(minSize * pow((end-start), 0.5-AAC_EPSILON/2) + 1e-5);
        aac_data->init(dataSize);
        
        AAC_DataCoord coord = AAC_DataCoord(0, dataSize);
        std::vector<BVHBuildNode*> *clusters = new std::vector<BVHBuildNode*>();
        
        recursiveBuildAAC(buildData, mortonCodes, start, end, totalNodes,
                          MORTON_CODE_START, clusters, aac_data, &coord);
        
        
        
        tt = new TaskTree();
        
        tt->bvhNodes = clusters;
        tt->aacData = aac_data;
        tt->aacDataC = coord;
        tt->numPrims = end-start;
        
//        Warning("hybridRecursiveBuild: BASE CASE END");
        return;
    }
    
    
    uint32_t splitIdx = makePartition(mortonCodes, start, end, partitionBit);
    int dim = getDimForMorton3D(partitionBit);
    int newPartionBit = partitionBit - 1;
    uint32_t rightTotalnodes = 0;
    
    tt = new TaskTree();
    
//    Warning("SPAWN LEFT");
    cilk_spawn hybridRecursiveBuild(buildData, mortonCodes, start, splitIdx, totalNodes, newPartionBit, (tt->leftNode));
//    Warning("SPAWN RIGHT");
    hybridRecursiveBuild(buildData, mortonCodes, splitIdx, end, &rightTotalnodes, newPartionBit, (tt->rightNode));
    cilk_sync;
    
//    Warning("hybridRecursiveBuild: DONE SPAWNING");
    
    (*totalNodes) += rightTotalnodes;
    
    Assert((tt->rightNode) || (tt->leftNode));
    
    if (!(tt->rightNode)) {
//        Warning("hybridRecursiveBuild: CASE 1");
        std::vector<BVHBuildNode*> leftC = *(tt->leftNode->bvhNodes);
        (tt->bvhNodes) = new std::vector<BVHBuildNode*> (leftC);
        
//        return;
    } else if (!(tt->leftNode)) {
//        Warning("hybridRecursiveBuild: CASE 2");
        std::vector<BVHBuildNode*> rightC = *(tt->rightNode->bvhNodes);
        (tt->bvhNodes) = new std::vector<BVHBuildNode*> (rightC);
//        return;
    } else {
//        Warning("hybridRecursiveBuild: CASE 3");
        std::vector<BVHBuildNode*> leftC = *(tt->leftNode->bvhNodes);
        std::vector<BVHBuildNode*> rightC = *(tt->rightNode->bvhNodes);
        (tt->bvhNodes) = new std::vector<BVHBuildNode*> (leftC);
        tt->bvhNodes->insert( tt->bvhNodes->end(), rightC.begin(), rightC.end() );
    }
    
//    Warning("hybridRecursiveBuild: COMBINE CLUSTER START");
    *(tt->bvhNodes) = combineCluster(*(tt->bvhNodes), AAC_F(end-start), totalNodes, dim);
//    Warning("hybridRecursiveBuild: COMBINE CLUSTER END");
}



uint32_t BVHAccel::bvhDfs(BVHBuildNode* node, vector<BVHPrimitiveInfo> &buildData, vector<Reference<Primitive> > &orderedPrims, uint32_t *offset) {
    
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    uint32_t myOffset = (*offset)++;
    
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        uint32_t firstPrimOffset = orderedPrims.size();
        uint32_t primNum = buildData[node->firstPrimOffset].primitiveNumber;
        orderedPrims.push_back(primitives[primNum]);
        node->firstPrimOffset = firstPrimOffset;
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    } else {
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        bvhDfs(node->children[0], buildData, orderedPrims, offset);
        linearNode->secondChildOffset =
            bvhDfs(node->children[1], buildData, orderedPrims, offset);
    }
    
    delete node;
    return myOffset;
}

void BVHAccel::buildAAC(vector<BVHPrimitiveInfo> &buildData, uint32_t start,
                                       uint32_t end, uint32_t *totalNodes,
                                       vector<Reference<Primitive> > &orderedPrims) {
    Timer rootTimer;
    rootTimer.Start();
    Warning("buildAAC numprims == %u", end-start);

    std::pair<uint32_t, uint32_t> *mortonCodes = new std::pair<uint32_t, uint32_t>[end-start];
    
    BBox centerBBox = BBox();
    for (uint32_t i = start; i < end; ++i) {
        centerBBox = Union(centerBBox, buildData[i].centroid);
    }
    
    for (uint32_t i = start; i < end; ++i) {
        Point c = buildData[i].centroid;
        float newX = (c.x - centerBBox.pMin.x)/(centerBBox.pMax.x - centerBBox.pMin.x);
        float newY = (c.y - centerBBox.pMin.y)/(centerBBox.pMax.y - centerBBox.pMin.y);
        float newZ = (c.z - centerBBox.pMin.z)/(centerBBox.pMax.z - centerBBox.pMin.z);
        uint32_t mc = morton3D(newX, newY, newZ);
        
        mortonCodes[i-start] = std::make_pair(mc, i);
    }
    
    
    integerSort(mortonCodes, (long) (end-start));
    rootTimer.Stop();
    Warning("Elapsed time to sort morton codes: %.2f seconds\n", rootTimer.Time());
    rootTimer.Reset();
    std::sort(mortonCodes, mortonCodes+(end-start));
    
    float minSize = float(AAC_DELTA/2);
    int dataSize = 4*int(minSize*pow((end-start), 0.5-AAC_EPSILON/2)+1e-5);
    
    TaskTree *tt;
    rootTimer.Start();
    spawnRecursiveBuild(buildData, mortonCodes, start, end, totalNodes, MORTON_CODE_START, tt);
    rootTimer.Stop();
    Warning("Elapsed time to spawnRecursiveBuild: %.2f seconds\n", rootTimer.Time());
    rootTimer.Reset();
    
    uint32_t numPrims = 0;
    AAC_Data aac_data;
    rootTimer.Start();
    aac_data.init(dataSize);
    rootTimer.Stop();
    Warning("Elapsed time to init aac_data: %.2f seconds\n", rootTimer.Time());
    rootTimer.Reset();
    AAC_DataCoord coord = AAC_DataCoord(0, dataSize);
    std::vector<BVHBuildNode*> clusters;
    
    
    rootTimer.Start();
    
    sequentialBuild(&numPrims, totalNodes, MORTON_CODE_START, &clusters, &aac_data, &coord, tt);
    
    rootTimer.Stop();
    Warning("Elapsed time to sequentialBuild: %.2f seconds\n", rootTimer.Time());
    rootTimer.Reset();
    
    //HYBRID
//    TaskTree *tt;
//    rootTimer.Start();
//    Warning("START HYBRID BUILD");
//    hybridRecursiveBuild(buildData, mortonCodes, start, end, totalNodes, MORTON_CODE_START, tt);
//    rootTimer.Stop();
//    Warning("Elapsed time to spawnRecursiveBuild: %.2f seconds\n", rootTimer.Time());
//    rootTimer.Reset();
//    
//    std::vector<BVHBuildNode*> clusters = *(tt->bvhNodes);
//    
//    BVHBuildNode* root = combineCluster(clusters, 1, totalNodes, 2)[0];
    
    //HYBRID
    
//    AAC_Data aac_data;
//    int dataSize = end-start;
    
//    rootTimer.Start();
//    aac_data.init(dataSize);
//    rootTimer.Stop();
//    Warning("Elapsed time to init aac_data: %.2f seconds\n", rootTimer.Time());
//    rootTimer.Reset();
//    
//    AAC_DataCoord coord = AAC_DataCoord(0, dataSize);
//    std::vector<BVHBuildNode*> clusters;
//    
//    rootTimer.Start();
//    recursiveBuildAAC(buildData, mortonCodes, start, end, totalNodes,
//                      MORTON_CODE_START, &clusters, &aac_data, &coord);
    
//    rootTimer.Stop();
//    Warning("Elapsed time to recursive build: %.2f seconds\n", rootTimer.Time());
//    Warning("Elapsed time to combine clusters: %.2f seconds\n", combineTimer.Time());
//    Warning("Elapsed time to make Leafs: %.2f seconds\n", leafsTimer.Time());
//    Warning("Elapsed time to merge data: %.2f seconds\n", dataMergeTimer.Time());
//    Warning("Elapsed time to append data: %.2f seconds\n", appendTimer.Time());
//    Warning("Elapsed time to base case: %.2f seconds\n", baseCaseTimer.Time());
//    rootTimer.Reset();
    
//    rootTimer.Start();
//    BVHBuildNode* root = combineCluster(clusters, 1, totalNodes, 2, &aac_data, &coord)[0];
//    *totalNodes = faketotalNodes;
    // use 2 for dim because start with z axis
    BVHBuildNode* root = combineCluster(clusters, 1, totalNodes, 2, &aac_data, &coord)[0];
//    rootTimer.Stop();
//    Warning("Elapsed time to combine to root cluster: %.2f seconds\n", rootTimer.Time());
//    rootTimer.Reset();
    
    delete[] mortonCodes;
    
    rootTimer.Start();
    // Compute representation of depth-first traversal of BVH tree
    nodes = AllocAligned<LinearBVHNode>(*totalNodes);
    for (uint32_t i = 0; i < (*totalNodes); ++i)
        new (&nodes[i]) LinearBVHNode;
    uint32_t offset = 0;
    
    bvhDfs(root, buildData, orderedPrims, &offset);
    rootTimer.Stop();
    Warning("Elapsed time to combine to do DFS: %.2f seconds\n", rootTimer.Time());
    
    primitives.swap(orderedPrims);
}


BVHBuildNode *BVHAccel::recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start,
        uint32_t end, uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims) {
    Assert(start != end);
    (*totalNodes)++;
    BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
    // Compute bounds of all primitives in BVH node
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
    uint32_t nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _BVHBuildNode_
        uint32_t firstPrimOffset = orderedPrims.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        uint32_t mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // If nPrimitives is no greater than maxPrimsInNode,
            // then all the nodes can be stored in a compact bvh node.
            if (nPrimitives <= maxPrimsInNode) {
                // Create leaf _BVHBuildNode_
                uint32_t firstPrimOffset = orderedPrims.size();
                for (uint32_t i = start; i < end; ++i) {
                    uint32_t primNum = buildData[i].primitiveNumber;
                    orderedPrims.push_back(primitives[primNum]);
                }
                node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                return node;
            }
            else {
                // else if nPrimitives is greater than maxPrimsInNode, we
                // need to split it further to guarantee each node contains
                // no more than maxPrimsInNode primitives.
                node->InitInterior(dim,
                                   recursiveBuild(buildArena, buildData, start, mid,
                                                  totalNodes, orderedPrims),
                                   recursiveBuild(buildArena, buildData, mid, end,
                                                  totalNodes, orderedPrims));
                return node;
            }
        }

        // Partition primitives based on _splitMethod_
        switch (splitMethod) {
        case SPLIT_MIDDLE: {
            // Partition primitives through node's midpoint
            float pmid = .5f * (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]);
            BVHPrimitiveInfo *midPtr = std::partition(&buildData[start],
                                                      &buildData[end-1]+1,
                                                      CompareToMid(dim, pmid));
            mid = midPtr - &buildData[0];
            if (mid != start && mid != end)
                // for lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall through
                // to SPLIT_EQUAL_COUNTS
                break;
        }
        case SPLIT_EQUAL_COUNTS: {
            // Partition primitives into equally-sized subsets
            mid = (start + end) / 2;
            std::nth_element(&buildData[start], &buildData[mid],
                             &buildData[end-1]+1, ComparePoints(dim));
            break;
        }
        case SPLIT_SAH: default: {
            // Partition primitives using approximate SAH
            if (nPrimitives <= 4) {
                // Partition primitives into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                                 &buildData[end-1]+1, ComparePoints(dim));
            }
            else {
                // Allocate _BucketInfo_ for SAH partition buckets
                const int nBuckets = 12;
                struct BucketInfo {
                    BucketInfo() { count = 0; }
                    int count;
                    BBox bounds;
                };
                BucketInfo buckets[nBuckets];

                // Initialize _BucketInfo_ for SAH partition buckets
                for (uint32_t i = start; i < end; ++i) {
                    int b = nBuckets *
                        ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                         (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                    if (b == nBuckets) b = nBuckets-1;
                    Assert(b >= 0 && b < nBuckets);
                    buckets[b].count++;
                    buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
                }

                // Compute costs for splitting after each bucket
                float cost[nBuckets-1];
                for (int i = 0; i < nBuckets-1; ++i) {
                    BBox b0, b1;
                    int count0 = 0, count1 = 0;
                    for (int j = 0; j <= i; ++j) {
                        b0 = Union(b0, buckets[j].bounds);
                        count0 += buckets[j].count;
                    }
                    for (int j = i+1; j < nBuckets; ++j) {
                        b1 = Union(b1, buckets[j].bounds);
                        count1 += buckets[j].count;
                    }
                    cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) /
                              bbox.SurfaceArea();
                }

                // Find bucket to split at that minimizes SAH metric
                float minCost = cost[0];
                uint32_t minCostSplit = 0;
                for (int i = 1; i < nBuckets-1; ++i) {
                    if (cost[i] < minCost) {
                        minCost = cost[i];
                        minCostSplit = i;
                    }
                }

                // Either create leaf or split primitives at selected SAH bucket
                if (nPrimitives > maxPrimsInNode ||
                    minCost < nPrimitives) {
                    BVHPrimitiveInfo *pmid = std::partition(&buildData[start],
                        &buildData[end-1]+1,
                        CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
                    mid = pmid - &buildData[0];
                }
                
                else {
                    // Create leaf _BVHBuildNode_
                    uint32_t firstPrimOffset = orderedPrims.size();
                    for (uint32_t i = start; i < end; ++i) {
                        uint32_t primNum = buildData[i].primitiveNumber;
                        orderedPrims.push_back(primitives[primNum]);
                    }
                    node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                    return node;
                }
            }
            break;
        }
        }
        node->InitInterior(dim,
                           recursiveBuild(buildArena, buildData, start, mid,
                                          totalNodes, orderedPrims),
                           recursiveBuild(buildArena, buildData, mid, end,
                                          totalNodes, orderedPrims));
    }
    return node;
}


uint32_t BVHAccel::flattenBVHTree(BVHBuildNode *node, uint32_t *offset) {
    LinearBVHNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
    uint32_t myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
        flattenBVHTree(node->children[0], offset);
        linearNode->secondChildOffset = flattenBVHTree(node->children[1],
                                                       offset);
    }
    
    if ((splitMethod == SPLIT_AAC) && node) {
        delete node;
    }
    
    return myOffset;
}


BVHAccel::~BVHAccel() {
    FreeAligned(nodes);
}


bool BVHAccel::Intersect(const Ray &ray, Intersection *isect) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTION_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    bool hit = false;
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    // Follow ray through BVH nodes to find primitive intersections
    uint32_t todoOffset = 0, nodeNum = 0;
    uint32_t todo[64];
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        // Check ray against BVH node
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            if (node->nPrimitives > 0) {
                // Intersect ray with primitives in leaf BVH node
                PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                for (uint32_t i = 0; i < node->nPrimitives; ++i)
                {
                    PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->Intersect(ray, isect))
                    {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        hit = true;
                    }
                    else {
                        PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                   }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                // Put far BVH node on _todo_ stack, advance to near node
                PBRT_BVH_INTERSECTION_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTION_FINISHED();
    return hit;
}


bool BVHAccel::IntersectP(const Ray &ray) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTIONP_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    uint32_t todo[64];
    uint32_t todoOffset = 0, nodeNum = 0;
    while (true) {
        const LinearBVHNode *node = &nodes[nodeNum];
        if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
            // Process BVH node _node_ for traversal
            if (node->nPrimitives > 0) {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearBVHNode *>(node));
                  for (uint32_t i = 0; i < node->nPrimitives; ++i) {
                    PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    if (primitives[node->primitivesOffset+i]->IntersectP(ray)) {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
                        return true;
                    }
                else {
                        PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
                    }
                }
                if (todoOffset == 0) break;
                nodeNum = todo[--todoOffset];
            }
            else {
                PBRT_BVH_INTERSECTIONP_TRAVERSED_INTERIOR_NODE(const_cast<LinearBVHNode *>(node));
                if (dirIsNeg[node->axis]) {
                   /// second child first
                   todo[todoOffset++] = nodeNum + 1;
                   nodeNum = node->secondChildOffset;
                }
                else {
                   todo[todoOffset++] = node->secondChildOffset;
                   nodeNum = nodeNum + 1;
                }
            }
        }
        else {
            if (todoOffset == 0) break;
            nodeNum = todo[--todoOffset];
        }
    }
    PBRT_BVH_INTERSECTIONP_FINISHED();
    return false;
}


BVHAccel *CreateBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    string splitMethod = ps.FindOneString("splitmethod", "sah");
    uint32_t maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return new BVHAccel(prims, maxPrimsInNode, splitMethod);
}


