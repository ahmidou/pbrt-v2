
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_ACCELERATORS_BVH_H
#define PBRT_ACCELERATORS_BVH_H

// accelerators/bvh.h*
#include "pbrt.h"
#include "primitive.h"
struct BVHBuildNode;

// BVHAccel Forward Declarations
struct BVHPrimitiveInfo;
struct LinearBVHNode;

struct AAC_Data {
    float **area, *minArea;
    //    int *label, *minPos, *minLabel;
    int *minPos;
    int size;
    int inited;
    
    void init(int size) {
        this->size = size;
        //        label = new int [size]; //idx in clusters
        minArea = new float [size]; // min area for each
        minPos = new int [size]; // min position of each
        //        minLabel = new int [size]; //
        area = new float* [size];
        for (int i=0; i<size; i++)
            area[i] = new float [size];
        
        inited = 1;
    }
    
    AAC_Data() {
        inited = 0;
    }
    
    ~AAC_Data() {
        if (inited) {
            //        delete[] label;
            delete[] minArea;
            delete[] minPos;
            //        delete[] minLabel;
            for (int i=0; i<size; i++)
                delete[] area[i];
            delete[] area;
            inited = 0;
        }
    }
};

struct AAC_DataCoord {
    int start;
    int end;
    
    AAC_DataCoord() {}
    AAC_DataCoord(int start, int end) : start(start), end(end) {}
    AAC_DataCoord(AAC_DataCoord const &c) : start(c.start), end(c.end) {}
};

struct TaskTree {
    std::vector<BVHBuildNode*> *bvhNodes;
    AAC_Data *aacData;
    AAC_DataCoord aacDataC;
    uint32_t numPrims;
    
    TaskTree *leftNode, *rightNode;
    
    TaskTree():leftNode(NULL), rightNode(NULL) {}
    
    TaskTree(std::vector<BVHBuildNode*> *bvhNodes, AAC_Data *aacData, AAC_DataCoord aacDataC, uint32_t numPrims):
    bvhNodes(bvhNodes), aacData(aacData), aacDataC(aacDataC), numPrims(numPrims), leftNode(NULL), rightNode(NULL) {}
    
    
    ~TaskTree() {
        if (aacData)
            delete aacData;
        if (bvhNodes)
            delete bvhNodes;
    }
} ;


// BVHAccel Declarations
class BVHAccel : public Aggregate {
public:
    // BVHAccel Public Methods
    BVHAccel(const vector<Reference<Primitive> > &p, uint32_t maxPrims = 1,
             const string &sm = "sah");
    BBox WorldBound() const;
    bool CanIntersect() const { return true; }
    ~BVHAccel();
    bool Intersect(const Ray &ray, Intersection *isect) const;
    bool IntersectP(const Ray &ray) const;
private:
    // BVHAccel Private Methods
    uint32_t bvhDfs(BVHBuildNode* node, vector<BVHPrimitiveInfo> &buildData, vector<Reference<Primitive> > &orderedPrims, uint32_t *offset);
    void buildAAC(vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end,
                                 uint32_t *totalNodes, vector<Reference<Primitive> > &orderedPrims);
    
    
    
    BVHBuildNode *recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end,
        uint32_t *totalNodes, vector<Reference<Primitive> > &orderedPrims);
    uint32_t flattenBVHTree(BVHBuildNode *node, uint32_t *offset);

    // BVHAccel Private Data
    uint32_t maxPrimsInNode;
    enum SplitMethod { SPLIT_MIDDLE, SPLIT_EQUAL_COUNTS, SPLIT_AAC, SPLIT_SAH };
    SplitMethod splitMethod;
    vector<Reference<Primitive> > primitives;
    LinearBVHNode *nodes;
};


BVHAccel *CreateBVHAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps);

#endif // PBRT_ACCELERATORS_BVH_H
