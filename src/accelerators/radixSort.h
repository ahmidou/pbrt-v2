// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef A_RADIX_INCLUDED
#define A_RADIX_INCLUDED

#include <iostream>
#include <algorithm>
#include <math.h>
#include <climits>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "transpose.h"
using namespace std;

typedef int intT;
typedef uint32_t uintT;
#define INT_T_MAX INT_MAX
#define UINT_T_MAX UINT_MAX

#define _BSIZE 2048
#define _SCAN_LOG_BSIZE 10
#define _SCAN_BSIZE (1 << _SCAN_LOG_BSIZE)

#define blocked_for(_i, _s, _e, _bsize, _body)  {	\
intT _ss = _s;					\
intT _ee = _e;					\
intT _n = _ee-_ss;					\
intT _l = nblocks(_n,_bsize);			\
cilk_for (intT _i = 0; _i < _l; _i++) {		\
intT _s = _ss + _i * (_bsize);			\
intT _e = min(_s + (_bsize), _ee);			\
_body						\
}						\
}

#define newA(__E,__n) (__E*) malloc((__n)*sizeof(__E))
#define nblocks(_n,_bsize) (1 + ((_n)-1)/(_bsize))

/************************** UTILS BEGIN **************************/

//namespace utils {

    // returns the log base 2 rounded up (works on ints or longs or unsigned versions)
    template <class T>
    static int log2Up(T i) {
        int a=0;
        T b=i-1;
        while (b > 0) {b = b >> 1; a++;}
        return a;
    }


    template <class E>
    struct identityF { E operator() (const E& x) {return x;}};


    template <class E>
    struct addF { E operator() (const E& a, const E& b) const {return a+b;}};


    template <class E>
    struct maxF { E operator() (const E& a, const E& b) const {return (a>b) ? a : b;}};


    template <class E>
    struct minF { E operator() (const E& a, const E& b) const {return (a<b) ? a : b;}};


    template <class E1, class E2>
    struct firstF {E1 operator() (std::pair<E1,E2> a) {return a.first;} };
    
//}

/************************** UTILS END **************************/

/************************** SEQUENCE OPS BEGIN **************************/

//namespace sequence {

    template <class ET, class intT>
    struct getA {
        ET* A;
        getA(ET* AA) : A(AA) {}
        ET operator() (intT i) {return A[i];}
    };

    template <class IT, class OT, class intT, class F>
    struct getAF {
        IT* A;
        F f;
        getAF(IT* AA, F ff) : A(AA), f(ff) {}
        OT operator () (intT i) {return f(A[i]);}
    };

    template <class OT, class intT, class F, class G>
    OT reduceSerial(intT s, intT e, F f, G g) {
        OT r = g(s);
        for (intT j=s+1; j < e; j++) r = f(r,g(j));
        return r;
    }

    template <class OT, class intT, class F, class G>
    OT reduce(intT s, intT e, F f, G g) {
        intT l = nblocks(e-s, _SCAN_BSIZE);
        if (l <= 1) return reduceSerial<OT>(s, e, f , g);
        OT *Sums = newA(OT,l);
        blocked_for (i, s, e, _SCAN_BSIZE,
                     Sums[i] = reduceSerial<OT>(s, e, f, g););
        OT r = reduce<OT>((intT) 0, l, f, getA<OT,intT>(Sums));
        free(Sums);
        return r;
    }

    template <class OT, class intT, class F>
    OT reduce(OT* A, intT n, F f) {
        return reduce<OT>((intT)0,n,f,getA<OT,intT>(A));
    }

    // g is the map function (applied to each element)
    // f is the reduce function
    // need to specify OT since it is not an argument
    template <class OT, class IT, class intT, class F, class G>
    OT mapReduce(IT* A, intT n, F f, G g) {
        return reduce<OT>((intT) 0,n,f,getAF<IT,OT,intT,G>(A,g));
    }

    template <class ET, class intT, class F, class G>
    ET scanSerial(ET* Out, intT s, intT e, F f, G g, ET zero, bool inclusive, bool back) {
        ET r = zero;
        
        if (inclusive) {
            if (back) for (long i = e-1; i >= s; i--) Out[i] = r = f(r,g(i));
            else for (intT i = s; i < e; i++) Out[i] = r = f(r,g(i));
        } else {
            if (back)
                for (long i = e-1; i >= s; i--) {
                    ET t = g(i);
                    Out[i] = r;
                    r = f(r,t);
                }
            else
                for (intT i = s; i < e; i++) {
                    ET t = g(i);
                    Out[i] = r;
                    r = f(r,t);
                }
        }
        return r;
    }

    template <class ET, class intT, class F>
    ET scanSerial(ET *In, ET* Out, intT n, F f, ET zero) {
        return scanSerial(Out, (intT) 0, n, f, getA<ET,intT>(In), zero, false, false);
    }

    // back indicates it runs in reverse direction
    template <class ET, class intT, class F, class G>
    ET scan(ET* Out, intT s, intT e, F f, G g,  ET zero, bool inclusive, bool back) {
        intT n = e-s;
        intT l = nblocks(n,_SCAN_BSIZE);
        if (l <= 2) return scanSerial(Out, s, e, f, g, zero, inclusive, back);
        ET *Sums = newA(ET,nblocks(n,_SCAN_BSIZE));
        blocked_for (i, s, e, _SCAN_BSIZE,
                     Sums[i] = reduceSerial<ET>(s, e, f, g););
        ET total = scan(Sums, (intT) 0, l, f, getA<ET,intT>(Sums), zero, false, back);
        blocked_for (i, s, e, _SCAN_BSIZE,
                     scanSerial(Out, s, e, f, g, Sums[i], inclusive, back););
        free(Sums);
        return total;
    }

    template <class ET, class intT, class F>
    ET scan(ET *In, ET* Out, intT n, F f, ET zero) {
        return scan(Out, (intT) 0, n, f, getA<ET,intT>(In), zero, false, false);}

    template <class ET, class intT, class F>
    ET scanIBack(ET *In, ET* Out, intT n, F f, ET zero) {
        return scan(Out, (intT) 0, n, f, getA<ET,intT>(In), zero, true, true);}

//}

/************************** SEQUENCE OPS END **************************/

namespace intSort {
    
    // Cannot be greater than 8 without changing definition of bIndexT
    //    from unsigned char to unsigned int (or unsigned short)
#define MAX_RADIX 8
#define BUCKETS 256    // 1 << MAX_RADIX
    
    
    // a type that must hold MAX_RADIX bits
    typedef unsigned char bIndexT;
    
    template <class E, class F, class bint>
    void radixBlock(E* A, E* B, bIndexT *Tmp,
                    bint counts[BUCKETS], bint offsets[BUCKETS],
                    bint Boffset, long n, long m, F extract) {
        
        for (long i = 0; i < m; i++)  counts[i] = 0;
        for (long j = 0; j < n; j++) {
            bint k = Tmp[j] = extract(A[j]);
            counts[k]++;
        }
        bint s = Boffset;
        for (long i = 0; i < m; i++) {
            s += counts[i];
            offsets[i] = s;
        }
        for (long j = n-1; j >= 0; j--) {
            bint x =  --offsets[Tmp[j]];
            B[x] = A[j];
        }
    }
    
    template <class E, class F, class bint>
    void radixStepSerial(E* A, E* B, bIndexT *Tmp, bint buckets[BUCKETS],
                         long n, long m, F extract) {
        radixBlock(A, B, Tmp, buckets, buckets, (bint)0, n, m, extract);
        for (long i=0; i < n; i++) A[i] = B[i];
        return;
    }
    
    // A is the input and sorted output (length = n)
    // B is temporary space for copying data (length = n)
    // Tmp is temporary space for extracting the bytes (length = n)
    // BK is an array of bucket sets, each set has BUCKETS integers
    //    it is used for temporary space for bucket counts and offsets
    // numBK is the length of BK (number of sets of buckets)
    // the first entry of BK is also used to return the offset of each bucket
    // m is the number of buckets per set (m <= BUCKETS)
    // extract is a function that extract the appropriate bits from A
    //  it must return a non-negative integer less than m
    template <class E, class F, class bint>
    void radixStep(E* A, E* B, bIndexT *Tmp, bint (*BK)[BUCKETS],
                   long numBK, long n, long m, bool top, F extract) {
        
        
        // need 3 bucket sets per block
        long expand = (sizeof(E)<=4) ? 64 : 32;
        long blocks = min(numBK/3,(1+n/(BUCKETS*expand)));
        
        if (blocks < 2) {
            radixStepSerial(A, B, Tmp, BK[0], n, m, extract);
            return;
        }
        long nn = (n+blocks-1)/blocks;
        bint* cnts = (bint*) BK;
        bint* oA = (bint*) (BK+blocks);
        bint* oB = (bint*) (BK+2*blocks);
        
        _Pragma("cilk grainsize = 1")
        cilk_for(long i=0; i < blocks; i++) {
            bint od = i*nn;
            long nni = min(max<long>(n-od,0),nn);
            radixBlock(A+od, B, Tmp+od, cnts + m*i, oB + m*i, od, nni, m, extract);
        }
        
        transpose<bint,bint>(cnts, oA).trans(blocks, m);
        
        long ss;
        if (top)
            ss = scan(oA, oA, blocks*m, addF<bint>(),(bint)0);
        else
            ss = scanSerial(oA, oA, blocks*m, addF<bint>(),(bint)0);
        //myAssert(ss == n, "radixStep: sizes don't match");
        
        blockTrans<E,bint>(B, A, oB, oA, cnts).trans(blocks, m);
        
        // put the offsets for each bucket in the first bucket set of BK
        for (long j = 0; j < m; j++) BK[0][j] = oA[j*blocks];
    }
    
    // a function to extract "bits" bits starting at bit location "offset"
    template <class E, class F>
    struct eBits {
        F _f;  long _mask;  long _offset;
        eBits(long bits, long offset, F f): _mask((1<<bits)-1),
        _offset(offset), _f(f) {}
        long operator() (E p) {return _mask&(_f(p)>>_offset);}
    };
    
    // Radix sort with low order bits first
    template <class E, class F, class bint>
    void radixLoopBottomUp(E *A, E *B, bIndexT *Tmp, bint (*BK)[BUCKETS],
                           long numBK, long n, long bits, bool top, F f) {
        long rounds = 1+(bits-1)/MAX_RADIX;
        long rbits = 1+(bits-1)/rounds;
        long bitOffset = 0;
        while (bitOffset < bits) {
            if (bitOffset+rbits > bits) rbits = bits-bitOffset;
            radixStep(A, B, Tmp, BK, numBK, n, 1 << rbits, top,
                      eBits<E,F>(rbits,bitOffset,f));
            bitOffset += rbits;
        }
    }
    
    // Radix sort with high order bits first
    template <class E, class F, class bint>
    void radixLoopTopDown(E *A, E *B, bIndexT *Tmp, bint (*BK)[BUCKETS],
                          long numBK, long n, long bits, F f) {
        if (n == 0) return;
        if (bits <= MAX_RADIX) {
            radixStep(A, B, Tmp, BK, numBK, n, ((long) 1) << bits, true,
                      eBits<E,F>(bits,0,f));
        } else if (numBK >= BUCKETS+1) {
            radixStep(A, B, Tmp, BK, numBK, n, (long) BUCKETS, true,
                      eBits<E,F>(MAX_RADIX,bits-MAX_RADIX,f));
            bint* offsets = BK[0];
            long remain = numBK - BUCKETS - 1;
            float y = remain / (float) n;
            cilk_for (int i = 0; i < BUCKETS; i++) {
                long segOffset = offsets[i];
                long segNextOffset = (i == BUCKETS-1) ? n : offsets[i+1];
                long segLen = segNextOffset - segOffset;
                long blocksOffset = ((long) floor(segOffset * y)) + i + 1;
                long blocksNextOffset = ((long) floor(segNextOffset * y)) + i + 2;
                long blockLen = blocksNextOffset - blocksOffset;
                radixLoopTopDown(A + segOffset, B + segOffset, Tmp + segOffset,
                                 BK + blocksOffset, blockLen, segLen,
                                 bits-MAX_RADIX, f);
            }
        } else {
            radixLoopBottomUp(A, B, Tmp, BK, numBK, n, bits, false, f);
        }
    }
    
    template <class E>
    long iSortSpace(long n) {
        long esize = (n >= INT_MAX) ? sizeof(long) : sizeof(int);
        long numBK = 1+n/(BUCKETS*8);
        return sizeof(E)*n + esize*n + esize*BUCKETS*numBK;
    }
    
    // Sorts the array A, which is of length n.
    // Function f maps each element into an integer in the range [0,m)
    // If bucketOffsets is not NULL then it should be an array of length m
    // The offset in A of each bucket i in [0,m) is placed in location i
    //   such that for i < m-1, offsets[i+1]-offsets[i] gives the number
    //   of keys=i.   For i = m-1, n-offsets[i] is the number.
    template <class bint, class E, class F, class oint>
    void iSortX(E *A, oint* bucketOffsets, long n, long m, bool bottomUp,
                char* tmpSpace, F f) {
        typedef bint bucketsT[BUCKETS];
        
        long bits = log2Up(m);
        long numBK = 1+n/(BUCKETS*8);
        
        // the temporary space is broken into 3 parts: B, Tmp and BK
        E *B = (E*) tmpSpace;
        long Bsize =sizeof(E)*n;
        bIndexT *Tmp = (bIndexT*) (tmpSpace+Bsize); // one byte per item
        long tmpSize = sizeof(bIndexT)*n;
        bucketsT *BK = (bucketsT*) (tmpSpace+Bsize+tmpSize);
        if (bits <= MAX_RADIX) {
            radixStep(A, B, Tmp, BK, numBK, n, (long) 1 << bits, true,
                      eBits<E,F>(bits,0,f));
            if (bucketOffsets != NULL) {
                cilk_for(long i = 0; i < m; i++)
                bucketOffsets[i] = BK[0][i];
            }
            return;
        } else if (bottomUp)
            radixLoopBottomUp(A, B, Tmp, BK, numBK, n, bits, true, f);
        else
            radixLoopTopDown(A, B, Tmp, BK, numBK, n, bits, f);
        if (bucketOffsets != NULL) {
            {cilk_for(long i=0; i < m; i++) bucketOffsets[i] = n;}
            {cilk_for(long i=0; i < n-1; i++) {
                long v = f(A[i]);
                long vn = f(A[i+1]);
                if (v != vn) bucketOffsets[vn] = i+1;
            }}
            bucketOffsets[f(A[0])] = 0;
            scanIBack(bucketOffsets, bucketOffsets, m,
                                minF<oint>(), (oint) n);
        }
    }
    
    template <class E, class F, class oint>
    void iSort(E *A, oint* bucketOffsets, long n, long m, bool bottomUp,
               char* tmpSpace, F f) {
        // if n fits in 32 bits then use unsigned ints for bucket counts
        // otherwise use unsigned longs
        // Doesn't make much difference in performance
        if (n < UINT_MAX)
            iSortX<unsigned int>(A, bucketOffsets, n, m, bottomUp, tmpSpace, f);
        else iSortX<unsigned long>(A, bucketOffsets, n, m, bottomUp, tmpSpace, f);
    }
    
    // THE REST ARE JUST SPECIAL CASES
    
    template <class E, class F, class oint>
    void iSort(E *A, oint* bucketOffsets, long n, long m, bool bottomUp, F f) {
        long x = iSortSpace<E>(n);
        char* s = (char*) malloc(x);
        iSort(A, bucketOffsets, n, m, bottomUp, s, f);
        free(s);
    }
    
    template <class E, class F, class oint>
    void iSort(E *A, oint* bucketOffsets, long n, long m, F f) {
        iSort(A, bucketOffsets, n, m, false, f);}
    
    // A version that uses a NULL bucketOffset
    template <class E, class Func>
    void iSort(E *A, long n, long m, Func f) { 
        iSort(A, (unsigned long*) NULL, n, m, false, f);
    }
    
    template <class E, class Func>
    void iSort(E *A, long n, long m, char* s, Func f) { 
        iSort(A, (unsigned long*) NULL, n, m, false, s, f);
    }
    
    template <class E, class F>
    void iSortBottomUp(E *A, long n, long m, F f) {
        iSort(A, (unsigned long*) NULL, n, m, true, f);
    }
};

//static void integerSort(uintT *A, long n) {
//    long maxV = sequence::reduce(A, n, maxF<uintT>());
//    intSort::iSort(A, n, maxV+1,  identityF<uintT>());
//}
//
//static void integerSort(uintT *A, long n, char* s) {
//    long maxV = reduce(A,n,maxF<uintT>());
//    intSort::iSort(A, n, maxV+1, s, identityF<uintT>());
//}

template <class T>
void integerSort(pair<uintT,T> *A, long n) {
    long maxV = mapReduce<uintT>(A,n,maxF<uintT>(),
                                           firstF<uintT,T>());
    intSort::iSort(A, n, maxV+1, firstF<uintT,T>());
}

template <class T>
void integerSort(pair<uintT,T> *A, long n, char* s) {
    long maxV = mapReduce<uintT>(A,n,maxF<uintT>(),
                                           firstF<uintT,T>());
    intSort::iSort(A, n, maxV+1, s, firstF<uintT,T>());
}

#endif

