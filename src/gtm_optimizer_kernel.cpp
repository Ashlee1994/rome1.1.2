#include "./time.h"


namespace GTM_Optimizer_Kernel {


#define __assume_aligned(P,A)


typedef __int64 INDEXSZ;

static const INDEXSZ SCALE = 10;
static const INDEXSZ K = 10000  /SCALE;
static const INDEXSZ M = 1000   /SCALE;
static const INDEXSZ N = 500000 /SCALE;
static const INDEXSZ D = 1;                                 // for our purposes, this is irrelevant


// Loops below are
// Computes
//      mkl_solve       MMM     Believe this is right       1e9
//      getWd
//                      KN                                  5e9
//                      KMM                                 1e9
//                      KMM
// Space
//                      KN                                  5e9
//                      KM                                  1e7
//


// This is from the DanaFaber gtm code

static size_t mallocAllocated = 0;
static double* mallocDoublesBase(size_t len) {
    mallocAllocated += sizeof(double)*len;
    //std::cerr << " Allocating " << sizeof(double)*len << " bytes," << "   total is now " << mallocAllocated << " bytes" << std::endl;
    if (mallocAllocated > 2e9) {
        std::cerr << " Allocating " << sizeof(double)*len << " bytes," << " total is now " << mallocAllocated << " bytes" << std::endl;
		EXIT_ABNORMALLY;
    }
    return (double*)aMalloc(sizeof(double)*len,64);
}
static void freeDoubles(double* ptr, size_t len) {
    aFree(ptr);
    mallocAllocated -= sizeof(double)*len;
    //std::cerr << " Freeing " << sizeof(double)*len << " bytes" << std::endl;;
}

#undef mallocDoubles
static double* mallocDoubles(size_t len) {
    auto p = mallocDoublesBase(len);
    for (size_t i = 0; i < len; i++) p[i] = double((i+1)*1000 + 666);       // catches uninitialized case
    return p;
}

static double* mallocSetDoubles(size_t len) {                                       // use when data is initialized in the original algorithm
    auto p = mallocDoublesBase(len);
    for (size_t i = 0; i < len; i++) p[i] = double(i % 16);
    return p;
}

//#define BOUNDSCHECK
#ifdef BOUNDSCHECK
class DoubleVector {
public:
    DoubleVector(size_t size) : v(size) {
        for (size_t i = 0; i < size; i++) v[i] = double(i);
    }
    double operator[](INDEXSZ i) const { return v[i]; }
    double* ptr() { return &v[0]; }
protected:
    std::vector<double> v;
};
double* ptr(DoubleVector& dv) { return dv.ptr(); }

class WritableDoubleVector : public DoubleVector {
public:
    WritableDoubleVector(size_t size) : DoubleVector(size) {}
    double& operator[](INDEXSZ i) { return v[i]; }
};

#else
typedef double* DoubleVector;
typedef double* WritableDoubleVector;
double* ptr(DoubleVector& dv) { return dv; }
#endif

#define GLOBALS \
    ELT(CTFvalue, D*N)  \
    ELT(R       , K*N)  \
    ELT(TReal   , D*N)  \
    ELT(TImag   , D*N)  \
    ELT(PHI     , K*M)  // end of macro

#ifdef BOUNDSCHECK
#define ELT(NAME,SIZE) const DoubleVector NAME(SIZE);
#else
#define ELT(NAME,SIZE) const DoubleVector NAME = mallocSetDoubles(SIZE);
#endif
    GLOBALS
#undef ELT

const double AverageAlpha = 3.0;
const double AverageBeta  = 2.0;

auto const ipiv = (int*)aMalloc(sizeof(int)*M,64);

double mkl_solve_elapsed;
void mkl_solve(
    DoubleVector& A,
    int N,              // Caller uses M, not N from the above numbers
    DoubleVector& b,
    int M) {            // Caller uses 2, not M from the above numbers

    double mkl_solve_elapsed_start = dtime();
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR,N,N,ptr(A),N,ipiv);
    LAPACKE_dgetrs(LAPACK_ROW_MAJOR,'N',N,M,ptr(A),N,ipiv,ptr(b),M);
    mkl_solve_elapsed = dtime() - mkl_solve_elapsed_start;
}

const INDEXSZ MM      = M*M;
const INDEXSZ Aoffset = MM + 4096;		// To put on different pages
const INDEXSZ boffset = M  + 4096;		// To put on different pages


void getWd0(
    const int d,
    WritableDoubleVector& A,
    WritableDoubleVector& b,
    int maxthreads,
    WritableDoubleVector& AA,
    WritableDoubleVector& bRealArray, WritableDoubleVector& bImagArray)
{
    __assume_aligned(WReal,64);
    __assume_aligned(WImag,64);
    __assume_aligned(CTFvalue,64);
    __assume_aligned(TReal,64);
    __assume_aligned(TImag,64);
    __assume_aligned(PHI,64);
    __assume_aligned(R,64);
    __assume_aligned(AA,64);
    __assume_aligned(A,64);
    __assume_aligned(b,64);
    __assume_aligned(bRealArray,64);
    __assume_aligned(bImagArray,64);

    const size_t usedCapacity = size_t(maxthreads);
    bool* used = vNew(bool, usedCapacity);
    for (int i = 0; i < usedCapacity; i++) used[i] = false;

    #pragma omp parallel for schedule(dynamic)
    for(int k = 0; k < K; k++)      // K is only about 1000
    {
        const int tid = omp_get_thread_num();
        assert(tid < maxthreads);
        used[tid] = true;

        double pSA = 0;
        double pSb1 = 0;
        double pSb2 = 0;
        double CTFVdNn;
        double RCTFV;

        // Each thread doing unit stride array accesses
        for(int n = 0; n < N; n++)   // N is massive
        {
            // HOT
            CTFVdNn = CTFvalue  [INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)];
            RCTFV   = R         [INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)]*CTFVdNn;
            pSb1  += RCTFV*TReal[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)]; //partialSumsb1[k] += RCTFV*TReal[d*N+n];//AverageBeta
            pSb2  += RCTFV*TImag[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)];
            pSA   += RCTFV*CTFVdNn;
        }

        int offsetA = tid*Aoffset;
        int offsetb = tid*boffset;

        double PHIM;
        for(int m = 0;m < M;m++)   // M is about 600
        {  //for each A's rows
            PHIM = PHI[k*M+m];

            // Sum this product separately for each thread iterating through
            // all k - we'll boil this down at the end
            bRealArray[offsetb + m] += pSb1*PHIM;
            bImagArray[offsetb + m] += pSb2*PHIM;
            for(int m2 = 0;m2 < M;m2++)
            { //for each A's col
                // HOT
                // thread-safe due to use of offset
               AA[offsetA + m*M + m2] += pSA*PHIM*PHI[k*M+m2];//AverageBeta
            } // m2
        } // m
    } // k

    // Now we reduce the results from AA into A.   Needs to be in a
    // separate loop for thread-safety
    #pragma omp parallel for schedule(dynamic)
    for(int m = 0;m < M;m++)
    {
        // Set A
        for(int m2 = 0;m2 < M;m2++)
        { //for each A's col
            double sumReal = 0;

            // maxthreads is relatively small, and each thread is striding by maxthreads
            // in PHI - hopefully won't hurt the caches too much
            // Iterating with k as the outside loop as above would cause a race
            for(int k = 0; k < maxthreads; k++) {
                // HOT
                sumReal += AA[k*Aoffset + m*M + m2];
                AA[k*Aoffset + m*M + m2] = 0;  // Zero for next time around
            }

            A[m*M + m2] = sumReal;
        }  // m2

        // Set b
        // Also reduce the results for b
        double bReal = 0;
        double bImag = 0;

        for(int t = 0;t < maxthreads;t++)
        {
            if (!used[t]) continue;     // BUG this is missing
            bReal += bRealArray[t*boffset + m];
            bRealArray[t*boffset + m] = 0; // Zero for next time around
        }
        for(int t = 0;t < maxthreads;t++)
        {
            if (!used[t]) continue;     // BUG this is missing
            bImag += bImagArray[t*boffset + m];
            bImagArray[t*boffset + m] = 0; // Zero for next time around
        }
        b[2*m  ] = bReal;
        b[2*m+1] = bImag;
    } // m

    for(int m = 0;m < M-1;m++)   //for each A's rows
        A[m*M + m] += AverageAlpha/AverageBeta;

    mkl_solve(A,M,b,2);
}

void getWd1(
    const int d,
    WritableDoubleVector& A,
    WritableDoubleVector& b,
    int maxthreads,
    WritableDoubleVector& AA,
    WritableDoubleVector& bRealArray, WritableDoubleVector& bImagArray)
{
    __assume_aligned(WReal,64);
    __assume_aligned(WImag,64);
    __assume_aligned(CTFvalue,64);
    __assume_aligned(TReal,64);
    __assume_aligned(TImag,64);
    __assume_aligned(PHI,64);
    __assume_aligned(R,64);
    __assume_aligned(AA,64);
    __assume_aligned(A,64);
    __assume_aligned(b,64);
    __assume_aligned(bRealArray,64);
    __assume_aligned(bImagArray,64);

    const size_t usedCapacity = size_t(maxthreads);
    bool* used = vNew(bool,usedCapacity);
    for (int i = 0; i < usedCapacity; i++) used[i] = false;

    struct S {
        double pSA ;
        double pSb1;
        double pSb2;
        S() { pSA = pSb1 = pSb2 = 0.0; }
    };
    std::vector<S> SForK(K);

    static const int nStep = 240;   // See below
    static const int kStep = 8;     // See below

    #pragma omp parallel for schedule(dynamic) collapse(2)
    for (int nLo = 0; nLo < N; nLo += nStep) {
        for (int kLo = 0; kLo < K; kLo += kStep) {      // K is only about 1000
            const int  nHi = std::min(int(N), nLo+nStep);

            // The following loop reads 3     doubles (32 bytes) indexed by [n] but not by k to produce its outputs
            // The following loop reads kStep doubles ( 8 bytes) indexed by [n] and     by k to produce its outputs
            // So each loop reads nStep*(32 + kStep*8)*nStep bytes
            // 32000 < nStep * (32 + 8*8)
            // 32000 / 96 < nStep
            // nStep should be less than 320.  Make it 256 to get plenty of room
            //
            #define HEAD(Q) double pSA##Q = 0, pSb1##Q=0, pSb2##Q=0;
                HEAD(0) HEAD(1) HEAD(2) HEAD(3) HEAD(4) HEAD(5) HEAD(6) HEAD(7)
            #undef HEAD

            #define BODY_PREFIX \
                double CTFVdNn  = CTFvalue[INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)];              \
                double TRealdNn = TReal   [INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)];              \
                double TImagdNn = TImag   [INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) + INDEXSZ(n)];              \
            // end of macro

            #define BODY(kOffset) \
            {                                                                                           \
                int k = kLo + kOffset;                                                                  \
                double RCTFV   = R [INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) + INDEXSZ(n)]*CTFVdNn;             \
                pSA ##kOffset += RCTFV*CTFVdNn;                                                         \
                pSb1##kOffset += RCTFV*TRealdNn;                                                        \
                pSb2##kOffset += RCTFV*TImagdNn;                                                        \
            };                                                                                          \
            // end of macro

            switch (K - kLo) {
            case 1:         for (int n = nLo; n < nHi; n++) { BODY_PREFIX BODY(0); } break;
            case 2:     for (int n = nLo; n < nHi; n++) { BODY_PREFIX BODY(0); BODY(1); } break;
            case 3:     for (int n = nLo; n < nHi; n++) { BODY_PREFIX BODY(0); BODY(1); BODY(2); } break;
            case 4:     for (int n = nLo; n < nHi; n++) { BODY_PREFIX BODY(0); BODY(1); BODY(2); BODY(3); } break;
            case 5:     for (int n = nLo; n < nHi; n++) { BODY_PREFIX BODY(0); BODY(1); BODY(2); BODY(3); BODY(4); } break;
            case 6:     for (int n = nLo; n < nHi; n++) { BODY_PREFIX BODY(0); BODY(1); BODY(2); BODY(3); BODY(4); BODY(5); } break;
            case 7:     for (int n = nLo; n < nHi; n++) { BODY_PREFIX BODY(0); BODY(1); BODY(2); BODY(3); BODY(4); BODY(5); BODY(6); } break;
            static_assert(kStep == 8, "Wrong number of kStep cases" );
            default:
                // For these to be vector aligned, the following must be true
                //      The whole arrays must be 64B aligned
                //      INDEXSZ(INDEXSZ(d) * INDEXSZ(N)) must be a multiple of 8
                //      INDEXSZ(INDEXSZ(k) * INDEXSZ(N)) must be a multiple of 8
                // This will happen when N is a multiple of 8 - which 1000
                //
                static_assert(N % 8 == 0, "N must be a multiple of eight");

                #pragma ivdep
                #pragma vector aligned
                for (int n = nLo; n < nHi; n++) {
                    BODY_PREFIX
                    BODY(0); BODY(1); BODY(2); BODY(3); BODY(4); BODY(5); BODY(6); BODY(7);
                    static_assert(kStep == 8, "Wrong number of kStep cases" );
                }
                break;
            }
            #undef BODY
            #define TAIL(Q) if (kLo+Q < K) { \
                S& s = SForK[kLo+Q];            \
                s.pSA  += pSA##Q;               \
                s.pSb1 += pSb1##Q;              \
                s.pSb2 += pSb2##Q;              \
            }                                   // end of macro
            TAIL(0) TAIL(1) TAIL(2) TAIL(3) TAIL(4) TAIL(5) TAIL(6) TAIL(7)
            #undef TAIL
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for(int k = 0; k < K; k++)      // K is only about 1000
    {
        S& s = SForK[k];

        const int tid = omp_get_thread_num();
        assert(tid < maxthreads);
        used[tid] = true;

        const double pSA  = s.pSA;
        const double pSb1 = s.pSb1;
        const double pSb2 = s.pSb2;

        const int offsetA = tid*Aoffset;
        const int offsetb = tid*boffset;

        static const int mStep = 8;
        int m = 0;
        for(; m+mStep <= M; m += mStep)   // M is about 600
        {   // for each A's rows
            // Sum this product separately for each thread iterating through
            // all k - we'll boil this down at the end
#define HEAD(Q) \
            double PHIM##Q = PHI[k*M+m+Q];                  \
            bRealArray[offsetb + m + Q] += pSb1*PHIM##Q;    \
            bImagArray[offsetb + m + Q] += pSb2*PHIM##Q;    // end of macro
            HEAD(0) HEAD(1) HEAD(2) HEAD(3) HEAD(4) HEAD(5) HEAD(6) HEAD(7)

            #pragma ivdep
            for(int m2 = 0;m2 < M;m2++)
            { //for each A's col
                // HOT
                // thread-safe due to use of offset
#define BODY(Q) AA[offsetA + (m+Q)*M + m2] += pSA*PHIM##Q*PHI[k*M+m2];//AverageBeta
                BODY(0) BODY(1) BODY(2) BODY(3) BODY(4) BODY(5) BODY(6) BODY(7)
            } // m2
        } // m
        for (; m < M; m += 1)   // M is about 600
        {
            HEAD(0)
            #pragma ivdep
            for (int m2 = 0; m2 < M; m2++) BODY(0)
        }
    } // k

    // Now we reduce the results from AA into A.   Needs to be in a
    // separate loop for thread-safety
    #pragma omp parallel for schedule(dynamic)
    for(int m = 0;m < M;m++)
    {
        // Set A
        for(int m2 = 0;m2 < M;m2++)
        { //for each A's col
            double sumReal = 0;

            // maxthreads is relatively small, and each thread is striding by maxthreads
            // in PHI - hopefully won't hurt the caches too much
            // Iterating with k as the outside loop as above would cause a race
            for(int k = 0; k < maxthreads; k++) {
                // HOT
                sumReal += AA[k*Aoffset + m*M + m2];
                AA[k*Aoffset + m*M + m2] = 0;  // Zero for next time around
            }

            A[m*M + m2] = sumReal;
        }  // m2

        // Set b
        // Also reduce the results for b
        double bReal = 0;
        double bImag = 0;

        for(int t = 0;t < maxthreads;t++)
        {
            if (!used[t]) continue;     // BUG this is missing
            bReal += bRealArray[t*boffset + m];
            bRealArray[t*boffset + m] = 0; // Zero for next time around
        }
        for(int t = 0;t < maxthreads;t++)
        {
            if (!used[t]) continue;     // BUG this is missing
            bImag += bImagArray[t*boffset + m];
            bImagArray[t*boffset + m] = 0; // Zero for next time around
        }
        b[2*m  ] = bReal;
        b[2*m+1] = bImag;
    } // m

    for(int m = 0;m < M-1;m++)   //for each A's rows
        A[m*M + m] += AverageAlpha/AverageBeta;

    mkl_solve(A,M,b,2);
}

struct S {
#define SELTS \
    ELT(A          , Aoffset)                           SEP \
    ELT(b          , boffset)                           SEP \
    ELT(AA         , K*boffset)                         SEP \
    ELT(bRealArray , omp_get_max_threads()*boffset)     SEP \
    ELT(bImagArray , omp_get_max_threads()*boffset)         // end of macro
                    /* BUG these had the wrong size */

#define SEP
#define ELT(NAME,SIZE) WritableDoubleVector NAME;
    SELTS
#undef ELT
#undef SEP

    S(int which) :
#define SEP ,
#ifdef BOUNDSCHECK
#define ELT(NAME,SIZE) NAME(SIZE)
#else
#define ELT(NAME,SIZE) NAME(mallocSetDoubles(SIZE))
#endif
        SELTS
#undef ELT
#undef SEP
    {
        for (int i = 0; i < omp_get_max_threads()*boffset; i++) bRealArray[i] = bImagArray[i] = 0;

        double start = dtime();

        int maxthreads = omp_get_max_threads();
        int d = 0;
        switch (which) {
        case 0: getWd0(d, A, b, maxthreads, AA, bRealArray, bImagArray);    break;
        case 1: getWd1(d, A, b, maxthreads, AA, bRealArray, bImagArray);    break;
        }

        double done = dtime() - start;
        std::cout << " total time was____" << done
            << " of which mkl_solve took__" << mkl_solve_elapsed 
			<< std::endl;
    }

    ~S() {
#define SEP
#ifdef BOUNDSCHECK
        #define ELT(NAME,SIZE)
#else
        #define ELT(NAME,SIZE) freeDoubles(NAME,SIZE);
#endif
        SELTS
#undef ELT
#undef SEP
    }

    bool nearEnough(double lhs, double rhs) {
        auto diff = std::abs(lhs-rhs);
        auto max = std::max(std::abs(lhs), std::abs(rhs));
        if ((diff < max*1e-6) || (diff < 1e-6)) return true;
        return false;
    }

    void compare(S & rhs) {
        const char* name = NULL;
        INDEXSZ error = 0;
        double lhsY, rhsY;
#define SEP
#define ELT(NAME,SIZE) \
        for (INDEXSZ x = 0; x < SIZE; x++) {                \
            if (nearEnough(NAME[x],rhs.NAME[x])) continue;  \
            name = #NAME;                                   \
            error = x;                                      \
            lhsY = NAME[x]; rhsY = rhs.NAME[x];             \
            break;                                          \
        }                                                   \
        if (name) goto Diff;                                // end of macro
        SELTS
#undef ELT
#undef SEP
        if (!name) return;
    Diff:
        std::cerr << "Different " << name << "[" << error << "]" << " lhs:" << lhsY << " rhs:" << rhsY << " Diff:" << lhsY - rhsY << std::endl;
    }
};


void unitTestCorrectness() {
	auto max_threads = omp_get_max_threads();
    std::cerr << "starting trials with " << max_threads << " threads" << std::endl;
    std::cout << "s0" << std::endl; S s0(0);
    std::cout << "s1" << std::endl; S s1(1);
    s0.compare(s1);
}

void unitTestPerformance() {
	auto max_threads = omp_get_max_threads();
    std::cerr << "starting trials with " << max_threads << " threads" << std::endl;
    std::cout << "s0" << std::endl; S s0(0);
    std::cout << "s1" << std::endl; S s1(1);
}

}	// namespace GTM_Optimizer_Kernel
