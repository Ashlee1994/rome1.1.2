/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include "util.h"		// used for building precompiled headers on Windows

#include "map3d_optimizer_kernel.h"

namespace Map3dOptimizer_kernel
{

// Note: FourierShellTranslation:: functions moved to exp_model.cpp

IntPerformanceCounter transformOneTileAndGetDiff_performanceCounter("transformOneTileAndGetDiff_performanceCounter");

void ShiftImageInFourierTransformImplemention
	::transformOneTileAndGetDiff(const int shell_n_start, const int shell_n_end,
                             	 const double* Fimag_real, const double* Fimag_imag,
                             	 const double* Fref_real, const double* Fref_imag,
                             	 const double* local_Fctfs,const double* Minvsigma2,
                             	 double myscale,double& diff2)
{
	transformOneTileAndGetDiff_performanceCounter.count.v++;

    assert(0 <= shell_n_start); assert(shell_n_end <= inout_Fsize2);
    assert(shell_n_start<shell_n_end);
    if (/*justMove*/false) {
        for (int n = shell_n_start; n < shell_n_end; n++) {
            double diff_real = Fref_real[n]*local_Fctfs[n]*myscale - Fimag_real[n];
            double diff_imag = Fref_imag[n]*local_Fctfs[n]*myscale - Fimag_imag[n];
            auto term = (diff_real * diff_real + diff_imag * diff_imag) * Minvsigma2[n];
            diff2 += term;
            //ERROR_REPORT("TODO ShiftImageInFourierTransformAndGetDiff..");
        }
        return;
    }
    
#define ALIGNED_DATA(INDEX)	\
	isVectorAligned(&aTable[INDEX])													\
    && isVectorAligned(&bTable[INDEX])												\
    && isVectorAligned(&Fimag_real[INDEX])											\
    && isVectorAligned(&Fimag_imag[INDEX])											\
    && isVectorAligned(&Fref_real[INDEX])												\
    && isVectorAligned(&Fref_imag[INDEX])												\
    && isVectorAligned(&local_Fctfs[INDEX])											\
    && isVectorAligned(&Minvsigma2[INDEX]									// end of macro
                     
#define LOOP \
    for (int n = shell_n_start; n < shell_n_end; n++) {								\
        const double a = aTable[n];													\
        const double b = bTable[n];													\
        const double c = Fimag_real[n];												\
        const double d = Fimag_imag[n];												\
        const double ac = a * c;													\
        const double bd = b * d;													\
        const double ab_cd = (a + b) * (c + d);										\
        double fout_shift_real = CHECK_NOT_IND(ac - bd        	        );			\
        double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd			);			\
        double diff_real = Fref_real[n]*local_Fctfs[n]*myscale - fout_shift_real;	\
        double diff_imag = Fref_imag[n]*local_Fctfs[n]*myscale - fout_shift_imag;	\
        auto term = (diff_real * diff_real + diff_imag * diff_imag) * Minvsigma2[n];\
        diff2 += term;																\
    }  // end of macro
	
    if (ALIGNED_DATA(shell_n_start))
        ) {
#pragma vector aligned
#pragma ivdep
        LOOP
    } else {
#pragma ivdep
        LOOP
    }
    //
#undef ALIGNED_DATA
#undef LOOP
}
    
//
IntPerformanceCounter transformOneTileAndGetSeveralDiff_performanceCounter("transformOneTileAndGetSeveralDiff_performanceCounter");

void ShiftImageInFourierTransformImplemention
    ::transformOneTileAndGetSeveralDiff(const int shell_n_start, const int shell_n_end,
                                 		const double* Fimag_real, const double* Fimag_imag,
                                 		const double* local_Fctfs, const double* Minvsigma2,
                                 		double myscale, Diff* diffs,  int diffs_size)
{
	transformOneTileAndGetSeveralDiff_performanceCounter.count.v += diffs_size;

    assert(0 <= shell_n_start); assert(shell_n_end <= inout_Fsize2);
    assert(shell_n_start < shell_n_end);
    if (/*justMove*/false) {
        for (int n = shell_n_start; n < shell_n_end; n++) {
            for (int iDiff = 0; iDiff < diffs_size; iDiff++)
            {
                auto &		task 		= diffs[iDiff];
                auto & 		diff2		= task.diff2;
                const auto	Fref_real 	= task.Fref_real;
                const auto	Fref_imag 	= task.Fref_imag;
                
                double diff_real = Fref_real[n]*local_Fctfs[n]*myscale - Fimag_real[n];
                double diff_imag = Fref_imag[n]*local_Fctfs[n]*myscale - Fimag_imag[n];
                auto term = (diff_real * diff_real + diff_imag * diff_imag) * Minvsigma2[n];
                diff2 += term;
            }
            // ERROR_REPORT("TODO ShiftImageInFourierTransformAndGetDiff..");
        }
        return;
    }
    
#define ALIGNED_DATA(INDEX)															\
	isVectorAligned(&aTable[INDEX])													\
    && isVectorAligned(&bTable[INDEX])												\
    && isVectorAligned(&Fimag_real[INDEX])											\
    && isVectorAligned(&Fimag_imag[INDEX])											\
    && isVectorAligned(&local_Fctfs[INDEX])											\
    && isVectorAligned(&Minvsigma2[INDEX])								// end of macro

	// Describe what the following accesses to the L2 cache model
	// even though the following is more efficient than this
	//
	for (int iDiff = 0; iDiff < diffs_size; iDiff++) {
		auto &		task	  = diffs[iDiff];
		const auto	Fref_real = task.Fref_real;
		const auto	Fref_imag = task.Fref_imag;
		L2CacheModel::seqAcc("Fref_real"  , -1, &Fref_real  [shell_n_start], shell_n_end - shell_n_start, sizeof(Fref_real  [shell_n_start]));
		L2CacheModel::seqAcc("Fref_imag"  , -1, &Fref_imag  [shell_n_start], shell_n_end - shell_n_start, sizeof(Fref_imag  [shell_n_start]));
		L2CacheModel::seqAcc("aTable"     , -1, &aTable     [shell_n_start], shell_n_end - shell_n_start, sizeof(aTable     [shell_n_start]));
		L2CacheModel::seqAcc("bTable"     , -1, &bTable     [shell_n_start], shell_n_end - shell_n_start, sizeof(bTable     [shell_n_start]));
		L2CacheModel::seqAcc("local_Fctfs", -1, &local_Fctfs[shell_n_start], shell_n_end - shell_n_start, sizeof(local_Fctfs[shell_n_start]));
		L2CacheModel::seqAcc("Minvsigma2" , -1, &local_Fctfs[shell_n_start], shell_n_end - shell_n_start, sizeof(Minvsigma2 [shell_n_start]));
	}

    //
#define HEAD_TASK(TASK) \
	double task##TASK##_diff2       = 0.0;											\
    auto &		task##TASK = diffs[iDiff+TASK];										\
    const auto	task##TASK##_Fref_real = task##TASK.Fref_real;						\
    const auto	task##TASK##_Fref_imag = task##TASK.Fref_imag;			// end of macro
    
#define TAIL_TASK(TASK)																\
    task##TASK.diff2 += task##TASK##_diff2;								 // end of macro

#define FREF_ALIGNED_TASKS(INDEX)													\
    isVectorAligned(&task0_Fref_real[INDEX])										\
	&& isVectorAligned(&task0_Fref_imag[INDEX])								// end of macro
    
#define LOOP_TASK(TASK) {															\
    	auto diff_real = task##TASK##_Fref_real[n]*local_Fctfs[n]*myscale - fout_shift_real;\
    	auto diff_imag = task##TASK##_Fref_imag[n]*local_Fctfs[n]*myscale - fout_shift_imag;\
    	auto term = (diff_real * diff_real + diff_imag * diff_imag) * Minvsigma2[n];\
    	task##TASK##_diff2 += term;													\
	}																		// end of macro
    
#define LOOP \
    for (int n = shell_n_start; n < shell_n_end; n++) {								\
        const double a = aTable[n];													\
        const double b = bTable[n];													\
        const double c = Fimag_real[n];												\
        const double d = Fimag_imag[n];												\
        const double ac = a * c;													\
        const double bd = b * d;													\
        const double ab_cd = (a + b) * (c + d);										\
        double fout_shift_real = CHECK_NOT_IND(ac - bd        	        );			\
        double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd			);			\
        LOOP_TASKS																	\
    }  // end of macro
    
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)	HEAD_TASK(3)	HEAD_TASK(4)	HEAD_TASK(5)	HEAD_TASK(6)	HEAD_TASK(7)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)	TAIL_TASK(3)	TAIL_TASK(4)	TAIL_TASK(5)	TAIL_TASK(6)	TAIL_TASK(7)
#define LOOP_TASKS	LOOP_TASK(0)	LOOP_TASK(1)	LOOP_TASK(2)	LOOP_TASK(3)	LOOP_TASK(4)	LOOP_TASK(5)	LOOP_TASK(6)	LOOP_TASK(7)

    int iDiff = 0;
    
    for (; iDiff+7 < diffs_size; iDiff += 8)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
#pragma vector aligned
#pragma ivdep // LOOP - expanded for performance analysis
            for (int n = shell_n_start; n < shell_n_end; n++) {
                const double a = aTable[n];
                const double b = bTable[n];
                const double c = Fimag_real[n];
                const double d = Fimag_imag[n];
                const double ac = a * c;
                const double bd = b * d;
                const double ab_cd = (a + b) * (c + d);
                double fout_shift_real = CHECK_NOT_IND(ac - bd        	    );
                double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd		);
                auto diff_real = task0_Fref_real[n]*local_Fctfs[n]*myscale - fout_shift_real;
                auto diff_imag = task0_Fref_imag[n]*local_Fctfs[n]*myscale - fout_shift_imag;
                auto term = (diff_real * diff_real + diff_imag * diff_imag) * Minvsigma2[n];
                task0_diff2 += term;
                LOOP_TASK(1)	LOOP_TASK(2)	LOOP_TASK(3)	LOOP_TASK(4)	LOOP_TASK(5)	LOOP_TASK(6)	LOOP_TASK(7)
            }
        } else {
#pragma ivdep
            LOOP
        }
        
        TAIL_TASKS
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)	HEAD_TASK(3)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)	TAIL_TASK(3)
#define LOOP_TASKS	LOOP_TASK(0)	LOOP_TASK(1)	LOOP_TASK(2)	LOOP_TASK(3)
    for (; iDiff+3 < diffs_size; iDiff += 4)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
    #pragma vector aligned
    #pragma ivdep
            LOOP
        } else {
    #pragma ivdep
            LOOP
        }
        
        TAIL_TASKS
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
//
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)
#define LOOP_TASKS	LOOP_TASK(0)	LOOP_TASK(1)	LOOP_TASK(2)
    for (; iDiff+2 < diffs_size; iDiff += 3)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
    #pragma vector aligned
    #pragma ivdep
            LOOP
        } else {
    #pragma ivdep
            LOOP
        }
        
        TAIL_TASKS
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
//
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)
#define LOOP_TASKS	LOOP_TASK(0)	LOOP_TASK(1)
    for (; iDiff+1 < diffs_size; iDiff += 2)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
    #pragma vector aligned
    #pragma ivdep
            LOOP
        } else {
    #pragma ivdep
            LOOP
        }
        
        TAIL_TASKS
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
//
#define HEAD_TASKS	HEAD_TASK(0)
#define TAIL_TASKS	TAIL_TASK(0)
#define LOOP_TASKS	LOOP_TASK(0)
    for (; iDiff+0 < diffs_size; iDiff += 1)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
    #pragma vector aligned
    #pragma ivdep
            LOOP
        } else {
    #pragma ivdep
            LOOP
        }
        
        TAIL_TASKS
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
	//
#undef HEAD_TASK
#undef TAIL_TASK
#undef FREF_ALIGNED_TASKS
#undef LOOP_TASK
#undef LOOP
    //
#undef ALIGNED_DATA
}

void ShiftImageInFourierTransformImplemention
	::transformOneTileAndGetSumWDiff(const int shell_n_start, const int shell_n_end,
                                 	 const int shell_norm_start, const int shell_norm_end,
                                 	 const int shell_scale_start, const int shell_scale_end,bool do_scale_correction,
                                 	 const double* Fimag_real, const double* Fimag_imag,
                                 	 const double* Fref_real, const double* Fref_imag,
                                 	 double myscale,const double* local_Fctfs,
                                 	 double* wsum_sigma2_noise,double* wsum_scale_correction_XA,
                                 	 double* wsum_scale_correction_AA,double weight,double& sum_wdiff2_out)
{
    assert(0 <= shell_n_start); assert(shell_n_start < shell_n_end);
    assert(shell_norm_start==0);assert(shell_norm_start==shell_scale_start);
    assert(shell_scale_end<=shell_norm_end);
    if (/*justMove*/false) {// TODO
        for (int n = shell_n_start; n < shell_n_end; n++) {
            ERROR_REPORT("TODO transformOneTileAndUpdateModel..");
        }
        return;
    }
    double sum_wdiff2 = 0;
    
#define NORM_ALIGNED(INDEX) \
    isVectorAligned(&aTable[INDEX])													\
    && isVectorAligned(&bTable[INDEX])												\
    && isVectorAligned(&Fimag_real[INDEX])											\
    && isVectorAligned(&Fimag_imag[INDEX])											\
    && isVectorAligned(&Fref_real[INDEX])												\
    && isVectorAligned(&Fref_imag[INDEX])												\
    && isVectorAligned(&local_Fctfs[INDEX])											\
    && isVectorAligned(&wsum_sigma2_noise[INDEX])	// end of macro
    
#define SCALE_ALIGNED(INDEX) \
    isVectorAligned(&wsum_scale_correction_XA[INDEX])									\
    && isVectorAligned(&wsum_scale_correction_AA[INDEX]) // end of macro

#define DO_NORM \
    const double a = aTable[n];														\
    const double b = bTable[n];														\
    const double c = Fimag_real[n];													\
    const double d = Fimag_imag[n];													\
    const double ac = a * c;														\
    const double bd = b * d;														\
    const double ab_cd = (a + b) * (c + d);											\
    double fout_shift_real = CHECK_NOT_IND(ac - bd          );						\
    double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd	);						\
    double frefctf_real = Fref_real[n]*local_Fctfs[n]*myscale;						\
    double frefctf_imag = Fref_imag[n]*local_Fctfs[n]*myscale;						\
    /* set sigma2_noise */															\
    /* Use FT of masked image for noise estimation!	*/								\
    double diff_real 	= frefctf_real - fout_shift_real;							\
    double diff_imag 	= frefctf_imag - fout_shift_imag;							\
    double wdiff2 		= weight * (diff_real * diff_real + diff_imag * diff_imag);	\
    /* group-wise sigma2_noise */													\
    wsum_sigma2_noise[n] += wdiff2;													\
    /* For norm_correction */														\
    sum_wdiff2 += wdiff2;															// end of macro
    
#define LOOP_SCALE \
    for (int n = scale_start; n < scale_end; n++) {									\
        DO_NORM																		\
        /* set scale_corr */														\
        double sumXA  = frefctf_real * fout_shift_real;								\
        sumXA +=	frefctf_imag * fout_shift_imag;									\
        double sumA2  = frefctf_real * frefctf_real;								\
        sumA2 += frefctf_imag * frefctf_imag;										\
        wsum_scale_correction_XA[n] += weight * sumXA;								\
        wsum_scale_correction_AA[n] += weight * sumA2;								\
    } // end of macro
    //
#define LOOP_NORM \
    for (int n = norm_start; n < norm_end; n++) {									\
        DO_NORM																		\
    } // end of macro
    
    if (do_scale_correction) // do scale correction
    {
        //
        int scale_start = std::max(shell_scale_start, shell_n_start);
        int scale_end = std::min(shell_scale_end, shell_n_end);
        if (NORM_ALIGNED(scale_start) &&
            SCALE_ALIGNED(scale_start)
            ){
#pragma vector aligned
#pragma ivdep
            LOOP_SCALE
        }
        else{
#pragma ivdep
            LOOP_SCALE
        }
        //
        int norm_start = std::max(shell_scale_end, shell_n_start);
        int norm_end = std::min(shell_norm_end, shell_n_end);
        if (NORM_ALIGNED(norm_start)
            ){
#pragma vector aligned
#pragma ivdep
            LOOP_NORM
        }
        else{
#pragma ivdep
            for (int n = norm_start; n < norm_end; n++) {
                const double a = aTable[n];
                const double b = bTable[n];
                const double c = Fimag_real[n];
                const double d = Fimag_imag[n];
                const double ac = a * c;
                const double bd = b * d;
                const double ab_cd = (a + b) * (c + d);
                double fout_shift_real = CHECK_NOT_IND(ac - bd          );
                double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd	);
                double frefctf_real = Fref_real[n]*local_Fctfs[n]*myscale;
                double frefctf_imag = Fref_imag[n]*local_Fctfs[n]*myscale;
                /* set sigma2_noise */
                /* Use FT of masked image for noise estimation!	*/
                double diff_real 	= frefctf_real - fout_shift_real;
                double diff_imag 	= frefctf_imag - fout_shift_imag;
                double wdiff2 		= weight * (diff_real * diff_real + diff_imag * diff_imag);
                /* group-wise sigma2_noise */
                wsum_sigma2_noise[n] += wdiff2;
                /* For norm_correction */
                sum_wdiff2 += wdiff2;															// end of macro
            } // end of macro
        }
    } // end scale correction
    else // do norm correction
    {
        //
        int norm_start = std::max(shell_norm_start, shell_n_start);
        int norm_end = std::min(shell_norm_end, shell_n_end);
        if (NORM_ALIGNED(norm_start)
            ){
#pragma vector aligned
#pragma ivdep
            LOOP_NORM
        }
        else{
#pragma ivdep
            LOOP_NORM
        }
    } // end norm correction
    //
    sum_wdiff2_out = sum_wdiff2;
    //
#undef NORM_ALIGNED
#undef SCALE_ALIGNED
#undef DO_NORM
#undef LOOP_SCALE
#undef LOOP_NORM
}

    
void ShiftImageInFourierTransformImplemention
    ::transformOneTileAndGetSeveralSumWDiff(const int shell_n_start, const int shell_n_end,
                                            const int shell_norm_start, const int shell_norm_end,
                                            const int shell_scale_start, const int shell_scale_end,bool do_scale_correction,
                                            const double* Fimag_real, const double* Fimag_imag,
                                            double myscale,const double* local_Fctfs,
                                            double* wsum_sigma2_noise,double* wsum_scale_correction_XA,
                                            double* wsum_scale_correction_AA,
                                            SumWDiff2* sumWdiff2s,int sumwdiff2_size)
{
    //
#define HEAD_TASK(TASK) \
    auto &		task##TASK = sumWdiff2s[iSumwdiff2+TASK];							\
    const auto  task##TASK##_weight	   = task##TASK.weight;							\
    const auto	task##TASK##_Fref_real = task##TASK.Fref_real;						\
    const auto	task##TASK##_Fref_imag = task##TASK.Fref_imag;				// end of macro

    //
#define NORM_ALIGNED(INDEX) \
    isVectorAligned(&aTable[INDEX])													\
    && isVectorAligned(&bTable[INDEX])												\
    && isVectorAligned(&Fimag_real[INDEX])											\
    && isVectorAligned(&Fimag_imag[INDEX])											\
    && isVectorAligned(&local_Fctfs[INDEX])											\
    && isVectorAligned(&wsum_sigma2_noise[INDEX])								// end of macro
    
#define SCALE_ALIGNED(INDEX) \
    isVectorAligned(&wsum_scale_correction_XA[INDEX])									\
    && isVectorAligned(&wsum_scale_correction_AA[INDEX]) 						// end of macro
    
#define FREF_ALIGNED_TASKS(INDEX)													\
    isVectorAligned(&task0_Fref_real[INDEX])											\
    && isVectorAligned(&task0_Fref_imag[INDEX])								// end of macro
    
#define SCALE_TASK(TASK){ 	\
        auto frefctf_real = task##TASK##_Fref_real[n]*local_Fctfs[n]*myscale;			\
        auto frefctf_imag = task##TASK##_Fref_imag[n]*local_Fctfs[n]*myscale;			\
        /* set sigma2_noise */															\
        /* Use FT of masked image for noise estimation!	*/								\
        auto diff_real = frefctf_real - fout_shift_real;								\
        auto diff_imag = frefctf_imag - fout_shift_imag;								\
        auto wdiff2 = task##TASK##_weight*(diff_real*diff_real + diff_imag*diff_imag);	\
        /* group-wise sigma2_noise */													\
        sigma2_noise += wdiff2;															\
        /* For norm_correction */														\
        sum_wdiff2 += wdiff2;															\
        /* set scale_corr */															\
        auto sumXA  = frefctf_real * fout_shift_real;									\
        sumXA +=	frefctf_imag * fout_shift_imag;										\
        auto sumA2  = frefctf_real * frefctf_real;										\
        sumA2 += frefctf_imag * frefctf_imag;											\
        scale_correction_XA += task##TASK##_weight * sumXA;								\
        scale_correction_AA += task##TASK##_weight * sumA2;								\
    }																		// end of macro
    
#define NORM_TASK(TASK){	\
        auto frefctf_real = task##TASK##_Fref_real[n]*local_Fctfs[n]*myscale;			\
        auto frefctf_imag = task##TASK##_Fref_imag[n]*local_Fctfs[n]*myscale;			\
        /* set sigma2_noise */															\
        /* Use FT of masked image for noise estimation!	*/								\
        auto diff_real = frefctf_real - fout_shift_real;								\
        auto diff_imag = frefctf_imag - fout_shift_imag;								\
        auto wdiff2 = task##TASK##_weight*(diff_real*diff_real + diff_imag*diff_imag);	\
        /* group-wise sigma2_noise */													\
        sigma2_noise += wdiff2;															\
        /* For norm_correction */														\
        sum_wdiff2 += wdiff2;															\
    }																		// end of macro
    
#define LOOP_MAIN \
    const double a = aTable[n];															\
    const double b = bTable[n];															\
    const double c = Fimag_real[n];														\
    const double d = Fimag_imag[n];														\
    const double ac = a * c;															\
    const double bd = b * d;															\
    const double ab_cd = (a + b) * (c + d);												\
    double fout_shift_real = CHECK_NOT_IND(ac - bd          );							\
    double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd	);				// end of macro
    
#define LOOP_SCALE \
	for(int n = scale_start;n < scale_end;n++){											\
    	LOOP_MAIN																		\
		double sigma2_noise = 0;														\
		double scale_correction_XA = 0;													\
    	double scale_correction_AA = 0;													\
		SCALE_TASKS																		\
		/* NOTE : write on the end to reduce store...*/									\
		wsum_sigma2_noise[n] += sigma2_noise;											\
		wsum_scale_correction_XA[n] += scale_correction_XA;								\
		wsum_scale_correction_AA[n] += scale_correction_AA;								\
    }																		// end of macro
    
#define LOOP_NORM \
    for(int n = norm_start;n < norm_end;n++){											\
        LOOP_MAIN																		\
        double sigma2_noise = 0;														\
        NORM_TASKS																		\
		wsum_sigma2_noise[n] += sigma2_noise;											\
    }																		// end of macro
    
    int iSumwdiff2 = 0;
    int scale_start,scale_end,norm_start,norm_end;
    if (do_scale_correction) {
        scale_start = std::max(shell_scale_start, shell_n_start);
        scale_end = std::min(shell_scale_end, shell_n_end);
        norm_start = std::max(shell_scale_end, shell_n_start);
        norm_end = std::min(shell_norm_end, shell_n_end);
    }
    else{
        norm_start = std::max(shell_norm_start, shell_n_start);
        norm_end = std::min(shell_norm_end, shell_n_end);
    }
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)	HEAD_TASK(3)	HEAD_TASK(4)	HEAD_TASK(5)	HEAD_TASK(6)	HEAD_TASK(7)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)    TAIL_TASK(3)	TAIL_TASK(4)	TAIL_TASK(5)	TAIL_TASK(6)	TAIL_TASK(7)
#define SCALE_TASKS	SCALE_TASK(0)	SCALE_TASK(1)	SCALE_TASK(2)	SCALE_TASK(3)	SCALE_TASK(4)	SCALE_TASK(5)	SCALE_TASK(6)	SCALE_TASK(7)
#define NORM_TASKS  NORM_TASK(0)    NORM_TASK(1)    NORM_TASK(2)    NORM_TASK(3)    NORM_TASK(4)    NORM_TASK(5)    NORM_TASK(6)    NORM_TASK(7)
    for (; iSumwdiff2+7 < sumwdiff2_size; iSumwdiff2 += 8)
    {
        double sum_wdiff2 = 0;
        HEAD_TASKS
        if (do_scale_correction){ // do scale and norm
            if (FREF_ALIGNED_TASKS(scale_start) &&
                NORM_ALIGNED(scale_start) &&
                SCALE_ALIGNED(scale_start)
                ){
                #pragma vector aligned
                #pragma ivdep
                LOOP_SCALE
            }
            else{
                #pragma ivdep
                LOOP_SCALE
            }
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
                #pragma vector aligned
                #pragma ivdep
                LOOP_NORM
            }
            else{
                #pragma ivdep // LOOP - expanded for performance analysis
                for(int n = norm_start;n < norm_end;n++){
                    const double a = aTable[n];
                    const double b = bTable[n];
                    const double c = Fimag_real[n];
                    const double d = Fimag_imag[n];
                    const double ac = a * c;
                    const double bd = b * d;
                    const double ab_cd = (a + b) * (c + d);
                    double fout_shift_real = CHECK_NOT_IND(ac - bd          );
                    double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd	);
                    double sigma2_noise = 0;
                    auto frefctf_real = task0_Fref_real[n]*local_Fctfs[n]*myscale;
                    auto frefctf_imag = task0_Fref_imag[n]*local_Fctfs[n]*myscale;
                    /* set sigma2_noise */
                    /* Use FT of masked image for noise estimation!	*/
                    auto diff_real = frefctf_real - fout_shift_real;
                    auto diff_imag = frefctf_imag - fout_shift_imag;
                    auto wdiff2 = task0_weight*(diff_real*diff_real + diff_imag*diff_imag);
                    /* group-wise sigma2_noise */
                    sigma2_noise += wdiff2;
                    /* For norm_correction */
                    sum_wdiff2 += wdiff2;
                    NORM_TASK(1)    NORM_TASK(2)    NORM_TASK(3)    NORM_TASK(4)    NORM_TASK(5)    NORM_TASK(6)    NORM_TASK(7)
                    wsum_sigma2_noise[n] += sigma2_noise;
                }
            }
        }// end do scale correction
        else{	// only do norm
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
                #pragma vector aligned
                #pragma ivdep
                LOOP_NORM
            }
            else{
                #pragma ivdep
                LOOP_NORM
            }
        }// end do norm correction
        task0.sum_wdiff2 += sum_wdiff2;
    }
#undef HEAD_TASKS
#undef TAIL_TASKS
#undef SCALE_TASKS
#undef NORM_TASKS
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)	HEAD_TASK(3)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)	TAIL_TASK(3)
#define SCALE_TASKS	SCALE_TASK(0)	SCALE_TASK(1)	SCALE_TASK(2)	SCALE_TASK(3)
#define NORM_TASKS  NORM_TASK(0)    NORM_TASK(1)    NORM_TASK(2)    NORM_TASK(3)
    for (; iSumwdiff2+3 < sumwdiff2_size; iSumwdiff2 += 4)
    {
        double sum_wdiff2 = 0;
        HEAD_TASKS
        if (do_scale_correction){ // do scale and norm
            if (FREF_ALIGNED_TASKS(scale_start) &&
                NORM_ALIGNED(scale_start) &&
                SCALE_ALIGNED(scale_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_SCALE
            }
            else{
#pragma ivdep
                LOOP_SCALE
            }
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_NORM
            }
            else{
#pragma ivdep
                LOOP_NORM
            }
        }// end do scale correction
        else{	// only do norm
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_NORM
            }
            else{
#pragma ivdep
                LOOP_NORM
            }
        }// end do norm correction
        task0.sum_wdiff2 += sum_wdiff2;
    }
#undef HEAD_TASKS
#undef TAIL_TASKS
#undef SCALE_TASKS
#undef NORM_TASKS
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)
#define SCALE_TASKS	SCALE_TASK(0)	SCALE_TASK(1)	SCALE_TASK(2)
#define NORM_TASKS  NORM_TASK(0)    NORM_TASK(1)    NORM_TASK(2)
    for (; iSumwdiff2+2 < sumwdiff2_size; iSumwdiff2 += 3)
    {
        double sum_wdiff2 = 0;
        HEAD_TASKS
        if (do_scale_correction){ // do scale and norm
            if (FREF_ALIGNED_TASKS(scale_start) &&
                NORM_ALIGNED(scale_start) &&
                SCALE_ALIGNED(scale_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_SCALE
            }
            else{
#pragma ivdep
                LOOP_SCALE
            }
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_NORM
            }
            else{
#pragma ivdep
                LOOP_NORM
            }
        }// end do scale correction
        else{	// only do norm
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_NORM
            }
            else{
#pragma ivdep
                LOOP_NORM
            }
        }// end do norm correction
        task0.sum_wdiff2 += sum_wdiff2;
    }
#undef HEAD_TASKS
#undef TAIL_TASKS
#undef SCALE_TASKS
#undef NORM_TASKS
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)
#define SCALE_TASKS	SCALE_TASK(0)	SCALE_TASK(1)
#define NORM_TASKS  NORM_TASK(0)    NORM_TASK(1)
    for (; iSumwdiff2+1 < sumwdiff2_size; iSumwdiff2 += 2)
    {
        double sum_wdiff2 = 0;
        HEAD_TASKS
        if (do_scale_correction){ // do scale and norm
            if (FREF_ALIGNED_TASKS(scale_start) &&
                NORM_ALIGNED(scale_start) &&
                SCALE_ALIGNED(scale_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_SCALE
            }
            else{
#pragma ivdep
                LOOP_SCALE
            }
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_NORM
            }
            else{
#pragma ivdep
                LOOP_NORM
            }
        }// end do scale correction
        else{	// only do norm
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_NORM
            }
            else{
#pragma ivdep
                LOOP_NORM
            }
        }// end do norm correction
        task0.sum_wdiff2 += sum_wdiff2;
    }
#undef HEAD_TASKS
#undef TAIL_TASKS
#undef SCALE_TASKS
#undef NORM_TASKS
    //
#define HEAD_TASKS	HEAD_TASK(0)
#define TAIL_TASKS	TAIL_TASK(0)
#define SCALE_TASKS	SCALE_TASK(0)
#define NORM_TASKS  NORM_TASK(0)
    for (; iSumwdiff2+0 < sumwdiff2_size; iSumwdiff2 += 1)
    {
        double sum_wdiff2 = 0;
        HEAD_TASKS
        if (do_scale_correction){ // do scale and norm
            if (FREF_ALIGNED_TASKS(scale_start) &&
                NORM_ALIGNED(scale_start) &&
                SCALE_ALIGNED(scale_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_SCALE
            }
            else{
#pragma ivdep
                LOOP_SCALE
            }
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_NORM
            }
            else{
#pragma ivdep
                LOOP_NORM
            }
        }// end do scale correction
        else{	// only do norm
            if (FREF_ALIGNED_TASKS(norm_start) &&
                NORM_ALIGNED(norm_start)
                ){
#pragma vector aligned
#pragma ivdep
                LOOP_NORM
            }
            else{
#pragma ivdep
                LOOP_NORM
            }
        }// end do norm correction
        task0.sum_wdiff2 += sum_wdiff2;
    }
#undef HEAD_TASKS
#undef TAIL_TASKS
#undef SCALE_TASKS
#undef NORM_TASKS
    //
#undef HEAD_TASK
#undef NORM_ALIGNED
#undef SCALE_ALIGNED
#undef FREF_ALIGNED_TASKS
#undef SCALE_TASK
#undef NORM_TASK
#undef LOOP_MAIN
#undef LOOP_SCALE
#undef LOOP_NORM
}

void ShiftImageInFourierTransformImplemention
    ::transformOneTileAndGetWeightSum(const int shell_n_start,const int shell_n_end,
                                      const double* Fimag_real,const double* Fimag_imag,const double* Minvsigma2s_X_Fctfs2,
                                      double* Frefctf_real,double* Frefctf_imag,double* Fweight,double weight)
{
    assert(0 <= shell_n_start); assert(shell_n_end <= inout_Fsize2);
    assert(shell_n_start<shell_n_end);
    if (/*justMove*/false) {
        for (int n = shell_n_start; n < shell_n_end; n++) {
            Frefctf_real[n] += Fimag_real[n]*weight;
            Frefctf_imag[n] += Fimag_imag[n]*weight;
            Fweight[n] += Minvsigma2s_X_Fctfs2[n]*weight;
            //ERROR_REPORT("TODO ShiftImageInFourierTransformAndGetDiff..");
        }
        return;
    }
    
#define ALIGNED_DATA(INDEX) 	\
    isVectorAligned(&aTable[INDEX])													\
    && isVectorAligned(&bTable[INDEX])												\
    && isVectorAligned(&Fimag_real[INDEX])											\
    && isVectorAligned(&Fimag_imag[INDEX])											\
    && isVectorAligned(&Frefctf_real[INDEX])											\
    && isVectorAligned(&Frefctf_imag[INDEX])											\
    && isVectorAligned(&Fweight[INDEX])												\
    && isVectorAligned(&Minvsigma2s_X_Fctfs2[INDEX]				// end of macro
                     
#define LOOP \
    for (int n = shell_n_start; n < shell_n_end; n++) {								\
        const double a = aTable[n];													\
        const double b = bTable[n];													\
        const double c = Fimag_real[n];												\
        const double d = Fimag_imag[n];												\
        const double ac = a * c;													\
        const double bd = b * d;													\
        const double ab_cd = (a + b) * (c + d);										\
        double fout_shift_real = CHECK_NOT_IND(ac - bd        	        );			\
        double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd			);			\
		Frefctf_real[n] += fout_shift_real*weight;									\
		Frefctf_imag[n] += fout_shift_imag*weight;									\
		Fweight[n] += Minvsigma2s_X_Fctfs2[n]*weight;								\
    }  															// end of macro
    
    if (ALIGNED_DATA(shell_n_start))
        ) {
#pragma vector aligned
#pragma ivdep
        LOOP
    } else {
#pragma ivdep
        LOOP
    }
    //
#undef ALIGNED_DATA
#undef LOOP
}

void ShiftImageInFourierTransformImplemention
    ::transformOneTileAndGetSeveralWeightSum(const int shell_n_start, const int shell_n_end,
                                             const double* Fimag_real, const double* Fimag_imag,
                                             const double* Minvsigma2s_X_Fctfs2,
                                             WeightSum* weightSum,int weightSum_size)
{
    
#define ALIGNED_DATA(INDEX)															\
    isVectorAligned(&aTable[INDEX])													\
    && isVectorAligned(&bTable[INDEX])												\
    && isVectorAligned(&Fimag_real[INDEX])											\
    && isVectorAligned(&Fimag_imag[INDEX])											\
    && isVectorAligned(&Minvsigma2s_X_Fctfs2[INDEX])						// end of macro
    
    //
#define HEAD_TASK(TASK) \
    auto &	task##TASK 				  = weightSum[iWeightSum+TASK];					\
	auto	task##TASK##_weight		  = task##TASK.weight;							\
    auto	task##TASK##_Frefctf_real = task##TASK.Frefctf_real;					\
    auto	task##TASK##_Frefctf_imag = task##TASK.Frefctf_imag;					\
    auto	task##TASK##_Fweight	  = task##TASK.Fweight;   			// end of macro

    
#define FREF_ALIGNED_TASKS(INDEX)													\
    isVectorAligned(&task0_Frefctf_real[INDEX])										\
    && isVectorAligned(&task0_Frefctf_imag[INDEX])									\
    && isVectorAligned(&task0_Fweight[INDEX])								// end of macro
    
#define LOOP_TASK(TASK) {															\
        task##TASK##_Frefctf_real[n] += fout_shift_real*task##TASK##_weight;		\
		task##TASK##_Frefctf_imag[n] += fout_shift_imag*task##TASK##_weight;		\
        task##TASK##_Fweight[n]   += Minvsigma2s_X_Fctfs2[n]*task##TASK##_weight;	\
    }																	// end of macro
    
#define LOOP \
for (int n = shell_n_start; n < shell_n_end; n++) {								\
        const double a = aTable[n];													\
        const double b = bTable[n];													\
        const double c = Fimag_real[n];												\
        const double d = Fimag_imag[n];												\
        const double ac = a * c;													\
        const double bd = b * d;													\
        const double ab_cd = (a + b) * (c + d);										\
        double fout_shift_real = CHECK_NOT_IND(ac - bd        	        );			\
        double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd			);			\
        LOOP_TASKS																	\
    }  // end of macro
    
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)	HEAD_TASK(3)	HEAD_TASK(4)	HEAD_TASK(5)	HEAD_TASK(6)	HEAD_TASK(7)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)	TAIL_TASK(3)	TAIL_TASK(4)	TAIL_TASK(5)	TAIL_TASK(6)	TAIL_TASK(7)
#define LOOP_TASKS	LOOP_TASK(0)	LOOP_TASK(1)	LOOP_TASK(2)	LOOP_TASK(3)	LOOP_TASK(4)	LOOP_TASK(5)	LOOP_TASK(6)	LOOP_TASK(7)
    
    int iWeightSum = 0;
    
    for (; iWeightSum+7 < weightSum_size; iWeightSum += 8)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
#pragma vector aligned
#pragma ivdep // LOOP - expanded for performance analysis
            for (int n = shell_n_start; n < shell_n_end; n++) {
                const double a = aTable[n];
                const double b = bTable[n];
                const double c = Fimag_real[n];
                const double d = Fimag_imag[n];
                const double ac = a * c;
                const double bd = b * d;
                const double ab_cd = (a + b) * (c + d);
                double fout_shift_real = CHECK_NOT_IND(ac - bd        	    );
                double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd		);
                task0_Frefctf_real[n] += fout_shift_real*task0_weight;
                task0_Frefctf_imag[n] += fout_shift_imag*task0_weight;
                task0_Fweight[n]   += Minvsigma2s_X_Fctfs2[n]*task0_weight;
                LOOP_TASK(1)	LOOP_TASK(2)	LOOP_TASK(3)	LOOP_TASK(4)	LOOP_TASK(5)	LOOP_TASK(6)	LOOP_TASK(7)
            }
        } else {
#pragma ivdep
            LOOP
        }
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
	//
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)	HEAD_TASK(3)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)	TAIL_TASK(3)
#define LOOP_TASKS	LOOP_TASK(0)	LOOP_TASK(1)	LOOP_TASK(2)	LOOP_TASK(3)
    
    for (; iWeightSum+3 < weightSum_size; iWeightSum += 4)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP
        } else {
#pragma ivdep
            LOOP
        }
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)	HEAD_TASK(2)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)	TAIL_TASK(2)
#define LOOP_TASKS	LOOP_TASK(0)	LOOP_TASK(1)	LOOP_TASK(2)
    
    for (; iWeightSum+2 < weightSum_size; iWeightSum += 3)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP
        } else {
#pragma ivdep
            LOOP
        }
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
    //
#define HEAD_TASKS	HEAD_TASK(0)	HEAD_TASK(1)
#define TAIL_TASKS	TAIL_TASK(0)	TAIL_TASK(1)
#define LOOP_TASKS	LOOP_TASK(0)	LOOP_TASK(1)
    
    for (; iWeightSum+1 < weightSum_size; iWeightSum += 2)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP
        } else {
#pragma ivdep
            LOOP
        }
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
    //
#define HEAD_TASKS	HEAD_TASK(0)
#define TAIL_TASKS	TAIL_TASK(0)
#define LOOP_TASKS	LOOP_TASK(0)
    
    for (; iWeightSum+0 < weightSum_size; iWeightSum += 1)
    {
        HEAD_TASKS
        
        if (ALIGNED_DATA(shell_n_start) &&
            FREF_ALIGNED_TASKS(shell_n_start)
            ) {
#pragma vector aligned
#pragma ivdep
            LOOP
        } else {
#pragma ivdep
            LOOP
        }
    }
#undef LOOP_TASKS
#undef TAIL_TASKS
#undef HEAD_TASKS
    //
#undef ALIGNED_DATA
#undef HEAD_TASK
#undef FREF_ALIGNED_TASKS
#undef LOOP_TASK
#undef LOOP
    
}
    
//
// --------------------  GetAllSquaredDifferencesFine_Kernel
void GetAllSquaredDifferencesFine_Kernel
	::acquireImage(const double* _Fimag_real,const double* _Fimag_imag,
                   const double* _local_Fctfs,const double* _Minvsigma2,
                   double _myscale,int _shell_n_start,int _shell_n_end)
{
    shell_n_start = _shell_n_start;shell_n_end = _shell_n_end;
    Fimag_real = _Fimag_real;Fimag_imag = _Fimag_imag;
    local_Fctfs = _local_Fctfs;Minvsigma2 = _Minvsigma2;
    myscale = _myscale;
}

void GetAllSquaredDifferencesFine_Kernel
    ::acquireTable(const double* _aTable,const double* _bTable,int _tableSize)
{
    shiftImageAndGetDiffImpl.setABTable(_aTable,_bTable,_tableSize);
}
    
void GetAllSquaredDifferencesFine_Kernel
    ::compute()
{
//     assert(diff2_size==capacity);
//     assert(Fref_size==capacity); // FALSE in local searching,ipsi_tile = 1
//#define DEBUG_GETDIFF
#ifdef DEBUG_GETDIFF
    for (int idiff = 0; idiff < diff2_size; idiff++) {
        shiftImageAndGetDiffImpl.transformOneTileAndGetDiff(shell_n_start, shell_n_end, Fimag_real, Fimag_imag,
                                                            diffs[idiff].Fref_real, diffs[idiff].Fref_imag,
                                                            local_Fctfs, Minvsigma2, myscale, diffs[idiff].diff2);
    }
#else
    shiftImageAndGetDiffImpl.transformOneTileAndGetSeveralDiff(shell_n_start, shell_n_end, Fimag_real, Fimag_imag,
                                                               local_Fctfs, Minvsigma2, myscale, diffs, diff2_size);
#endif
}

void GetAllSquaredDifferencesFine_Kernel
    ::release(double* diff2, int count)
{
	if (count != diff2_size) EXIT_ABNORMALLY;
    for (int i = 0; i < diff2_size; i++)
        diff2[i] = diffs[i].diff2;
    diff2_size = 0;
}
    
//
// -------------   UpdateModel_Kernel
void UpdateModel_Kernel
	::acquireImage(const double* _Fimag_real,const double* _Fimag_imag,const double* _local_Fctfs,
                   double* _wsum_sigma2_noise,double* _wsum_scale_correction_XA,double* _wsum_scale_correction_AA,
                   double _myscale)
{
    Fimag_real = _Fimag_real; Fimag_imag = _Fimag_imag; local_Fctfs = _local_Fctfs;
    wsum_sigma2_noise = _wsum_sigma2_noise; wsum_scale_correction_XA = _wsum_scale_correction_XA;
    wsum_scale_correction_AA = _wsum_scale_correction_AA;myscale = _myscale;
}

void UpdateModel_Kernel
    ::acquireTable(const double* _aTable,const double* _bTable,int _tableSize)
{
    shiftImageAndUpdateModelImpl.setABTable(_aTable, _bTable, _tableSize);
}

void UpdateModel_Kernel
    ::compute()
{
    assert(sumwdiff2_size<=capacity);
    
#if 1
    shiftImageAndUpdateModelImpl
    .transformOneTileAndGetSeveralSumWDiff(shell_n_start, shell_n_end, shell_norm_start, shell_norm_end,
                                           shell_scale_start, shell_scale_end, do_scale_correction,
                                           Fimag_real, Fimag_imag, myscale, local_Fctfs, wsum_sigma2_noise,
                                           wsum_scale_correction_XA, wsum_scale_correction_AA, sumWdiff2s, sumwdiff2_size);
#else
    for (int i = 0; i < sumwdiff2_size; i++){
        shiftImageAndUpdateModelImpl
        .transformOneTileAndGetSumWDiff(shell_n_start, shell_n_end,shell_norm_start,shell_norm_end,
                                        shell_scale_start, shell_scale_end,do_scale_correction,
                                        Fimag_real, Fimag_imag, sumWdiff2s[i].Fref_real, sumWdiff2s[i].Fref_imag,
                                        myscale, local_Fctfs, wsum_sigma2_noise, wsum_scale_correction_XA, wsum_scale_correction_AA,
                                        sumWdiff2s[i].weight, sumWdiff2s[i].sum_wdiff2);
    }
#endif
}

void UpdateModel_Kernel
    ::release(double& sum_wdiff2)
{
    for (int i = 0; i < sumwdiff2_size; i++){
         sum_wdiff2 += sumWdiff2s[i].sum_wdiff2;
    }
    sumwdiff2_size = 0;
}
    
void BackProjection_Kernel
	::acquireImage(const double* _Fimag_real,const double* _Fimag_imag,
                   const double* _Minvsigma2s_X_Fctfs2)
{
    Fimag_real = _Fimag_real;Fimag_imag = _Fimag_imag;
    Minvsigma2s_X_Fctfs2 = _Minvsigma2s_X_Fctfs2;
}

void BackProjection_Kernel
    ::acquireTable(const double* _aTable,const double* _bTable,int _tableSize)
{
    shiftImageAndBackprojectionImpl.setABTable(_aTable, _bTable, _tableSize);
}
    
void BackProjection_Kernel
    ::compute()
{
#if 1
    shiftImageAndBackprojectionImpl
    .transformOneTileAndGetSeveralWeightSum(shell_n_start, shell_n_end, Fimag_real, Fimag_imag,
                                            Minvsigma2s_X_Fctfs2, weightSum, Frefctf_size);
#else
    for (int i = 0; i < Frefctf_size; i++){
        shiftImageAndBackprojectionImpl
        .transformOneTileAndGetWeightSum(shell_n_start, shell_n_end, Fimag_real, Fimag_imag, Minvsigma2s_X_Fctfs2,
                                         weightSum[i].Frefctf_real, weightSum[i].Frefctf_imag, weightSum[i].Fweight, weightSum[i].weight);
    }
#endif
}
    
void BackProjection_Kernel
    ::release()
{
    Frefctf_size = 0;
}
    
    
    
    
    
    
    
    
    
}// end Map3dOptimizer_kernel