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

#ifndef MAP3D_OPTIMIZER_KERNEL_H_
#define MAP3D_OPTIMIZER_KERNEL_H_

#include <omp.h>
#include <vector>

#include "./mpi.h"
#include "./sampling.h"
#include "./map_model.h"
#include "./exp_model.h"
#include "./array_vector.h"
#include "./map_optimizer_base_old.h"

namespace Map3dOptimizer_kernel
{
struct Diff {
    double const * Fref_real;
    double const * Fref_imag;
    double         diff2;
    double		   pleaceholder;// pleaceholder to make four double in struct
    Diff() : diff2(0), Fref_real(nullptr), Fref_imag(nullptr) {}
};
struct SumWDiff2 {
    double const * Fref_real;
    double const * Fref_imag;
    double		   weight;
    double         sum_wdiff2;// TODO
    SumWDiff2() : sum_wdiff2(0), weight(0), Fref_real(nullptr), Fref_imag(nullptr) {}
};
struct WeightSum{
    double * Frefctf_real;
    double * Frefctf_imag;
    double * Fweight;
    double weight;
    WeightSum() : weight(0),Frefctf_real(nullptr), Frefctf_imag(nullptr), Fweight(nullptr) {}
};
    
class ShiftImageInFourierTransformImplemention /*: public ShiftImageInFourierTransform<double, double>*/{
private:
    const double* aTable;
    const double* bTable;
    int inout_size;
    int	 inout_Fsize;
    int	 inout_Fsize2;
    /*bool justMove;*/
    bool isVectorAligned(void const* p) { return (size_t(p) % 64 == 0); }
public:
    ShiftImageInFourierTransformImplemention(int _inout_size){
        inout_size		=	_inout_size;
        inout_Fsize 	= 	inout_size/2+1;
        inout_Fsize2 	= 	inout_size*inout_Fsize;
    }
    ~ShiftImageInFourierTransformImplemention(){}
    //
    inline void setABTable(const double* _aTable,const double* _bTable,int _tableSize){
		if (inout_Fsize2 != _tableSize) {
			std::cerr << __FILE__ << ":" << __LINE__
				<< " inout_Fsize2:" << inout_Fsize2 << " != _tableSize:" << _tableSize
				<< std::endl;
			EXIT_ABNORMALLY;
		}
        aTable = _aTable;
        bTable = _bTable;
    }
    // shift image in fourier space and get diff
    void transformOneTileAndGetDiff(const int shell_n_start, const int shell_n_end,
                                    const double* Fimag_real, const double* Fimag_imag,
                                    const double* Fref_real, const double* Fref_imag,
                                    const double* local_Fctfs, const double* Minvsigma2,
                                    double myscale, double& diff2);
    //
    void transformOneTileAndGetSeveralDiff(const int shell_n_start, const int shell_n_end,
                                           const double* Fimag_real, const double* Fimag_imag,
                                           const double* local_Fctfs, const double* Minvsigma2,
                                           double myscale, Diff* diffs, int diffs_size);
    //
    void transformOneTileAndGetSumWDiff(const int shell_n_start, const int shell_n_end,
                                        const int shell_norm_start, const int shell_norm_end,
                                        const int shell_scale_start, const int shell_scale_end,bool do_scale_correction,
                                        const double* Fimag_real, const double* Fimag_imag,
                                        const double* Fref_real, const double* Fref_imag,
                                        double myscale,const double* local_Fctfs,
                                        double* wsum_sigma2_noise,double* wsum_scale_correction_XA,
                                        double* wsum_scale_correction_AA,double weight,double& sum_wdiff2_out);
    //
    void transformOneTileAndGetSeveralSumWDiff(const int shell_n_start, const int shell_n_end,
                                               const int shell_norm_start, const int shell_norm_end,
                                               const int shell_scale_start, const int shell_scale_end,bool do_scale_correction,
                                               const double* Fimag_real, const double* Fimag_imag,
                                               double myscale,const double* local_Fctfs,
                                               double* wsum_sigma2_noise,double* wsum_scale_correction_XA,
                                               double* wsum_scale_correction_AA,
                                               SumWDiff2* sumWdiff2s,int sumwdiff2_size);
    //
    void transformOneTileAndGetWeightSum(const int shell_n_start,const int shell_n_end,
                                         const double* Fimag_real,const double* Fimag_imag,const double* Minvsigma2s_X_Fctfs2,
                                         double* Frefctf_real,double* Frefctf_imag,double* Fweight,double weight);
    //
    void transformOneTileAndGetSeveralWeightSum(const int shell_n_start, const int shell_n_end,
                                         		const double* Fimag_real, const double* Fimag_imag,
                                                const double* Minvsigma2s_X_Fctfs2,
                                                WeightSum* weightSum,int weightSum_size);
};

class GetAllSquaredDifferencesFine_Kernel{
private:
    ShiftImageInFourierTransformImplemention shiftImageAndGetDiffImpl;
    int shell_n_start,shell_n_end;
    const double *Fimag_real,*Fimag_imag;
    const double *local_Fctfs;
    const double *Minvsigma2;
    double myscale;
private:
    static const int capacity = 8;
    Diff   diffs[capacity];
    int diff2_size,Fref_size;
public:
    GetAllSquaredDifferencesFine_Kernel(int exp_current_size):
    shiftImageAndGetDiffImpl(exp_current_size),
    shell_n_start(0),shell_n_end(0),Fimag_real(nullptr),Fimag_imag(nullptr),
    local_Fctfs(nullptr),Minvsigma2(nullptr),myscale(0),
    diff2_size(0),Fref_size(0){}
    ~GetAllSquaredDifferencesFine_Kernel(){/*TODO*/}
    void acquireImage(const double* _Fimag_real,const double* _Fimag_imag,
                      const double* _local_Fctfs,const double* _Minvsigma2,
                      double _myscale,int _shell_n_start,int _shell_n_end);
    void acquireTable(const double* _aTable,const double* _bTable,int _tableSize);
    inline void appendFref(const double* _Fref_real, const double* _Fref_imag){
        assert(Fref_size < capacity);
        auto & d = diffs[Fref_size++];
        d.Fref_real	=	_Fref_real;
        d.Fref_imag	=	_Fref_imag;
    }
    inline void appendDiff2(int count){
		assert(diff2_size == 0);
        assert(diff2_size + count <= capacity);
        for (int i = 0; i < count; i++) diffs[diff2_size+i].diff2 = 0;
        diff2_size += count;
    }
    void compute();
    void release(double* diff2, int count);
};

class UpdateModel_Kernel{
private:
    ShiftImageInFourierTransformImplemention shiftImageAndUpdateModelImpl;
    int shell_n_start,shell_n_end;
    int shell_scale_start,shell_scale_end;
    int shell_norm_start,shell_norm_end;
    bool do_scale_correction;
    const double *Fimag_real,*Fimag_imag;
    double* wsum_sigma2_noise;
    double* wsum_scale_correction_XA;
    double* wsum_scale_correction_AA;
    const double *local_Fctfs;
    double myscale;
private:
    static const int capacity = 8;
    SumWDiff2   sumWdiff2s[capacity];
    int sumwdiff2_size;
public:
    UpdateModel_Kernel(int _exp_current_size, int _shell_n_start, int _shell_n_end,
                       int _shell_norm_start, int _shell_norm_end, int _shell_scale_start,
                       int _shell_scale_end, bool _do_scale_correction):
    shiftImageAndUpdateModelImpl(_exp_current_size),do_scale_correction(_do_scale_correction),
    shell_n_start(_shell_n_start),shell_n_end(_shell_n_end),shell_norm_start(_shell_norm_start),
    shell_norm_end(_shell_norm_end),shell_scale_start(_shell_scale_start),shell_scale_end(_shell_scale_end),
    Fimag_real(nullptr),Fimag_imag(nullptr),wsum_sigma2_noise(nullptr),wsum_scale_correction_XA(nullptr),
    wsum_scale_correction_AA(nullptr),local_Fctfs(nullptr),myscale(0),sumwdiff2_size(0){}
    ~UpdateModel_Kernel(){/*TODO*/}
    void acquireImage(const double* _Fimag_real,const double* _Fimag_imag,const double* _local_Fctfs,
                      double* _wsum_sigma2_noise,double* _wsum_scale_correction_XA,double* _wsum_scale_correction_AA,
                      double _myscale);
    void acquireTable(const double* _aTable,const double* _bTable,int _tableSize);
    inline void appendSumWdiff2s(const double* _Fref_real, const double* _Fref_imag,double weight){
        assert(sumwdiff2_size < capacity);
        auto & s = sumWdiff2s[sumwdiff2_size++];
        s.Fref_real	= _Fref_real;
        s.Fref_imag	= _Fref_imag;
        s.sum_wdiff2= 0;
        s.weight	= weight;
    }
    void compute();
    void release(double& sum_wdiff2);
};

class BackProjection_Kernel{
private:
    ShiftImageInFourierTransformImplemention shiftImageAndBackprojectionImpl;
    int shell_n_start,shell_n_end;
    const double *Fimag_real,*Fimag_imag;
    const double *Minvsigma2s_X_Fctfs2;
private:
    static const int capacity = 8;
    WeightSum  weightSum[capacity];
    int Frefctf_size;
public:
    BackProjection_Kernel(int _exp_current_size,int _shell_n_start,int _shell_n_end):
    shiftImageAndBackprojectionImpl(_exp_current_size),shell_n_start(_shell_n_start),
    shell_n_end(_shell_n_end),Fimag_real(nullptr),Fimag_imag(nullptr),
    Minvsigma2s_X_Fctfs2(nullptr),Frefctf_size(0){}
    ~BackProjection_Kernel(){/*TODO*/}
    void acquireImage(const double* _Fimag_real,const double* _Fimag_imag,
                      const double* _Minvsigma2s_X_Fctfs2);
    void acquireTable(const double* _aTable,const double* _bTable,int _tableSize);
    inline void appendFref(double* _Frefctf_real,double* _Frefctf_imag,double* _Fweight,double _weight){
        assert(Frefctf_size<capacity);
        auto & w = weightSum[Frefctf_size++];
        w.Frefctf_real	= _Frefctf_real;
        w.Frefctf_imag	= _Frefctf_imag;
        w.Fweight		= _Fweight;
        w.weight		= _weight;
    }
    void compute();
    void release();
};
    
};


#endif /* defined(MAP3D_OPTIMIZER_KERNEL_H_) */
