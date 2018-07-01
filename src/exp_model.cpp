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

#include "exp_model.h"

//
void ParticleModel::initialize(int _ori_size,double _pixel_size,double _particle_diameter,int _width_mask_edge,
                               double _sigma2_fudge,int _random_seed,bool _do_norm_correction,bool _do_zero_mask,
                               bool _do_shifts_onthefly,int _nr_threads,DataStream* _global_data_stream)
{
    ori_size = _ori_size;ori_size2 = _ori_size*ori_size;
    ori_Fsize = (_ori_size/2+1);ori_Fsize2 = (_ori_size/2+1)*_ori_size;
    pixel_size = _pixel_size;
    particle_diameter = _particle_diameter;
    width_mask_edge = _width_mask_edge;
    do_norm_correction = _do_norm_correction;
    do_zero_mask = _do_zero_mask;
    random_seed = _random_seed;
    sigma2_fudge = _sigma2_fudge;
    // Prepare thread data
    nr_threads = _nr_threads;//omp_get_max_threads();
    threadFimages_real	.init(nr_threads, ori_size2);	threadFimages_real	.fill_with_first_touch(0.);
    threadFimages_imag	.init(nr_threads, ori_size2);	threadFimages_imag	.fill_with_first_touch(0.);
    threadImages		.init(nr_threads, ori_size2);	threadImages		.fill_with_first_touch(0.);
    threadFFTtransformer.resize(nr_threads);
    threadCTFer			.resize(nr_threads);
    //
#if defined(FLOAT_PRECISION)
    for(auto &ffter : threadFFTtransformer) ffter = new FFTWFTransformer(ori_size,ori_size);
#else
    for(auto &ffter : threadFFTtransformer) ffter = new FFTWTransformer(ori_size,ori_size);
#endif
    for(auto &ctfer : threadCTFer) ctfer = new CTF;
    //
    do_shifts_onthefly = _do_shifts_onthefly;
    global_data_stream = _global_data_stream;
}

//
void ParticleModel::finalize()
{
    threadFimages_real	.fini();
    threadFimages_imag	.fini();
    threadImages		.fini();
    exp_metadata		.fini();
    
    for(auto &ffter : threadFFTtransformer	) delete ffter;
    for(auto &ctfer : threadCTFer			) delete ctfer;
    threadFFTtransformer.resize(0);
    threadCTFer			.resize(0);
}

// Set up the Images Transformer
void ParticleModel::setup(int _nr_pool,int _current_size,int _coarse_size,int exp_nr_trans/* = 0*/,int exp_nr_over_trans/* = 0*/)
{
    if (exp_nr_trans==0) assert(exp_nr_over_trans==0 && do_shifts_onthefly==false);
    if (exp_nr_trans!=0) assert(exp_nr_over_trans!=0 && do_shifts_onthefly==true);
    exp_nr_images = _nr_pool;
    fine_size = _current_size;fine_size2 = _current_size*_current_size;
    fine_Fsize = (_current_size/2+1);fine_Fsize2 = (_current_size/2+1)*_current_size;
    coarse_size = _coarse_size;coarse_size2 = _coarse_size*_coarse_size;
    coarse_Fsize = (_coarse_size/2+1);coarse_Fsize2 = (_coarse_size/2+1)*_coarse_size;
    
    // Prepare the images data
    exp_images	.init(exp_nr_images, ori_size2);	exp_images	.fill_with_first_touch(0.);
    Fctfs		.init(exp_nr_images, fine_Fsize2);	Fctfs		.fill_with_first_touch(0.);
    exp_metadata.init(exp_nr_images);
    
    // Prepare the fourier images data
    Fimages_mask_fine_real	.init(exp_nr_images, fine_Fsize2);	Fimages_mask_fine_real	.fill_with_first_touch(0.);
    Fimages_mask_fine_imag	.init(exp_nr_images, fine_Fsize2);	Fimages_mask_fine_imag	.fill_with_first_touch(0.);
    Fimages_mask_coarse_real.init(exp_nr_images, coarse_Fsize2);Fimages_mask_coarse_real.fill_with_first_touch(0.);
    Fimages_mask_coarse_imag.init(exp_nr_images, coarse_Fsize2);Fimages_mask_coarse_imag.fill_with_first_touch(0.);
    Fimages_nomask_real		.init(exp_nr_images, fine_Fsize2);	Fimages_nomask_real		.fill_with_first_touch(0.);
    Fimages_nomask_imag		.init(exp_nr_images, fine_Fsize2);	Fimages_nomask_imag		.fill_with_first_touch(0.);
}
    
void ParticleModel::destory()
{
    exp_images				.fini();
    Fctfs					.fini();
    Fimages_mask_fine_real	.fini();
    Fimages_mask_fine_imag	.fini();
    Fimages_mask_coarse_real.fini();
    Fimages_mask_coarse_imag.fini();
    Fimages_nomask_real		.fini();
    Fimages_nomask_imag		.fini();
}

// TODO,use optimized one
void ParticleModel::preShiftedImagesCtfsAndInvSigma2s(Aligned3dArray<FDOUBLE>& exp_Fimgs_shifted_real,
                                                      Aligned3dArray<FDOUBLE>& exp_Fimgs_shifted_imag,
                                                      VectorOfArray2d<FDOUBLE>& exp_local_Minvsigma2s,
                                                      VectorOfArray1d<FDOUBLE>& exp_local_sqrtXi2,bool do_cc,
                                                      VectorOfArray2d<FDOUBLE>& exp_local_Fctfs,int exp_current_size,
                                                      SamplingGrid& samplingGrid,MLModel& mlModel,
                                                      VectorOfInt& Mresol_coarse,VectorOfInt& Mresol_fine)
{
#ifdef DATA_STREAM
    global_data_stream->foutDouble(0, "##################start_precalculateShiftedImagesCtfsAndInvSigma2s#####################", __FILE__, __LINE__);
#endif
    
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    bool do_coarse_search = (exp_current_size == coarse_size);
    int exp_nr_trans = samplingGrid.exp_nr_trans;
    int exp_nr_over_trans = samplingGrid.exp_nr_over_trans;
    //
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
#ifdef DATA_STREAM
        global_data_stream->foutInt(iimage, "precalculateShiftedImagesCtfsAndInvSigma2s()_start------", __FILE__, __LINE__);
#endif
        int tid = omp_get_thread_num();
        
        int igroup = exp_metadata[iimage].GROUP_NO-1;
        
        auto tempFimage_real = threadFimages_real[tid].wptr(ori_size2);
        auto tempFimage_imag = threadFimages_imag[tid].wptr(ori_size2);
        
        // Downsize Fimg and Fctf to exp_current_image_size, also initialise Fref and Fimg_shift to the right size
        // In the second pass of the adaptive approach this will have no effect,
        // since then exp_current_image_size will be the same as the size of exp_Fctfs
        windowFourierTransform(Fimages_mask_fine_real[iimage].wptr(fine_Fsize2),
                               Fimages_mask_fine_imag[iimage].wptr(fine_Fsize2), fine_size,
                               tempFimage_real, tempFimage_imag, exp_current_size);
        
        if (do_cc)
        {
            double sumxi2 = 0.;
            for (int n = 0; n < exp_current_Fsize2; n++) {
                sumxi2 += tempFimage_real[n]*tempFimage_real[n] \
                        + tempFimage_imag[n]*tempFimage_imag[n];
            }
            // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
            exp_local_sqrtXi2.wptrAll()[iimage] = sqrt(sumxi2);
        }
        
#ifdef DATA_STREAM
        global_data_stream->foutDouble(tempFimage_real, exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_Fimg_real", __FILE__, __LINE__);
        global_data_stream->foutDouble(tempFimage_imag, exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_Fimg_imag", __FILE__, __LINE__);
#endif
        
        windowTransform(Fctfs[iimage].wptr(fine_Fsize2), fine_size, exp_local_Fctfs[iimage].wptr(exp_current_Fsize2), exp_current_size);
        
#ifdef DATA_STREAM
        global_data_stream->foutDouble(exp_local_Fctfs[iimage].wptr(exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_local_Fctfs", __FILE__, __LINE__);
        global_data_stream->check();global_data_stream->flush();
#endif
        for (int itrans = 0; itrans < exp_nr_trans; itrans++)
        {
            for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
            {
                // Shift through phase-shifts in the Fourier transform
                // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
                FDOUBLE shiftx,shifty;
                samplingGrid.getShiftxy(shiftx, shifty, itrans, iover_trans);
                
                int itrans_over_trans = itrans*exp_nr_over_trans+iover_trans;
                auto exp_Fimgs_shifted_real_aux = exp_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                auto exp_Fimgs_shifted_imag_aux = exp_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                
                shiftImageInFourierTransform(tempFimage_real, tempFimage_imag,
                                             exp_Fimgs_shifted_real_aux,exp_Fimgs_shifted_imag_aux,
                                             exp_current_size,shiftx,shifty,ori_size);
            }
        }
        
        auto Minvsigma2 = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
        auto myMresol = do_coarse_search ? Mresol_coarse.rptrAll() : Mresol_fine.rptrAll();
        auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].wptr(ori_Fsize);
        // With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
        for (int n = 0; n < exp_current_Fsize2; n++)
        {
            int ires = myMresol[n];
            // Exclude origin (ires==0) from the Probability-calculation
            // This way we are invariant to additive factors
            if (ires > 0){
                Minvsigma2[n] = CHECK_NOT_NAN_OR_INF(1. / (sigma2_fudge * sigma2_noise_igroup[ires]));
            }
            else
                Minvsigma2[n] = 0;
        }
#ifdef DATA_STREAM
        global_data_stream->foutDouble(Minvsigma2, exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_local_Minvsigma2s", __FILE__, __LINE__);
        global_data_stream->foutDouble(exp_Fimgs_shifted_real.wptr(iimage, 0, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_real1", __FILE__, __LINE__);
        global_data_stream->foutDouble(exp_Fimgs_shifted_imag.wptr(iimage, 0, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_imagn", __FILE__, __LINE__);
        global_data_stream->foutDouble(exp_Fimgs_shifted_real.wptr(iimage, exp_nr_trans*exp_nr_over_trans-1, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_realN", __FILE__, __LINE__);
        global_data_stream->foutDouble(exp_Fimgs_shifted_imag.wptr(iimage, exp_nr_trans*exp_nr_over_trans-1, exp_current_Fsize2), exp_current_Fsize2, "precalculateShiftedImagesCtfsAndInvSigma2s()_exp_Fimgs_shifted_imagN", __FILE__, __LINE__);
        global_data_stream->check();global_data_stream->flush();
#endif
    }
}


void ParticleModel::pregetShiftedImagesNomask(Aligned3dArray<FDOUBLE>& exp_nomask_Fimgs_shifted_real,
                                              Aligned3dArray<FDOUBLE>& exp_nomask_Fimgs_shifted_imag,
                                              int exp_current_size,SamplingGrid& samplingGrid)
{
    assert(exp_current_size==fine_size);
    int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
    int exp_nr_trans = samplingGrid.exp_nr_trans;
    int exp_nr_over_trans = samplingGrid.exp_nr_over_trans;
    //
#pragma omp parallel for collapse(3)
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        for (int itrans = 0; itrans < exp_nr_trans; itrans++)
        {
            for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
            {
                int itrans_over_trans = itrans*exp_nr_over_trans + iover_trans;
                
                int tid = omp_get_thread_num();
                
                // auto tempFimage_real = threadFimages_real[tid];
                // auto tempFimage_imag = threadFimages_imag[tid];
                
                // exp_current_size is same as current_size
                // windowFourierTransform(Fimages_nomask_real[iimage], Fimages_nomask_imag[iimage], fine_size,
                //                        tempFimage_real, tempFimage_imag, exp_current_size);
                
                // Shift through phase-shifts in the Fourier transform
                // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
                FDOUBLE shiftx,shifty;
                samplingGrid.getShiftxy(shiftx, shifty, itrans, iover_trans);
                
                // NOTE : using exp_Fimgs_shifted as exp_Fimg_shifted_nomask
                auto exp_Fimgs_shifted_nomask_real_aux = exp_nomask_Fimgs_shifted_real.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                auto exp_Fimgs_shifted_nomask_imag_aux = exp_nomask_Fimgs_shifted_imag.wptr(iimage, itrans_over_trans, exp_current_Fsize2);
                
                //shiftImageInFourierTransform(Fimg_nomask_aux, exp_Fimgs_shifted_nomask_aux,exp_current_size,shiftx,shifty,ori_size);
                shiftImageInFourierTransform(Fimages_nomask_real[iimage].wptr(fine_Fsize2),
                                             Fimages_nomask_imag[iimage].wptr(fine_Fsize2),
                                             exp_Fimgs_shifted_nomask_real_aux,exp_Fimgs_shifted_nomask_imag_aux,
                                             exp_current_size,shiftx,shifty,ori_size);
            }
        }
    }
}
    
void ParticleModel::unitTest()
{
    // Vector2D test;
}
