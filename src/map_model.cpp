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

#include "map_model.h"

void MAPModel::initialize(int _nr_classes,int _ori_size,double _particle_diameter,
                          double _pixel_size,double _ini_high,int _nr_threads,
                          double _width_mask_edge/* = 5*/,double _width_fmask_edge/* = 2*/,
                          int _minres_map/* = 5*/,bool _do_gridding/* = true*/,bool _do_map/* = true*/,
                          int _ref_dim/* = 3*/,std::string _fn_sym/* = "C1"*/)
{
    //
    nr_classes = _nr_classes;
    particle_diameter = _particle_diameter;
    ini_high = _ini_high;
    pixel_size = _pixel_size;
    width_mask_edge = _width_mask_edge;
    width_fmask_edge = _width_fmask_edge;
    minres_map = _minres_map;
    nr_threads = _nr_threads;
    do_map = _do_map;do_gridding = _do_gridding;
    ori_size = _ori_size;ori_Fsize = _ori_size/2+1;
    //
    ref_dim = _ref_dim;
    fn_sym = _fn_sym;
    // initialize reference
    Irefs.resize(nr_classes);
    if(ref_dim==3)
        for(auto& Iref : Irefs) Iref.init(ori_size, ori_size, ori_size, true);
    else
        for(auto& Iref : Irefs) Iref.init(1, ori_size, ori_size, true);
    //
    projector.resize(nr_classes);
    backprojector.resize(nr_classes);
}

void MAPModel::initializeRef(std::string ref_fn)
{
    assert(nr_classes!=0);
    assert(ref_dim==3);
    if (ref_fn == "RANDOM") {
        std::cout<<"ref_fn is NULL,use random data"<<std::endl;
        // copy data
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            for (int i = 0; i < Irefs[iclass].dimzyx; i++) {
                Irefs[iclass](0,0,i) = (FDOUBLE)rand()/(FDOUBLE)RAND_MAX;
            }
        }
    }
    else {
        // read head
        Mrcs::readMrcHead(ref_fn, refHead, ori_size);
        // read volume
        ERROR_CHECK(ref_dim!=3, "-ref *.mrc only used in Map3D.");
        Mrcs::MrcVolume volume(ori_size,0);
        volume.read(ref_fn);
        auto volume_data = volume.ptr();
        // copy data
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            for (int i = 0; i < Irefs[iclass].dimzyx; i++) {
                Irefs[iclass](0,0,i) = volume_data[i];
            }
        }
    }
}

void MAPModel::initializeRef(const std::vector<std::string>& ref_fns)
{
    assert(nr_classes!=0);
    assert(ref_dim==3);
    assert(nr_classes==ref_fns.size());
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        // read volume
        ERROR_CHECK(ref_dim!=3, "-ref *.mrc only used in Map3D.");
        Mrcs::MrcVolume volume(ori_size,0);
        volume.read(ref_fns[iclass]);
        auto volume_data = volume.ptr();
        // copy data
        for (int i = 0; i < Irefs[iclass].dimzyx; i++) {
        	ACCESS(Irefs[iclass], 0, 0, i) = volume_data[i];
        }
    }
}

void MAPModel::initializeRef(Image& averageImage)
{
    assert(nr_classes!=0);
    assert(ref_dim==2);
    for (int iclass = 0; iclass < nr_classes; iclass++){
        memcpy(Irefs[iclass].wptr(), averageImage.rptr(ori_size*ori_size), sizeof(FDOUBLE)*ori_size*ori_size);
    }
    /*
    if (ref_fn != "NULL") {
        FILE* refFile = fopen((ref_fn+".mrcs").c_str(),"rb");
        
        int ref_head[256];
        fread((char*)ref_head,256*sizeof(float),1,refFile);
        float *model_Iref_aux = (float*)_mm_malloc(sizeof(float)*ref_head[2]*ori_size*ori_size,64);
        
        if (ref_head[0] != ori_size || ref_head[2] != nr_classes) {
            std::cerr<<"reference file wrong. "<<ref_head[0]<<" "<<ref_head[2]<<std::endl;
            exit(EXIT_FAILURE);
        }
        fread((char*)model_Iref_aux,nr_classes*ori_size*ori_size*sizeof(float),1,refFile);
        
        for (int i = 0; i < ori_size*ori_size; i++) {
            model_Irefs[i] = model_Iref_aux[i];
        }
        
        _mm_free(model_Iref_aux);
        fclose(refFile);
    }
     */
}

void MAPModel::setFourierTransformMaps(std::vector<VectorOfFDOUBLE>& power_spectrum,int current_size)
{
#pragma omp parallel for
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        projector[iclass].computeFourierTransformMap(Irefs[iclass], ori_size, current_size, power_spectrum[iclass].wptr(ori_Fsize),do_gridding);
        for (int n = 0; n < projector[iclass].data.dimzyx; n++) {
            checkFrefPadValue(projector[iclass].data(0,0,n).real);
            checkFrefPadValue(projector[iclass].data(0,0,n).imag);
        }
        backprojector[iclass].initialize(ori_size, current_size, ref_dim, fn_sym);
    }
    
    if(ref_dim==3)
    {
        int pad_size = backprojector[0].pad_size;
        int pad_Fsize = backprojector[0].pad_Fsize;
        thread_data.resize(0);thread_weight.resize(0);
        thread_data.resize(nr_threads*nr_classes);thread_weight.resize(nr_threads*nr_classes);
#pragma omp parallel for
        for (int thread = 0; thread < nr_threads; thread++) {
            for (int iclass = 0; iclass < nr_classes; iclass++) {
                int thread_iclass = thread*nr_classes+iclass;
                thread_data[thread_iclass].init(pad_size, pad_size, pad_Fsize, true);
                thread_weight[thread_iclass].init(pad_size, pad_size, pad_Fsize, true);
            }
        }
    }
}

// get 2D slice form 3D map
void MAPModel::get2DFourierTransform(int iclass,FDOUBLE* img_out_real, FDOUBLE* img_out_imag,int img_out_size,FDOUBLE A[][3], bool inv)
{
    if (ref_dim == 3)
        projector[iclass].project(img_out_real, img_out_imag, img_out_size, A, inv);
    else
        projector[iclass].rotate2D(img_out_real, img_out_imag, img_out_size, A, inv);
}

// get 2D slice form 3D map
void MAPModel::get2DFourierTransformOneTile(int iclass,FDOUBLE* img_out_real, FDOUBLE* img_out_imag,
                                  			int n_start,int n_end,int img_out_size,const FDOUBLE A[][3], bool inv)
{
    if (ref_dim == 3)
        projector[iclass].projectOneTile(img_out_real, img_out_imag, n_start, n_end, img_out_size, A, inv);
    else
        ERROR_REPORT("get2DFourierTransformOneTile....rotate2D..");
}

//
void MAPModel::set2DFourierTransform(int thread,int iclass,const FDOUBLE* img_in_real,const FDOUBLE* img_in_imag,int img_in_size,
                                     const FDOUBLE A[][3], bool inv,const FDOUBLE* Mweight/* = NULL*/)
{
    int thread_iclass = thread*nr_classes+iclass;
    if (ref_dim == 3)
        backprojector[iclass].backproject(img_in_real,img_in_imag, img_in_size, A, inv,
                                          thread_data[thread_iclass],thread_weight[thread_iclass], Mweight);
    else
        backprojector[iclass].backrotate2D(img_in_real,img_in_imag, img_in_size, A, inv,
                                           backprojector[iclass].data,backprojector[iclass].weight, Mweight);
}

//
void MAPModel::reduceThreadMap()
{
    for (int thread = 0; thread < nr_threads; thread++)
    {
        for (int iclass = 0; iclass < nr_classes; iclass++)
        {
            int thread_iclass = thread*nr_classes+iclass;
            // add data
            auto backprojector_data_iclass = (FDOUBLE*)backprojector[iclass].data.wptr();
            auto backprojector_thread_data_iclass = (FDOUBLE*)thread_data[thread_iclass].wptr();
            size_t real_image_size = 2*backprojector[iclass].data.dimzyx;
            #pragma ivdep
            for (int i = 0; i < real_image_size; i++)
                backprojector_data_iclass[i] += backprojector_thread_data_iclass[i];
            // add weight
            auto backprojector_weight_iclass = backprojector[iclass].weight.wptr();
            auto backprojector_thread_weight_iclass = thread_weight[thread_iclass].wptr();
            size_t size = backprojector[iclass].weight.dimzyx;
            #pragma ivdep
            for (int i = 0; i < size; i++)
                backprojector_weight_iclass[i] += backprojector_thread_weight_iclass[i];
            //
            thread_data[thread_iclass].zero();thread_weight[thread_iclass].fill(0);
        }
    }
}

//
void MAPModel::reconstruction(int iclass,int max_iter_preweight,bool do_map,double tau2_fudge,FDOUBLE* tau2,FDOUBLE* sigma2,FDOUBLE* data_vs_prior,
                              const FDOUBLE* fsc /* only input*/, double normalise/* = 1.*/, bool update_tau2_with_fsc /* = false*/,
                              bool is_whole_instead_of_half /* = false*/,int _nr_threads /* = 1*/)
{
    backprojector[iclass].reconstruct(Irefs[iclass], max_iter_preweight, do_map, tau2_fudge, tau2, sigma2, data_vs_prior,
                                      fsc, normalise, update_tau2_with_fsc, is_whole_instead_of_half, minres_map, _nr_threads);
}


//
void MAPModel::writeResult(std::string filename)
{
    if (ref_dim == 3) {
        for (int iclass = 0; iclass < nr_classes; iclass++)
        {
            std::string fn_mrcs = pathRemoveSuffix(filename) + num2str(iclass+1) + ".mrc";
            Mrcs::MrcVolume oneImage(Irefs[iclass].wptr(),ori_size,pixel_size);
            oneImage.write(fn_mrcs,refHead);
        }
    }
    else{
        auto all_Iref = (FDOUBLE*)_mm_malloc(sizeof(FDOUBLE)*nr_classes*ori_size*ori_size,64);
        for (int iclass = 0; iclass < nr_classes; iclass++)
            memcpy(all_Iref+iclass*ori_size*ori_size, Irefs[iclass].wptr(), sizeof(FDOUBLE)*ori_size*ori_size);
        Mrcs::MrcsImages listOfImages(all_Iref,ori_size,nr_classes);
        listOfImages.write(filename);
        _mm_free(all_Iref);
    }
}

//
void MAPModel::applyLowPassFilter()
{
    if(ini_high <= 0.) return;

    // Make a soft (raised cosine) filter in Fourier space to prevent artefacts in real-space
    // The raised cosine goes through 0.5 at the filter frequency and has a width of width_mask_edge fourier pixels
    FDOUBLE radius = ori_size * pixel_size / ini_high;
    radius -= width_fmask_edge / 2.;
    FDOUBLE radius_p = radius + width_fmask_edge;
    FourierTransformer *transformer;
    Vol<FDOUBLE> Faux_real,Faux_imag;
    if (ref_dim == 2){
        transformer = new FourierTransformer(ori_size,ori_size);
        Faux_real.init(1, ori_size, ori_Fsize);
        Faux_imag.init(1, ori_size, ori_Fsize);
    }
    else{
        transformer = new FourierTransformer(ori_size,ori_size,ori_size);
        Faux_real.init(ori_size, ori_size, ori_Fsize);
        Faux_imag.init(ori_size, ori_size, ori_Fsize);
    }
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        transformer->FourierTransform(Irefs[iclass].wptr(), Faux_real.wptr(), Faux_imag.wptr());
        //FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
        for (int k = 0,kp = 0; k < Faux_real.dimz; k++,kp = k < Faux_real.dimx?k:(k-Faux_real.dimz))
            for (int i = 0,ip = 0;i < Faux_real.dimy; i++,ip = i < Faux_real.dimx?i:(i-Faux_real.dimy))
                for (int j = 0,jp = 0; j < Faux_real.dimx; j++,jp = j)
                {
                    FDOUBLE r = sqrt((FDOUBLE)(kp*kp + ip*ip + jp*jp));
                    if (r < radius)
                        continue;
                    else if (r > radius_p){
                        Faux_real(k, i, j) = 0.;
                        Faux_imag(k, i, j) = 0.;
                    }
                    else
                    {
                        Faux_real(k, i, j) *= 0.5 - 0.5 * cos(rome_pi * (radius_p - r) / width_fmask_edge);
                        Faux_imag(k, i, j) *= 0.5 - 0.5 * cos(rome_pi * (radius_p - r) / width_fmask_edge);
                    }
                }
        transformer->inverseFourierTransform(Faux_real.wptr(), Faux_imag.wptr(), Irefs[iclass].wptr());
    }
    delete transformer;
    Faux_real.fini();Faux_imag.fini();
}


// Apply a solvent flattening to a map
void MAPModel::applySolventFlatten()
{
    // First read solvent mask from disc, or pre-calculate it
    // NOTE : donot read solvent
    Vol<FDOUBLE> Isolvent;
    if (ref_dim == 2)
        Isolvent.init(1, ori_size, ori_size);
    else
        Isolvent.init(ori_size, ori_size, ori_size);
    
    FDOUBLE radius = particle_diameter / (2. * pixel_size);
    FDOUBLE radius_p = radius + width_mask_edge;
    int Isolvent_origin_z = XMIPP_ORIGIN(Isolvent.dimz);
    int Isolvent_origin_y = XMIPP_ORIGIN(Isolvent.dimy);
    int Isolvent_origin_x = XMIPP_ORIGIN(Isolvent.dimx);
    // FOR_ALL_ELEMENTS_IN_ARRAY3D(Isolvent())
    for (int k = 0; k < Isolvent.dimz; k++)
        for (int i = 0; i < Isolvent.dimy; i++)
            for (int j = 0; j < Isolvent.dimx; j++)
            {
                int kp = k + Isolvent_origin_z;
                int ip = i + Isolvent_origin_y;
                int jp = j + Isolvent_origin_x;
                FDOUBLE r = sqrt((FDOUBLE)(kp*kp + ip*ip + jp*jp));
                if (r < radius)
                    Isolvent(k, i, j) = 1.;
                else if (r > radius_p)
                    Isolvent(k, i, j) = 0.;
                else
                    Isolvent(k, i, j) = 0.5 - 0.5 * cos(rome_pi * (radius_p - r) / width_mask_edge );
            }
    
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        
        // Then apply the expanded solvent mask to the map
        for (int i = 0; i < Irefs[iclass].dimzyx; i++)
            (Irefs[iclass])(0,0,i) *= Isolvent(0,0,i);
        
        // Apply a second solvent mask if necessary
        // This may for example be useful to set the interior of icosahedral viruses to a constant density value that is higher than the solvent
        // Invert the solvent mask, so that an input mask can be given where 1 is the masked area and 0 is protein....
        // if (!fn_mask2.contains("None"))
        //     softMaskOutsideMap(mymodel.Iref[iclass], Isolvent2(), true);
    } // end for iclass
    Isolvent.fini();
}

// set_by_ini_high = ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc) )
void MAPModel::updateCurrentResolution(std::vector<VectorOfFDOUBLE>& data_vs_prior_class,bool set_by_ini_high)
{
    static FDOUBLE best_resol_thus_far = 1./999.;
    static int nr_iter_wo_resol_gain = 0;
    
    int maxres = 0;
    if (do_map)
    {
        // Set current resolution
        if (set_by_ini_high)
        {
            maxres = round(ori_size * pixel_size / ini_high);
        }
        else
        {
            // Calculate at which resolution shell the data_vs_prior drops below 1
            int ires;
            for (int iclass = 0; iclass < nr_classes; iclass++)
            {
                for (ires = 1; ires < ori_size/2; ires++)
                {
                    if (data_vs_prior_class[iclass][ires] < 1.)
                        break;
                }
                // Subtract one shell to be back on the safe side
                ires--;
                if (ires > maxres)
                    maxres = ires;
            }
            
            // Never allow smaller maxres than minres_map
            maxres = std::max(maxres, minres_map);
        }
    }
    else
    {
        // If we are not doing MAP-estimation, set maxres to Nyquist
        maxres = ori_size/2;
    }
    FDOUBLE newres = getResolution(maxres, pixel_size, ori_size);
    
    // Check whether resolution improved, if not increase nr_iter_wo_resol_gain
    if (newres <= current_resolution+0.0001) // Add 0.0001 to avoid problems due to rounding error
        nr_iter_wo_resol_gain++;
    else
        nr_iter_wo_resol_gain = 0;
    
    // Store best resolution thus far (but no longer do anything with it anymore...)
    if (newres > best_resol_thus_far)
        best_resol_thus_far = newres;
    
    current_resolution = newres;
}

//
void MAPModel::printResolution(std::vector<VectorOfFDOUBLE>& data_vs_prior_class,bool set_by_ini_high)
{
    if (!do_map || set_by_ini_high) return;
    
    for (int iclass = 0; iclass < nr_classes; iclass++)
    {
        int ires = 0;
        // Calculate at which resolution shell the data_vs_prior drops below 1
        for (ires = 1; ires < ori_size/2; ires++)
        {
            if (data_vs_prior_class[iclass][ires] < 1.)
                break;
        }
        // Subtract one shell to be back on the safe side
        int maxres = ires - 1;
        // Never allow smaller maxres than minres_map
        maxres = std::max(maxres, minres_map);
        FDOUBLE newres = getResolution(maxres, pixel_size, ori_size);
        std::cout<<" class = "<<iclass+1<<" resolution = "<<(float)1./newres<<std::endl;
    }
}

//
void MAPModel::updateImageSizeAndResolutionPointers(VectorOfInt& Npix_per_shell,VectorOfInt& Mresol_coarse,VectorOfInt& Mresol_fine,
                                                    int& coarse_size,int& current_size,
                                                    int adaptive_oversampling,double angularSampling,
                                                    FDOUBLE model_ave_Pmax,bool has_high_fsc_at_limit /*= false*/)
{
    int max_coarse_size = -1;
    // Default max_coarse_size is original size
    if (max_coarse_size < 0)
        max_coarse_size = ori_size;
    
    // Increment the current_size
    // If we are far from convergence (in the initial stages of refinement) take steps of 25% the image size
    // Do this whenever the FSC at the current_size is larger than 0.2, but NOT when this is in combination with very low Pmax values,
    // in the latter case, over-marginalisation may lead to spuriously high FSCs (2 smoothed maps may look very similar at high-res: all zero!)
    //
    int maxres = getPixelFromResolution(current_resolution,pixel_size,ori_size);
    if (model_ave_Pmax > 0.1 && has_high_fsc_at_limit)
    {
        maxres += round(0.25 * ori_size / 2);
    }
    else
    {
        // If we are near our resolution limit, use incr_size (by default 10 shells)
        //		maxres += incr_size;
        maxres += 10;
    }
    
    // Go back from resolution shells (i.e. radius) to image size, which are BTW always even...
    current_size = maxres * 2;
    
    // current_size can never be larger than ori_size:
    current_size = std::min(current_size, ori_size);
    // The current size is also used in wsum_model (in unpacking)
    
    if (adaptive_oversampling > 0.)
    {
        // Dependency of coarse_size on the angular sampling used in the first pass
        FDOUBLE rotated_distance = (angularSampling / 360.) * PI * particle_diameter;
        FDOUBLE keepsafe_factor = (ref_dim == 3) ? 1.2 : 1.5;
        
        FDOUBLE coarse_resolution = rotated_distance / keepsafe_factor;
        // Note coarse_size should be even-valued!
        coarse_size = 2 * ceil(pixel_size * ori_size / coarse_resolution);
        // Coarse size can never be larger than max_coarse_size
        coarse_size = std::min(max_coarse_size, coarse_size);
    }
    else
        coarse_size = current_size;
    
    coarse_size = std::min(current_size, coarse_size);
    
    // Also update the resolution pointers here
    // Calculate number of pixels per resolution shell
    int ori_Fsize = (ori_size/2+1);
    int ori_Fsize2 = ori_size*(ori_size/2+1);
    Npix_per_shell.init(ori_Fsize);
    Npix_per_shell.zero();
    for (int i = 0 ; i<ori_size; i++){
        for (int j = 0 ; j<ori_Fsize; j++)
        {
            int ip = (i < ori_Fsize) ? i : i - ori_size;
            int jp = j;
            int ires = round(sqrt((FDOUBLE)(ip*ip + jp*jp)));
            // TODO: better check for volume_refine, but the same still seems to hold... Half of the yz plane (either ip<0 or kp<0 is redundant at jp==0)
            // Exclude points beyond XSIZE(Npix_per_shell), and exclude half of the x=0 column that is stored twice in FFTW
            if (ires < ori_Fsize && !(jp==0 && ip < 0))
                Npix_per_shell[ires] += 1;
        }
    }
    
    int current_Fsize = (current_size/2+1);
    int current_Fsize2 = current_size*(current_size/2+1);
    Mresol_fine.init(current_Fsize2);
    for (int i = 0 ; i < current_size; i++){
        for (int j = 0 ; j < current_Fsize; j++)
        {
            int ip = (i < current_Fsize) ? i : i - current_size;
            int jp = j;
            int ires = round(sqrt((FDOUBLE)(ip*ip + jp*jp)));
            // TODO: better check for volume_refine, but the same still seems to hold... Half of the yz plane (either ip<0 or kp<0 is redundant at jp==0)
            // Exclude points beyond ires, and exclude and half (y<0) of the x=0 column that is stored twice in FFTW
            if (ires < current_Fsize  && !(jp==0 && ip < 0))
                Mresol_fine[i*current_Fsize+j] = ires;
            else
                Mresol_fine[i*current_Fsize+j] = -1.;
        }
    }
    
    int coarse_Fsize = (coarse_size/2+1);
    int coarse_Fsize2 = coarse_size*(coarse_size/2+1);
    Mresol_coarse.init(coarse_Fsize2);
    for (int i = 0 ; i < coarse_size; i++){
        for (int j = 0 ; j < coarse_Fsize; j++)
        {
            int ip = (i < coarse_Fsize) ? i : i - coarse_size;
            int jp = j;
            int ires = round(sqrt((FDOUBLE)(ip*ip + jp*jp)));
            // Exclude points beyond ires, and exclude and half (y<0) of the x=0 column that is stored twice in FFTW
            // exclude lowest-resolution points
            if (ires < coarse_Fsize && !(jp==0 && ip < 0))
                Mresol_coarse[i*coarse_Fsize+j] = ires;
            else
                Mresol_coarse[i*coarse_Fsize+j] = -1.;
        }
    }
}

//
void MAPModel::printSpaceInfo(std::ostream &os)
{
#define PRINT(NAME,SIZE) os<<std::setw(40)<<NAME<<SIZE<<" GB."<<std::endl;totalSpace+=SIZE;
    double perGB = 1024.*1024.*1024.;
    double totalSpace = 0.;
    os<<"Memory space for MAPModel : "<<std::endl;
    PRINT("Iref : ", sizeof(double)*Irefs[0].dimzyx*Irefs.size()/perGB);
    PRINT("projector.data : ", sizeof(MKL_Complex16)*projector[0].data.dimzyx*projector.size()/perGB);
    PRINT("backprojector.data : ", sizeof(MKL_Complex16)*backprojector[0].data.dimzyx*backprojector.size()/perGB);
    PRINT("backprojector.weight : ", sizeof(double)*backprojector[0].weight.dimzyx*backprojector.size()/perGB);
    PRINT("thread data : ", sizeof(MKL_Complex16)*thread_data[0].dimzyx*thread_data.size()/perGB);
    PRINT("thread weight : ", sizeof(double)*thread_weight[0].dimzyx*thread_weight.size()/perGB);
    os<<std::setw(40)<<"MAPMODEL TotalSpace"<<": "<<totalSpace<<" GB."<<std::endl;
    std::cout<<"----------------------------------------------"<<std::endl;
}

#ifdef USEMPI
void addContiguous(const void *invec, void *inoutvec, int *len,const MPI::Datatype& datatype)
{
    int length = datatype.Get_size()/sizeof(FDOUBLE);
    for (int i = 0; i < length; i++ )
        ((FDOUBLE*)inoutvec)[i] += ((FDOUBLE*)invec)[i];
}
#endif

// reduce Fourier-space sampling data to each node
void MAPModel::reduceData(bool use_allreduce)
{
#ifdef USEMPI
    if (use_allreduce) {
        int nodes = MPI::COMM_WORLD.Get_size();
        int local_temp_size = backprojector[0].data.dimzyx*2;
        FDOUBLE* local_temp = (FDOUBLE*)_mm_malloc(sizeof(FDOUBLE)*local_temp_size,64);
        // TOOD ,use MPI_Reduce_Scatter
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            memcpy(local_temp, (FDOUBLE*)backprojector[iclass].data.wptr(), sizeof(FDOUBLE)*local_temp_size);
            MPI::COMM_WORLD.Allreduce(local_temp,(FDOUBLE*)backprojector[iclass].data.wptr(),local_temp_size,MPI_FDOUBLE,MPI::SUM);
            memcpy(local_temp, backprojector[iclass].weight.wptr(), sizeof(FDOUBLE)*local_temp_size/2);
            MPI::COMM_WORLD.Allreduce(local_temp,backprojector[iclass].weight.wptr(),local_temp_size/2,MPI_FDOUBLE,MPI::SUM);
        }
        _mm_free(local_temp);
        return;
    }
    // 2x large than current_size
    // if image size is 300,this data maybe get to 600*600*600 ~= 1.6 GB
    // but much of time the size will be much smaller than 300,like 50~100..
    int fVolumeSize = backprojector[0].data.dimzyx;
    // size shoud be large enough
    int nodes = MPI::COMM_WORLD.Get_size();
    int node = MPI::COMM_WORLD.Get_rank();
    int ranks = nodes;
    int first_local_class,last_local_class;
    int nr_local_classes = divide_equally(nr_classes, nodes, node, first_local_class, last_local_class);
    // TODO : if senBuf is too large,send 2~3 classes each time...
    MASTERNODE std::cout<<"Try to reduce-scatter backprojector.data...size : "<<(nr_classes*fVolumeSize*2./1024./1024./1024.)<<" GB."<<std::endl;
    // pack data for backprojector.data
    MPI::Datatype fVolumeComplex = MPI_FDOUBLE.Create_contiguous( 2*fVolumeSize );
    fVolumeComplex.Commit();
    MPI::Op MPIContiguousSum;
    MPIContiguousSum.Init((MPI::User_function *)addContiguous, 1);
    FDOUBLE* sendBuf = (FDOUBLE*)_mm_malloc(sizeof(FDOUBLE)*nr_classes*fVolumeSize*2,64); // for complex datatype
    FDOUBLE* recvBuf = (FDOUBLE*)_mm_malloc(sizeof(FDOUBLE)*nr_local_classes*fVolumeSize*2,64);
    int* recvcounts = (int*)_mm_malloc(sizeof(int)*ranks,64);
    for (int iclass = 0; iclass < nr_classes; iclass++)
        memcpy(sendBuf+iclass*fVolumeSize*2, (FDOUBLE*)backprojector[iclass].data.wptr(), sizeof(FDOUBLE)*fVolumeSize*2);
    for	(int rank = 0; rank < ranks; rank++){
        int first_local_class,last_local_class;
        int nr_local_classes = divide_equally(nr_classes, ranks, rank, first_local_class, last_local_class);
        recvcounts[rank] = nr_local_classes;
    }
    //
    MPI::COMM_WORLD.Reduce_scatter(sendBuf, recvBuf, recvcounts, fVolumeComplex, MPIContiguousSum);
    // unpack data for backprojector.data
    for (int iclass = first_local_class; iclass <= last_local_class; iclass++)
        memcpy((FDOUBLE*)backprojector[iclass].data.wptr(), recvBuf+(iclass-first_local_class)*fVolumeSize*2, sizeof(FDOUBLE)*fVolumeSize*2);
    fVolumeComplex.Free();
    // pack data for backprojector.weight
    MPI::Datatype fVolume = MPI_FDOUBLE.Create_contiguous( fVolumeSize );
    fVolume.Commit();
    for (int iclass = 0; iclass < nr_classes; iclass++)
        memcpy(sendBuf+iclass*fVolumeSize, backprojector[iclass].weight.wptr(), sizeof(FDOUBLE)*fVolumeSize);
    //
    MPI::COMM_WORLD.Reduce_scatter(sendBuf, recvBuf, recvcounts, fVolume, MPIContiguousSum);
    // unpack data for backprojector.weight
    for (int iclass = first_local_class; iclass <= last_local_class; iclass++)
        memcpy(backprojector[iclass].weight.wptr(), recvBuf+(iclass-first_local_class)*fVolumeSize, sizeof(FDOUBLE)*fVolumeSize);
    fVolume.Free();
    //
    MPIContiguousSum.Free();_mm_free(sendBuf);_mm_free(recvBuf);_mm_free(recvcounts);
#endif
}

// broadcast Volume data to each nodes
void MAPModel::broadcastData()
{
//#ifdef USEMPI
//    // TODO,use MPI_Allgather
//    int nodes = MPI::COMM_WORLD.Get_size();
//    int node = MPI::COMM_WORLD.Get_rank();
//    int ranks = nodes;
//    int first_local_class,last_local_class;
//    int nr_local_classes = divide_equally(nr_classes, nodes, node, first_local_class, last_local_class);
//    int IrefSize = Irefs[0].dimzyx;
//    MASTERNODE std::cout<<"Try to Allgather Iref...size : "<<(nr_classes*IrefSize/1024./1024./1024.)<<" GB."<<std::endl;
//    MPI::Datatype IrefType = MPI::DOUBLE.Create_contiguous( IrefSize );
//    IrefType.Commit();
//    MPI::Op MPIContiguousSum;
//    MPIContiguousSum.Init((MPI::User_function *)addContiguous, 1);
//    double* sendBuf = (double*)_mm_malloc(sizeof(double)*nr_local_classes*IrefSize,64);
//    double* recvBuf = (double*)_mm_malloc(sizeof(double)*nr_classes*IrefSize,64);
//    int* recvcounts = (int*)_mm_malloc(sizeof(int)*ranks,64);
//    int* displs = (int*)_mm_malloc(sizeof(int)*ranks,64);
//    NODE0ONLY std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
//    // pack the data
//    for (int iclass = first_local_class; iclass <= last_local_class; iclass++)
//        memcpy(sendBuf+(iclass-first_local_class)*IrefSize, Irefs[iclass].wptr(), sizeof(double)*IrefSize);
//    for	(int rank = 0; rank < ranks; rank++){
//        int first_local_class,last_local_class;
//        int nr_local_classes = divide_equally(nr_classes, ranks, rank, first_local_class, last_local_class);
//        recvcounts[rank] = nr_local_classes;
//        if (rank==0) displs[rank] = 0;
//        else displs[rank] = displs[rank-1] + recvcounts[rank-1];
//    }
//    NODE0ONLY std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
//    //
//    MPI::COMM_WORLD.Allgatherv(sendBuf, nr_local_classes, IrefType,recvBuf, recvcounts, displs, IrefType);
//    // unpack the data
//    NODE0ONLY std::cout<<__FILE__<<" "<<__LINE__<<std::endl;
//    for (int iclass = 0; iclass < nr_classes; iclass++)
//        memcpy(Irefs[iclass].wptr(), recvBuf+iclass*IrefSize,sizeof(double)*IrefSize);
//    _mm_free(sendBuf);_mm_free(recvBuf);_mm_free(recvcounts);_mm_free(displs);
//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#if USE_DISCARD
//    int nodes = MPI::COMM_WORLD.Get_size();
//    for (int iclass = 0; iclass < nr_classes; iclass++) {
//        // Only reduce this iclass data to toNode...
//        int fromNode = divide_equally_which_group(nr_classes, nodes, iclass);
//        MPI::COMM_WORLD.Bcast(Iref[iclass].wptr(),Iref[iclass].dimzyx,MPI::DOUBLE,fromNode);
//    }
//#endif
//    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#endif
}

