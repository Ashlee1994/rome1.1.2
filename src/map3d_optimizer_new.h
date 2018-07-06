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
#ifndef MAP3D_OPTIMIZER_ORIGINAL_H_
#define MAP3D_OPTIMIZER_ORIGINAL_H_

// #include "tbb/parallel_sort.h"
#include "./map_optimizer_base_old.h"
#include "./map3d_optimizer_kernel.h"

#define EXP_MWEIGHT_NEW

namespace Map3dOptimizer_new
{
    using namespace MapOptimizerBase_old;
    using namespace Map3dOptimizer_kernel;
    //
    extern DataStream global_data_stream;
    // model
    extern HealpixSampler sampler3d;
    extern MAPModel mapModel;
    extern MLModel mlModel;
    extern ParticleModel particleModel;
    extern HiddenVariableMonitor hiddenVarMonitor;
    
    // image data
    extern Images images;
    extern MetaDataTable metadata;
    
    // sampling
    extern double offset_step;
    extern double offset_range;
    extern int sampler3d_healpix_order;
    extern std::string sampler3d_fn_sym;
    extern SamplingGrid samplingGrid;
    
    // ------------ link  some variable ------------- //
    static auto& particle_diameter = particleModel.particle_diameter;
    static auto& only_flip_phases = particleModel.only_flip_phases;
    static auto& ctf_phase_flipped = particleModel.ctf_phase_flipped;
    static auto& intact_ctf_first_peak = particleModel.intact_ctf_first_peak;
    
    // ------------ variable for expectation step  ------------- //
#define SEP
#define ELTONE(T,N,S1,S2) extern T N;
#define ELTVE1(T,N,S1,S2) extern VectorOfArray1d<T> N;
#define ELTVE2(T,N,S1,S2) extern VectorOfArray2d<T> N;
#define ELTVE3(T,N,S1,S2) extern Aligned3dArray <T> N;
    MAPOPTIMIZER_OLD_EXP_VARS
#undef SEP
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef ELTVE3
    //
    extern VectorOfStruct<MetaDataElem> exp_metadata;
    
    // ---------------   thread variable   ------------ //
#define SEP
#define ELTONE(T,N,S1,S2) extern T N;
#define ELTVE1(T,N,S1,S2) extern VectorOfArray1d<T> N;
#define ELTVE2(T,N,S1,S2) extern VectorOfArray2d<T> N;
#define ELTVE3(T,N,S1,S2) extern VectorOfArray3d<T> N;
    MAPOPTIMIZER_OLD_THREAD_VARS
#undef SEP
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef ELTVE3
    //
    extern VectorOfArray2d<double> thread_exp_Mweight_sub;
    extern VectorOfArray3d<FDOUBLE> thread_Frefctf_real;
    extern VectorOfArray3d<FDOUBLE> thread_Frefctf_imag;
    extern VectorOfArray3d<FDOUBLE> thread_Fimg_real;
    extern VectorOfArray3d<FDOUBLE> thread_Fimg_imag;
    extern VectorOfArray3d<FDOUBLE> thread_Fweight;
    extern VectorOfArray3d<FDOUBLE> thread_wsum_pdf_direction;
    extern VectorOfArray2d< char > threadfake_do_scale_norm_class;
    extern std::vector<std::vector<GridIndex>> thread_exp_max_weight_index;
    //
    extern bool do_local_searching;
    extern Exp_Mcoarse_Rot_significant_new exp_Mcoarse_Rot_significant;
    //  ---------------   setup      ---------------------- //
    void setupMLoptimizer();
    
    void destroyMLoptimizer();
    
    // Interpret command line for the initial start of a run
    void prepare();
    
    // Perform expectation-maximization iterations
    void iterate();
    
    // ------------   EM-Iteration     ----------------- //
    
    void expectation();
    
    // expectation nr_pool image each time
    void expectationSomeParticles();
    
    // (de)allocate memory space for each expectation step
    void prepareExpMap();
    void prepareExpData();
    void endExpData();
    
    // ------------ some function in expectationsomeparticles functon    --------- //
    
    // get all rotated reference  and the significant rotation
    void getReferenceAllOrientations();
    
    // get all reference and images 's squared differences
    void getAllSquaredDifferences(bool do_coarse_search);
    
    // calculates exp_sum_weight and, for adaptive approach, also exp_significant_weight
    void findAllSignificantPoints(
#ifdef EXP_MWEIGHT_NEW
    Exp_Mweight_new& exp_Mweight
#else
    Exp_Mweight_old& exp_Mweight
#endif
    );
    
    // convert all squared difference to weight(P = exp(-x))
    void convertSquaredDifferencesToWeights(
#ifdef EXP_MWEIGHT_NEW
    Exp_Mweight_new& exp_Mweight
#else
    Exp_Mweight_old& exp_Mweight
#endif
    );
    
    // calculate norm_correction,dLL and Pmax
    void storeWeightedSums();
    
    // update mlModel
    void updateOtherParams();
    
    // add all shifted and rotated images back
    void backProjection();
    
    // ------------  Maximization step    ------------ //
    
    void maximization(bool update_tau2_with_fsc = false,
                      bool is_whole_instead_of_half = false);
    
    // Perform the actual reconstructions
    void maximizationReconstructClass(int iclass);
    
    // Updates all other model parameters (besides the reconstructions)
    void maximizationOtherParameters();
    
    // ------------  read and Write files     ------------- //
    
    void readResult();
    void writeResult(bool do_write_sampling = true, bool do_write_data = true, bool do_write_optimiser = true,
                     bool do_write_model = true, bool do_write_mrc = true, int random_subset = -1);
    void checkResult();
    void printMem(int set_nr_pool);
    
    // ------------  some help function     ------------- //
    //
#ifdef USEMPI
    inline MPI::Intracomm getNewIntracomm(){
        static MPI::Intracomm split_half_world;
        static bool initialized = false;
        if (do_split_random_halves) {
            if (!initialized) {
                int color;
                if (node == nodes-1) color = 0;// final node as master node
                else if (node < (nodes-1)/2) color = 1;
                else color = 2;
                split_half_world = MPI::COMM_WORLD.Split(color, node);
                initialized = true;
            }
            if (node==nodes-1) return split_half_world;
            //
            int split_half_nodes = split_half_world.Get_size();
            assert(split_half_nodes==(nodes-1)/2);
        }
        else{
            split_half_world = MPI::COMM_WORLD;
        }
        return split_half_world;
    }
#endif
    //
    // do maximumu likehood in fourier shell increased space
    // for getAllSquaredDifferences() and updateOtherParams()
    // which means the frequency outside the maximum frequency (size/2) will be ignored
    // this transform function is used to rearrange the normal FFT array to fourier shell increased array
    void transformFourierArrayToShellIncreasedCoarse();
    void transformFourierArrayToShellIncreasedFine();

    //
    // use StatusTracer to write data to disk
    // to compare the different.....
    
};


#endif
