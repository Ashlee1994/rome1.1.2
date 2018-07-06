/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Bevin Brett
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

#include "./map3d_optimizer_new.h"
#include "./primitives_for_each_os.h"
#include "./findBestPoints.h"
#include "./checker.h"
#include "./peakFinder.h"
#include "./imagePermutor.h"

#define USE_PROJECTION_CACHE
//#define USE_BEFORE_PROJECTION_CACHE


namespace Map3dOptimizer_new {

    using namespace MapOptimizerBase_old;
    using namespace Map3dOptimizer_kernel;

    DataStream global_data_stream;
    int data_stream_node = 0;
    
    // model
    HealpixSampler sampler3d;
    MAPModel mapModel;
    MLModel mlModel;
    ParticleModel particleModel;
    HiddenVariableMonitor hiddenVarMonitor;
    FourierShellTranslation fourierShellTrans;
    // image data
    Images images;
    MetaDataTable metadata;

    // sampling
    double offset_step;
    double offset_range;
    int sampler3d_healpix_order;
    std::string sampler3d_fn_sym;
    SamplingGrid samplingGrid;

    // ------------ variable for expectation step -------------- //
#define SEP
#define ELTONE(T,N,S1,S2) T N;
#define ELTVE1(T,N,S1,S2) VectorOfArray1d<T> N;
#define ELTVE2(T,N,S1,S2) VectorOfArray2d<T> N;
#define ELTVE3(T,N,S1,S2) Aligned3dArray <T> N;
    MAPOPTIMIZER_OLD_EXP_VARS
#undef SEP
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef ELTVE3
    //
    VectorOfStruct<MetaDataElem> exp_metadata;
    //
	//
    //
#ifdef EXP_MWEIGHT_NEW
    Exp_Mweight_new exp_Mweight_coarse;
    Exp_Mweight_new exp_Mweight_fine;
#else
    Exp_Mweight_old exp_Mweight_coarse;
    Exp_Mweight_old exp_Mweight_fine;
#endif
    //
    bool do_local_searching;
    Exp_Mcoarse_Rot_significant_new exp_Mcoarse_Rot_significant;
    // ------------ thread variable ----------------- //
#define SEP
#define ELTONE(T,N,S1,S2) T N;
#define ELTVE1(T,N,S1,S2) VectorOfArray1d<T> N;
#define ELTVE2(T,N,S1,S2) VectorOfArray2d<T> N;
#define ELTVE3(T,N,S1,S2) VectorOfArray3d<T> N;
    MAPOPTIMIZER_OLD_THREAD_VARS
#undef SEP
#undef ELTONE
#undef ELTVE1
#undef ELTVE2
#undef ELTVE3
    VectorOfArray2d<double> thread_exp_Mweight_sub;
    VectorOfArray3d<FDOUBLE> thread_Frefctf_real;
    VectorOfArray3d<FDOUBLE> thread_Frefctf_imag;
    VectorOfArray3d<FDOUBLE> thread_Fimg_real;
    VectorOfArray3d<FDOUBLE> thread_Fimg_imag;
    VectorOfArray3d<FDOUBLE> thread_Fweight;
    VectorOfArray3d<FDOUBLE> thread_wsum_pdf_direction;
    VectorOfArray2d< char > threadfake_do_scale_norm_class;
    std::vector<std::vector<GridIndex>> thread_exp_max_weight_index;
    //
    static const int const_fix_iover_rot_tile = 8;
	static const int const_max_ipsi_tile      = 8;
	int const_max_iimage_tile;
	//
    class TileCfg{
    public:
        TileCfg(){}
        ~TileCfg(){}
        //
        void choiceTileForBDW()
        {
            //
        }
        //
        void choiceTileForKNL(int oversampling)
        {
			//
        }
        void getTileForFineCC(int N,int& _N_tile,int& _ipsi_tile,int& _iover_rot_tile,int& _iimage_tile,int& _iimage_sub_tile,int& _itrans_sub_tile)
        {
            _N_tile = N;
            _ipsi_tile       = 1;
            _iover_rot_tile  = exp_nr_over_rot;	// other values lead to changes in output		std::min(4,exp_nr_over_rot);
            _iimage_tile     = 1;
            _iimage_sub_tile = 2;
            _itrans_sub_tile = 2;
        }
        // tile size for getAllSquaredDifferencesCoarse.
        void getTileForCoarseSearch(int N, bool can_N_tile, 
			int& _N_tile,int& _ipsi_tile,int& _iover_rot_tile,int& _iimage_tile,int& _iimage_sub_tile,int& _itrans_sub_tile)
        {
//             _N_tile			= std::min(128,N);			// was 512
//             _ipsi_tile		= const_max_ipsi_tile;		// fix
//             _iover_rot_tile = -1;						// not used on coarse searching
//             _iimage_tile	= const_max_iimage_tile;	// fix equal to maxthreads or maxthreads*2
// 
//             const double L2_kb_one_thread	= 256/2*0.9;						// no Hyper-Threading, two thread share one L2 on boardwell
//             const double rot_image_kb		= 2*_ipsi_tile*_N_tile*8/1024.;		// Fref_real,Fref_imag
//             const double per_trans_image_kb	= 4*_N_tile*8/1024.;				// Fimage_real,Fimage_imag,Minvsigma2s,local_ctf
// 			int squared_tile_size;
// #if defined(TRIPLE_TRANSLATION)
//             double itrans_abTable_kb = 2*offset_range/offset_step*_N_tile*8/1024.;	// all positive shift
//             double iover_trans_kb	 = exp_nr_over_trans*_N_tile*8/1024.;			// all over shift
//             // may no need for sub tile
//             squared_tile_size = const_max_iimage_tile;
// #else
//             double per_abTable_kb = 2*_N_tile*8/1024.;// aTable,bTable
//             computeSquaredTile(squared_tile_size, per_abTable_kb+per_trans_image_kb, L2_kb_one_thread-rot_image_kb,"getTileForCoarseSearch",false);
// #endif
//             _iimage_sub_tile = squared_tile_size;
//             _itrans_sub_tile = squared_tile_size;

											//       N etc by experimentation
											// 192 etc by experimentation
#ifdef TRIPLE_TRANSLATION
			_N_tile				= can_N_tile ? 192 : N;
			_ipsi_tile			= const_max_ipsi_tile;
			_iimage_tile		= can_N_tile ?   8 : 5;
			_iimage_sub_tile	= can_N_tile ?   9 : 5;
			_itrans_sub_tile	= can_N_tile ?   1 : 1;
#else
			_N_tile				= can_N_tile ? 192 : N;
			_ipsi_tile			= const_max_ipsi_tile;
			_iimage_tile		= can_N_tile ?  12 : 5;
			_iimage_sub_tile	= can_N_tile ?  15 : 5;
			_itrans_sub_tile	= can_N_tile ?   2 : 1;
#endif
        }
        void getTileForFineSearch(int N,int& _N_tile,int& _ipsi_tile,int& _iover_rot_tile,int& _iimage_tile,int& _iimage_sub_tile,int& _itrans_sub_tile)
        {
            _N_tile = std::min(512,N);
            _ipsi_tile = 1;
            _iover_rot_tile = const_fix_iover_rot_tile;// fix
            _iimage_tile = const_max_iimage_tile;// fix equal to maxthreads or maxthreads*2
            int squared_tile_size;
            double L2_kb_one_thread = 256/2*0.9;// no Hyper-Threading,two thread share one L2 on boardwell
            double rot_image_kb = 2*_iover_rot_tile*_N_tile*8/1024.;// Fref_real,Fref_imag
            double per_trans_image_kb = 4*_N_tile*8/1024.;// Fimage_real,Fimage_imag,Minvsigma2s,local_ctf
#if defined(TRIPLE_TRANSLATION)
            double itrans_abTable_kb = 2*offset_range/offset_step*_N_tile*8/1024.;
            double iover_trans_kb = exp_nr_over_trans*_N_tile*8/1024.;
            // may no need for sub tile
            squared_tile_size = const_max_iimage_tile;
#elif defined(DOUBLE_TRANSLATION)
            double over_trans_kb = exp_nr_over_trans*2*_N_tile*8/1024.;// iover_trans aTable,bTable
            double per_abTable_kb = 2*_N_tile*8/1024.;// aTable,bTable
            computeSquaredTile(squared_tile_size, per_trans_image_kb+per_abTable_kb, L2_kb_one_thread-rot_image_kb-over_trans_kb,"getTileForFineSearch",false);
#else
            double per_abTable_kb = exp_nr_over_trans*2*_N_tile*8/1024.;// iover_trans aTable,bTable
            computeSquaredTile(squared_tile_size, per_trans_image_kb+per_abTable_kb, L2_kb_one_thread-rot_image_kb,"getTileForFineSearch",false);
#endif
            _iimage_sub_tile = squared_tile_size;
            _itrans_sub_tile = squared_tile_size;
        }
        void getTileForUpdateModel(int N,int& _N_tile,int& _ipsi_tile,int& _iover_rot_tile,int& _iimage_tile,int& _iimage_sub_tile,int& _itrans_sub_tile)
        {
            _N_tile = std::min(512,N);
            _ipsi_tile = 1;
            _iover_rot_tile = const_fix_iover_rot_tile;
            _iimage_tile = const_max_iimage_tile;
            int squared_tile_size;
            double L2_kb_one_thread = 256/2*0.9;// no Hyper-Threading,two thread share one L2 on boardwell
            double rot_image_kb = 2*_iover_rot_tile*_N_tile*8/1024.;// Fref_real,Fref_imag
            double per_trans_image_kb = (2+0.8+0.3+0.3)*_N_tile*8/1024.;// Fimag_real,Fimage_imag,part of sigma2_noise/correction_XA/correction_XX
#if defined(TRIPLE_TRANSLATION)
            double itrans_abTable_kb = 2*offset_range/offset_step*_N_tile*8/1024.;
            double iover_trans_kb = exp_nr_over_trans*_N_tile*8/1024.;
            // may no need for sub tile
            squared_tile_size = const_max_iimage_tile;
#elif defined(DOUBLE_TRANSLATION)
            double over_trans_kb = exp_nr_over_trans*2*_N_tile*8/1024.;// iover_trans aTable,bTable
            double per_abTable_kb = 2*_N_tile*8/1024.;// aTable bTable
            computeSquaredTile(squared_tile_size, per_trans_image_kb+per_abTable_kb, L2_kb_one_thread-rot_image_kb-over_trans_kb,"getTileForUpdateModel",false);
#else
            double per_abTable_kb = 2*_N_tile*8/1024.;// aTable,bTable
            computeSquaredTile(squared_tile_size, per_trans_image_kb+per_abTable_kb, L2_kb_one_thread-rot_image_kb,"getTileForUpdateModel",false);
#endif
            _iimage_sub_tile = squared_tile_size;
            _itrans_sub_tile = squared_tile_size;
        }
        void getTileForBackproject(int N,int& _N_tile,int& _ipsi_tile,int& _iover_rot_tile,int& _iimage_tile,int& _iimage_sub_tile,int& _itrans_sub_tile)
        {
            _N_tile = std::min(512,N);
            _ipsi_tile = 1;
            _iover_rot_tile = const_fix_iover_rot_tile;
            _iimage_tile = const_max_iimage_tile;
            int squared_tile_size;
            double L2_kb_one_thread = 256/2*0.9;// no Hyper-Threading,two thread share one L2 on boardwell
            double rot_image_kb = 3*_iover_rot_tile*_N_tile*8/1024.;// Fref_real,Fref_imag,Fweight
            double per_trans_image_kb = 3*_N_tile*8/1024.;// Fimag_real,Fimage_imag,Minvsigma2s_X_Fctfs2
#if defined(TRIPLE_TRANSLATION)
            double itrans_abTable_kb = 2*offset_range/offset_step*_N_tile*8/1024.;
            double iover_trans_kb = exp_nr_over_trans*_N_tile*8/1024.;
            // may no need for sub tile
            squared_tile_size = const_max_iimage_tile;
#elif defined(DOUBLE_TRANSLATION)
            double over_trans_kb = exp_nr_over_trans*2*_N_tile*8/1024.;// iover_trans aTable,bTable
            double per_abTable_kb = 2*_N_tile*8/1024.;// aTable bTable
            computeSquaredTile(squared_tile_size, per_trans_image_kb+per_abTable_kb, L2_kb_one_thread-rot_image_kb-over_trans_kb,"getTileForBackproject",false);
#else
            double per_abTable_kb = 2*_N_tile*8/1024.;// aTable,bTable
            computeSquaredTile(squared_tile_size, per_trans_image_kb+per_abTable_kb, L2_kb_one_thread-rot_image_kb,"getTileForBackproject",false);
#endif
            _iimage_sub_tile = squared_tile_size;
            _itrans_sub_tile = squared_tile_size;
        }
    private:
        void computeSquaredTile(int &tile_size,double one_cell_kb,double avaiable_L2_kb,const char* name=nullptr,bool show=false)
        {
            tile_size = 1;
            for (; tile_size < const_max_iimage_tile; tile_size++) {
                double cell_kb = one_cell_kb*tile_size;
                if (cell_kb > avaiable_L2_kb){
                    if (show) {
                        std::cout<<name<<",sub_tile : "<<tile_size<<",cell_kb : "<<cell_kb
                                <<" KB, L2 avaiable : "<<avaiable_L2_kb<<" KB."<<std::endl;}
                    break;
                }
            }
            tile_size = tile_size==1?1:tile_size-1;
        }
    };
    //
void setupMLoptimizer()
{
#ifdef USEMPI
    // MPI::Init();//init mpi outside
    nodes = MPI::COMM_WORLD.Get_size();
    node = MPI::COMM_WORLD.Get_rank();
    MPI::Get_processor_name(nodeName,nodeNameLen);
    // std::cout<<nodeName<<" "<<node<<" "<<nodes<<std::flush<<std::endl;
#else
    nodes = 1;
    node = 0;
#endif
    //
    ERROR_CHECK(do_split_random_halves && (nodes==1 || nodes%2==0), "do_split_random_halves,Please set MPI Rank a odd number(>=3)");
    
    if (continue_fn!="NULL") {
        if (continue_fn.find("_optimiser.star")!=std::string::npos)
            star_fn = continue_fn.substr(0,continue_fn.find("_optimiser"))+"_data.star";
        else if(continue_fn.find("_backup.back")!=std::string::npos)
            star_fn = continue_fn.substr(0,continue_fn.find("_backup"))+"_data.star";
        else
            ERROR_REPORT("Wrong continue file name "+continue_fn+",use *_optimiser.star or *_backup.back");
        iter = atoi( continue_fn.substr(continue_fn.find("_it")+3,continue_fn.find_last_of('_')-continue_fn.find("_it")-3).c_str() );
        NODE0ONLY std::cout<<"###### Continue starting from iter "<<iter<<std::endl;
    }

#ifdef DATA_STREAM
    global_data_stream.init(data_stream_out_fn, data_stream_in_fn);
    if (do_split_random_halves) {
        if (node==data_stream_node) global_data_stream.doDataStream=true;
        else global_data_stream.doDataStream=false;
    }
    global_data_stream.foutInt(&nr_iter, 1, "nr_iter", __FILE__, __LINE__);
    global_data_stream.foutInt(&nr_classes, 1, "nr_classes", __FILE__, __LINE__);
    global_data_stream.foutDouble(&pixel_size, 1, "pixel_size", __FILE__, __LINE__);
    global_data_stream.foutInt(&random_seed, 1, "random_seed", __FILE__, __LINE__);
    global_data_stream.foutDouble(&ini_high, 1, "ini_high", __FILE__, __LINE__);
    global_data_stream.foutDouble(&tau2_fudge_factor, 1, "tau2_fudge_factor", __FILE__, __LINE__);
    global_data_stream.foutDouble(&particle_diameter, 1, "particle_diameter", __FILE__, __LINE__);
    global_data_stream.foutInt(&adaptive_oversampling, 1, "adaptive_oversampling", __FILE__, __LINE__);
    global_data_stream.foutInt(&sampler3d_healpix_order, 1, "sampling.healpix_order", __FILE__, __LINE__);
    global_data_stream.foutDouble(&sampler3d.psi_step, 1, "sampling.psi_step", __FILE__, __LINE__);
    global_data_stream.foutDouble(-91, "sampling.limit_tilt", __FILE__, __LINE__);
    global_data_stream.foutDouble(&offset_range, 1, "sampling.offset_range", __FILE__, __LINE__);
    global_data_stream.foutDouble(&offset_step, 1, "sampling.offset_step", __FILE__, __LINE__);
    global_data_stream.foutDouble(0.5, "sampling.perturbation_factor", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
    //
    do_local_searching = (sigma2_angle!=0);
    // ------------ initialize sampling ---------- //
    sampler3d.initialize(offset_step,offset_range,-1,sampler3d_healpix_order,sampler3d_fn_sym,do_local_searching?PRIOR_ROTTILT_PSI:NOPRIOR);

    maxthreads = omp_get_max_threads_or_one_if_no_openmp();
    
    const_max_iimage_tile = maxthreads*2;
    NODE0ONLY std::cout<<"Threads number is "<<maxthreads<<",max image tile size is "<<const_max_iimage_tile<<"?"<<std::endl;
    
    //omp_set_num_threads(8);
    if(debug_flag) NODE0ONLY std::cout<<"Number of threads = "<<maxthreads<<std::endl;
    
    metadata.readFromStar(star_fn);

    Mrcs::MrcsHead mrcsHead;
    Mrcs::readMrcsHead(metadata[0].IMAGE.NAME, mrcsHead);
    NODE0ONLY std::cout<<"The image size is "<<mrcsHead.NC<<"x"<<mrcsHead.NC<<",it store as mode "<<mrcsHead.MODE<<"(2:float) in *.mrcs file."<<std::endl;

    ori_size = mrcsHead.NC;

	metadata.shuffle(random_seed,do_split_random_halves);
    setImagesDistribution(nr_global_images, nr_local_images, nodes, node, first_local_image, last_local_image, metadata, do_split_random_halves);
    if (do_split_random_halves) metadata.selectHalfMetadata(node<(nodes-1)/2?1:2);
    
    NODE0ONLY std::cout<<"There are "<<nr_local_images<<" images process on node 0."<<std::endl;

    // read the local images data
    images.resize(nr_local_images);
    for (auto& image : images) image.init(ori_size*ori_size);

    // All nodes read the file at once rather than broadcasting
    float *buffer = (float*)aMalloc(sizeof(float)*ori_size*ori_size,64);

    std::string preFileName = "";
	FILE* mrcsFile = nullptr;
	
    NODE0ONLY std::cout<<"Start reading image data, it my take a while....."<<std::endl;
    
    for (int iimage = 0; iimage < nr_local_images; iimage++) {

		if (metadata[iimage + first_local_image].IMAGE.NAME != preFileName)
		{
			if (mrcsFile) fclose(mrcsFile);
			preFileName = metadata[iimage + first_local_image].IMAGE.NAME;
			mrcsFile = fopen(preFileName.c_str(), "rb");
			if (NULL == mrcsFile)
				ERROR_REPORT("make sure you have " + preFileName + " accessible from " + currentDirectory());
		}

        //change type to int may cause bug,avoid offset out of range
        long image_id = metadata[iimage+first_local_image].IMAGE.INDEX;

        long offset = (256+(image_id-1)*ori_size*ori_size)*sizeof(float);

        //std::cout<<image_id<<" "<<offset<<std::endl;

        fseek(mrcsFile,offset,SEEK_SET);

        if(fread((char*)buffer,ori_size*ori_size*sizeof(float),1,mrcsFile) == NULL){
            std::cerr<<"read file failed."<<std::endl;
            EXIT_ABNORMALLY;
        }

        auto image_data = images[iimage].wptr(ori_size*ori_size);
        for (int i = 0; i < ori_size*ori_size; i++) {
            image_data[i] = buffer[i];
        }

    }

    aFree(buffer);

    if(mrcsFile) fclose(mrcsFile);
    
    // Also randomize random-number-generator for perturbations on the angles
    dontShare_Random_generator.init(random_seed);
}
    
void prepare()
{
#ifdef USEMPI //
    int split_half_world_node = -1;
    const auto split_half_world = getNewIntracomm();
    split_half_world_node = split_half_world.Get_rank();
    if (do_split_random_halves && node == nodes-1) {
        assert(split_half_world_node == 0);
        return;
    }
#endif
    // --------------- check pixel size ------------------ //
    if (metadata[0].MAGNIFICATION!=0 && metadata[0].DETECTOR_PIXEL_SIZE!=0) {
        double new_pixel_size = 10000. * metadata[0].DETECTOR_PIXEL_SIZE / metadata[0].MAGNIFICATION;
        if (fabs(new_pixel_size-pixel_size)>0.01) {
            NODE0ONLY std::cout<<"MODIFYING pixel size from "<<pixel_size<<" to "<<new_pixel_size
            			<<" based on magnification information in the input STAR file"<<std::endl;
            pixel_size = new_pixel_size;
        }
        for (int iimage = 1; iimage < nr_local_images; iimage++) {
            double my_pixel_size = 10000. * metadata[iimage].DETECTOR_PIXEL_SIZE / metadata[iimage].MAGNIFICATION;
            ERROR_CHECK(my_pixel_size != new_pixel_size, "ERROR inconsistent magnification and detector pixel sizes in images in input STAR file");
        }
    }
    // --------------- initialize MAP Model ---------------- //
    if (particle_diameter < 0.)
        particle_diameter = (ori_size - width_mask_edge) * pixel_size;
    
    // --------------- initialize Particle Model ----------------- //
    particleModel.initialize(ori_size, pixel_size, particle_diameter, width_mask_edge,
                             sigma2_fudge, random_seed, do_norm_correction, do_zero_mask,
                             do_shifts_onthefly,maxthreads,&global_data_stream);

    // ini_high ,user set
    mapModel.initialize(nr_classes, ori_size, particle_diameter, pixel_size, ini_high, maxthreads,
                        width_mask_edge, width_fmask_edge, 5, true, do_map, 3, sampler3d_fn_sym);
    mapModel.initializeRef(ref_fn);

    // ------------ initialize model and wsum --------------- //
    mlModel.initialize(ori_size, nr_classes, metadata.numberOfGroups(), sampler3d.NrDir(), sigma2_angle);

    if (mlModel.nr_groups == 1) do_scale_correction = false;
    NODE0ONLY {
        std::cout<<"The data will ";
        if (do_scale_correction) std::cout<<"do";else std::cout<<"not do";
        std::cout<<" scale correction."<<std::endl;
    }
    // ------------ initialize model ----------- //

    // Calculate initial sigma noise model from power_class spectra of the individual images
    Image Mavg; Mavg.init(ori_size*ori_size); Mavg.zero();
    // check the data whether normalized
    // checkNormalize(images_data, ori_size, nr_local_images, particle_diameter, pixel_size);

    // initialize Mavg,wsum_sigma2_noise,wsum_sumw_group,
    mlModel.calculateSumOfPowerSpectraAndAverageImage(Mavg, images, do_zero_mask, metadata, first_local_image, mapModel);

#ifdef USEMPI
    // TODO : some diff after 3 iteration when use "--scale"
    // reduce data
    int local_temp_size = std::max(mlModel.nr_groups,ori_size*ori_size);
    double *local_temp = (double*)aMalloc(sizeof(double)*local_temp_size,64);
    
    auto MPI_Intracomm_Self_Reduce = [&](void *v2, int v3, const MPI::Datatype &v4, const MPI::Op &v5, int v6){
        assert(v4==MPI::DOUBLE);
        memcpy(local_temp, v2, sizeof(double)*v3);
        if (do_split_random_halves) split_half_world.Reduce(local_temp, v2, v3, v4, v5, v6);
        else MPI::COMM_WORLD.Reduce(local_temp, v2, v3, v4, v5, v6);
    };
    //
    MPI_Intracomm_Self_Reduce(Mavg.wptrAll(), ori_size*ori_size, MPI::DOUBLE, MPI::SUM, 0);
    // NOTE : because we only get average sigma2_noise,do not need to reduce all nr_groups's sigma2 noise
    MPI_Intracomm_Self_Reduce(mlModel.wsum_sigma2_noise[0].wptrAll(), ori_size/2+1, MPI::DOUBLE, MPI::SUM, 0);
	//
    MPI_Intracomm_Self_Reduce(mlModel.wsum_sumw_group.wptrAll(), mlModel.nr_groups, MPI::DOUBLE, MPI::SUM, 0);
    
    aFree(local_temp);

    int* local_temp2 = (int*)aMalloc(sizeof(int)*mlModel.nr_groups,64);
    
    auto MPI_Intracomm_Self_AllReduce = [&](void *v2, int v3, const MPI::Datatype &v4, const MPI::Op &v5){
        assert(v4==MPI::INT);
        memcpy(local_temp2, v2, sizeof(int)*v3);
        if (do_split_random_halves) split_half_world.Allreduce(local_temp2, v2, v3, v4, v5);
        else MPI::COMM_WORLD.Allreduce(local_temp2, v2, v3, v4, v5);
    };
    //
    MPI_Intracomm_Self_AllReduce(mlModel.nr_particles_group.wptrAll(), mlModel.nr_groups, MPI::INT, MPI::SUM);

    aFree(local_temp2);
    if (node == 1) {
        for (int i = 0; i < mlModel.nr_groups; i++) {
            std::cerr<<mlModel.nr_particles_group.rptrAll()[i]<<" ";
        }
        std::cerr<<std::endl;
    }
#endif

#ifdef DATA_STREAM
    //TODO : check why different with relion.
    global_data_stream.foutInt(metadata.numberOfGroups(), "numberOfGroups", __FILE__, __LINE__);
    global_data_stream.foutInt(metadata.numberOfMicrographs(), "numberOfMicrographs", __FILE__, __LINE__);
    global_data_stream.foutInt(metadata.numberOfParticles(), "numberOfOriginalParticles", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), ori_size*ori_size*ori_size, "ref_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), ori_size*ori_size*ori_size, "ref_1", __FILE__, __LINE__);
    global_data_stream.foutInt(mlModel.orientational_prior_mode, "mymodel.orientational_prior_mode!!!", __FILE__, __LINE__);
    global_data_stream.foutInt(mapModel.ref_dim, "mymodel.ref_dim", __FILE__, __LINE__);
    global_data_stream.foutInt(sampler3d.NrDir(), "sampling.NrDirections()", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[0].wptr(sampler3d.NrDir()), sampler3d.NrDir(), "mymodel.pdf_direction[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[nr_classes-1].wptr(sampler3d.NrDir()), sampler3d.NrDir(), "mymodel.pdf_direction_nr_class", __FILE__, __LINE__);

    global_data_stream.foutDouble(Mavg.wptrAll(), ori_size*ori_size, "Mavg", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_noise[0].wptr(ori_size/2+1), ori_size/2+1, "wsum_model_sigma2_noise", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[0], "wsum_model_sumw_group1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[mlModel.nr_groups-1], "wsum_model_sumw_groupN", __FILE__, __LINE__);
    global_data_stream.foutInt(mlModel.nr_particles_group.wptrAll()[0], "mymodel_sumw_nr_particles_group1", __FILE__, __LINE__);
    global_data_stream.foutInt(mlModel.nr_particles_group.wptrAll()[mlModel.nr_groups-1], "mymodel_sumw_nr_particles_groupN", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif

    for (int i = 0; i < ori_size*ori_size; i++)
        Mavg[i] /= mlModel.wsum_sumw_group.rptrAll()[0];

    // Set model_sigma2_noise and model_Iref from averaged poser spectra and Mavg
    NODEONLY(split_half_world_node,0) mlModel.setSigmaNoiseEstimates(Mavg);

#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.sigma2_noise[0].wptr(ori_size/2+1), ori_size/2+1, "sigma2_noise1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.sigma2_noise[mlModel.nr_groups-1].wptr(ori_size/2+1), ori_size/2+1, "sigma2_noiseN", __FILE__, __LINE__);
    global_data_stream.foutInt(ori_size/2+1, "sigma2_size", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif

    // First low-pass filter the initial references
    mapModel.applyLowPassFilter();

#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), ori_size*ori_size*ori_size, "ref_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), ori_size*ori_size*ori_size, "ref_N", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif

    // Initialise the model_data_versus_prior ratio to get the initial current_size right
    // model_tau2_class,model_data_vs_prior_class
    NODEONLY(split_half_world_node,0) mlModel.initialiseDataVersusPrior(mapModel,tau2_fudge_factor); // fix_tau was set in initialiseGeneral

#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[0].wptr(ori_size/2+1), ori_size/2+1, "data_vs_prior_class_1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "data_vs_prior_class_N", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif

#ifdef USEMPI
    auto MPI_Intracomm_Bcast = [&](void *v1, int v2, const MPI::Datatype &v3, int v4){
        if (do_split_random_halves) split_half_world.Bcast(v1, v2, v3, v4);
        else MPI::COMM_WORLD.Bcast(v1, v2, v3, v4);
    };

    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        MPI_Intracomm_Bcast(mlModel.sigma2_noise[igroup].wptrAll(),mlModel.ori_Fsize,MPI::DOUBLE,0);

    for (int iclass = 0; iclass < nr_classes; iclass++){
        // in 3D case,only read reference from mrc file
        // MPI::COMM_WORLD.Bcast(mapModel.Irefs[iclass].wptr(),mapModel.Irefs[iclass].dimzyx,MPI::DOUBLE,0);
        MPI_Intracomm_Bcast(mlModel.tau2_class[iclass].wptrAll(),mlModel.ori_Fsize,MPI::DOUBLE,0);
        MPI_Intracomm_Bcast(mlModel.data_vs_prior_class[iclass].wptrAll(),mlModel.ori_Fsize,MPI::DOUBLE,0);
    }
#endif

    Mavg.fini();
}

void destroyMLoptimizer()
{
    images			.resize(0);
    mapModel		.finalize();
    mlModel			.finalize();
    particleModel	.finalize();
}

StatusTracer statusTracer;

void setupStatusTracer()
{
    statusTracer.clear();
    for (int i = 0; i < mapModel.Irefs.size(); i++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mapModel.Irefs[i].wptr(), mapModel.Irefs[i].dimzyx, "mapmodel_Irefs", true);
#else
        statusTracer.appendDoublePtr(mapModel.Irefs[i].wptr(), mapModel.Irefs[i].dimzyx, "mapmodel_Irefs");
#endif
    }
    statusTracer.appendDoublePtr(&mapModel.current_resolution, 1, "mapModel_current_resolution");
    //
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_sigma2_noise", true);
        statusTracer.appendFloatPtr(mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_signal_product_spectra", true);
        statusTracer.appendFloatPtr(mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_reference_power_spectra", true);
#else
        statusTracer.appendDoublePtr(mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_sigma2_noise");
        statusTracer.appendDoublePtr(mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_signal_product_spectra");
        statusTracer.appendDoublePtr(mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_reference_power_spectra");
#endif
    }
    for (int iclass = 0; iclass < mlModel.nr_classes; iclass++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mlModel.tau2_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_tau2_class", true);
        statusTracer.appendFloatPtr(mlModel.sigma2_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_sigma2_class", true);
        statusTracer.appendFloatPtr(mlModel.data_vs_prior_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_data_vs_prior_class", true);
        statusTracer.appendFloatPtr(mlModel.fsc_halves_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_fsc_halves_class", true);
        statusTracer.appendFloatPtr(mlModel.pdf_direction[iclass].wptr(mlModel.nr_directions), mlModel.nr_directions, "mlmodel_pdf_direction", true);
#else
        statusTracer.appendDoublePtr(mlModel.tau2_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_tau2_class");
        statusTracer.appendDoublePtr(mlModel.sigma2_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_sigma2_class");
        statusTracer.appendDoublePtr(mlModel.data_vs_prior_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_data_vs_prior_class");
        statusTracer.appendDoublePtr(mlModel.fsc_halves_class[iclass].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_fsc_halves_class");
        statusTracer.appendDoublePtr(mlModel.pdf_direction[iclass].wptr(mlModel.nr_directions), mlModel.nr_directions, "mlmodel_pdf_direction");
#endif
    }

#if defined(FLOAT_PRECISION)
    statusTracer.appendFloatPtr(mlModel.scale_correction.wptr(mlModel.nr_groups), mlModel.nr_groups, "mlmodel_scale_correction", true);
    statusTracer.appendFloatPtr(mlModel.prior_offsetx_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_prior_offsetx_class", true);
    statusTracer.appendFloatPtr(mlModel.prior_offsety_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_prior_offsety_class", true);
    statusTracer.appendFloatPtr(mlModel.pdf_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_pdf_class", true);
    statusTracer.appendFloatPtr(&mlModel.avg_norm_correction, 1, "mlmodel_avg_norm_correction", true);
    statusTracer.appendFloatPtr(&mlModel.ave_Pmax, 1, "mlmodel_ave_Pmax", true);
    statusTracer.appendFloatPtr(&mlModel.sigma2_offset, 1, "mlmodel_sigma2_offset", true);
#else
    statusTracer.appendDoublePtr(mlModel.scale_correction.wptr(mlModel.nr_groups), mlModel.nr_groups, "mlmodel_scale_correction");
    statusTracer.appendDoublePtr(mlModel.prior_offsetx_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_prior_offsetx_class");
    statusTracer.appendDoublePtr(mlModel.prior_offsety_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_prior_offsety_class");
    statusTracer.appendDoublePtr(mlModel.pdf_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_pdf_class");
    statusTracer.appendDoublePtr(&mlModel.avg_norm_correction, 1, "mlmodel_avg_norm_correction");
    statusTracer.appendDoublePtr(&mlModel.ave_Pmax, 1, "mlmodel_ave_Pmax");
    statusTracer.appendDoublePtr(&mlModel.sigma2_offset, 1, "mlmodel_sigma2_offset");
#endif

    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_sigma2_noise", true);
#else
        statusTracer.appendDoublePtr(mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "mlmodel_wsum_sigma2_noise");
#endif
    }
    
#if defined(FLOAT_PRECISION)
    statusTracer.appendFloatPtr(mlModel.wsum_sumw_group.wptr(mlModel.nr_groups), mlModel.nr_groups, "mlmodel_wsum_sumw_group", true);
    statusTracer.appendFloatPtr(mlModel.wsum_pdf_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_pdf_class", true);
    statusTracer.appendFloatPtr(mlModel.wsum_prior_offsetx_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_prior_offsetx_class", true);
    statusTracer.appendFloatPtr(mlModel.wsum_prior_offsety_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_prior_offsety_class", true);
#else
    statusTracer.appendDoublePtr(mlModel.wsum_sumw_group.wptr(mlModel.nr_groups), mlModel.nr_groups, "mlmodel_wsum_sumw_group");
    statusTracer.appendDoublePtr(mlModel.wsum_pdf_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_pdf_class");
    statusTracer.appendDoublePtr(mlModel.wsum_prior_offsetx_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_prior_offsetx_class");
    statusTracer.appendDoublePtr(mlModel.wsum_prior_offsety_class.wptr(mlModel.nr_classes), mlModel.nr_classes, "mlmodel_wsum_prior_offsety_class");
#endif
    
    for (int iclass = 0; iclass < mlModel.nr_classes; iclass++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(mlModel.wsum_pdf_direction[iclass].wptr(mlModel.nr_directions), mlModel.nr_directions, "mlmodel_wsum_pdf_direction", true);
#else
        statusTracer.appendDoublePtr(mlModel.wsum_pdf_direction[iclass].wptr(mlModel.nr_directions), mlModel.nr_directions, "mlmodel_wsum_pdf_direction");
#endif
    }
    //
#if defined(FLOAT_PRECISION)
    statusTracer.appendFloatPtr(&sampler3d.random_perturbation,1,"sampling3d.random_perturbation",true);
#else
    statusTracer.appendDoublePtr(&sampler3d.random_perturbation,1,"sampling3d.random_perturbation");
#endif
    //
    for (int iimage = 0; iimage < nr_global_images; iimage++) {
#if defined(FLOAT_PRECISION)
        statusTracer.appendFloatPtr(&metadata[iimage].NORM, 1, "metadata["+num2str(iimage,6)+"].NORM", true);
        statusTracer.appendFloatPtr(&metadata[iimage].XOFF, 1, "metadata["+num2str(iimage,6)+"].XOFF", true);
        statusTracer.appendFloatPtr(&metadata[iimage].YOFF, 1, "metadata["+num2str(iimage,6)+"].YOFF", true);
        statusTracer.appendFloatPtr(&metadata[iimage].PSI, 1, "metadata["+num2str(iimage,6)+"].PSI", true);
#else
        statusTracer.appendDoublePtr(&metadata[iimage].NORM, 1, "metadata["+num2str(iimage,6)+"].NORM");
        statusTracer.appendDoublePtr(&metadata[iimage].XOFF, 1, "metadata["+num2str(iimage,6)+"].XOFF");
        statusTracer.appendDoublePtr(&metadata[iimage].YOFF, 1, "metadata["+num2str(iimage,6)+"].YOFF");
        statusTracer.appendDoublePtr(&metadata[iimage].PSI, 1, "metadata["+num2str(iimage,6)+"].PSI");
#endif
        statusTracer.appendIntPtr(&metadata[iimage].CLASS, 1, "metadata["+num2str(iimage,6)+"].CLASS");
    }
}

// ------------------------- EM-Iteration  ------------------------- //

void iterate()
{

	// Update the current resolution and image sizes, and precalculate resolution pointers
	// The rest of the time this will be done after maximization and before writing output files,
	// so that current resolution is in the output files of the current iteration
    bool set_by_ini_high = ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc) );
    static int nr_iter_wo_resol_gain = 0;
	mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
#endif

    setupStatusTracer();

    // continue
    if (continue_fn!="NULL") {
        readResult();
        // After the first iteration the references are always CTF-corrected
        if (do_ctf_correction)
            refs_are_ctf_corrected = true;
    }

	bool has_already_reached_convergence = false;
	for (iter = iter + 1; iter <= nr_iter; iter++)
    {
        NODE0ONLY std::cout<<"Start iteration "<<iter<<",Current resolution is "<<mapModel.current_resolution<<"."<<std::endl;
		if(debug_flag) NODE0ONLY showPopulation();

#ifdef DATA_STREAM
        global_data_stream.foutInt(iter, "iterate()_iter", __FILE__, __LINE__);
#endif
        
        const double starttime = dtime();
        // update coarse_size,current_size,Npix_per_shell,Mresol_coarse,Mresol_fine
        double angularSampler = sampler3d.getAngularSampling();
        mapModel.updateImageSizeAndResolutionPointers(	Npix_per_shell, Mresol_coarse, Mresol_fine,coarse_size, current_size,
                                                    	adaptive_oversampling, angularSampler,mlModel.ave_Pmax,
                                                      	false, false, debug_flag);
        
#ifdef DATA_STREAM
        global_data_stream.foutInt(current_size, "updateImageSizeAndResolutionPointers()_current_size", __FILE__, __LINE__);
        global_data_stream.foutInt(coarse_size, "updateImageSizeAndResolutionPointers()_coarse_size", __FILE__, __LINE__);
        global_data_stream.foutInt(Npix_per_shell.wptrAll(), (ori_size/2+1), "updateImageSizeAndResolutionPointers()_Npix_per_shell", __FILE__, __LINE__);
        global_data_stream.foutInt(Mresol_fine.wptrAll(), current_size*(current_size/2+1), "updateImageSizeAndResolutionPointers()_Mresol_fine", __FILE__, __LINE__);
        global_data_stream.foutInt(Mresol_coarse.wptrAll(), coarse_size*(coarse_size/2+1), "updateImageSizeAndResolutionPointers()_Mresol_coarse", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
        
        NODE0ONLY std::cout<<"The image size on coarse and fine searching step is "<<coarse_size<<" and "<<current_size<<"."<<std::endl;

        // set the exp_metadata and exp_image_data and some other data
        // allocate the maximum data
        prepareExpMap();
        prepareExpData();

		expectation();


#ifdef USEMPI

        mlModel.reduceData(MPI::COMM_WORLD);

//#define DO_RECONSTRUCT_EACH_NODE
#ifdef DO_RECONSTRUCT_EACH_NODE
        mapModel.reduceData(MPI::COMM_WORLD,false);
#else
        mapModel.reduceData(MPI::COMM_WORLD);
#endif
        gatherMetaDataToMaster(metadata);

        MPI::COMM_WORLD.Barrier();

#endif

		maximization();

		// Apply masks to the reference images
		// At the last iteration, do not mask the map for validation purposes
        if(do_solvent && !has_converged)
            mapModel.applySolventFlatten(mask_fn);

#ifdef DATA_STREAM
        global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), mapModel.Irefs[0].dimzyx, "iterate()_Iref[0]", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), mapModel.Irefs[nr_classes-1].dimzyx, "iterate()_Iref[nr_classes-1]", __FILE__, __LINE__);
#endif

		// Re-calculate the current resolution, do this before writing to get the correct values in the output files
        bool set_by_ini_high = ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc) );
		mapModel.updateCurrentResolution(mlModel.data_vs_prior_class,set_by_ini_high,nr_iter_wo_resol_gain);
        mapModel.printResolution(mlModel.data_vs_prior_class,set_by_ini_high,debug_flag);

#ifdef DATA_STREAM
        global_data_stream.foutDouble(mapModel.current_resolution, "updateCurrentResolution()_current_resolution", __FILE__, __LINE__);
#endif

        // Write output files
        NODE0ONLY
        {
            checkResult();
            writeResult();
        }

		// Show more stats
		if(debug_flag) NODE0ONLY  PerformanceCounter::showAll(std::cout, false);
		if(debug_flag) NODE0ONLY  PerformanceCounter::showAll(std::cerr, true);

        // update the metadata,free exp_metadata and exp_image_data
        endExpData();
        const double endtime = dtime();
        NODE0ONLY std::cout<<"**** iteration "<<iter<<" completed in "<<endtime-starttime<<" seconds ****"<<std::endl<<std::flush;
    }

	if(debug_flag) NODE0ONLY showPopulation();
}

void prepareExpMap()
{
#ifdef DATA_STREAM
    // ------------ initialize model and wsum ---------------- //
    global_data_stream.foutDouble(mlModel.tau2_class[0].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.tau2_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif

    int pad_size = 2*(2*std::min(current_size/2, ori_size/2)+1)+1;
    int pad_Fsize = pad_size/2+1;
    if (debug_flag) NODE0ONLY {
        std::cout<<"setFourierTransformMaps, "<<__FILE__<<" , "<<__LINE__<<std::endl;
        std::cout<<"Create temporary threads data for 3D Volume."<<std::endl;
        std::cout<<"Memory needed for threads 3D Volume data ~= "<<maxthreads*nr_classes*(pad_size*pad_size*pad_Fsize/1024./1024./1024.)<<" GB."<<std::endl;
    }
    // Waiting to complete computeFourierTransformMap in Reference class
    mapModel.setFourierTransformMaps(mlModel.tau2_class, current_size);
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.tau2_class[0].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.tau2_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "expectationSetup()_tau2_class", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}

void prepareExpData()
{
    // Initialise all weighted sums to zero
    mlModel.resetZero();

    dontShare_Random_generator.init(random_seed + iter);
    // Reset the random perturbation for this sampling
    sampler3d.resetRandomlyPerturbedSampling();

    // the first time we use small nr_pool to initialize model_Iref,
    // in order to peak much performance,set it larger after iter 1.
    if(debug_flag) NODE0ONLY printMem(nr_pool);

    int current_Fsize2 = current_size*(current_size/2 + 1);
    int nr_images = nr_pool;
    exp_highres_Xi2_imgs		.init(nr_images);
    exp_min_diff2				.init(nr_images);
    exp_sum_weight				.init(nr_images);
    exp_significant_weight		.init(nr_images);
    exp_old_offsetx				.init(nr_images);
    exp_old_offsety				.init(nr_images);
    exp_wsum_norm_correction	.init(nr_images);
    exp_local_sqrtXi2			.init(nr_images);
    exp_local_Fctfs				.init(nr_images, current_Fsize2);
    exp_local_Minvsigma2s		.init(nr_images, current_Fsize2);
    exp_wsum_scale_correction_XA.init(nr_images, ori_size/2+1);
    exp_wsum_scale_correction_AA.init(nr_images, ori_size/2+1);
    exp_imgs					.init(nr_images, ori_size*ori_size);
    exp_power_imgs				.init(nr_images, ori_size/2+1);
    //
    exp_metadata				.init(nr_images);

	// initialize sampling
    samplingGrid.initialize(sampler3d, adaptive_oversampling);
    exp_nr_trans = samplingGrid.exp_nr_trans();
    exp_nr_psi = samplingGrid.exp_nr_psi();
    exp_nr_dir = samplingGrid.exp_nr_dir();
    int max_nr_over_rot   = sampler3d.oversamplingFactorOrientations(adaptive_oversampling);
    int max_nr_over_trans = sampler3d.oversamplingFactorTranslations(adaptive_oversampling);

    // intialzie some shift images and its ctf
    int exp_nr_all_trans    		= sampler3d.NrTrans(adaptive_oversampling);
    exp_Fimgs_shifted_real			.init(nr_images, exp_nr_all_trans, current_Fsize2);
    exp_Fimgs_shifted_imag			.init(nr_images, exp_nr_all_trans, current_Fsize2);
    exp_local_Fctfs					.init(nr_images, current_Fsize2);
    exp_local_Minvsigma2s			.init(nr_images, current_Fsize2);

    exp_Mcoarse_Rot_significant		.init(nr_classes,exp_nr_dir,exp_nr_psi,nr_images,exp_nr_trans,do_local_searching);
    
    // base on the density of this matrix
    // for high-density case(in coarse easrch step or previous iterations ) use c-array
    // for low-density case(in fine search step or later iterations ) use stack data-structure
    // for low-density case,also need to improve the nr_pool to increase the image-processing each time
#ifdef EXP_MWEIGHT_NEW
    exp_Mweight_coarse	.init(	nr_images, nr_classes, exp_nr_dir, exp_nr_psi,
								1, exp_nr_trans, 1);
    exp_Mweight_fine	.init(	nr_images, nr_classes, exp_nr_dir, exp_nr_psi,
                          		max_nr_over_rot, exp_nr_trans, max_nr_over_trans);
#else
    exp_Mweight_coarse		.init(	nr_images, nr_classes, exp_nr_dir, exp_nr_psi, exp_nr_trans);
    exp_Mweight_fine		.init(	nr_images, nr_classes, exp_nr_dir, exp_nr_psi,
                                  max_nr_over_rot*exp_nr_trans*max_nr_over_trans);
#endif
    // thread data
    // -------------   thread data     --------------- //
    thread_exp_max_weight_index		.resize(maxthreads);
    for (auto &exp_max_weight_index : thread_exp_max_weight_index)
        exp_max_weight_index.resize(nr_images);
    thread_Frefctf_real				.init(maxthreads, exp_nr_psi*max_nr_over_rot, current_Fsize2);
    thread_Frefctf_real				.fill_with_first_touch(0.);
    thread_Frefctf_imag				.init(maxthreads, exp_nr_psi*max_nr_over_rot, current_Fsize2);
    thread_Frefctf_imag				.fill_with_first_touch(0.);
    thread_Fimg_real				.init(maxthreads, max_nr_over_trans, ori_size*ori_size);
    thread_Fimg_real				.fill_with_first_touch(0.);
    thread_Fimg_imag				.init(maxthreads, max_nr_over_trans, ori_size*ori_size);
    thread_Fimg_imag				.fill_with_first_touch(0.);
    thread_Fweight					.init(maxthreads, exp_nr_psi, ori_size*ori_size);
    thread_Fweight					.fill_with_first_touch(0.);
    // ------- //
    thread_exp_sum_weight			.init(maxthreads, nr_images);
    thread_exp_min_diff2			.init(maxthreads, nr_images);
    thread_exp_max_weight			.init(maxthreads, nr_images);
    thread_wsum_norm_correction		.init(maxthreads, nr_images);
    thread_wsum_pdf_class			.init(maxthreads, nr_classes);
    thread_sumw_group				.init(maxthreads, nr_images);
    threadfake_do_scale_norm_class	.init(nr_classes, current_Fsize2);
    thread_wsum_sigma2_offset		.init(maxthreads, 1);
    thread_wsum_pdf_direction		.init(maxthreads, nr_classes, exp_nr_dir);
    thread_wsum_sigma2_noise		.init(maxthreads, nr_images, current_Fsize2);
    thread_wsum_scale_correction_XA	.init(maxthreads, nr_images, current_Fsize2);
    thread_wsum_scale_correction_AA	.init(maxthreads, nr_images, current_Fsize2);
    thread_exp_Mweight_sub			.init(maxthreads, const_max_iimage_tile*exp_nr_trans*std::max(const_max_ipsi_tile,max_nr_over_trans*max_nr_over_rot));

    //
    assert(do_shifts_onthefly==true);
    particleModel.setup(nr_images, current_size, coarse_size, exp_nr_trans, max_nr_over_trans);

    //
    fourierShellTrans.setCoarseImageBoundary(Mresol_coarse.rptrAll(), coarse_size*(coarse_size/2+1));
    fourierShellTrans.setNormCorrectionBoundary(Mresol_fine.rptrAll(),current_Fsize2);
    fourierShellTrans.setScaleCorrectionBoundary(mlModel.data_vs_prior_class);
    fourierShellTrans.setTransformFineAndCoarse(maxthreads);
}

void endExpData()
{
    //
    Npix_per_shell					.fini();
    Mresol_coarse					.fini();
    Mresol_fine						.fini();
    //
    exp_highres_Xi2_imgs			.fini();
    exp_min_diff2					.fini();
    exp_sum_weight					.fini();
    exp_significant_weight			.fini();
    exp_old_offsetx					.fini();
    exp_old_offsety					.fini();
    exp_wsum_norm_correction		.fini();
    exp_local_sqrtXi2				.fini();
    exp_local_Fctfs					.fini();
    exp_local_Minvsigma2s			.fini();
    exp_wsum_scale_correction_XA	.fini();
    exp_wsum_scale_correction_AA	.fini();
    exp_imgs						.fini();
    exp_power_imgs					.fini();
    //
    exp_metadata					.fini();
    //
    exp_Fimgs_shifted_real			.fini();
    exp_Fimgs_shifted_imag			.fini();
    //
    exp_Mcoarse_Rot_significant		.fini();
    exp_Mweight_coarse				.fini();
    exp_Mweight_fine				.fini();
    //
    samplingGrid					.finalize();
    // some thread variable
    thread_exp_max_weight_index		.resize(0);
    //
    thread_Frefctf_real				.fini();
    thread_Frefctf_imag				.fini();
    thread_Fimg_real				.fini();
    thread_Fimg_imag				.fini();
    thread_Fweight					.fini();
    //
    thread_exp_sum_weight			.fini();
    thread_exp_min_diff2			.fini();
    thread_exp_max_weight			.fini();
    thread_wsum_norm_correction		.fini();
    thread_wsum_pdf_class			.fini();
    thread_sumw_group				.fini();
    threadfake_do_scale_norm_class	.fini();
    thread_wsum_sigma2_offset		.fini();
    //
    thread_wsum_pdf_direction		.fini();
    thread_wsum_sigma2_noise		.fini();
    thread_wsum_scale_correction_XA	.fini();
    thread_wsum_scale_correction_AA	.fini();
    thread_exp_Mweight_sub			.fini();
    //
    particleModel					.destroy();
    //
    fourierShellTrans				.clear();
}

void expectation()
{
	maybeGoSerial("expectation");

    NODE0ONLY std::cerr << " Expectation iteration " << iter<< " of " << nr_iter<<std::endl;

    int my_first_image,my_last_image,nr_images_done = 0;

    NODE0ONLY showProgressBar(0, nr_local_images);

    while (nr_images_done < nr_local_images)
    {
        my_first_image = first_local_image + nr_images_done;

        if ((!do_firstiter_cc && iter == 1) || (do_firstiter_cc && iter == 2)){
            // equally divided like nr_pool=1 case
            int iclass = divide_equally_which_group(nr_global_images, nr_classes, my_first_image);
            int first,last;
            divide_equally(nr_global_images, nr_classes, iclass, first, last);
            int suitable_pool = std::min(nr_pool, last-my_first_image+1);
            my_last_image = std::min(first_local_image+nr_local_images - 1, my_first_image + suitable_pool - 1);
        }
        else
            my_last_image = std::min(first_local_image+nr_local_images - 1, my_first_image + nr_pool - 1);

        exp_first_image = my_first_image;
        exp_last_image = my_last_image;
        exp_nr_images = my_last_image - my_first_image + 1;

#ifdef DATA_STREAM
        global_data_stream.foutInt(exp_first_image, "expectation()_my_first_ori_particle", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_last_image, "expectation()_my_last_ori_particle", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_nr_images, "expectation()_nr_pool", __FILE__, __LINE__);
#endif

        // first time divided the images by nr_classes peiece,make the initialized reference different
        exp_iclass_min = 0;
        exp_iclass_max = nr_classes - 1;
        // low-pass filter again and generate the seeds
        if (true/*do_generate_seeds*/)
        {
            if (do_firstiter_cc && iter == 1)
            {
                // In first (CC) iter, use a single reference (and CC)
                exp_iclass_min = exp_iclass_max = 0;
            }
            else if ( (do_firstiter_cc && iter == 2) || (!do_firstiter_cc && iter == 1))
            {
                // In second CC iter, or first iter without CC: generate the seeds
                // Now select a single random class
                // exp_part_id is already in randomized order (controlled by -seed)
                // WARNING: USING SAME iclass_min AND iclass_max FOR SomeParticles!!
                exp_iclass_min = exp_iclass_max = divide_equally_which_group(nr_global_images, nr_classes, exp_first_image);
            }
        }

#ifdef DATA_STREAM
        global_data_stream.foutInt(do_firstiter_cc, "expectationOneParticle()_do_firstiter_cc", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_iclass_min, "expectationOneParticle()_exp_iclass_min", __FILE__, __LINE__);
        global_data_stream.foutInt(exp_iclass_max, "expectationOneParticle()_exp_iclass_max", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif

        // prepare exp_image_data and exp_metadata
        // each node keep nr_local_images's images data and all(nr_global_images) metadata
        for (int iimage = 0; iimage < exp_nr_images; iimage++) {
            exp_metadata[iimage] = metadata[iimage+my_first_image];
            ::copy(exp_imgs[iimage].wptrAll(), ori_size*ori_size, images[my_first_image+iimage-first_local_image].rptrAll(), ori_size*ori_size);
            // set the scale correction to prevent NUMA-ISSUA
            if (do_scale_correction){exp_metadata[iimage].SCALE = mlModel.scale_correction.wptrAll()[exp_metadata[iimage].GROUP_NO-1];}
            else{exp_metadata[iimage].SCALE = 1.0;}
        }

        // get Image(with mask and nomask) and CTF in Fourier Space
        // it will downsize the image from ori_size to current_size
        particleModel.prepare(exp_imgs, exp_metadata, exp_first_image, exp_nr_images);
        particleModel.setFourierTransforms(	exp_power_imgs, exp_highres_Xi2_imgs,
                                            exp_old_offsetx, exp_old_offsety,
                                            mlModel,do_ctf_correction);

        // If do local searching
        if (do_local_searching)
        {
            assert(mlModel.orientational_prior_mode == PRIOR_ROTTILT_PSI);
            assert(sampler3d.orientational_prior_mode == PRIOR_ROTTILT_PSI);
            assert(sigma2_angle!=0);
            samplingGrid.selectOrientationsForLocalSearch(sampler3d, mlModel.sigma2_rot, mlModel.sigma2_tilt,
                                                          mlModel.sigma2_psi, &exp_metadata[0], exp_nr_images);
        }
        // perform the actual expectation step on several particles
        expectationSomeParticles();

        // Also monitor the changes in the optimal orientations and classes
        hiddenVarMonitor.monitorHiddenVariableChanges(sampler3d, metadata, my_first_image, &exp_metadata[0], exp_nr_images);

        // update the metadata
        for (int iimage = 0; iimage < exp_nr_images; iimage++) metadata[iimage+my_first_image] = exp_metadata[iimage];

        //
        nr_images_done += my_last_image - my_first_image + 1;

        NODE0ONLY showProgressBar(nr_images_done, nr_local_images);
    }
    NODE0ONLY std::cerr<<std::endl;
}



void expectationSomeParticles()
{
	maybeGoSerial("expectationSomeParticles");

	// Only perform a second pass when using adaptive oversampling
	const int  nr_sampling_passes   = (adaptive_oversampling > 0) ? 2 : 1;
    const bool do_cross_correlation = (iter == 1 && do_firstiter_cc) || do_always_cc;
	// Pass twice through the sampling of the entire space of rot, tilt and psi
	// The first pass uses a coarser angular sampling and possibly smaller FFTs than the second pass.
	// Only those sampling points that contribute to the highest x% of the weights in the first pass are oversampled in the second pass
	// Only those sampling points will contribute to the weighted sums in the third loop below
	for (exp_ipass = 0; exp_ipass < nr_sampling_passes; exp_ipass++)
	{
        const bool do_coarse_search = (exp_ipass == 0);
		
        // Use smaller images in the first pass, larger ones in the second pass
        exp_current_size = (do_coarse_search && adaptive_oversampling > 0) ? coarse_size : current_size;

		// Use coarse sampling in the first pass, oversampled one the second pass
        // and initialize sampling data
		exp_current_oversampling = do_coarse_search ? 0 : adaptive_oversampling;
		if (0) std::cerr << __FILE__ << ":" << __LINE__ << " calling samplingGrid.computeGrid(sampler3d, exp_current_oversampling:" << exp_current_oversampling << ");" << std::endl;
        samplingGrid.computeGrid3D(sampler3d, exp_current_oversampling);
        samplingGrid.testGetShift();
        exp_nr_over_rot = samplingGrid.exp_nr_over_rot();
        exp_nr_over_trans = samplingGrid.exp_nr_over_trans();
		particleModel.shiftAssistorsSetCurrDims(exp_nr_trans, exp_nr_over_trans);
        //
        particleModel.preShiftedImagesCtfsAndInvSigma2s(exp_Fimgs_shifted_real, exp_Fimgs_shifted_imag,
                                                        exp_local_Minvsigma2s, exp_local_sqrtXi2, exp_local_Fctfs,
                                                        exp_current_size, do_cross_correlation, do_coarse_search,
                                                        samplingGrid, mlModel, Mresol_coarse, Mresol_fine);

        // get all rotated reference  and the significant rotation
        if(!do_coarse_search) getReferenceAllOrientations();

        // get all reference and images 's squared differences
        getAllSquaredDifferences(do_coarse_search);

		// convert all squared differences to weight,and find significant(maximum) weight for each image
        if (do_coarse_search){
            convertSquaredDifferencesToWeights(exp_Mweight_coarse);
            findAllSignificantPoints(exp_Mweight_coarse);
        }
        else{
            convertSquaredDifferencesToWeights(exp_Mweight_fine);
            findAllSignificantPoints(exp_Mweight_fine);
        }

        NODE0ONLY
        {
            std::string iterStr = num2str(iter);
            if (do_coarse_search) {
                std::string fn_exp_mweight = write_path+write_fn+"_iter"+iterStr+"_exp_Mweight_coarse";
                std::string note = "coarse_search_"+std::to_string((long long)exp_first_image);
                exp_Mweight_coarse.analysis(fn_exp_mweight,note);
            }
            else {
                std::string fn_exp_mweight = write_path+write_fn+"_iter"+iterStr+"_exp_Mweight_fine";
                std::string note = "fine_search_"+std::to_string((long long)exp_first_image);
                exp_Mweight_fine.analysis(fn_exp_mweight,note);
            }
        }
	}// end loop over 2 exp_ipass iterations

    // update some parameter for maximization
	maybeGoSerial("updateOtherParams");
    updateOtherParams();

    //
	maybeGoSerial("backProjection");
    backProjection();

    //
	maybeGoSerial("storeWeightedSums");
    storeWeightedSums();
}


void getReferenceAllOrientations()
{
    // !!!! exp_Rot_significant will be not call on coarse step !!!!!
    assert(exp_ipass == 1);
    // In the first pass, always proceed
    // In the second pass, check whether one of the translations for this orientation of any of
    // the particles had a significant weight in the first pass
    // if so, proceed with projecting the reference in that direction
    exp_Mcoarse_Rot_significant.resetRotSignificant();
    for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++){
        if (mlModel.pdf_class.rptrAll()[iclass] <= 0.)
            continue;
        for (int idir = 0; idir < exp_nr_dir; idir++){
            for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++){
                exp_Mcoarse_Rot_significant.allocAndAssignRotSignificant(iclass, idir, ipsi, true);
            }// end loop of ipsi
        }// end loop of idir
    }// end loop of iclass
}

// TODO,better one
inline void decodeRotOverTrans(int ihidden,int& itrans,int& iover_trans,int& iover_rot)
{
    // decode the real position,psi->trans->over_rot->over_trans
    int denominator = exp_nr_over_trans*exp_nr_over_rot;
    itrans = ihidden / denominator;

    ihidden = ihidden % denominator;
    denominator = exp_nr_over_rot;
    iover_trans = ihidden / denominator;
    iover_rot = ihidden % denominator;
}

static void check_scale(double& myscale,int& igroup){
    static bool have_warned_small_scale = false;
    if (myscale > 10000.)
    {
        std::cerr << " rlnMicrographScaleCorrection= " << myscale << " group= " << igroup + 1 << std::endl;
        ERROR_REPORT("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
    }
    else if (myscale < 0.001)
    {
        if (!have_warned_small_scale)
        {
            std::cout << " WARNING: ignoring group " << igroup + 1 << " with very small or negative scale (" << myscale <<
            "); Use larger groups for more stable scale estimates." << std::endl;
            have_warned_small_scale = true;
        }
        myscale = 0.001;
    }
};

namespace CoarseAndFineSearch
{
	int exp_current_Fsize2;
	bool do_cross_correlation;
	bool do_coarse_search;
	// --------------- get projected class
	inline void getProjectedClass(int iclass, FDOUBLE over_rot, FDOUBLE over_tilt, FDOUBLE over_psi,
                                  FDOUBLE* Fref_real, FDOUBLE* Fref_imag, int n_start, int n_end)
	{
		FDOUBLE A[3][3];
        Euler_angles2matrix(over_rot, over_tilt, over_psi, A);
		mapModel.get2DFourierTransformOneTile(iclass, Fref_real, Fref_imag, n_start, n_end, exp_current_size, A);
	}
	inline void getProjectedClassByShellCoarse(int iclass, FDOUBLE over_rot, FDOUBLE over_tilt, FDOUBLE over_psi,
                                               FDOUBLE* Fref_real, FDOUBLE* Fref_imag, int shell_n_start, int shell_n_end)
	{
		FDOUBLE A[3][3];
		Euler_angles2matrix(over_rot, over_tilt, over_psi, A);
		mapModel.get2DFourierTransformOneTile(iclass, Fref_real, Fref_imag, shell_n_start, shell_n_end,
			exp_current_size, A, fourierShellTrans.rptr_nIndexCoarse());
	}
	inline void getProjectedClassByShellFine(int iclass, FDOUBLE over_rot, FDOUBLE over_tilt, FDOUBLE over_psi,
                                             FDOUBLE* Fref_real, FDOUBLE* Fref_imag,int shell_n_start, int shell_n_end)
	{
		FDOUBLE A[3][3];
		Euler_angles2matrix(over_rot, over_tilt, over_psi, A);
		mapModel.get2DFourierTransformOneTile(iclass, Fref_real, Fref_imag, shell_n_start, shell_n_end,
			exp_current_size, A, fourierShellTrans.rptr_nIndexFine());
	}
	//--------------------

	//
	void global_data_stream_before()
	{
#ifdef DATA_STREAM
		int tid = 0;
		for (int iimage = 0; iimage < exp_nr_images; iimage++) {
			auto Minvsigma2_iimage = exp_local_Minvsigma2s[iimage].wptr(exp_current_Fsize2);
			for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++) {
				for (int idir = 0; idir < exp_nr_dir; idir++) {
					if (mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
					for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
						// exp_iclass loop does not always go from 0 to nr_classes!
						// int iClassOrient = iclass*exp_nr_dir + idir;
						// iClassOrient = iClassOrient*exp_nr_psi + ipsi;
						if (exp_ipass==1 && !exp_Mcoarse_Rot_significant.isRotSignificant(iclass, idir, ipsi)) continue;
                        // If do local searching
                        double pdf_orientation;
                        if (do_local_searching) {
                            if (!samplingGrid.isSignificantOrientation(idir, ipsi)) continue;
                            pdf_orientation = (*samplingGrid.pdfOrientationForNonzeroOrient(idir, ipsi))[iimage];
                            if (pdf_orientation <= 0) continue;
                        }
                        else{
                            pdf_orientation = mlModel.pdf_direction[iclass].rptrAll()[idir];
                            if (pdf_orientation <= 0.) continue;
                        }

                        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
                        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                                    exp_over_rot, exp_over_tilt, exp_over_psi);
						for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++)
						{
                            // std::cout<<exp_over_rot[iover_rot]<<" "<<exp_over_tilt[iover_rot]<<" "<<exp_over_psi[iover_rot]<<std::endl;
							auto Frefctf_real = thread_Frefctf_real[tid][ipsi].wptr(exp_current_Fsize2);
							auto Frefctf_imag = thread_Frefctf_imag[tid][ipsi].wptr(exp_current_Fsize2);
							// get projected image
							getProjectedClass(iclass, exp_over_rot[iover_rot], exp_over_tilt[iover_rot], exp_over_psi[iover_rot],
                                              Frefctf_real, Frefctf_imag, 0, exp_current_Fsize2);
							// write out data stream
							global_data_stream.foutInt(iclass, "getAllSquaredDifferences()_iclass", __FILE__, __LINE__);
							global_data_stream.foutInt(iimage, "getAllSquaredDifferences()_iimage", __FILE__, __LINE__);
							global_data_stream.foutInt(idir, "getAllSquaredDifferences()_idir", __FILE__, __LINE__);
							global_data_stream.foutInt(ipsi, "getAllSquaredDifferences()_ipsi", __FILE__, __LINE__);
							global_data_stream.foutInt(iover_rot, "getAllSquaredDifferences()_iover_rot", __FILE__, __LINE__);
							global_data_stream.foutDouble(mlModel.pdf_class.wptrAll()[iclass], "getAllSquaredDifferences()_pdf_class[iclass]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(pdf_orientation, "getAllSquaredDifferences()_pdf_direction[iclass][idir]", __FILE__, __LINE__);
							global_data_stream.foutDouble(exp_over_rot[iover_rot], "getAllSquaredDifferences()_exp_over_rot[iorient][iover_rot]", __FILE__, __LINE__);
							global_data_stream.foutDouble(exp_over_tilt[iover_rot], "getAllSquaredDifferences()_exp_over_tilt[iorient][iover_rot]", __FILE__, __LINE__);
							global_data_stream.foutDouble(exp_over_psi[iover_rot], "getAllSquaredDifferences()_exp_over_psi[iorient][iover_rot]", __FILE__, __LINE__);
							global_data_stream.foutInt(do_firstiter_cc, "getAllSquaredDifferences()_do_firstiter_cc", __FILE__, __LINE__);
							global_data_stream.foutInt(do_always_cc, "getAllSquaredDifferences()_do_always_cc", __FILE__, __LINE__);
							global_data_stream.foutInt(do_scale_correction, "getAllSquaredDifferences()_do_scale_correction", __FILE__, __LINE__);
							global_data_stream.foutInt(do_ctf_correction, "getAllSquaredDifferences()_do_ctf_correction", __FILE__, __LINE__);
							global_data_stream.foutInt(refs_are_ctf_corrected, "getAllSquaredDifferences()_refs_are_ctf_corrected", __FILE__, __LINE__);
							global_data_stream.foutInt(exp_metadata[iimage].GROUP_NO - 1, "getAllSquaredDifferences()_GROUP_NO", __FILE__, __LINE__);
							global_data_stream.foutDouble(exp_metadata[iimage].SCALE, "getAllSquaredDifferences()_scale_correction", __FILE__, __LINE__);
							global_data_stream.foutDouble(Frefctf_real, exp_current_Fsize2, "getAllSquaredDifferences()_Fref_real", __FILE__, __LINE__);
							global_data_stream.foutDouble(Frefctf_imag, exp_current_Fsize2, "getAllSquaredDifferences_Fref_imag", __FILE__, __LINE__);
							global_data_stream.foutDouble(Minvsigma2_iimage, exp_current_Fsize2, "getAllSquaredDifferences_Fref_imag()_Minvsigma2", __FILE__, __LINE__);
							global_data_stream.check(); global_data_stream.flush();
						}
					}
				}
			}
		}
#endif
	}
	//
	void global_data_stream_after()
	{
#ifdef DATA_STREAM
		for (int iimage = 0; iimage < exp_nr_images; iimage++)
		{
			if (do_cross_correlation)
				global_data_stream.foutDouble(exp_local_sqrtXi2.wptrAll()[iimage], "getAllSquaredDifferences()_exp_local_sqrtXi2", __FILE__, __LINE__);
			else
				global_data_stream.foutDouble(exp_highres_Xi2_imgs.wptrAll()[iimage], "getAllSquaredDifferences()_exp_highres_Xi2_imgs", __FILE__, __LINE__);
		}
		// NOTE,exp_Mweight is checked in convertSquaredDifferencesToWeights()....
		global_data_stream.foutDouble(exp_min_diff2.wptrAll()[0], "getAllSquaredDifferences()_exp_min_diff2[0]", __FILE__, __LINE__);
		global_data_stream.check(); global_data_stream.flush();
#endif
	}

	// do cross-correlation
	IntPerformanceCounter transformAndGetDiffCC_performanceCounter("transformAndGetDiffCC_performanceCounter");
		// incremented by the callers of transformAndGetDiffCC because they can do it outside a loop

	void transformAndGetDiffCC(int iimage, double& diff2, int n_start, int n_end,
		const FDOUBLE* Fimag_real, const FDOUBLE* Fimag_imag,
		const FDOUBLE* aTable, const FDOUBLE* bTable,
		const FDOUBLE* exp_local_Fctfs, double myscale,
		const FDOUBLE* Fref_real, const FDOUBLE* Fref_imag)
	{
		L2CacheModel::seqAcc("transformAndGetDiffCC aTable",     -1, &aTable    [n_start], n_end - n_start, sizeof(aTable    [n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC bTable",     -1, &bTable    [n_start], n_end - n_start, sizeof(bTable    [n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fimag_real", -1, &Fimag_real[n_start], n_end - n_start, sizeof(Fimag_real[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fimag_imag", -1, &Fimag_imag[n_start], n_end - n_start, sizeof(Fimag_imag[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_real",  -1, &Fref_real [n_start], n_end - n_start, sizeof(Fref_real [n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_imag",  -1, &Fref_imag [n_start], n_end - n_start, sizeof(Fref_imag [n_start]));

		assert(n_start == 0); assert(n_end == exp_current_Fsize2);
		diff2 = 0.;
		// Do not calculate squared-differences, but signal product
		// Negative values because smaller is worse in this case
		double suma2 = 0.;

		if (!(true
			&& isVectorAligned(aTable          + n_start)
			&& isVectorAligned(bTable          + n_start)
			&& isVectorAligned(Fimag_real      + n_start)
			&& isVectorAligned(Fimag_imag      + n_start)
			&& isVectorAligned(Fref_real       + n_start)
			&& isVectorAligned(Fref_imag       + n_start)
			&& isVectorAligned(exp_local_Fctfs + n_start)
			&& isVectorAligned(Fimag_imag      + n_start))) {
			EXIT_ABNORMALLY;
		}

#pragma vector aligned
#pragma ivdep
		for (int n = n_start; n < n_end; n++)
		{
			const double a = aTable[n];
			const double b = bTable[n];
			const double c = Fimag_real[n];
			const double d = Fimag_imag[n];
			const double ac = a * c;
			const double bd = b * d;
			const double ab_cd = (a + b) * (c + d);
			double fout_shift_real = CHECK_NOT_IND(ac - bd);
			double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd);
            double frefctf_real = Fref_real[n] * exp_local_Fctfs[n];// * myscale;
            double frefctf_imag = Fref_imag[n] * exp_local_Fctfs[n];// * myscale;
			diff2 -= frefctf_real * fout_shift_real;
			diff2 -= frefctf_imag * fout_shift_imag;
			suma2 += frefctf_real * frefctf_real + frefctf_imag * frefctf_imag;
		}
		// Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
		diff2 /= sqrt(suma2) * exp_local_sqrtXi2.rptrAll()[iimage];
	}

	typedef double Diff2_4[4];
	void transformAndGetDiffCC(int iimage, Diff2_4& diff2, int n_start, int n_end,
		const FDOUBLE* Fimag_real, const FDOUBLE* Fimag_imag,
		const FDOUBLE* aTable,     const FDOUBLE* bTable,
		const FDOUBLE* exp_local_Fctfs, double myscale,
		const FDOUBLE* Fref_real_0, const FDOUBLE* Fref_imag_0,
		const FDOUBLE* Fref_real_1, const FDOUBLE* Fref_imag_1,
		const FDOUBLE* Fref_real_2, const FDOUBLE* Fref_imag_2,
		const FDOUBLE* Fref_real_3, const FDOUBLE* Fref_imag_3)
	{
		L2CacheModel::seqAcc("transformAndGetDiffCC aTable",      -1, &aTable     [n_start], n_end - n_start, sizeof(aTable     [n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC bTable",      -1, &bTable     [n_start], n_end - n_start, sizeof(bTable     [n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fimag_real",  -1, &Fimag_real [n_start], n_end - n_start, sizeof(Fimag_real [n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fimag_imag",  -1, &Fimag_imag [n_start], n_end - n_start, sizeof(Fimag_imag [n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_real_0", -1, &Fref_real_0[n_start], n_end - n_start, sizeof(Fref_real_0[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_imag_0", -1, &Fref_imag_0[n_start], n_end - n_start, sizeof(Fref_imag_0[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_real_1", -1, &Fref_real_1[n_start], n_end - n_start, sizeof(Fref_real_1[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_imag_1", -1, &Fref_imag_1[n_start], n_end - n_start, sizeof(Fref_imag_1[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_real_2", -1, &Fref_real_2[n_start], n_end - n_start, sizeof(Fref_real_2[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_imag_2", -1, &Fref_imag_2[n_start], n_end - n_start, sizeof(Fref_imag_2[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_real_3", -1, &Fref_real_3[n_start], n_end - n_start, sizeof(Fref_real_3[n_start]));
		L2CacheModel::seqAcc("transformAndGetDiffCC Fref_imag_3", -1, &Fref_imag_3[n_start], n_end - n_start, sizeof(Fref_imag_3[n_start]));

		assert(n_start == 0); assert(n_end == exp_current_Fsize2);
		double diff2_0(0.0), diff2_1(0.0), diff2_2(0.0), diff2_3(0.0);
		// Do not calculate squared-differences, but signal product
		// Negative values because smaller is worse in this case
		double suma2_0(0.0), suma2_1(0.0), suma2_2(0.0), suma2_3(0.0);

		assert(true
				&& isVectorAligned(aTable			+ n_start)
				&& isVectorAligned(bTable			+ n_start)
				&& isVectorAligned(Fimag_real		+ n_start)
				&& isVectorAligned(Fimag_imag		+ n_start)
				&& isVectorAligned(Fref_real_0		+ n_start)
				&& isVectorAligned(Fref_real_1		+ n_start)
				&& isVectorAligned(Fref_real_2		+ n_start)
				&& isVectorAligned(Fref_real_3		+ n_start)
				&& isVectorAligned(Fref_imag_0		+ n_start)
				&& isVectorAligned(Fref_imag_1		+ n_start)
				&& isVectorAligned(Fref_imag_2		+ n_start)
				&& isVectorAligned(Fref_imag_3		+ n_start)
				&& isVectorAligned(exp_local_Fctfs	+ n_start)
				&& isVectorAligned(Fimag_imag		+ n_start));

#pragma vector aligned
#pragma ivdep
		for (int n = n_start; n < n_end; n++)
		{
			const double a = aTable[n];
			const double b = bTable[n];
			const double c = Fimag_real[n];
			const double d = Fimag_imag[n];
			const double ac = a * c;
			const double bd = b * d;
			const double ab_cd = (a + b) * (c + d);
			double fout_shift_real = CHECK_NOT_IND(ac - bd);
			double fout_shift_imag = CHECK_NOT_IND(ab_cd - ac - bd);
			double frefctf_real;
			double frefctf_imag;
#define ONE_FREF(N)																	\
			frefctf_real = Fref_real_##N[n] * exp_local_Fctfs[n]; /* * myscale; */	\
			frefctf_imag = Fref_imag_##N[n] * exp_local_Fctfs[n]; /* * myscale; */	\
			diff2_##N -= frefctf_real * fout_shift_real;							\
			diff2_##N -= frefctf_imag * fout_shift_imag;							\
			suma2_##N += frefctf_real * frefctf_real + frefctf_imag * frefctf_imag;	\
			// end of macro
			ONE_FREF(0) ONE_FREF(1) ONE_FREF(2) ONE_FREF(3)
#undef ONE_FREF
		}
		// Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
		auto f = exp_local_sqrtXi2.rptrAll()[iimage];
		diff2[0] = diff2_0 / (sqrt(suma2_0) * f);
		diff2[1] = diff2_1 / (sqrt(suma2_1) * f);
		diff2[2] = diff2_2 / (sqrt(suma2_2) * f);
		diff2[3] = diff2_3 / (sqrt(suma2_3) * f);
	}

	void getAllSquaredDifferencesCoarse_loop2Body(TileCfg& tileCfg,
		const int N, 
		const int n_start, 
		const int N_tile,
		const int iclass, 
		const int idir,
		const int ipsi_tile_start,
		const int ipsi_tile,
		const int iimage_tile,
		const int iimage_sub_tile,
		const int itrans_sub_tile) {

		const int  tid = omp_get_thread_num();
		const int  iover_rot = 0;
		const int  iover_trans = 0;
		const int  n_end                       = std::min(n_start + N_tile, N);
		const int  ipsi_tile_end               = std::min(ipsi_tile_start + ipsi_tile, exp_nr_psi);
		const auto thread_exp_min_diff2_tid    = thread_exp_min_diff2[tid].wptr(exp_nr_images);

		// Get these (an expensive process) and save them for use across many images and slides.
		// Since there are many images and many slides and they are all being done, this cost is spread across much work
		// and if even if these images don't stay in the cache, it takes a lot of more memory traffic to compute them later
		// and that traffic won't be cached either.
		//
		GetAllSquaredDifferencesFine_Kernel kernel(exp_current_size);
        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
		for (int ipsi = ipsi_tile_start; ipsi < ipsi_tile_end; ipsi++)
		{
            samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                        exp_over_rot, exp_over_tilt, exp_over_psi);
			auto Fref_real = thread_Frefctf_real[tid][ipsi - ipsi_tile_start].wptrAll();
			auto Fref_imag = thread_Frefctf_imag[tid][ipsi - ipsi_tile_start].wptrAll();
			// get projected image
			if (do_cross_correlation) {
				getProjectedClass(iclass, exp_over_rot[0], exp_over_tilt[0], exp_over_psi[0], Fref_real, Fref_imag, n_start, n_end);
			} else {
				getProjectedClassByShellCoarse(iclass, exp_over_rot[0], exp_over_tilt[0], exp_over_psi[0], Fref_real, Fref_imag, n_start, n_end);
				kernel.appendFref(Fref_real, Fref_imag);
			}
		}
		
        // If do local search
        const std::vector<FDOUBLE>* pdf_orientation_iimage;
        if (do_local_searching) {
            assert(1==ipsi_tile);
            pdf_orientation_iimage = samplingGrid.pdfOrientationForNonzeroOrient(idir, ipsi_tile_start);
        }
		// Do the images and the slides in the right order to minimize the memory traffic.
		// 
		// Since the size of the slide-specific info greatly exceeds the image-specific info
		// we load as many shifters as will fit, and then bring in an image, 
		// The inner loop
		//
		for (int iimage_tile_start = 0; iimage_tile_start < exp_nr_images; iimage_tile_start += iimage_tile) {
			int  const iimage_tile_end = std::min(iimage_tile_start + iimage_tile, exp_nr_images);

			// The following uses 
			//		auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid + (iimage - iimage_tile_start) * exp_nr_trans * ipsi_tile;
			//		     thread_exp_Mweight_sub_sub_tid                                                                    [itrans * ipsi_tile + (ipsi - ipsi_tile_start)]
			//	as its indexing scheme within thread_exp_Mweight_sub[tid]
			//
			auto &     thread_exp_Mweight_sub_tid_vector = thread_exp_Mweight_sub[tid];
			const int  thread_exp_Mweight_sub_size = (iimage_tile_end - iimage_tile_start)*exp_nr_trans*ipsi_tile;
			const auto thread_exp_Mweight_sub_tid  = thread_exp_Mweight_sub_tid_vector.wptr(thread_exp_Mweight_sub_size);

			for (int iimage_sub_tile_start = iimage_tile_start; iimage_sub_tile_start < iimage_tile_end; iimage_sub_tile_start += iimage_sub_tile)
			{
				for (int itrans_sub_tile_start = 0; itrans_sub_tile_start < exp_nr_trans; itrans_sub_tile_start += itrans_sub_tile)
				{// 100(offset_range=10,step=2) or 400(offset_range=10,step=1)
					int iimage_sub_tile_end = std::min(iimage_sub_tile_start + iimage_sub_tile, iimage_tile_end);
					int itrans_sub_tile_end = std::min(itrans_sub_tile_start + itrans_sub_tile, exp_nr_trans);
					for (int iimage = iimage_sub_tile_start; iimage < iimage_sub_tile_end; iimage++)
                    {
                        // If do local search
                        if (do_local_searching) {
                            assert(1==ipsi_tile);
                            if((*pdf_orientation_iimage)[iimage] <= 0) continue;
                        }
						auto Minvsigma2_iimage = exp_local_Minvsigma2s[iimage].wptrAll();
						auto local_ctf = refs_are_ctf_corrected ? exp_local_Fctfs[iimage].wptrAll() : particleModel.FctfsOne[iimage].wptrAll();
						//int igroup = exp_metadata[iimage].GROUP_NO-1;
						// double myscale = mlModel.scale_correction.rptrAll()[igroup];
						//check_scale(myscale,igroup);
						const double myscale = exp_metadata[iimage].SCALE;// TODO : NUMA ISSUA ???
						const auto thread_exp_Mweight_sub_sub_tid_offset = (iimage - iimage_tile_start) * exp_nr_trans * ipsi_tile;
						const auto thread_exp_Mweight_sub_sub_tid_size   = thread_exp_Mweight_sub_size - thread_exp_Mweight_sub_sub_tid_offset;
						if (thread_exp_Mweight_sub_sub_tid_size < (itrans_sub_tile_end-1)*ipsi_tile + (ipsi_tile_end-1 - ipsi_tile_start)) {
							EXIT_ABNORMALLY;
						}
						const auto thread_exp_Mweight_sub_sub_tid        = thread_exp_Mweight_sub_tid  + thread_exp_Mweight_sub_sub_tid_offset;

						if (do_cross_correlation)
						{
							auto &     thread_Frefctf_real_tid  = thread_Frefctf_real[tid];
							auto &     thread_Frefctf_imag_tid  = thread_Frefctf_imag[tid];
							auto const Fimages_mask_coarse_real = particleModel.Fimages_mask_coarse_real[iimage].rptrAll();
							auto const Fimages_mask_coarse_imag = particleModel.Fimages_mask_coarse_imag[iimage].rptrAll();

							for (int itrans = itrans_sub_tile_start; itrans < itrans_sub_tile_end; itrans++)
							{
								auto shifter = particleModel.shiftImageAssistor.acquireShifter(itrans, 0, itrans*exp_nr_over_trans + 0);
								auto aTable = shifter->aTable_rptr();	// There are too many of these
								auto bTable = shifter->bTable_rptr();

								transformAndGetDiffCC_performanceCounter.count.v += ipsi_tile_end - ipsi_tile_start;

								int ipsi = ipsi_tile_start;

								double diff2v[4];
								for (; ipsi+3 < ipsi_tile_end; ipsi+=4)
								{
									auto Fref_real_0 = thread_Frefctf_real_tid[ipsi+0 - ipsi_tile_start].rptrAll();
									auto Fref_imag_0 = thread_Frefctf_imag_tid[ipsi+0 - ipsi_tile_start].rptrAll();
									auto Fref_real_1 = thread_Frefctf_real_tid[ipsi+1 - ipsi_tile_start].rptrAll();
									auto Fref_imag_1 = thread_Frefctf_imag_tid[ipsi+1 - ipsi_tile_start].rptrAll();
									auto Fref_real_2 = thread_Frefctf_real_tid[ipsi+2 - ipsi_tile_start].rptrAll();
									auto Fref_imag_2 = thread_Frefctf_imag_tid[ipsi+2 - ipsi_tile_start].rptrAll();
									auto Fref_real_3 = thread_Frefctf_real_tid[ipsi+3 - ipsi_tile_start].rptrAll();
									auto Fref_imag_3 = thread_Frefctf_imag_tid[ipsi+3 - ipsi_tile_start].rptrAll();
									// compute difference
									transformAndGetDiffCC(iimage, diff2v, n_start, n_end,
										Fimages_mask_coarse_real,
										Fimages_mask_coarse_imag,
										aTable,
										bTable,
										local_ctf, myscale, 
										Fref_real_0, Fref_imag_0,
										Fref_real_1, Fref_imag_1,
										Fref_real_2, Fref_imag_2,
										Fref_real_3, Fref_imag_3);
									
									// NOTE : this inorder array access may cause big slow???
									for (int i = 0; i < 4; i++) {
										thread_exp_Mweight_sub_sub_tid[itrans*ipsi_tile + (ipsi+i - ipsi_tile_start)] = diff2v[i];
									}
								} // end loop ipsi

								for (; ipsi < ipsi_tile_end; ipsi++)
								{
									auto Fref_real = thread_Frefctf_real_tid[ipsi - ipsi_tile_start].rptrAll();
									auto Fref_imag = thread_Frefctf_imag_tid[ipsi - ipsi_tile_start].rptrAll();
									// compute difference
									double diff2 = 0;
									transformAndGetDiffCC(iimage, diff2, n_start, n_end,
										Fimages_mask_coarse_real,
										Fimages_mask_coarse_imag,
										aTable,
										bTable,
										local_ctf, myscale, Fref_real, Fref_imag);

									// NOTE : this inorder array access may cause big slow???
									thread_exp_Mweight_sub_sub_tid[itrans*ipsi_tile + (ipsi - ipsi_tile_start)] = diff2;
								} // end loop ipsi

								particleModel.shiftImageAssistor.releaseShifter(shifter);
							}
						} else {
							for (int itrans = itrans_sub_tile_start; itrans < itrans_sub_tile_end; itrans++)
							{
								{
									kernel.acquireImage(particleModel.Fimages_mask_coarse_real[iimage].rptrAll(),
										particleModel.Fimages_mask_coarse_imag[iimage].rptrAll(),
										local_ctf, Minvsigma2_iimage, myscale, n_start, n_end);
									auto shifter = particleModel.shiftImageAssistor.acquireShifter(itrans, 0, itrans*exp_nr_over_trans + 0);

									kernel.acquireTable(
										shifter->aTable_rptr(),
										shifter->bTable_rptr(),
										shifter->tableSize());

									const auto offset_within_sub_sub = itrans*ipsi_tile;
									const auto size_within_sub_sub = ipsi_tile_end - ipsi_tile_start;
									if (offset_within_sub_sub + size_within_sub_sub > thread_exp_Mweight_sub_sub_tid_size) {
										assert(false);
										EXIT_ABNORMALLY;
									}

									kernel.appendDiff2(size_within_sub_sub);

									kernel.compute();

									// Store all diff2 in exp_Mweight, TODO....
									kernel.release(thread_exp_Mweight_sub_sub_tid + offset_within_sub_sub, size_within_sub_sub);
									// WHY IS THIS PASSING A BAD PTR IN?	ptr + 36*24 + 16

									particleModel.shiftImageAssistor.releaseShifter(shifter);
								}
							} // end loop itrans
						}
					}// end loop iimage
				}// end loop itrans_sub_tile
			}// end loop iimage_sub_tile

			for (int iimage = iimage_tile_start; iimage < iimage_tile_end; iimage++) {
                // If do local search
                if (do_local_searching) {
                    assert(1==ipsi_tile);
                    if((*pdf_orientation_iimage)[iimage] <=0 ) continue;
                }
				auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid + (iimage - iimage_tile_start) * exp_nr_trans * ipsi_tile;
				for (int ipsi = ipsi_tile_start; ipsi < ipsi_tile_end; ipsi++) // 24(healpix_order=1)~48(healpix_order=2)
				{
#ifdef EXP_MWEIGHT_NEW
					auto& exp_Mweight_coarse_sub = exp_Mweight_coarse.mweightsForSomeSpinsAndSlides(iimage, iclass, idir, ipsi);
#else
					auto& exp_Mweight_coarse_sub = exp_Mweight_coarse.wptr_sparse(iimage, iclass, idir, ipsi);
#endif
					for (int itrans = 0; itrans < exp_nr_trans; itrans++)
					{
#ifdef EXP_MWEIGHT_NEW
						auto diff2 = exp_Mweight_coarse_sub.value(itrans);
#else
						int ihidden = itrans*exp_nr_over_rot + iover_rot;
						ihidden = ihidden*exp_nr_over_trans + iover_trans;
						exp_Mweight_coarse_sub[itrans].first = ihidden;
						double& diff2 = exp_Mweight_coarse_sub[itrans].second;
#endif
						diff2 += thread_exp_Mweight_sub_sub_tid[itrans*ipsi_tile + (ipsi - ipsi_tile_start)];
						if (n_end == N) {
							if (!do_cross_correlation) {
								// Calculate the actual squared difference term of the Gaussian probability function
								// If current_size < ori_size diff2 is initialised to the sum of
								// all |Xij|2 terms that lie between current_size and ori_size
								// Factor two because of factor 2 in division below, NOT because of 2-dimensionality of the complex plane!
								diff2 = (diff2 + exp_highres_Xi2_imgs.rptrAll()[iimage]) / 2;

							}
							if (diff2 < thread_exp_min_diff2_tid[iimage]) {
								thread_exp_min_diff2_tid[iimage] = diff2;
							}
						}
#ifdef EXP_MWEIGHT_NEW
						exp_Mweight_coarse_sub.setValue(itrans, diff2);
#endif
					}// end loop itrans
				}// end loop ipsi
			}// end loop iimage
		}
	}

	class TaskManagerForExpMweight {
		struct Task { int idir; int iclass; int ipsi; int iimage; };
		int maxthreads;
		std::vector<Task> AllTasks;
		std::vector<Task> scheduledTasks[256];
	public:
        TaskManagerForExpMweight() : maxthreads(omp_get_max_threads_or_one_if_no_openmp()) {
        }
		~TaskManagerForExpMweight() { AllTasks.clear(); for (auto&v : scheduledTasks) v.clear(); }
		void insertTask(int idir,int iclass,int ipsi,int iimage) {
			scheduledTasks[omp_get_thread_num()].push_back(Task{ idir,iclass,ipsi,iimage });
		}
		// rearrangment,set each thread same number of tasks
		void rearrangeTasks() {
			for (int thread = 0; thread < maxthreads; thread++) {
				for (auto& task : scheduledTasks[thread])
					AllTasks.push_back(task);
				scheduledTasks[thread].clear();
			}
			for (int i = 0; i < AllTasks.size(); i++)
				scheduledTasks[i % maxthreads].push_back(AllTasks[i]);
		}
		std::vector<Task>& getTasks() {
			return scheduledTasks[omp_get_thread_num()];
		}
	};

	void getAllSquaredDifferencesCoarse(TileCfg& tileCfg)
	{
		// set tile
		int N_tile, ipsi_tile, iimage_tile;
		int itrans_sub_tile, iimage_sub_tile;
		int N, shell_size2 = -1;
		int dummy;
		if (do_cross_correlation) {
			N = exp_current_Fsize2;
			tileCfg.getTileForCoarseSearch(N, false, N_tile, ipsi_tile, dummy, iimage_tile, iimage_sub_tile, itrans_sub_tile);
			assert(N_tile == N);	// N_tile not supported yet
		} else {
			N = fourierShellTrans.coarseFsize2();
			tileCfg.getTileForCoarseSearch(N, true,  N_tile, ipsi_tile, dummy, iimage_tile, iimage_sub_tile, itrans_sub_tile);
		}
		// If do local searching
        if (do_local_searching) {ipsi_tile = 1;}
        // ipsi_tile = 1;// DEBUG
        
		{
			MAJOR_SCOPE(getAllSquaredDifferencesCoarse_loop1)

			TaskManagerForExpMweight taskManagerForExpMweight;

			auto setMweightsForSomeSpinsAndSlides = [&](int iimage, int iclass, int idir, int ipsi) {
				const int iover_rot = 0;
				const int iover_trans = 0;
				const double diff2 = 0;
#ifdef EXP_MWEIGHT_NEW
				auto& exp_Mweight_coarse_sub = exp_Mweight_coarse.mweightsForSomeSpinsAndSlides(iimage, iclass, idir, ipsi);
#else
				auto& exp_Mweight_coarse_sub = exp_Mweight_coarse.wptr_sparse(iimage, iclass, idir, ipsi);
#endif
				for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
#ifdef EXP_MWEIGHT_NEW
					exp_Mweight_coarse_sub.insert(iover_rot, itrans, iover_trans, diff2);
#else
					int ihidden = itrans*exp_nr_over_rot + iover_rot;
					ihidden = ihidden*exp_nr_over_trans + iover_trans;
					exp_Mweight_coarse_sub.push_back(std::pair<int, double>(ihidden, diff2));
#endif
				}// end loop itrans
			};

			CHECKER_PARALLEL_COUNTED_FOR4(
				int, idir,		0,				exp_nr_dir,			/* 192(healpix_order=1)~768(healpix_order=2) */
				int, iclass,	exp_iclass_min,	exp_iclass_max+1,	/* 4~10 */
				int, ipsi,		0,				exp_nr_psi,
				int, iimage,	0,				exp_nr_images) {

				if (mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
                
                // If do local searching
                if (do_local_searching) {
                    if (!samplingGrid.isSignificantOrientation(idir, ipsi)) continue;
                    auto pdf_orientation_iimage = samplingGrid.pdfOrientationForNonzeroOrient(idir, ipsi);
                    if ((*pdf_orientation_iimage)[iimage] <= 0) continue;
                }
                else{
                    if (mlModel.pdf_direction[iclass].rptrAll()[idir] <= 0.) continue;
                }
                
				MAJOR_SCOPE2(getAllSquaredDifferencesCoarse_loop1_iter, life)

				//L2CacheModel::Interval l2CacheModelInterval;
				//l2CacheModelInterval.begin();
				taskManagerForExpMweight.insertTask(idir, iclass, ipsi, iimage);
				//setMweightsForSomeSpinsAndSlides(iimage, iclass, idir, ipsi);
				//l2CacheModelInterval.end();
				//l2CacheModelInterval.showL2CacheIntervalLocked(&life);
			} CHECKER_PARALLEL_FOR_END

			taskManagerForExpMweight.rearrangeTasks();
			//
			//std::cout << "before : " << std::endl;
			//exp_Mweight_coarse.statisticHeap();
			//
			#pragma omp parallel
			{
				for (auto &task : taskManagerForExpMweight.getTasks())
					setMweightsForSomeSpinsAndSlides(task.iimage, task.iclass, task.idir, task.ipsi);
			}

			//
			//std::cout << "after : " << std::endl;
			//exp_Mweight_coarse.statisticHeap();
		}

		{
			MAJOR_SCOPE(getAllSquaredDifferencesCoarse_loop2)

			struct Task { int idir; int iclass; };
			std::vector<Task> toDo;
			for (int iclass = exp_iclass_min; iclass < exp_iclass_max + 1; iclass++) {		// 4~10
				if (mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
				for (int idir = 0; idir < exp_nr_dir; idir++) {								// 192(healpix_order=1)~768(healpix_order=2)
                    if (!do_local_searching){
                        if (mlModel.pdf_direction[iclass].rptrAll()[idir] <= 0.) continue;
                    }
					toDo.push_back({ idir, iclass });
				}
			}

			const int toDoSize = int(toDo.size());

			// There are several choices that will affect performance in the following loop
			// They do extra computes to decrease L2 misses
			//															Lo		Hi
			//		ShiftImageInFourierTransformNew_useSinCosTable		false	true
			//		iimage_tile											1		exp_nr_images
			//		iimage_sub_tile										1		iimage_tile
			//		itrans_sub_tile										1		exp_nr_trans
			//		N_tile												1		N-n_start
			//
			//	While in theory it may be possible to determine the best values in advance, 
			//	in practise they are probably found by experimentation.
			//
			//	The following should do that experimentation.
			//	Meanwhile it is just doing L2 miss rate measurement.
			//
			bool		useRecommendedTiling		= true;
			float		recommendedTilingHitRate	= 1.0;		// until measured otherwise
			const bool	doingPeakFinding			=
#ifdef L2_CACHE_MODELING
				false && !do_cross_correlation;
#else
				false;
#endif

#ifdef L2_CACHE_MODELING
			PeakFinder::Bounds peakFinderBounds;
			auto pushBound = [&](const char* name, int default, int lo, int hi) {
				peakFinderBounds.push_back({ lo,  hi+1 });
				std::cerr << "Bound " << name << " [" << lo << ".." << hi << "] default:" << default << std::endl;
				std::cout << "Bound " << name << " [" << lo << ".." << hi << "] default:" << default << std::endl;
			};

			if (doingPeakFinding) {
				if (do_cross_correlation)	// N_tile not supported yet
					pushBound("N_tile",			N,	                N,							N                   );
				else
					pushBound("N_tile_div32",	(N_tile + 31) / 32, 4,							divRoundUp(N,32)    );

				pushBound("iimage_tile",		iimage_tile,		(iimage_tile     + 1) / 2,	2 * iimage_tile     );
				pushBound("iimage_sub_tile",	iimage_sub_tile,	(iimage_sub_tile + 1) / 2,	2 * iimage_sub_tile );
				pushBound("itrans_sub_tile",	itrans_sub_tile,	1,							8                   );
			}

			PeakFinder peakFinder(peakFinderBounds);
				// TODO learns across iterations also
#endif

#ifdef L2_CACHE_MODELING
			do 
#endif
			{
#ifdef L2_CACHE_MODELING
				if (doingPeakFinding) {
					std::cerr << (do_cross_correlation ? "GASDC, do_cc" : "GASDC, not do_cc") << ". Next trial N:" << N << std::endl;
					std::cout << (do_cross_correlation ? "GASDC, do_cc" : "GASDC, not do_cc") << ". Next trial N:" << N << std::endl;
				}
#endif

				for (int n_start = 0; n_start < N; n_start += N_tile) {
					CHECKER_PARALLEL_COUNTED_FOR2_PLUS(
						int, toDoI, 0, toDoSize, 1,
						int, ipsi_tile_start, 0, exp_nr_psi, ipsi_tile) {

						MAJOR_SCOPE2(getAllSquaredDifferencesCoarse_loop2_iter, life);

#ifdef L2_CACHE_MODELING
						L2CacheModel::Interval l2CacheModelInterval;
						if (doingPeakFinding) l2CacheModelInterval.begin();

						class Assignment : public PeakFinder::Assignment {
							const int original_N_tile;
						public:
							Assignment(int original_N_tile) : original_N_tile(original_N_tile) {}
							int N_tile, iimage_tile, iimage_sub_tile, itrans_sub_tile;
							PeakFinder::ResultIndex resultIndex;
							virtual void coords(
								PeakFinder::ResultIndex resultIndex,
								PeakFinder::Coordinates const & coordinates) {
								N_tile          = coordinates[0] * (do_cross_correlation ? 1 : 32);
								iimage_tile     = coordinates[1];
								iimage_sub_tile = coordinates[2];
								itrans_sub_tile = coordinates[3];
								this->resultIndex = resultIndex;
							}
						} assignment(N_tile);

						const auto resultIndex = useRecommendedTiling ? PeakFinder::noAssignment : peakFinder.askForAssignment(assignment);
						if (resultIndex == PeakFinder::noAssignment) {
							if (!useRecommendedTiling) {
								assignment.coords(PeakFinder::noAssignment, peakFinder.bestCoordinates());
							} else {
								assignment.N_tile          = N_tile;
								assignment.iimage_tile     = iimage_tile;
								assignment.iimage_sub_tile = iimage_sub_tile;
								assignment.itrans_sub_tile = itrans_sub_tile;
								assignment.resultIndex     = resultIndex;
							}
						}
#define ASSIGNMENT_PREFIX assignment.
#else
#define ASSIGNMENT_PREFIX
#endif

                        
                        // If do local searching
                        if (do_local_searching) {
                            int ipsi = ipsi_tile_start;
                            assert(ipsi_tile==1);
                            if (!samplingGrid.isSignificantOrientation(toDo[toDoI].idir, ipsi)) continue;
                        }
                        
						getAllSquaredDifferencesCoarse_loop2Body(tileCfg,
							N, n_start, ASSIGNMENT_PREFIX N_tile,
							toDo[toDoI].iclass, toDo[toDoI].idir,
							ipsi_tile_start, ipsi_tile,
							ASSIGNMENT_PREFIX iimage_tile,
							ASSIGNMENT_PREFIX iimage_sub_tile,
							ASSIGNMENT_PREFIX itrans_sub_tile);
#undef ASSIGNMENT_PREFIX
						if (!doingPeakFinding) continue;

#ifdef L2_CACHE_MODELING
						l2CacheModelInterval.end();
						auto l2HitRate = l2CacheModelInterval.hitRate();
						// l2CacheModelInterval.showL2CacheIntervalLocked(&life);

						if (resultIndex != PeakFinder::noAssignment) {
							peakFinder.reportResult(resultIndex, l2HitRate);
						}

						#pragma omp critical
						{
							if (useRecommendedTiling) {
								recommendedTilingHitRate = std::min(recommendedTilingHitRate, l2HitRate);
							}

							auto msg = [&](std::ostream&os) {
								os << "GASDC l2HitRate:" << std::setprecision(3) << l2HitRate << " tid:" << omp_get_thread_num() << " ri:" << resultIndex
									<< " n_tile:"          << assignment.N_tile
									<< " iimage_tile:"	   << assignment.iimage_tile
									<< " iimage_sub_tile:" << assignment.iimage_sub_tile
									<< " itrans_sub_tile:" << assignment.itrans_sub_tile
									<< ", best yet:"       << peakFinder.bestResult()
									<< "  default:"		   << recommendedTilingHitRate
									<< std::endl;
							};
							msg(std::cerr);
							msg(std::cout);
						}
#endif
					} CHECKER_PARALLEL_FOR_END
				}	// end loop n_tile
				
#ifdef L2_CACHE_MODELING
				if (!doingPeakFinding) break;

				useRecommendedTiling = false;
#endif
			} 
#ifdef L2_CACHE_MODELING
			while (!peakFinder.hasFoundPeak());

			if (doingPeakFinding) {
				auto msg = [&](std::ostream& os) {
					os << "GASDC PeakFinder::bestResult:" << peakFinder.bestResult();
					auto & bestCoords = peakFinder.bestCoordinates();
					const char* sep = " found at ";
					for (auto c : bestCoords) {
						os << sep << c;
						sep = ",";
					}
					os << " in a volume " << peakFinder.volume() << " searched with " << peakFinder.numberAssignmentsReported() << " trials" << std::endl;
				};
				msg(std::cout);
			}
#endif
		}
	}

	struct ProjectedClassCache_Key {
		// in most-likely-to-differ-first order
#define PROJECTEDCLASSCACHE_KEY_ELTS \
            ELT(int		, iclass	)  SEP \
            ELT(int		, idir		)  SEP \
			ELT(int		, ipsi		)  SEP \
            ELT(int		, iover_rot	)  SEP \
			ELT(FDOUBLE	, over_rot 	)  SEP \
			ELT(FDOUBLE	, over_tilt )  SEP \
            ELT(FDOUBLE	, over_psi 	)  SEP \
			ELT(int		, n_start	)  SEP \
			ELT(int		, n_end		)	   \
			// end of macro
#define SEP
#define ELT(T,N) T N;
		PROJECTEDCLASSCACHE_KEY_ELTS
#undef ELT
#undef SEP
			ProjectedClassCache_Key(
			) :
#define SEP ,
#define ELT(T,N) N(-1)
			PROJECTEDCLASSCACHE_KEY_ELTS
#undef ELT
#undef SEP
		{}
			ProjectedClassCache_Key(
#define SEP ,
#define ELT(T,N) T N
				PROJECTEDCLASSCACHE_KEY_ELTS
#undef ELT
#undef SEP
			) :
#define SEP ,
#define ELT(T,N) N(N)
			PROJECTEDCLASSCACHE_KEY_ELTS
#undef ELT
#undef SEP
		{}
			bool operator==(ProjectedClassCache_Key const & rhs) const {
			return
#define SEP &&
#define ELT(T,N) (N == rhs.N)
				PROJECTEDCLASSCACHE_KEY_ELTS
#undef ELT
#undef SEP
				;
		}
	};

	static std::ostream& operator<< (std::ostream & os, ProjectedClassCache_Key const & key) {
		os << "{"
#define SEP << ", "
#define ELT(T,N) << #N << ":" << std::setw(4) << key.N
			PROJECTEDCLASSCACHE_KEY_ELTS
#undef ELT
#undef SEP
			<< "}";
		return os;
	}

	class ProjectedClassCache {
	public:
		typedef ProjectedClassCache_Key Key;

		static const size_t capacity  = 16;
		static const bool   showStats = false;
		//#define UNDERSTAND_MISSES

	private:

		struct Entry : public NoCopy {
			FDOUBLE * real;
			FDOUBLE * imag;
			size_t    inUse;
			Entry() : real(nullptr), imag(nullptr), inUse(0) {}
			~Entry() {
				if (real) {
					Heap::freeScalars<FDOUBLE>(real);
					Heap::freeScalars<FDOUBLE>(imag);
				}
			}
		};

		char  keyIndexs[capacity];	// avoids having to move the keys and entries
		Key   keys     [capacity];
		Entry entries  [capacity];

#ifndef  UNDERSTAND_MISSES
		void maybeShowMissLog(bool found, Key const & key) {}
		void maybeShowMissLog() {}
#else
		InitZeroSingleton<size_t> missLogSize, missLogShows, missCount, longDistanceHits;
		Key                       missLog[capacity];
		void maybeShowMissLog(bool found, Key const & key) {
			if (omp_get_thread_num() != 0) return;
			if (found) {
				if (missLogSize.v == capacity) {
					if (missLogShows.v++ < 10) showMissLog();
				}
				missLogSize.v = 0;
				missCount  .v = 0;
			} else {
				if ((missLogSize.v > 1) && (key == missLog[0])) {
					longDistanceHits.v++;
					std::cerr << "Found " << key << " after " << missCount.v << " lookups <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
					missLogSize.v = 0;
					missCount  .v = 0;
				} else {
					missCount.v++;
				}
			}
			if (missLogSize.v < capacity) {
				missLog[missLogSize.v++] = key;
			}
		}
		void maybeShowMissLog() {
			std::cerr << "~ProjectedClassCache longDistanceHits:" << longDistanceHits.v << std::endl;
			if (missLogShows.v == 0) showMissLog();
		}
		void showMissLog() {
			for (size_t i = 0; i < missLogSize.v; i++) {
				auto& key = missLog[i];
				std::cerr << "missLog[" << std::setw(2) << i << "]:" << key << std::endl;
			}
		}
#endif

		size_t lookups;
		size_t hits    [capacity];
		size_t misses;
	public:
		ProjectedClassCache() : lookups(0), misses(0) {
			for (size_t row = 0; row < capacity; row++) {
				keyIndexs[row] = char(row); hits[row] = 0;
			}
		}
		~ProjectedClassCache() {
			if (showStats) {
				bool headingDone = false;
				auto heading = [&]() {
					if (headingDone) return;
					headingDone = true;
					std::cerr << "ProjectedClassCache stats" << std::endl;
				};
				size_t sum(0);
				for (size_t row = 0; row < capacity; row++) {
					if (!hits[row]) continue;
					sum += hits[row];
					heading(); std::cerr << "row:" << row << " hits:" << hits[row] << std::endl;
				}
				sum += misses;
				if (sum) { heading(); std::cerr << "Misses:" << misses << "/" << sum << std::endl; }
			}
			maybeShowMissLog();
		}

		class Projection {
			const FDOUBLE * _real;
			const FDOUBLE * _imag;
			int             _keyIndex;
		public:
			Projection() : _real(nullptr), _imag(nullptr), _keyIndex(-1) {}
			Projection(const FDOUBLE * const real, const FDOUBLE * const imag, int keyIndex) : _real(real), _imag(imag), _keyIndex(keyIndex) {}
			Projection& operator=(Projection const & rhs) {
				_real = rhs._real; _imag = rhs._imag; _keyIndex = rhs._keyIndex;
				return *this;
			}
			const FDOUBLE * real    () const { return _real    ;}
			const FDOUBLE * imag    () const { return _imag    ;}
			int             keyIndex() const { return _keyIndex;}
		};

		Projection cachedProjectedClass(Key const & key) {

			// search for this one
			auto prev = int(keyIndexs[capacity - 1]);		// assumes don't find it, assumes this one is not in use
			assert(entries[prev].inUse == 0);
			bool found = false;
			for (size_t row = 0; row < capacity; row++) {
				auto curr = keyIndexs[row];
				keyIndexs[row] = prev;						// rotates the misses deeper into the cache
				auto& rowKey = keys[curr];
				found = (rowKey == key);
				if (found) {								// hit
					hits[row]++;
					keyIndexs[0] = curr;					// move it to the front
					break;
				}
				prev = curr;
			}

			lookups++;
			if (showStats && (lookups < 2*capacity) && (omp_get_thread_num() == 0))
#pragma omp critical
			{
				std::cerr << "cachedProjectedClass(" << key << ")" << (found ? "    hit" : "    miss") << std::endl;
			}

			// make one if not found
			if (found) {
				L2CacheModel::seqAcc("cachedProjectedClass hit", -1, (FDOUBLE*)nullptr, 0, 0);
				prev = keyIndexs[0];
				maybeShowMissLog(true, key);
			} else {
				maybeShowMissLog(false, key);
				misses++;
				assert(prev == keyIndexs[0]);
				auto & entry = entries[prev];
				assert(entry.inUse == 0);
				auto& rowKey = keys[prev];
				rowKey = key;
				if (!entry.real) {
					entry.real = Heap::allocScalars<FDOUBLE>(exp_current_Fsize2, __FILE__, __LINE__);
					entry.imag = Heap::allocScalars<FDOUBLE>(exp_current_Fsize2, __FILE__, __LINE__);
				}
				if (do_cross_correlation) {
					L2CacheModel::seqAcc("cachedProjectedClass miss CC", -1, (FDOUBLE*)nullptr, 0, 0);
					getProjectedClass(key.iclass, key.over_rot, key.over_tilt, key.over_psi, entry.real, entry.imag, key.n_start, key.n_end);
				} else if (do_coarse_search) {
					L2CacheModel::seqAcc("cachedProjectedClass miss coarse", -1, (FDOUBLE*)nullptr, 0, 0);
					getProjectedClassByShellCoarse(key.iclass, key.over_rot, key.over_tilt, key.over_psi, entry.real, entry.imag, key.n_start, key.n_end);
				} else {
					L2CacheModel::seqAcc("cachedProjectedClass miss fine", -1, (FDOUBLE*)nullptr, 0, 0);
					getProjectedClassByShellFine(key.iclass, key.over_rot, key.over_tilt, key.over_psi, entry.real, entry.imag, key.n_start, key.n_end);
				}
			}
			// increase the count
			assert(prev == keyIndexs[0]);
			auto& rowKey = keys[prev];
			assert(rowKey == key);
			auto& entry = entries[prev];
			entry.inUse++;
			// return it
			return Projection(entry.real, entry.imag, prev);
		}

		void release(Projection const & projection) {
			auto& entry = entries[projection.keyIndex()];
			assert(entry.inUse > 0);
			entry.inUse--;
		}
	};

	IntPerformanceCounter gasdf_loadingKernel       ("gasdf_loadingKernel");
	IntPerformanceCounter gasdf_avoidedLoadingKernel("gasdf_avoidedLoadingKernel");

	size_t getAllSquaredDifferencesFineCC(					// returns count of significant
		ProjectedClassCache& projectedClassCache,
		double* thread_exp_Mweight_sub_tid,
		const char* exp_Mcoarse_significant_rptr,
		int tid, int iclass, int idir, int ipsi,
		int n_start, int n_end,
		int iover_rot_tile_start, int iover_rot_tile_end,
		int iimage_tile_start, int iimage_tile_end)
	{
#ifndef USE_BEFORE_PROJECTION_CACHE
		// poison these so not accidently used
		typedef void thread_Frefctf_real;
		typedef void thread_Frefctf_imag;
#endif

		assert(n_start == 0); assert(n_end == exp_current_Fsize2);
        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                    exp_over_rot, exp_over_tilt, exp_over_psi);
        
        // If do local search
        const std::vector<FDOUBLE>* pdf_orientation_iimage;
        if (do_local_searching) {
            pdf_orientation_iimage = samplingGrid.pdfOrientationForNonzeroOrient(idir, ipsi);
        }
		// NOTE : no need for tile,only one (iimage,itrans) is significant
		//
		size_t significantCount(0);
		for (int iimage = iimage_tile_start; iimage < iimage_tile_end; iimage++)
		{
            // If do local search
            if (do_local_searching) if(!(*pdf_orientation_iimage)[iimage]) continue;
			for (int itrans = 0; itrans < exp_nr_trans; itrans++) // 100(offset_range=10,step=2) or 400(offset_range=10,step=1)
			{
				if (!exp_Mcoarse_significant_rptr[iimage*exp_nr_trans + itrans])
					continue;

				significantCount += 1;
				transformAndGetDiffCC_performanceCounter.count.v += exp_nr_over_trans*(iover_rot_tile_end - iover_rot_tile_start);

				auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid +
					((iimage - iimage_tile_start) * exp_nr_trans + itrans) *
					exp_nr_over_trans * exp_nr_over_rot;

				for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)//4(oversampling=1)
				{
					for (int iover_rot = iover_rot_tile_start; iover_rot < iover_rot_tile_end; iover_rot++) //8(oversampling=1)
					{
						// Store all diff2 in exp_Mweight
						int ihidden_tid = iover_trans*exp_nr_over_rot + iover_rot;
						thread_exp_Mweight_sub_sub_tid[ihidden_tid] = 0;
						assert(thread_exp_Mweight_sub_sub_tid[ihidden_tid] == 0);
						auto& diff2 = thread_exp_Mweight_sub_sub_tid[ihidden_tid];
                        auto over_rot = exp_over_rot[iover_rot];
                        auto over_tilt = exp_over_tilt[iover_rot];
                        auto over_psi = exp_over_psi[iover_rot];
						// projected class
						auto projection = projectedClassCache.cachedProjectedClass(ProjectedClassCache::Key(
#define SEP ,
#define ELT(T,N) N
							PROJECTEDCLASSCACHE_KEY_ELTS
#undef ELT
#undef SEP
						));


#ifdef USE_PROJECTION_CACHE
						auto Fref_real = projection.real();
						auto Fref_imag = projection.imag();
#endif

#ifdef USE_BEFORE_PROJECTION_CACHE
						auto Fref_real_old = thread_Frefctf_real[tid][iover_rot - iover_rot_tile_start].wptr(exp_current_Fsize2);
						auto Fref_imag_old = thread_Frefctf_imag[tid][iover_rot - iover_rot_tile_start].wptr(exp_current_Fsize2);
                        getProjectedClass(iclass, idir, ipsi, iover_rot, n_start, n_end, Fref_real_old, Fref_imag_old);
						                

#ifdef USE_PROJECTION_CACHE
						// TODO - compare them
#else
						auto Fref_real = Fref_real_old;
						auto Fref_imag = Fref_imag_old;
#endif

#endif

						// compute difference
						auto local_ctf = refs_are_ctf_corrected ? exp_local_Fctfs[iimage].rptrAll() : particleModel.FctfsOne[iimage].rptrAll();
						double myscale = exp_metadata[iimage].SCALE;
						auto shifter = particleModel.shiftImageAssistor.acquireShifter(itrans, iover_trans, itrans*exp_nr_over_trans + iover_trans);
						transformAndGetDiffCC(iimage, diff2, n_start, n_end,
							particleModel.Fimages_mask_fine_real[iimage].rptrAll(),
							particleModel.Fimages_mask_fine_imag[iimage].rptrAll(),
							shifter->aTable_rptr(),
							shifter->bTable_rptr(),
							local_ctf, myscale, Fref_real, Fref_imag);
						particleModel.shiftImageAssistor.releaseShifter(shifter);

						projectedClassCache.release(projection);

					} // end loop iover_rot
				} // end loop iover_trans
			} // end loop itrans
		}// end loop iimage

		return significantCount;
	}

	// do maximum likelihood
	size_t getAllSquaredDifferencesFineML(					// returns count of significant
		ProjectedClassCache& projectedClassCache,
		double* thread_exp_Mweight_sub_tid,
		const char* exp_Mcoarse_significant_rptr,
		int tid, int iclass, int idir, int ipsi,
		int shell_n_start, int shell_n_end,
		int iover_rot_tile_start, int iover_rot_tile_end,
		int iimage_tile_start, int iimage_tile_end,
		int iimage_sub_tile, int itrans_sub_tile)
	{
		GetAllSquaredDifferencesFine_Kernel kernel(exp_current_size);
        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                    exp_over_rot, exp_over_tilt, exp_over_psi);
		// get projected class
		size_t significantCount(0);
        // If do local search
        const std::vector<FDOUBLE>* pdf_orientation_iimage;
        if (do_local_searching) pdf_orientation_iimage = samplingGrid.pdfOrientationForNonzeroOrient(idir, ipsi);
        
		std::vector<ProjectedClassCache::Projection> projections;

		auto loadProjectionsIntoKernel = [&]() {
			projections.resize(iover_rot_tile_end - iover_rot_tile_start);

			for (int iover_rot = iover_rot_tile_start; iover_rot < iover_rot_tile_end; iover_rot++)
			{
#ifdef USE_PROJECTION_CACHE
                auto over_rot = exp_over_rot[iover_rot];
                auto over_tilt = exp_over_tilt[iover_rot];
                auto over_psi = exp_over_psi[iover_rot];
				auto n_start = shell_n_start;
				auto n_end = shell_n_end;
				auto & projection = projections[iover_rot - iover_rot_tile_start];
				projection = projectedClassCache.cachedProjectedClass(ProjectedClassCache::Key(
#define SEP ,
#define ELT(T,N) N
					PROJECTEDCLASSCACHE_KEY_ELTS
#undef ELT
#undef SEP
				));

				auto Fref_real = projection.real();
				auto Fref_imag = projection.imag();
#endif

#ifdef USE_BEFORE_PROJECTION_CACHE
				auto Fref_real_old = thread_Frefctf_real[tid][iover_rot - iover_rot_tile_start].wptrAll();
				auto Fref_imag_old = thread_Frefctf_imag[tid][iover_rot - iover_rot_tile_start].wptrAll();
				// get projected image
				getProjectedClassByShellFine(iclass, idir, ipsi, iover_rot, shell_n_start, shell_n_end, Fref_real_old, Fref_imag_old);

#ifdef USE_PROJECTION_CACHE
				// TODO - compare
#else
				auto Fref_real = Fref_real_old;
				auto Fref_imag = Fref_imag_old;
#endif

#endif
				kernel.appendFref(Fref_real, Fref_imag);
			}
		};

		// the compare matrix is like :
		// 1) the significant image is less,0%~10%
		// 2) if the image is significant,then most of the translation is significant
		// so put the translation loop outside will be a good choice,and also the double-translation is not work very well
		// next try the triple-translation,it will reduce the ABTable to only ten.
		for (int iimage_sub_tile_start = iimage_tile_start; iimage_sub_tile_start < iimage_tile_end; iimage_sub_tile_start += iimage_sub_tile)
		{
			for (int itrans_sub_tile_start = 0; itrans_sub_tile_start < exp_nr_trans; itrans_sub_tile_start += itrans_sub_tile)
			{// 100(offset_range=10,step=2) or 400(offset_range=10,step=1)
				int iimage_sub_tile_end = std::min(iimage_sub_tile_start + iimage_sub_tile, iimage_tile_end);
				int itrans_sub_tile_end = std::min(itrans_sub_tile_start + itrans_sub_tile, exp_nr_trans);
				for (int iimage = iimage_sub_tile_start; iimage < iimage_sub_tile_end; iimage++) {
                    // If do local search
                    if (do_local_searching && !(*pdf_orientation_iimage)[iimage]) continue;
                    //
					for (int itrans = itrans_sub_tile_start; itrans < itrans_sub_tile_end; itrans++) {

						if (!exp_Mcoarse_significant_rptr[iimage*exp_nr_trans + itrans])
							continue;

						significantCount++;

						if (projections.size() == 0) loadProjectionsIntoKernel();

						//
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
						auto Fimg_real = thread_Fimg_real[tid][0].wptrAll();
						auto Fimg_imag = thread_Fimg_imag[tid][0].wptrAll();
#endif
#if defined(DOUBLE_TRANSLATION)
						particleModel.getLargeShiftedMaskImageOneTile(iimage, Fimg_real, Fimg_imag, shell_n_start, shell_n_end, itrans);
#elif defined(TRIPLE_TRANSLATION)
                        particleModel.getLargeShiftedMaskImageDecompOneTile(iimage, Fimg_real, Fimg_imag, shell_n_start, shell_n_end, itrans, samplingGrid);
#endif
						//
						auto Minvsigma2_iimage = exp_local_Minvsigma2s[iimage].wptrAll();
						double myscale = exp_metadata[iimage].SCALE;
                        auto local_ctf = refs_are_ctf_corrected ? exp_local_Fctfs[iimage].rptrAll() : particleModel.FctfsOne[iimage].rptrAll();
						auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid +
							((iimage - iimage_tile_start) * exp_nr_trans + itrans) *
							exp_nr_over_trans * exp_nr_over_rot;
						//
						kernel.acquireImage(
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
							Fimg_real, Fimg_imag,
#else
							particleModel.Fimages_mask_fine_real[iimage].rptrAll(),
							particleModel.Fimages_mask_fine_imag[iimage].rptrAll(),
#endif
							local_ctf,
							Minvsigma2_iimage, myscale, shell_n_start, shell_n_end
						);

						for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)//4(oversampling=1)
						{
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
							auto shifter = particleModel.smallShiftedABTable.acquireShifter(0,iover_trans,iover_trans);
#else
							auto shifter = particleModel.shiftImageAssistor.acquireShifter(itrans, iover_trans, itrans*exp_nr_over_trans + iover_trans);
#endif
							kernel.acquireTable(
								shifter->aTable_rptr(),
								shifter->bTable_rptr(),
								shifter->tableSize());

							kernel.appendDiff2(
								iover_rot_tile_end - iover_rot_tile_start);

							kernel.compute();
							// Store all diff2 in exp_Mweight,TODO....
							kernel.release(
								thread_exp_Mweight_sub_sub_tid + iover_trans*exp_nr_over_rot + iover_rot_tile_start,
								iover_rot_tile_end - iover_rot_tile_start);
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
							particleModel.smallShiftedABTable.releaseShifter(shifter);
#else
							particleModel.shiftImageAssistor.releaseShifter(shifter);
#endif

						} // end loop iover_trans
					} // end loop itrans
				} // end loop iimage
			} // end loop itrans sub tile
		}// end loop iimage sub tile

		((significantCount == 0) ? gasdf_avoidedLoadingKernel : gasdf_loadingKernel).count.v++;

		for (auto & projection : projections) {
			projectedClassCache.release(projection);
		}

		return significantCount;
	}

	typedef std::vector<int> SignificantCountHistogram;

	class GetAllSquaredDifferencesFineMain_Task {
	public:
#define GASDFM_INNER_ELTS \
		ELT(int, iclass				 )  SEP \
		ELT(int, idir				 )	SEP \
		ELT(int, n_start			 )	SEP \
		ELT(int, n_end			     )	SEP \
		ELT(int, N					 )	SEP \
		ELT(int, ipsi_tile_start	 )	SEP \
		ELT(int, ipsi_tile_end		 )	SEP \
		ELT(int, iover_rot_tile_start)	SEP \
		ELT(int, iover_rot_tile_end	 )	SEP \
		ELT(int, iimage_tile_start	 )	SEP \
		ELT(int, iimage_tile_end	 )	SEP \
		ELT(int, iimage_sub_tile	 )	SEP \
		ELT(int, itrans_sub_tile	 )		\
		// end of macro
#define SEP
#define ELT(T,N) T N;
		GASDFM_INNER_ELTS
#undef ELT
#undef SEP
		bool performWillDoAtLeastOne();
		int perform(
				ProjectedClassCache& projectedClassCache,
				SignificantCountHistogram & significantCountHistogram);
	};

	IntPerformanceCounter gasdf_inserted("GetAllSquaredDifferencesFineMain_Task_Buffer_inserted");
	IntPerformanceCounter gasdf_not_inserted("GetAllSquaredDifferencesFineMain_Task_Buffer_not_inserted");

	class GetAllSquaredDifferencesFineMain_Task_Buffer {
		SignificantCountHistogram & _significantCountHistogram;
	public:
		GetAllSquaredDifferencesFineMain_Task_Buffer(SignificantCountHistogram & significantCountHistogram) 
			: _capacity(4096), 
			_buffer(vNew(GetAllSquaredDifferencesFineMain_Task, _capacity)), 
			_size(0), 
			_projectedClassCache(omp_get_max_threads()), 
			_significantCountHistogram(significantCountHistogram),
			_significantCount(0) {
		}
		~GetAllSquaredDifferencesFineMain_Task_Buffer() { flush(); vDelete(_buffer); }
		bool maybeInsert(
#define SEP ,
#define ELT(T,N) T N
			GASDFM_INNER_ELTS
#undef ELT
#undef SEP
		) {
			if (_size == _capacity) flush();
			auto& t = _buffer[_size];
#define SEP
#define ELT(T,N) t.N = N;
			GASDFM_INNER_ELTS
#undef ELT
#undef SEP
			if (!t.performWillDoAtLeastOne()) {
				gasdf_not_inserted.count.v++;
				return false;
			}

			gasdf_inserted.count.v++;
			_size++;
			return true;
		}
	private:
		const int							    _capacity;
		GetAllSquaredDifferencesFineMain_Task*  _buffer;
		int										_size;
		std::vector<ProjectedClassCache>		_projectedClassCache;
		std::atomic<int>						_significantCount;
	public:
		void flush() {
			#pragma omp parallel for schedule(dynamic)
			for (int i = 0; i < _size; i++) {
				_significantCount += 
				_buffer[i].perform(
					_projectedClassCache[omp_get_thread_num()],
					_significantCountHistogram);
			}
			_size = 0;
		}
		size_t significantCount() const { return _significantCount; }
		void zeroSignificantCount() { _significantCount = 0; }
	};

    void getAllSquaredDifferencesFineMain(TileCfg tileCfg)
    {
        assert(exp_ipass==1);
        // set tile
        int N_tile,ipsi_tile,iover_rot_tile,iimage_tile;
        int itrans_sub_tile,iimage_sub_tile;
        int N,shell_size2 = -1;
        if (do_cross_correlation) {
            N	= exp_current_Fsize2;
            tileCfg.getTileForFineCC(N, N_tile, ipsi_tile, iover_rot_tile, iimage_tile, iimage_sub_tile, itrans_sub_tile);
			assert(N == N_tile); // because cross_correlation doesn't support N tiling yet
		} else {
            N 	= fourierShellTrans.fineFsize2();
            tileCfg.getTileForFineSearch(N, N_tile, ipsi_tile, iover_rot_tile, iimage_tile, iimage_sub_tile, itrans_sub_tile);
        }

        // If do local searching
        if (do_local_searching) {ipsi_tile = 1;}

        // particleModel.testDoubleOrTripleTranslation(samplingGrid, exp_current_size);
        // EXIT_ABNORMALLY;
        // set exp_Mweight
		{
			MAJOR_SCOPE(getAllSquaredDifferencesFineMain_loop1)

			TaskManagerForExpMweight taskManagerForExpMweight;

			CHECKER_PARALLEL_COUNTED_FOR3(
				int, iimage, 0, exp_nr_images,
				int, iclass, exp_iclass_min, exp_iclass_max + 1,
				int, idir, 0, exp_nr_dir) {

				if (mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;

				for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
				{
					// If do local searching
					if (do_local_searching) {
						if (!samplingGrid.isSignificantOrientation(idir, ipsi)) continue;// may not need
						auto pdf_orientation_iimage = samplingGrid.pdfOrientationForNonzeroOrient(idir, ipsi);
						if ((*pdf_orientation_iimage)[iimage] <= 0) continue;
					}
					else {
						if (mlModel.pdf_direction[iclass].rptrAll()[idir] <= 0.) continue;
					}
					// int iClassOrient = iclass*exp_nr_dir+idir;
					// iClassOrient = iClassOrient*exp_nr_psi+ipsi;
					if (!exp_Mcoarse_Rot_significant.isRotSignificant(iclass, idir, ipsi))
						continue;
					// exp_iclass loop does not always go from 0 to nr_classes!
					auto exp_Mcoarse_significant_rptr = exp_Mcoarse_Rot_significant.isMcoarseSignificantRptrAll(iclass, idir, ipsi);
					for (int itrans = 0; itrans < exp_nr_trans; itrans++)
					{
						if (exp_Mcoarse_significant_rptr[iimage*exp_nr_trans + itrans]) {
							taskManagerForExpMweight.insertTask(idir, iclass, ipsi, iimage);
							break;
						}
					}// end loop itrans
				}// end loop ipsi
			} CHECKER_PARALLEL_FOR_END

			auto setMweightsForSomeSpinsAndSlidesOrNullptr = [&](int iimage, int iclass, int idir, int ipsi)
			{
#ifdef EXP_MWEIGHT_NEW
				auto exp_Mweight_fine_sub_ptr = exp_Mweight_fine.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi, true);
				exp_Mweight_fine_sub_ptr->zeroSomeMakeSureRestUnused();
#else
				auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(iimage, iclass, idir, ipsi);
#endif
				auto exp_Mcoarse_significant_rptr = exp_Mcoarse_Rot_significant.isMcoarseSignificantRptrAll(iclass, idir, ipsi);
				for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
					if (!exp_Mcoarse_significant_rptr[iimage*exp_nr_trans + itrans])
						continue;
					// The set of entries inserted here is used later to
					// pull the values out of the thread_exp_Mweight_sub_sub_tid
					// Search for [[which values come from the thread_exp_Mweight_sub_sub_tid]]
					//
					for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
						for (int iover_rot = 0; iover_rot < exp_nr_over_rot; iover_rot++) {
#ifdef EXP_MWEIGHT_NEW
							exp_Mweight_fine_sub_ptr->insert(
								iover_rot, itrans, iover_trans, 0.0);
#else
							int ihidden = itrans*exp_nr_over_trans + iover_trans;
							ihidden = ihidden*exp_nr_over_rot + iover_rot;
							exp_Mweight_fine_sub.push_back(std::pair<int, double>(ihidden, 0));
#endif
						}// end loop iover_trans
					}// end loop iover_rot
				}// end loop itrans
			};

			taskManagerForExpMweight.rearrangeTasks();

			//std::cout << "before" << std::endl;
			//exp_Mweight_fine.statisticHeap();

			#pragma omp parallel
			{
				for (auto &task : taskManagerForExpMweight.getTasks())
					setMweightsForSomeSpinsAndSlidesOrNullptr(task.iimage, task.iclass, task.idir, task.ipsi);
			}

			//std::cout << "after" << std::endl;
			//exp_Mweight_fine.statisticHeap();
		}

		{
			MAJOR_SCOPE(getAllSquaredDifferencesFineMain_loop2)

			bool		useRecommendedTiling = true;		// The first time through do this, to get a base line
			float		recommendedTilingHitRate = 1.0;		// Assume perfect until measured otherwise
			const bool	doingPeakFinding =
#ifdef L2_CACHE_MODELING
				false && do_cross_correlation;				// Get the CC one tiled correctly first
#else
				false;
#endif
#ifdef L2_CACHE_MODELING
			PeakFinder::Bounds peakFinderBounds;
			auto pushBound = [&](const char* name, int default, int lo, int hi) {
				peakFinderBounds.push_back({ lo,  hi + 1 });
				std::cerr << "Bound " << name << " [" << lo << ".." << hi << "] default:" << default << std::endl;
				std::cout << "Bound " << name << " [" << lo << ".." << hi << "] default:" << default << std::endl;
			};
			if (doingPeakFinding) {
				if (do_cross_correlation)	// N_tile not supported yet
					pushBound("N_tile", N, N, N);
				else
					pushBound("N_tile_div32", (N_tile + 31) / 32, 4, divRoundUp(N, 32));
				pushBound("ipsi_tile", ipsi_tile, (ipsi_tile + 1) / 2, 2 * ipsi_tile);
				pushBound("iover_rot_tile", iover_rot_tile, (iover_rot_tile + 1) / 2, 2 * iover_rot_tile);
				pushBound("iimage_tile", iimage_tile, (iimage_tile + 1) / 2, 2 * iimage_tile);
				pushBound("iimage_sub_tile", iimage_sub_tile, (iimage_sub_tile + 1) / 2, 2 * iimage_sub_tile);
				pushBound("itrans_sub_tile", itrans_sub_tile, 1, 8);
			}
			PeakFinder peakFinder(peakFinderBounds);
#endif
#ifdef L2_CACHE_MODELING
			do
#endif
			{
#ifdef L2_CACHE_MODELING
				bool reportedNextTrial(false);
#endif
				if (do_cross_correlation) if (N != N_tile) EXIT_ABNORMALLY; // because cross_correlation doesn't support N tiling yet

				SignificantCountHistogram significantCountHistogram(100);
				for (auto & i : significantCountHistogram) i = 0;

				// in first iteration does not do Ntile because of the cross_correlation
				for (int n_start = 0; n_start < N; n_start += N_tile)
				{
					const int n_end = std::min(n_start + N_tile, N);

					GetAllSquaredDifferencesFineMain_Task_Buffer buffer(significantCountHistogram);

					for (int idir = 0; idir < exp_nr_dir; idir++) // 192(healpix_order=2)~768(healpix_order=3)
					{
						// There is a brick, the significant rays along the N axis in its volume have to be visited
						//
						//		One face of the brick is	indexed by idir, iclass, ipsi, iover_rot
						//			each of these is an expensive projection, and so you want to do a few of them into the L2 cache and reuse them from there as many times as possible
						//			before doing any more of them.
						//
						//			The projection requires a 2D slice of the 3D class, and some tables that are specific to the ipsi,iover_rot
						//			I believe the rot acceleration table data is probably larger than the class data, but this needs to be verified
						//			so it is more important to reuse the rot acceleration table data than to reuse the class data,
						//			but it would be even better to redo both by tiling
						//
						//		Another face of the brick is indexed by image, itrans, iovertrans
						//			The curious thing about this face is that translating by 2 units is the same as translating by 1 unit and then doing it again
						//			So if you need to visit cells (iimage, itrans, iovertrans==0), (iimage, itrans, iovertrans==1) (iimage, itrans, iovertrans==2) ...
						//			you can do it by doing iimage,itrans,0  and then (without even storing the data) compute iovertrans==1, and then iovertrans==2 (again without storing the data)
						//			which drastically reduces the number of tables that are needed and the amount of memory traffic
						//
						//		The depth of the brick is the N-length rays.
						//		This dimension can be chunked into portions at a slight cost to reduce the sizes of the above items to get them to fit into cache
						//
						//		The following code is attempting to subdivide this brick into subbricks, each of which can be efficiently processed within about 1/3rd of the L2 cache
						//		and get lots of L2 cache hits while doing so.
						//
						for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++) // 4~10
						{	// 24(healpix_order=2)~48(healpix_order=3) x 8(oversampling=1)

							if (mlModel.pdf_class.rptrAll()[iclass] <= 0.) continue;
                            if (!do_local_searching) {
                                if (mlModel.pdf_direction[iclass].rptrAll()[idir] <= 0.) continue;
                            }
							// Get the tiling to use for this group
							// Note: We explore different N_tile to get times, but the outputs are bogus because the result might not do all the 
							//		n values because we use the default N_tile in the above outer loop
							//
#ifdef L2_CACHE_MODELING
							class Assignment : public PeakFinder::Assignment {
								const int original_N_tile;
							public:
								Assignment(int original_N_tile) : original_N_tile(original_N_tile) {}
								int N_tile, ipsi_tile, iover_rot_tile, iimage_tile, iimage_sub_tile, itrans_sub_tile;
								PeakFinder::ResultIndex resultIndex;
								virtual void coords(
									PeakFinder::ResultIndex resultIndex,
									PeakFinder::Coordinates const & coordinates) {
									N_tile				= coordinates[0] * (do_cross_correlation ? 1 : 32);
									ipsi_tile			= coordinates[1];
									iover_rot_tile		= coordinates[2];
									iimage_tile			= coordinates[3];
									iimage_sub_tile		= coordinates[4];
									itrans_sub_tile		= coordinates[5];
									this->resultIndex	= resultIndex;
								}
							} assignment(N_tile);

							const auto resultIndex = useRecommendedTiling ? PeakFinder::noAssignment : peakFinder.askForAssignment(assignment);
							if (resultIndex == PeakFinder::noAssignment) {
								if (!useRecommendedTiling) {
									assignment.coords(PeakFinder::noAssignment, peakFinder.bestCoordinates());
								} else {
									assignment.N_tile          = N_tile;
									assignment.ipsi_tile       = ipsi_tile;
									assignment.iover_rot_tile  = iover_rot_tile;
									assignment.iimage_tile     = iimage_tile;
									assignment.iimage_sub_tile = iimage_sub_tile;
									assignment.itrans_sub_tile = itrans_sub_tile;
									assignment.resultIndex     = resultIndex;
								}
							}

							auto const N_tile		   = assignment.N_tile;
							auto const ipsi_tile       = assignment.ipsi_tile;
							auto const iover_rot_tile  = assignment.iover_rot_tile;
							auto const iimage_tile	   = assignment.iimage_tile;
							auto const iimage_sub_tile = assignment.iimage_sub_tile;
							auto const itrans_sub_tile = assignment.itrans_sub_tile;
#endif
							for (int ipsi_tile_start = 0; ipsi_tile_start < exp_nr_psi; ipsi_tile_start += ipsi_tile)
							{
								const int ipsi_tile_end = std::min(ipsi_tile_start + ipsi_tile, exp_nr_psi);

                                // If do local searching
                                if (do_local_searching) {
                                    assert(ipsi_tile==1);
                                    assert(ipsi_tile_start+1==ipsi_tile_end);
                                    if (!samplingGrid.isSignificantOrientation(idir, ipsi_tile_start)) continue;
                                }
                                
								for (int iover_rot_tile_start = 0; iover_rot_tile_start < exp_nr_over_rot; iover_rot_tile_start += iover_rot_tile)// 8
								{
									// IT LOOKS IMPLEMENTED SO TRY IT OUT	assert(iover_rot_tile == exp_nr_over_rot);	// TODO
									const int iover_rot_tile_end = std::min(iover_rot_tile_start + iover_rot_tile, exp_nr_over_rot);

									for (int iimage_tile_start = 0; iimage_tile_start < exp_nr_images; iimage_tile_start += iimage_tile) // =(nr_pool)
									{
										const int iimage_tile_end = std::min(iimage_tile_start + iimage_tile, exp_nr_images);

										if (buffer.maybeInsert(
#define SEP ,
#define ELT(T,N) N
											GASDFM_INNER_ELTS
#undef ELT
#undef SEP
										)) {
#ifdef L2_CACHE_MODELING
											if (doingPeakFinding && !reportedNextTrial) {
												reportedNextTrial = true;
												std::cerr << (do_cross_correlation ? "GASDF, do_cc" : "GASDF, not do_cc") << ". Next trial N:" << N << std::endl;
												std::cout << (do_cross_correlation ? "GASDF, do_cc" : "GASDF, not do_cc") << ". Next trial N:" << N << std::endl;
											}
#endif
										}
									} // end loop iimage_tile
								} // end loop iover_rot_tile
							} // end loop ipsi_tile
#ifdef L2_CACHE_MODELING
							if (!doingPeakFinding) continue;

							L2CacheModel::Interval l2CacheModelInterval;
							l2CacheModelInterval.begin();

							buffer.flush();
							int significantCount = buffer.significantCount();
							buffer.zeroSignificantCount();

							l2CacheModelInterval.end();
							auto l2HitRate = l2CacheModelInterval.hitRate();
							if (l2HitRate == 0) {
								significantCount = 0;	// This happens when none of the buffered tasks gets assigned to thread 0!
							}

							// l2CacheModelInterval.showL2CacheIntervalLocked(&life);

							if (resultIndex != PeakFinder::noAssignment) {
								if (significantCount == 0) {
									peakFinder.assignmentMustBeRedone(resultIndex);
								} else {
									if (l2HitRate < 0.1) std::cout << " significantCount:" << significantCount << " but surprisingly l2HitRate:" << l2HitRate << std::endl;
									peakFinder.reportResult(resultIndex, l2HitRate);
								}
							}

							if (significantCount > 0)
#pragma omp critical
							{
								if (useRecommendedTiling) {
									recommendedTilingHitRate = std::min(recommendedTilingHitRate, l2HitRate);
								}

								auto msg = [&](std::ostream&os) {
									os << "GASDF l2HitRate:" << std::setprecision(3) << l2HitRate << " tid:" << omp_get_thread_num() << " ri:" << resultIndex
										<< " n_tile:" << assignment.N_tile
										<< " ipsi_tile:" << assignment.ipsi_tile
										<< " iover_rot_tile:" << assignment.iover_rot_tile
										<< " iimage_tile:" << assignment.iimage_tile
										<< " iimage_sub_tile:" << assignment.iimage_sub_tile
										<< " itrans_sub_tile:" << assignment.itrans_sub_tile
										<< ", best yet:" << peakFinder.bestResult()
										<< "  default:" << recommendedTilingHitRate
										<< std::endl;
								};
								msg(std::cerr);
								msg(std::cout);
							}
#endif
						} // end loop idir
					} // end loop iclass
				} // end loop n_tile

#ifdef L2_CACHE_MODELING
				if (!doingPeakFinding) {
					auto showHistogram = [&](std::ostream& os) {
						os << "significantCountHistogram" << std::endl;
						for (size_t i = 0; i < significantCountHistogram.size(); i++) {
							auto count = significantCountHistogram[i];
							if (count == 0) continue;
							os << " " << i << ":" << count << std::endl;
						}
						os << "" << std::endl;
					};
					showHistogram(std::cerr);
					showHistogram(std::cout);
				}
				if (!doingPeakFinding) break;
#endif

				useRecommendedTiling = false;
			}
#ifdef L2_CACHE_MODELING
			while (!peakFinder.hasFoundPeak());

			if (doingPeakFinding) {
				auto msg = [&](std::ostream& os) {
					os << "GASDF PeakFinder::bestResult:" << peakFinder.bestResult();
					auto & bestCoords = peakFinder.bestCoordinates();
					const char* sep = " found at ";
					for (auto c : bestCoords) {
						os << sep << c;
						sep = ",";
					}
					os << " in a volume " << peakFinder.volume() << " searched with " << peakFinder.numberAssignmentsReported() << " trials" << std::endl;
				};
				msg(std::cout);
			}
#endif
		}
	}

    // TODO,YONGBEI
	bool GetAllSquaredDifferencesFineMain_Task::performWillDoAtLeastOne() {
		for (int ipsi = ipsi_tile_start; ipsi < ipsi_tile_end; ipsi++) {
			//int iClassOrient = iclass*exp_nr_dir + idir;
			//iClassOrient = iClassOrient*exp_nr_psi + ipsi;
			if (!exp_Mcoarse_Rot_significant.isRotSignificant(iclass, idir, ipsi))
				continue;
            auto exp_Mcoarse_significant_rptr = exp_Mcoarse_Rot_significant.isMcoarseSignificantRptrAll(iclass, idir, ipsi);
			for (int iimage = iimage_tile_start; iimage < iimage_tile_end; iimage++) {
				for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
					if (exp_Mcoarse_significant_rptr[iimage*exp_nr_trans + itrans]) {
						return true;
					}
				}
			}
		}
		return false;
	}
	int GetAllSquaredDifferencesFineMain_Task::perform(
		ProjectedClassCache& projectedClassCache, 
		SignificantCountHistogram & significantCountHistogram) {

		MAJOR_SCOPE2(GetAllSquaredDifferencesFineMain_Task::perform, life)

		assert(!(mlModel.pdf_class.rptrAll()[iclass] <= 0.));
        if (!do_local_searching) {
            assert(!(mlModel.pdf_direction[iclass].rptrAll()[idir] <= 0.));
        }
		L2CacheModel::Interval l2CacheModelInterval;
		l2CacheModelInterval.begin();

		// thread data
		const int tid = omp_get_thread_num();
		const int thread_exp_Mweight_sub_size = const_max_iimage_tile*exp_nr_trans*exp_nr_over_trans*exp_nr_over_rot;
		const auto thread_exp_Mweight_sub_tid = thread_exp_Mweight_sub[tid].wptr(thread_exp_Mweight_sub_size);

#ifdef L2_CACHE_MODELING
		auto showConfig = [&](std::ostream& os) {
			static int prev_classes     = 0;
			static int prev_ipsi_tile   = 0;
			static int prev_iimage_tile = 0;

			int classes     = exp_iclass_max + 1 - exp_iclass_min;
			int ipsi_tile   = ipsi_tile_end      - ipsi_tile_start;
			int iimage_tile = iimage_tile_end    - iimage_tile_start;

			if (	prev_classes == classes
				&&	prev_ipsi_tile == ipsi_tile
				&&	prev_iimage_tile == iimage_tile
				) return;

			prev_classes     = classes;
			prev_ipsi_tile   = ipsi_tile;
			prev_iimage_tile = iimage_tile;

			os << "getAllSquaredDifferencesFineMain - tiling configuration"
#define P(N) << " " << #N << ":" << N
				P(N) P(n_end) << std::endl
				P(exp_iclass_min) P(exp_iclass_max) << std::endl
				P(exp_nr_psi) P(ipsi_tile_start) P(ipsi_tile_end) << std::endl
				P(exp_nr_images) P(iimage_tile_start) P(iimage_tile_end) << std::endl
				<< std::endl;
		};

		if (tid == 0) 
		#pragma omp critical
		{
			showConfig(std::cerr);
		}
#endif

		int significantCount = 0;
		for (int ipsi = ipsi_tile_start; ipsi < ipsi_tile_end; ipsi++)
		{
			// exp_iclass loop does not always go from 0 to nr_classes!
			if (!exp_Mcoarse_Rot_significant.isRotSignificant(iclass, idir, ipsi)) continue;

            auto exp_Mcoarse_significant_rptr = exp_Mcoarse_Rot_significant.isMcoarseSignificantRptrAll(iclass, idir, ipsi);

			auto cacheProbeCount_begin = l2CacheModelInterval.cacheProbeCount();

			int partialSignificantCount = 0;
			if (do_cross_correlation) {
				partialSignificantCount =
					getAllSquaredDifferencesFineCC(
						projectedClassCache,
						thread_exp_Mweight_sub_tid,
						exp_Mcoarse_significant_rptr,
						tid, iclass, idir, ipsi, n_start, n_end,
						iover_rot_tile_start, iover_rot_tile_end,
						iimage_tile_start, iimage_tile_end);
			} else {
				partialSignificantCount =
					getAllSquaredDifferencesFineML(
						projectedClassCache,
						thread_exp_Mweight_sub_tid,
						exp_Mcoarse_significant_rptr,
						tid, iclass, idir, ipsi, n_start, n_end,
						iover_rot_tile_start, iover_rot_tile_end,
						iimage_tile_start, iimage_tile_end,
						iimage_sub_tile, itrans_sub_tile);
			}
			significantCount += partialSignificantCount;

#ifdef L2_CACHE_MODELING
			auto cacheProbeCount_end = l2CacheModelInterval.cacheProbeCount();
			if ((partialSignificantCount > 0) && (cacheProbeCount_begin == cacheProbeCount_end)) {
				std::cerr << "cacheProbeCount unchanged during "
					<< (do_cross_correlation ? "getAllSquaredDifferencesFineCC" : "getAllSquaredDifferencesFineML")
					<< std::endl;
			}
#endif
            // If do local search
            const std::vector<FDOUBLE>* pdf_orientation_iimage;
            if (do_local_searching) {
                pdf_orientation_iimage = samplingGrid.pdfOrientationForNonzeroOrient(idir, ipsi);
            }
			for (int iimage = iimage_tile_start; iimage < iimage_tile_end; iimage++)
			{
                // If do local search
                if (do_local_searching) {
                    if((*pdf_orientation_iimage)[iimage] <= 0) continue;
                }
#ifdef EXP_MWEIGHT_NEW
				auto  exp_Mweight_fine_sub_ptr = exp_Mweight_fine.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi);
				if (!exp_Mweight_fine_sub_ptr) continue;
				auto& exp_Mweight_fine_sub = *exp_Mweight_fine_sub_ptr;
#else
				auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(iimage, iclass, idir, ipsi);
#endif
				auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid +
					(iimage - iimage_tile_start) * exp_nr_trans *
					exp_nr_over_trans * exp_nr_over_rot;
				auto thread_exp_min_diff2_tid = thread_exp_min_diff2[tid].wptr(exp_nr_images);
				// Search for [[which values come from the thread_exp_Mweight_sub_sub_tid]]
				// to understand how this works
				//
				for (int i = 0; i < exp_Mweight_fine_sub.size(); i++)
				{
#ifdef EXP_MWEIGHT_NEW
					int itrans, iover_rot, iover_trans;
					exp_Mweight_fine_sub.key(i).decode(iover_rot, itrans, iover_trans);

					auto diff2 = exp_Mweight_fine_sub.value(i);
					int ihidden_tid = itrans*exp_nr_over_trans + iover_trans;// difference with i??
					ihidden_tid = ihidden_tid*exp_nr_over_rot + iover_rot;
#else
					auto ihidden_tid = exp_Mweight_fine_sub[i].first;
					auto& diff2 = exp_Mweight_fine_sub[i].second;
#endif
					//
					assert(thread_exp_Mweight_sub_sub_tid[ihidden_tid] != 0);
					const auto diff2_tid = thread_exp_Mweight_sub_sub_tid[ihidden_tid];
					diff2 += diff2_tid;
					//
					if (!do_cross_correlation && n_end == N)
					{
						// Calculate the actual squared difference term of the Gaussian probability function
						// If current_size < ori_size diff2 is initialised to the sum of
						// all |Xij|2 terms that lie between current_size and ori_size
						// Factor two because of factor 2 in division below, NOT because of 2-dimensionality of the complex plane!
						diff2 = (diff2 + exp_highres_Xi2_imgs.rptrAll()[iimage]) / 2;
					}
					if (n_end == N) {
						if (diff2 < thread_exp_min_diff2_tid[iimage])
							thread_exp_min_diff2_tid[iimage] = diff2;
					}
#ifdef EXP_MWEIGHT_NEW
					exp_Mweight_fine_sub.setValue(i, diff2);
#endif
				}
			}// end loop iimage
		}// end loop ipsi

#ifdef L2_CACHE_MODELING
		if (0)
		#pragma omp critical
		{
			significantCountHistogram[std::min(size_t(significantCount), significantCountHistogram.size() - 1)]++;

			if (significantCount > 0) {
				auto showSignificant = [&](std::ostream& os) {
					os << "GetAllSquaredDifferencesFineMain_Task::perform"
						<< " significantCount:" << significantCount 
						<< " tid:"				<< tid
						<< " accessesRecorded:" << l2CacheModelInterval.accessesRecordedCount()
						<< " cacheProbes:"      << l2CacheModelInterval.cacheProbeCount()
						<< std::endl;
				};
				showSignificant(std::cerr);
				showSignificant(std::cout);

				bool goodHitRate = l2CacheModelInterval.showL2CacheIntervalUnlocked(&life);

				if (!goodHitRate || (significantCount == 1)) {
					showConfig(std::cerr);
					showConfig(std::cout);
				}
			}
		}
#endif
		return significantCount;
	}

	//
    void vtuneAnalysis(){
#define VTUNE_ANA 1

#ifdef VTUNE_ANA
            //
            // printExp_Mcoarse_significant(write_path+write_fn+"_it"+num2str(iter)+"_significant");
            //
            TileCfg tilecfg;
            tilecfg.choiceTileForKNL(1);
            for (int i = 0; i < 1; i++)
            {
                double time_start = dtime();
#ifdef EXP_MWEIGHT_NEW
                exp_Mweight_fine.clearAllImages();
#else
                exp_Mweight_fine.reset();
#endif
                getAllSquaredDifferencesFineMain(tilecfg);
                double time_end = dtime();
                std::cout	<<"count : "<<i<<"time : "<<(time_end-time_start)<<std::endl;
            }
#else
            // detect best tile
#endif
            EXIT_ABNORMALLY;
    }
    //
    void getAllSquaredDifferences(bool _do_coarse_search)
    {
		MAJOR_SCOPE(getAllSquaredDifferences)

        exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
        do_cross_correlation = (iter == 1 && do_firstiter_cc) || do_always_cc;
        do_coarse_search = _do_coarse_search;
        TileCfg tilecfg;
        //
        thread_exp_min_diff2.fill_with_first_touch((std::numeric_limits<double>::max)());
#ifdef DATA_STREAM
        // turn off data stream inside openmp for loop
        // if(data_stream_node == node) global_data_stream.turnOff();
        global_data_stream_before();
#endif
        if (do_coarse_search)
        {
#ifdef EXP_MWEIGHT_NEW
            exp_Mweight_coarse.clearAllImages();
#else
            exp_Mweight_coarse.reset();
#endif
            // NOTE : transform fft array to shell increased format to reduce computing ...
            // only when do maximumu-likehood...
            if (!do_cross_correlation) {transformFourierArrayToShellIncreasedCoarse();}
            //
            getAllSquaredDifferencesCoarse(tilecfg);
        }
        else
        {
#ifdef EXP_MWEIGHT_NEW
            exp_Mweight_fine.clearAllImages();
#else
            exp_Mweight_fine.reset();
#endif
            //
            if (!do_cross_correlation) { transformFourierArrayToShellIncreasedFine();}
            //
            if (false) {vtuneAnalysis();}
            //
            getAllSquaredDifferencesFineMain(tilecfg);
            // updateModel will use these
            if (do_cross_correlation) { transformFourierArrayToShellIncreasedFine();}
        }
		//
        for (int iimage = 0; iimage < exp_nr_images; iimage++) {
            exp_min_diff2.wptrAll()[iimage] = (std::numeric_limits<double>::max)();
            for (int thread = 0; thread < maxthreads; thread++) {
                auto thread_exp_min_diff2_tid = thread_exp_min_diff2[thread].wptr(exp_nr_images);
                if (thread_exp_min_diff2_tid[iimage] < exp_min_diff2.rptrAll()[iimage]) {
                    exp_min_diff2.wptrAll()[iimage] = thread_exp_min_diff2_tid[iimage];
                }
            }
        }
#ifdef DATA_STREAM
        //
        if(data_stream_node == node) global_data_stream.turnOn();
        //
        global_data_stream_after();
#endif
    }

//
}

void getAllSquaredDifferences(bool do_coarse_search){
    CoarseAndFineSearch::getAllSquaredDifferences(do_coarse_search);
}

void convertSquaredDifferencesToWeights(
#ifdef EXP_MWEIGHT_NEW
    Exp_Mweight_new& exp_Mweight
#else
    Exp_Mweight_old& exp_Mweight
#endif
    )
{
    thread_exp_max_weight.fill_with_first_touch((std::numeric_limits<double>::min)());
    thread_exp_sum_weight.fill_with_first_touch(0.);
#pragma omp parallel for
    for (int thread = 0; thread < maxthreads; thread++)
    {
        auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[thread].data();
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
            thread_exp_max_weight_index_tid[iimage] = {0,0,0,0,0,0,0};
    }
#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    // if(data_stream_node == node) global_data_stream.turnOff();
#endif
    
    DISABLE_DENORMALS
    
    // do cross-correlation
    bool do_cross_correlation = (iter == 1 && do_firstiter_cc) || do_always_cc;
    if (do_cross_correlation)
    {
		CHECKER_PARALLEL_COUNTED_FOR1(int, iimage, 0, exp_nr_images)
//      #pragma omp parallel for
//      for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            int tid = omp_get_thread_num();
            auto thread_exp_max_weight_tid = thread_exp_max_weight[tid].wptr(exp_nr_images);
            auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[tid].data();

            double const mymin_diff2_init = (std::numeric_limits<double>::max)();
            double       mymin_diff2      = mymin_diff2_init;

#ifdef	EXP_MWEIGHT_NEW
			PackedOverRotTransOverTrans mymin_key;
            // Binarize the squared differences array to skip marginalisation
            // Find the smallest element in this row of exp_Mweight
            auto set_cc = [&](int iimage,int iclass,int idir,int ipsi, PackedOverRotTransOverTrans key, double cc){
                // ignore non-determined cc
                if (cc == -999.123456789)
                    return;

                // just search for the maximum
                if (cc < mymin_diff2)
                {
                    mymin_diff2 = cc;
                    mymin_key = key;
                    thread_exp_max_weight_tid[iimage] = 1.;
                    int iover_rot,itrans,iover_trans;
                    key.decode(iover_rot, itrans, iover_trans);
                    thread_exp_max_weight_index_tid[iimage] = {iimage,iclass,idir,ipsi,iover_rot,itrans,iover_trans};
                }
            };
            for (int iclass = 0; iclass < nr_classes; iclass++) {
                for (int idir = 0; idir < exp_nr_dir; idir++) {
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
						auto someWeights_ptr = exp_Mweight.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi);
						if (!someWeights_ptr) continue;
						auto& someWeights = *someWeights_ptr;
						for (size_t i = 0; i < someWeights.size(); i++) {
                            set_cc(iimage,iclass,idir,ipsi, someWeights.key(i),someWeights.value(i));
                        }
                    }
                }
            }

#else
            int mymin_ihidden;
            // Binarize the squared differences array to skip marginalisation
            // Find the smallest element in this row of exp_Mweight
            auto set_cc = [&](int iimage,int iclass,int idir,int ipsi,int ihidden,double cc){
                // ignore non-determined cc
                if (cc == -999.123456789)
                    return;

                // just search for the maximum
                if (cc < mymin_diff2)
                {
                    mymin_diff2 = cc;
                    mymin_ihidden = ihidden;
                    thread_exp_max_weight_tid[iimage] = 1.;
                    int iover_rot,itrans,iover_trans;
                    decodeRotOverTrans(ihidden, itrans, iover_trans, iover_rot);
                    thread_exp_max_weight_index_tid[iimage] = {iimage,iclass,idir,ipsi,iover_rot,itrans,iover_trans};
                }
            };
            for (int iclass = 0; iclass < nr_classes; iclass++) {
                for (int idir = 0; idir < exp_nr_dir; idir++) {
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
                        for (auto const&ihidden_weight : exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi)){
                            set_cc(iimage,iclass,idir,ipsi,ihidden_weight.first,ihidden_weight.second);
                        }
                    }
                }
            }
#endif
            //
            auto thread_exp_sum_weight_tid = thread_exp_sum_weight[tid].wptr(exp_nr_images);
            thread_exp_sum_weight_tid[iimage] += thread_exp_max_weight_tid[iimage];
            //
            ERROR_CHECK(mymin_diff2==mymin_diff2_init, "none significant point for image : "+std::to_string((long long)iimage));
            int mymin_iclass = thread_exp_max_weight_index_tid[iimage].iclass;
            int mymin_idir   = thread_exp_max_weight_index_tid[iimage].idir;
            int mymin_ipsi   = thread_exp_max_weight_index_tid[iimage].ipsi;
            // Set all except for the best hidden variable to zero and the smallest element to 1
#ifdef	EXP_MWEIGHT_NEW
            exp_Mweight.clear(iimage);
            auto& exp_Mweight_sub = exp_Mweight.mweightsForSomeSpinsAndSlides(iimage, mymin_iclass, mymin_idir, mymin_ipsi);
            exp_Mweight_sub.insert(mymin_key, 1.0);
#else
            exp_Mweight.clear_image(iimage);
            auto& exp_Mweight_sub = exp_Mweight.wptr_sparse(iimage, mymin_iclass, mymin_idir, mymin_ipsi);
            exp_Mweight_sub.push_back( std::pair<int, double>(mymin_ihidden,1) );
#endif

#ifdef DATA_STREAM
            global_data_stream.foutDouble(mymin_diff2, "convertSquaredDifferencesToWeights()_mymindiff2", __FILE__, __LINE__);
#endif
        } CHECKER_PARALLEL_FOR_END
    } // end do cross-correlation
    else // do maximum-likelihood
    {
		// BEVIN
		//	QUESTION WHY 3 NOT 4?
		//
		CHECKER_PARALLEL_COUNTED_FOR3(
			int, iimage, 0,				 exp_nr_images,
			int, iclass, exp_iclass_min, exp_iclass_max+1,
			int, idir,   0,              exp_nr_dir) {
//		#pragma omp parallel for collapse(3) schedule(dynamic)
//      for (int iimage = 0; iimage < exp_nr_images;iimage++)
//      {
//          for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
//          {
//              for (int idir = 0; idir < exp_nr_dir; idir++)
//              {
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                    {
#ifdef EXP_MWEIGHT_NEW
						auto const exp_Mweight_sub_ptr = exp_Mweight.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi);
						if (!exp_Mweight_sub_ptr) continue;
#endif
						int tid = omp_get_thread_num();
                        auto thread_exp_sum_weight_tid = thread_exp_sum_weight[tid].wptr(exp_nr_images);
                        auto thread_exp_max_weight_tid = thread_exp_max_weight[tid].wptr(exp_nr_images);
                        auto thread_exp_max_weight_index_tid =  thread_exp_max_weight_index[tid].data();
                        double pdf_orientation;
                        if (do_local_searching) {
                            pdf_orientation = (*samplingGrid.pdfOrientationForNonzeroOrient(idir, ipsi))[iimage];
                        }
                        else{
                            pdf_orientation = mlModel.pdf_direction[iclass][idir];
                        }
                        

#ifdef EXP_MWEIGHT_NEW
                        auto& exp_Mweight_sub = *exp_Mweight_sub_ptr;
#else
                        auto& exp_Mweight_sub = exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi);
#endif
                        FDOUBLE offset_x,offset_y;
                        double pdf_offset = (std::numeric_limits<double>::max)();
                        int pre_itrans = (std::numeric_limits<int>::max)();

                        // for itrans,iover_rot,iover_trans
						for (size_t i = 0; i < exp_Mweight_sub.size(); i++)
                        {
                            int itrans,iover_rot,iover_trans;
#ifdef EXP_MWEIGHT_NEW
							exp_Mweight_sub.key(i, iover_rot, itrans, iover_trans);
#else
                            decodeRotOverTrans(exp_Mweight_sub[i].first, itrans, iover_trans, iover_rot);
#endif
                            if (itrans != pre_itrans) // recompute pdf_offset
                            {
                                pre_itrans = itrans;
                                // NOTE : To speed things up, only calculate pdf_offset at the coarse sampling.
                                // That should not matter much, and that way one does not need to calculate all the OversampledTranslations
                                sampler3d.getTranslation(itrans, offset_x, offset_y);
                                offset_x += exp_old_offsetx.rptrAll()[iimage];
                                offset_y += exp_old_offsety.rptrAll()[iimage];
                                // Convert offsets back to Angstroms to calculate PDF!
                                pdf_offset = mlModel.calculatePdfOffset(offset_x,offset_y,
                                                                        mlModel.prior_offsetx_class.rptrAll()[iclass],
                                                                        mlModel.prior_offsety_class.rptrAll()[iclass]);
                            }

                            double weight = pdf_orientation * pdf_offset;

                            // Sacrifice some performance to get same result as relion
#ifdef 	EXP_MWEIGHT_NEW
                            double diff2 = exp_Mweight_sub.value(i) - exp_min_diff2.rptrAll()[iimage];
#else
                            double diff2 = exp_Mweight_sub[i].second - exp_min_diff2.rptrAll()[iimage];
#endif

                            // next line because of numerical precision of exp-function
                            if (diff2 > 700.) weight = 0.;
                            // TODO: use tabulated exp function?
                            else weight *= exp(-diff2);

                            // Store the weight
#ifdef EXP_MWEIGHT_NEW
							exp_Mweight_sub.setValue(i, weight);
#else
                            exp_Mweight_sub[i].second = weight;
#endif

                            // Use the weight
                            thread_exp_sum_weight_tid[iimage] += weight;
                            if (weight > thread_exp_max_weight_tid[iimage])
                            {
                                thread_exp_max_weight_tid[iimage] = weight;
                                thread_exp_max_weight_index_tid[iimage] = {iimage,iclass,idir,ipsi,iover_rot,itrans,iover_trans};
                            }
#ifdef DATA_STREAM
                            global_data_stream.foutInt(iimage, "convertSquaredDifferencesToWeights()_iimage", __FILE__, __LINE__);
                            global_data_stream.foutInt(iclass, "convertSquaredDifferencesToWeights()_iclass", __FILE__, __LINE__);
                            global_data_stream.foutInt(idir, "convertSquaredDifferencesToWeights()_idir", __FILE__, __LINE__);
                            global_data_stream.foutInt(ipsi, "convertSquaredDifferencesToWeights()_ipsi", __FILE__, __LINE__);
                            global_data_stream.foutInt(itrans, "convertSquaredDifferencesToWeights()_itrans", __FILE__, __LINE__);
                            global_data_stream.foutInt(iover_rot, "convertSquaredDifferencesToWeights()_iover_rot", __FILE__, __LINE__);
                            global_data_stream.foutInt(iover_trans, "convertSquaredDifferencesToWeights()_iover_trans", __FILE__, __LINE__);
                            global_data_stream.foutDouble(offset_x, "convertSquaredDifferencesToWeights()_offset_x", __FILE__, __LINE__);
                            global_data_stream.foutDouble(offset_y, "convertSquaredDifferencesToWeights()_offset_y", __FILE__, __LINE__);
                            global_data_stream.foutDouble(pdf_offset, "convertSquaredDifferencesToWeights()_pdf_offset", __FILE__, __LINE__);
                            global_data_stream.foutDouble(pdf_orientation, "convertSquaredDifferencesToWeights()_pdf_orientation", __FILE__, __LINE__);
                            global_data_stream.foutDouble(diff2, "convertSquaredDifferencesToWeights()_diff2", __FILE__, __LINE__);
                            global_data_stream.foutDouble(weight, "convertSquaredDifferencesToWeights()_weight", __FILE__, __LINE__);
                            global_data_stream.foutDouble(mlModel.prior_offsetx_class.wptrAll()[iclass], "convertSquaredDifferencesToWeights()_prior_offsetx_class[iclass]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(mlModel.prior_offsety_class.wptrAll()[iclass], "convertSquaredDifferencesToWeights()_prior_offsety_class[iclass]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(mlModel.sigma2_offset, "convertSquaredDifferencesToWeights()_sigma2_offset", __FILE__, __LINE__);
                            global_data_stream.foutDouble(exp_old_offsetx.wptrAll()[iimage], "convertSquaredDifferencesToWeights()_exp_old_offsetx[iimage]", __FILE__, __LINE__);
                            global_data_stream.foutDouble(exp_old_offsety.wptrAll()[iimage], "convertSquaredDifferencesToWeights()_exp_old_offsety[iimage]", __FILE__, __LINE__);
                            global_data_stream.check();global_data_stream.flush();
#endif
                        }// end loop ihidden(itrans,iover_rot,iover_trans)
                    }// end loop ipsi
//                }					// end loop idir
//            }						// end loop iclass
        } CHECKER_PARALLEL_FOR_END	// end loop iimage
    } // end do maximum-likelihood
#ifdef DATA_STREAM
    if(data_stream_node == node) global_data_stream.turnOn();
#endif
    //
    ENABLE_DENORMALS
    //
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        double exp_max_weight = 0;
        GridIndex exp_max_weight_index;
        exp_sum_weight.wptrAll()[iimage] = 0.;
        for (int thread = 0; thread < maxthreads; thread++)
        {
            auto thread_exp_sum_weight_tid = thread_exp_sum_weight[thread].wptr(exp_nr_images);
            exp_sum_weight.wptrAll()[iimage] += thread_exp_sum_weight_tid[iimage];
            auto thread_exp_max_weight_tid = thread_exp_max_weight[thread].wptr(exp_nr_images);
            auto thread_exp_max_weight_index_tid = thread_exp_max_weight_index[thread].data();
            if(thread_exp_max_weight_tid[iimage] > exp_max_weight){
                exp_max_weight = thread_exp_max_weight_tid[iimage];
                exp_max_weight_index = thread_exp_max_weight_index_tid[iimage];
            }
        }
        //
        int iimage = exp_max_weight_index.iimage;
        int iclass = exp_max_weight_index.iclass;
        int idir = exp_max_weight_index.idir;
        int ipsi = exp_max_weight_index.ipsi;
        int iover_rot = exp_max_weight_index.iover_rot;
        int itrans = exp_max_weight_index.itrans;
        int iover_trans = exp_max_weight_index.iover_trans;
        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                    exp_over_rot, exp_over_tilt, exp_over_psi);
        //
        exp_metadata[iimage].PMAX = exp_max_weight/exp_sum_weight.wptrAll()[iimage];
        exp_metadata[iimage].CLASS = iclass+1;
        exp_metadata[iimage].ROT = exp_over_rot[iover_rot];
        exp_metadata[iimage].TILT = exp_over_tilt[iover_rot];
        exp_metadata[iimage].PSI = exp_over_psi[iover_rot];
        exp_metadata[iimage].XOFF = exp_old_offsetx.wptrAll()[iimage] + samplingGrid.exp_over_trans_x[itrans][iover_trans];
        exp_metadata[iimage].YOFF = exp_old_offsety.wptrAll()[iimage] + samplingGrid.exp_over_trans_y[itrans][iover_trans];
    }
#ifdef DATA_STREAM
    global_data_stream.foutDouble(exp_sum_weight.wptrAll()[0], "convertSquaredDifferencesToWeights()_exp_sum_weight[0]", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}


class FindAllSignificantPoints : public FindBestPoints {
public:
	double	frac_weight;
	double	significant_weight_limit;
	int		significant_weight_count;

	static void estimateWork(
#ifdef EXP_MWEIGHT_NEW
		Exp_Mweight_new& exp_Mweight,
#else
		Exp_Mweight_old& exp_Mweight,
#endif
		int const	iimage,
		int &		weight_count,
		float&		work
	) {
		weight_count = calculate_weight_count(exp_Mweight,iimage);
		work = weight_count > 0 ? float(weight_count)*std::log(weight_count) : 0.0;
	}

	static int calculate_weight_count(
#ifdef EXP_MWEIGHT_NEW
		Exp_Mweight_new& exp_Mweight,
#else
		Exp_Mweight_old& exp_Mweight,
#endif
		int const			 iimage) {

		int weight_count = 0;

		// get weight count
		// Only select non-zero probabilities to speed up sorting
#ifdef EXP_MWEIGHT_NEW
		weight_count =
			exp_Mweight.sum_sizes_non_zero(
				iimage,
				0, nr_classes,
				0, exp_nr_dir,
				0, exp_nr_psi);
#else
		for (int iclass = 0; iclass < nr_classes; iclass++)
			for (int idir = 0; idir < exp_nr_dir; idir++)
				for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
                    for (auto const &ihidden_weight : exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi)){
                        auto value = ihidden_weight.second;
                        if(value > 0) weight_count++;
                    }
				}
#endif

		if (weight_count == 0) {
			std::cerr << "sorted_weight_count == 0!" << std::endl;
			ERROR_REPORT("sorted_weight_count <= 0");
		}

		return weight_count;
	}

	void find(
#ifdef EXP_MWEIGHT_NEW
		Exp_Mweight_new& exp_Mweight,
#else
		Exp_Mweight_old& exp_Mweight,
#endif
		int const			 iimage,
		int const			 weight_count) {

		auto const sum_weight    = exp_sum_weight.rptrAll()[iimage];
		auto const target_weight = adaptive_fraction * exp_sum_weight.rptrAll()[iimage];
		setWeightCount(weight_count, target_weight, sum_weight);

		// insert the weights
		for (int iclass = 0; iclass < nr_classes; iclass++)
			for (int idir = 0; idir < exp_nr_dir; idir++)
				for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
                    // Only select non-zero probabilities to speed up sorting
#ifdef EXP_MWEIGHT_NEW
					auto ptr = exp_Mweight.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi);
					if (!ptr) continue;
					for (size_t i = 0; i < ptr->size(); i++) {
                        auto value = ptr->value(i);
						if(value > 0) insert(value);
					}
#else
                    for (auto const &ihidden_weight : exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi)){
                        auto value = ihidden_weight.second;
						if(value > 0) insert(value);
                    }
#endif
				}

		// I suspect a load balancing problem because of the serialization hotspot on the sort

		prepareToPop();

		frac_weight = 0.;
		significant_weight_count = 0;

		for (int i = 0; i < weight_count; i++)
		{
			significant_weight_count++;
			significant_weight_limit = pop();
			frac_weight += significant_weight_limit;
			if (frac_weight > target_weight)
				break;
		}

		if (0)
		#pragma omp critical
		{
			std::cerr
				<< "tid:" << omp_get_thread_num()
				<< "FindAllSignificantPoints()"
				<< " weight_count: "			<< weight_count
				<< " significant_weight_count:" << significant_weight_count
				<< std::endl;
		}

		if (exp_ipass==0 && significant_weight_count == 0)
#pragma omp critical
		{
			std::cerr << " ipart= " << iimage << " adaptive_fraction= " << adaptive_fraction << std::endl;
			std::cerr << " frac_weight= " << frac_weight << std::endl;
			std::cerr << " exp_sum_weight[iimage] = " << exp_sum_weight.rptrAll()[iimage] << std::endl;
			std::cerr << " sorted_weight_count = " << weight_count << std::endl;
			ERROR_REPORT("significant_weight_count == 0");
		}

		if (exp_ipass==0)
		{
			exp_metadata[iimage].NR_SIGN = (double)significant_weight_count;

			// Keep track of which coarse samplings were significant were significant for this particle
			// in ipass = 0,exp_Mcoarse_xsize eqaul to exp_Mweight_xsize
			// TODO better data structure for exp_Mcoarse_significant....
			// transform exp_Mweight coordinate to exp_Mcoarse_significant
#ifdef 	EXP_MWEIGHT_NEW
			//
			for (int iclass = 0; iclass < nr_classes; iclass++) {
				for (int idir = 0; idir < exp_nr_dir; idir++) {
					for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
						//
						auto someMweights_ptr = exp_Mweight.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi);
						if (!someMweights_ptr) continue;
						auto & someMweights = *someMweights_ptr;
                        std::vector<int> iimage_trans;
						for (auto i = 0; i < someMweights.size(); i++)
						{
                            if (someMweights.value(i) >= significant_weight_limit){
                                int iover_rot,itrans,iover_trans;
                                someMweights.key(i).decode(iover_rot, itrans, iover_trans);
                                iimage_trans.push_back(iimage*exp_nr_trans+itrans);
                            }
						}
                        if (iimage_trans.size() > 0) {
                            // false sharing???
                            exp_Mcoarse_Rot_significant.assignMcoarseSignificant(iclass, idir, ipsi, iimage_trans, true);
                        }
					}
				}
			}
#else
			//
			for (int iclass = 0; iclass < nr_classes; iclass++) {
				for (int idir = 0; idir < exp_nr_dir; idir++) {
					for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
                        if (exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi).size() == 0) continue;
                        std::vector<int> iimage_trans;
						for (auto const & ihidden_weight : exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi))
						{
                            if (ihidden_weight.second >= significant_weight_limit){
                                int iover_rot, itrans, iover_trans;
                                decodeRotOverTrans(ihidden_weight.first, itrans, iover_trans, iover_rot);
                                iimage_trans.push_back(iimage*exp_nr_trans+itrans);
                            }
						}
                        if (iimage_trans.size() > 0) {
                            // false sharing???
                            exp_Mcoarse_Rot_significant.assignMcoarseSignificant(iclass, idir, ipsi, iimage_trans, true);
                        }
					}
				}
			}
#endif
		}
	}
};


void findAllSignificantPoints(
#ifdef EXP_MWEIGHT_NEW
    Exp_Mweight_new& exp_Mweight
#else
    Exp_Mweight_old& exp_Mweight
#endif
	)
{
	static double reported_adaptive_fraction = -1.0;
	if (reported_adaptive_fraction != adaptive_fraction) {
		reported_adaptive_fraction = adaptive_fraction;
		NODE0ONLY std::cerr << __FILE__ << ":" << __LINE__ << " findAllSignificantPoints adaptive_fraction= " << adaptive_fraction << std::endl;
	}

#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    // if(data_stream_node == node) global_data_stream.turnOff();
#endif

	// Decide what has to be done
	//
//#define LOAD_BALANCE_FASP
#ifdef LOAD_BALANCE_FASP
	struct Task {
		// The work to be done
		int				weight_count;
		float			work;
		// How it was done
		Microseconds	elapsed_getting_weight_count;
		int				tid;
		Microseconds	elapsed_getting_points;
	};
	std::vector<Task> tasks(exp_nr_images);

	#pragma omp parallel for
	for (int iimage = 0; iimage < exp_nr_images; iimage++) {
		auto & task = tasks[iimage];
		auto start = timeInMicroseconds();
		FindAllSignificantPoints::estimateWork(exp_Mweight, iimage, task.weight_count, task.work);
		task.elapsed_getting_weight_count = timeInMicroseconds() - start;
	}
#endif

    //
    if(exp_ipass==0)
    {
        exp_Mcoarse_Rot_significant.resetMcoarseSignificant();
        if (exp_Mcoarse_Rot_significant.use_sparse_data/*do_local_searching*/) {
#ifdef 	EXP_MWEIGHT_NEW
        for (int iclass = 0; iclass < nr_classes; iclass++)
            for (int idir = 0; idir < exp_nr_dir; idir++)
                for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                    for (int iimage = 0; iimage < exp_nr_images; iimage++)
                    {
                        auto someMweights_ptr = exp_Mweight.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi);
                        if (someMweights_ptr) {
                            exp_Mcoarse_Rot_significant.allocMcoarseSignificant(iclass, idir, ipsi);
                            break;
                        }
                    }
#else
        //
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
            for (int iclass = 0; iclass < nr_classes; iclass++)
                for (int idir = 0; idir < exp_nr_dir; idir++)
                    for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                    {
                        if (exp_Mweight.wptr_sparse(iimage, iclass, idir, ipsi).size() > 0)
                            exp_Mcoarse_Rot_significant.allocMcoarseSignificant(iclass, idir, ipsi);
                    }
#endif
        }
    }
	// Now, for each image,find the exp_significant_weight that encompasses adaptive_fraction of exp_sum_weight
	std::vector<FindAllSignificantPoints> FindAllSignificantPoints_for_each_thread(omp_get_max_threads());
		// Reuse the vector rather than resizing for each image

	#pragma omp parallel for
	for (int iimage = 0; iimage < exp_nr_images; iimage++) {
#ifdef LOAD_BALANCE_FASP
		auto & task = tasks[iimage];
#endif
		auto tid = omp_get_thread_num();
		auto& points = FindAllSignificantPoints_for_each_thread[tid];

		auto start = timeInMicroseconds();

#ifdef LOAD_BALANCE_FASP
		points.find(exp_Mweight, iimage, task.weight_count);
		task.tid = tid;
		task.elapsed_getting_points = timeInMicroseconds() - start;
#else
		{
			int weight_count;
			float work;
			FindAllSignificantPoints::estimateWork(exp_Mweight, iimage, weight_count, work);
			points.find(exp_Mweight, iimage, weight_count);
		}
#endif

		exp_significant_weight.wptrAll()[iimage] = points.significant_weight_limit;

		if (
#ifdef DATA_STREAM
			true
#else
			false
#endif
			) {
			global_data_stream.foutDouble(points.get_weight_count(), "findAllSignificantPoints()_np", __FILE__, __LINE__);
			if (auto ptr = points.get_sorted_weights_ptr_if_avail()) {
				global_data_stream.foutDouble(ptr, points.get_weight_count(), "findAllSignificantPoints()_sorted_weight", __FILE__, __LINE__);
			}
			global_data_stream.foutDouble(points.significant_weight_count, "findAllSignificantPoints()_my_nr_significant_coarse_samples", __FILE__, __LINE__);
			global_data_stream.foutDouble(points.significant_weight_limit, "findAllSignificantPoints()_my_significant_weight", __FILE__, __LINE__);
			global_data_stream.foutDouble(points.frac_weight, "findAllSignificantPoints()_frac_weight", __FILE__, __LINE__);
			global_data_stream.foutDouble(adaptive_fraction, "findAllSignificantPoints()_adaptive_fraction", __FILE__, __LINE__);
			global_data_stream.check(); global_data_stream.flush();
		}
    }

	// Describe what was done in the loops
	//
#ifdef LOAD_BALANCE_FASP
	if (false) {
		std::cerr << std::endl;
		struct Done {
			Microseconds elapsed_getting_weight_count;
			Microseconds elapsed_getting_points;
			int weight_count;
			float work;
			Done() : elapsed_getting_weight_count(0.0), elapsed_getting_points(0.0), weight_count(0), work(0.0) {}
		};
		std::vector<Done> perTid(omp_get_max_threads());
		std::cerr << "Task stats for exp_nr_images:" << exp_nr_images << std::endl;
		std::cerr << "weight_count, work, elapsed_getting_weight_count, elapsed_getting_points" << std::endl;
		for (auto i = 0; i < tasks.size(); i++) {
			auto& task = tasks[i];
			auto& t = perTid[task.tid];
			t.elapsed_getting_weight_count	+= task.elapsed_getting_weight_count;
			t.elapsed_getting_points		+= task.elapsed_getting_points;
			t.weight_count					+= task.weight_count;
			t.work							+= task.work;
			if (i < 50) std::cerr << task.weight_count << "," << task.work << ", " << task.elapsed_getting_weight_count << "," << task.elapsed_getting_points << std::endl;
		}
		std::cerr << "Tid stats for exp_nr_images:" << exp_nr_images << std::endl;
		std::cerr << "weight_count, work, elapsed_getting_weight_count, elapsed_getting_points" << std::endl;
		for (auto i = 0; i < perTid.size(); i++) {
			auto& t = perTid[i];
			std::cerr << t.weight_count << ", " << t.work << ", " << t.elapsed_getting_weight_count << "," << t.elapsed_getting_points << std::endl;
		}
		std::cerr << "-------------" << std::endl;
	}
#endif

#ifdef DATA_STREAM
    if(data_stream_node == node) global_data_stream.turnOn();
#endif
}

void storeWeightedSums()
{
    int current_Fsize = current_size/2+1;
    int current_Fsize2 = current_size*(current_size/2+1);
    int ori_Fsize = ori_size/2+1;
#ifdef DATA_STREAM
    // turn off data stream inside openmp for loop
    if(data_stream_node == node) global_data_stream.turnOff();
#endif
    // Extend norm_correction and sigma2_noise estimation to higher resolutions
    // NOTE : some memory conflict for group-data,like wsum_sigma2_noise,wsum_signal_product_spectra...
    // so removing openmp
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        int igroup = exp_metadata[iimage].GROUP_NO-1;
        auto exp_power_imgs_iimage = exp_power_imgs[iimage].wptr(ori_Fsize);
        auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
        // If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
        for (int ires = current_Fsize; ires < ori_Fsize; ires++)
        {
            // memory conflict
             wsum_sigma2_noise_igroup[ires] += exp_power_imgs_iimage[ires];
            // Also extend the weighted sum of the norm_correction

            // to remove outside????
            exp_wsum_norm_correction.wptrAll()[iimage] += exp_power_imgs_iimage[ires];
        }

        // Store norm_correction
        // Multiply by old value because the old norm_correction term was already applied to the image
        if (do_norm_correction)
        {
            double old_norm_correction = exp_metadata[iimage].NORM / mlModel.avg_norm_correction;
            // Now set the new norm_correction in the relevant position of exp_metadata22
            // The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
            // The variance of the total image (on which one normalizes) is twice this value!
            exp_metadata[iimage].NORM = old_norm_correction * sqrt(exp_wsum_norm_correction.rptrAll()[iimage] * 2.);
            if (!(iter == 1 && do_firstiter_cc) && exp_metadata[iimage].NORM > 10.)
            {
#pragma omp critical
                std::cerr<<"Warning in storeWeightedSums(),please debug this function."<<std::endl;
            }
#pragma omp atomic
            mlModel.wsum_avg_norm_correction += old_norm_correction * sqrt(exp_wsum_norm_correction.rptrAll()[iimage] * 2.);
        }
        // Store weighted sums for scale_correction
        if (do_scale_correction)
        {
            // Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
            for (int n = 0; n < mlModel.ori_Fsize; n++) {
                exp_wsum_scale_correction_XA[iimage].wptrAll()[n] /= mlModel.scale_correction.rptrAll()[igroup];
                exp_wsum_scale_correction_AA[iimage].wptrAll()[n] /= mlModel.scale_correction.rptrAll()[igroup] * mlModel.scale_correction.rptrAll()[igroup];
            }
            auto wsum_signal_product_spectra_group = mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize);
            auto wsum_reference_power_spectra_group = mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize);
            for (int n = 0; n < mlModel.ori_Fsize; n++) {
                wsum_signal_product_spectra_group[n] += exp_wsum_scale_correction_XA[iimage].rptrAll()[n];
                wsum_reference_power_spectra_group[n] += exp_wsum_scale_correction_AA[iimage].rptrAll()[n];
            }
        }

        // Some analytics...
        // Calculate normalization constant for dLL
        // loop over all particles inside this ori_particle
        double logsigma2 = 0.;

        for (int n = 0; n < current_Fsize2; n++)
        {
            int ires = Mresol_fine.rptrAll()[n];
            // Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
            // Also exclude origin from logsigma2, as this will not be considered in the P-calculations
            // use previous step's sigma2_noise???
            auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].rptr(mlModel.ori_Fsize);
            if (ires > 0)
                logsigma2 += log( 2. * PI * sigma2_noise_igroup[ires]);
        }


        if (exp_sum_weight.rptrAll()[iimage]==0)
        {
            ERROR_REPORT("ERROR: exp_sum_weight[ipart]==0");
        }

        double dLL = 0.;
        if ((iter==1 && do_firstiter_cc) || do_always_cc)
            dLL = -exp_min_diff2.rptrAll()[iimage];
        else
            dLL = log(exp_sum_weight.rptrAll()[iimage]) - exp_min_diff2.rptrAll()[iimage] - logsigma2;

        // Also store dLL of each image in the output array
        exp_metadata[iimage].DLL = dLL;

#pragma omp atomic
        mlModel.wsum_LL += dLL;
#pragma omp atomic
        mlModel.wsum_ave_Pmax += exp_metadata[iimage].PMAX;
#ifdef DATA_STREAM
        global_data_stream.foutInt(iimage, "storeWeightedSums()_iimage", __FILE__, __LINE__);
        global_data_stream.foutDouble(exp_wsum_norm_correction.wptrAll()[iimage], "storeWeightedSums()_exp_wsum_norm_correction", __FILE__, __LINE__);
        global_data_stream.foutDouble(exp_metadata[iimage].NORM, "storeWeightedSums()_exp_metadata[iimage].NORM", __FILE__, __LINE__);
        global_data_stream.foutDouble(exp_metadata[iimage].DLL, "storeWeightedSums()_exp_metadata[iimage].DLL", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
    } // end loop iimage
#ifdef DATA_STREAM
    if(data_stream_node == node) global_data_stream.turnOn();
    global_data_stream.foutDouble(mlModel.wsum_LL, "storeWeightedSums()_wsum_LL", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_ave_Pmax, "storeWeightedSums()_wsum_ave_Pmax", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_avg_norm_correction, "storeWeightedSums()_wsum_avg_norm_correction", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_signal_product_spectra[0].wptr(ori_Fsize), ori_Fsize, "storeWeightedSums()_wsum_signal_product_spectra1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_signal_product_spectra[mlModel.nr_groups-1].wptr(ori_Fsize), ori_Fsize, "toreWeightedSums()_wsum_signal_product_spectraN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_reference_power_spectra[0].wptr(ori_Fsize),ori_Fsize , "storeWeightedSums()_wsum_reference_power_spectra1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_reference_power_spectra[mlModel.nr_groups-1].wptr(ori_Fsize), ori_Fsize, "storeWeightedSums()_wsum_reference_power_spectraN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_noise[0].wptr(ori_Fsize), ori_Fsize, "storeWeightedSums()_wsum_sigma2_noise1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_noise[mlModel.nr_groups-1].wptr(ori_Fsize), ori_Fsize, "storeWeightedSums()_wsum_sigma2_noiseN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[0], "storeWeightedSums()_wsum_sumw_group1", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sumw_group.wptrAll()[mlModel.nr_groups-1], "storeWeightedSums()_wsum_sumw_groupN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_sigma2_offset, "storeWeightedSums()_wsum_sigma2_offset", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_pdf_class.wptrAll()[0], "storeWeightedSums()_wsum_pdf_class0", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_pdf_class.wptrAll()[nr_classes-1], "storeWeightedSums()_wsum_pdf_classN", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_pdf_direction[0].wptr(exp_nr_dir), exp_nr_dir, "storeWeightedSums()_wsum_pdf_direction0", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.wsum_pdf_direction[nr_classes-1].wptr(exp_nr_dir), exp_nr_dir, "storeWeightedSums()_wsum_pdf_directionN", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}

namespace DoUpdateOtherParams
{
    int exp_current_Fsize2;
    int shell_size2;
    FDOUBLE COMPUTE_FLAG = (std::numeric_limits<FDOUBLE>::max)();

    //
    inline void applyCtfToClass(int iimage,FDOUBLE* &Frefctf_real,FDOUBLE* &Frefctf_imag,int n_start,int n_end)
    {
        // Apply CTF to reference projection
        // after first iteration refs_are_ctf_corrected = true
        if (do_ctf_correction && refs_are_ctf_corrected)
        {
            auto exp_local_Fctfs_iimage = exp_local_Fctfs[iimage].wptr(exp_current_Fsize2);
#pragma vector aligned
#pragma ivdep
            for (int n = n_start; n < n_end; n++) {
                Frefctf_real[n] = Frefctf_real[n]*exp_local_Fctfs_iimage[n];
                Frefctf_imag[n] = Frefctf_imag[n]*exp_local_Fctfs_iimage[n];
            }
        }
        if (do_scale_correction)
        {
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            // double myscale = mlModel.scale_correction.rptrAll()[igroup];
            double myscale = exp_metadata[iimage].SCALE;
            check_scale(myscale,igroup);
#pragma vector aligned
#pragma ivdep
            for (int n = n_start; n < n_end; n++)
            {
                Frefctf_real[n] = Frefctf_real[n]*myscale;
                Frefctf_imag[n] = Frefctf_imag[n]*myscale;
            }
        }
    }

    //
    void updateModel(const double* thread_exp_Mweight_sub_tid,
                     const char* exp_Mcoarse_significant_rptr,
                     double* thread_wsum_norm_correction_tid,
                     int tid,int iclass, int idir, int ipsi,
                     int shell_n_start, int shell_n_end,
                     int shell_norm_start, int shell_norm_end,
                     int shell_scale_start, int shell_scale_end,
                     int iover_rot_tile_start, int iover_rot_tile_end,
					 ImagePermutor const & imagePermutor, int permuted_iimage_tile_start, int permuted_iimage_tile_end,
                     int iimage_sub_tile, int itrans_sub_tile)
    {
        UpdateModel_Kernel kernel(exp_current_size, shell_n_start, shell_n_end, shell_norm_start, shell_norm_end,
                                  shell_scale_start, shell_scale_end, do_scale_correction);
        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                    exp_over_rot, exp_over_tilt, exp_over_psi);
        // TODO : some iover_rot projected image maybe not need but still compute..
        for (int iover_rot = iover_rot_tile_start; iover_rot < iover_rot_tile_end; iover_rot++)
        {
            // get projected class and shifted image
            auto Frefctf_real = thread_Frefctf_real[tid][iover_rot-iover_rot_tile_start].wptr(exp_current_Fsize2);
            auto Frefctf_imag = thread_Frefctf_imag[tid][iover_rot-iover_rot_tile_start].wptr(exp_current_Fsize2);
            FDOUBLE A[3][3];
            Euler_angles2matrix(exp_over_rot[iover_rot], exp_over_tilt[iover_rot], exp_over_psi[iover_rot], A);
            mapModel.get2DFourierTransformOneTile(iclass, Frefctf_real, Frefctf_imag, shell_n_start, shell_n_end,
                                                  exp_current_size, A, fourierShellTrans.rptr_nIndexFine());
        }

        for (int permuted_iimage_sub_tile_start = permuted_iimage_tile_start; 
			permuted_iimage_sub_tile_start < permuted_iimage_tile_end; 
			permuted_iimage_sub_tile_start += iimage_sub_tile)
        {
            for (int itrans_sub_tile_start = 0; itrans_sub_tile_start < exp_nr_trans; itrans_sub_tile_start+=itrans_sub_tile)
            {// 100(offset_range=10,step=2) or 400(offset_range=10,step=1)
                int permuted_iimage_sub_tile_end = std::min(permuted_iimage_sub_tile_start+iimage_sub_tile, permuted_iimage_tile_end);
                int itrans_sub_tile_end = std::min(itrans_sub_tile_start+itrans_sub_tile, exp_nr_trans);
                for (int permuted_iimage = permuted_iimage_sub_tile_start; permuted_iimage < permuted_iimage_sub_tile_end; permuted_iimage++){
					const auto unpermuted_iimage = imagePermutor.image(permuted_iimage);
                    for (int itrans = itrans_sub_tile_start; itrans < itrans_sub_tile_end; itrans++)
                    {
                        // TOOD.. TODO...TODO,need to reset because it use new significant point
                        if (!exp_Mcoarse_significant_rptr[unpermuted_iimage*exp_nr_trans+itrans]) continue;
                        //
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
                        auto Fimg_real 	=	thread_Fimg_real[tid][0].wptrAll();
                        auto Fimg_imag 	=	thread_Fimg_imag[tid][0].wptrAll();
#endif
#if defined(DOUBLE_TRANSLATION)
                            particleModel.getLargeShiftedMaskImageOneTile(unpermuted_iimage, Fimg_real, Fimg_imag, shell_n_start, shell_n_end, itrans);
#elif defined(TRIPLE_TRANSLATION)
                            particleModel.getLargeShiftedMaskImageDecompOneTile(unpermuted_iimage, Fimg_real, Fimg_imag, shell_n_start, shell_n_end, itrans, samplingGrid);
#endif
                        //
                        auto thread_wsum_sigma2_noise_tid_iimage 		= thread_wsum_sigma2_noise       [tid][unpermuted_iimage].wptrAll();
                        auto thread_wsum_scale_correction_XA_tid_iimage = thread_wsum_scale_correction_XA[tid][unpermuted_iimage].wptrAll();
                        auto thread_wsum_scale_correction_AA_tid_iimage = thread_wsum_scale_correction_AA[tid][unpermuted_iimage].wptrAll();
                        //
                        auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid +
                                                                (permuted_iimage - permuted_iimage_tile_start) * exp_nr_trans *
                                                                exp_nr_over_trans * exp_nr_over_rot;
                        //
                        double myscale = exp_metadata[unpermuted_iimage].SCALE;
                        // NOTE : in first cc,refs_are_ctf_corrected = false
                        auto local_ctf = refs_are_ctf_corrected?exp_local_Fctfs[unpermuted_iimage].rptrAll():particleModel.FctfsOne[unpermuted_iimage].rptrAll();
                        kernel.acquireImage(
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
                                            Fimg_real, Fimg_imag,
#else
                                            particleModel.Fimages_mask_fine_real[unpermuted_iimage].rptrAll(),
                                            particleModel.Fimages_mask_fine_imag[unpermuted_iimage].rptrAll(),
#endif
                                            local_ctf, thread_wsum_sigma2_noise_tid_iimage,
                                            thread_wsum_scale_correction_XA_tid_iimage,
                                            thread_wsum_scale_correction_AA_tid_iimage, myscale);

                        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++)
                        {
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
							auto shifter = particleModel.smallShiftedABTable.acquireShifter(0,iover_trans, iover_trans);
#else
							auto shifter = particleModel.shiftImageAssistor.acquireShifter(itrans, iover_trans, itrans*exp_nr_over_trans+iover_trans);
#endif
							kernel.acquireTable(
                                 shifter->aTable_rptr(),
                                 shifter->bTable_rptr(),
                                 shifter->tableSize());

                            for (int iover_rot = iover_rot_tile_start; iover_rot < iover_rot_tile_end; iover_rot++)
                            {
                                double weight = thread_exp_Mweight_sub_sub_tid[(itrans*exp_nr_over_trans+iover_trans)*exp_nr_over_rot+iover_rot];
                                if (weight==0) continue;

                                auto Frefctf_real = thread_Frefctf_real[tid][iover_rot-iover_rot_tile_start].rptr(exp_current_Fsize2);
                                auto Frefctf_imag = thread_Frefctf_imag[tid][iover_rot-iover_rot_tile_start].rptr(exp_current_Fsize2);
                                kernel.appendSumWdiff2s(Frefctf_real, Frefctf_imag, weight);
                            }// end loop iover_rot
                            kernel.compute();
                            kernel.release(thread_wsum_norm_correction_tid[unpermuted_iimage]);
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
							particleModel.smallShiftedABTable.releaseShifter(shifter);
#else
							particleModel.shiftImageAssistor.releaseShifter(shifter);
#endif
                        }// end loop iover_trans
                    }// end loop itrans
                }// end loop iimage
        	}// end loop itrans_sub_tile
        }// end loop iimage_sub_tile
    }

    void updateModelMain(TileCfg tileCfg)
    {
        //
        exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
        //
        shell_size2 = fourierShellTrans.fineFsize2();
        // set tile
        int shell_N_tile,ipsi_tile,iover_rot_tile,iimage_tile,itrans_sub_tile,iimage_sub_tile;
        tileCfg.getTileForUpdateModel(shell_size2, shell_N_tile, ipsi_tile, iover_rot_tile, iimage_tile, iimage_sub_tile, itrans_sub_tile);
        //
        int shell_norm_start = fourierShellTrans.getNormCorrLo();
        int shell_norm_end   = fourierShellTrans.getNormCorrHi();
        //
		ImagePermutor imagePermutor;
		imagePermutor.setImageRange(0,exp_nr_images);
		// For each image, enter a (class,dir,ipsi) for which it is significant
		{
			for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++) {
				for (int idir = 0; idir < exp_nr_dir; idir++) {
					for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++) {
						int iClassOrient = iclass*exp_nr_dir + idir;
						iClassOrient = iClassOrient*exp_nr_psi + ipsi;

                        if (!exp_Mcoarse_Rot_significant.isRotSignificant(iclass, idir, ipsi));
							continue;

                        auto exp_Mcoarse_significant_rptr = exp_Mcoarse_Rot_significant.isMcoarseSignificantRptrAll(iclass, idir, ipsi);
						for (int iimage = 0; iimage < exp_nr_images; iimage++) {
							if (imagePermutor.keyAlreadySet(iimage)) 
								continue;
							for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
								if (!exp_Mcoarse_significant_rptr[iimage*exp_nr_trans + itrans]) continue;
								imagePermutor.setKey(iimage, iClassOrient);
								break;
							}
						}
					}
				}
			}
		}
		imagePermutor.permute();
		
        for (int shell_n_start = 0; shell_n_start < shell_size2; shell_n_start += shell_N_tile)
        {
#pragma omp parallel for collapse(4) schedule(dynamic)
            for (int idir = 0; idir < exp_nr_dir; idir++)
            {
                for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
                {
                    for (int ipsi_tile_start = 0; ipsi_tile_start < exp_nr_psi; ipsi_tile_start+=ipsi_tile)
                    {
                        for (int iover_rot_tile_start = 0; iover_rot_tile_start < exp_nr_over_rot; iover_rot_tile_start += iover_rot_tile)// 8
                        {
                            for (int permuted_iimage_tile_start = 0; permuted_iimage_tile_start < exp_nr_images; permuted_iimage_tile_start += iimage_tile)
                            {
                                // some thread data
                                const int tid = omp_get_thread_num();
                                const auto thread_wsum_norm_correction_tid = thread_wsum_norm_correction[tid].wptr(exp_nr_images);
                                const int shell_scale_start           = fourierShellTrans.getScaleCorrLo(iclass);
                                const int shell_scale_end             = fourierShellTrans.getScaleCorrHi(iclass);
                                const int shell_n_end                 = std::min(shell_n_start+shell_N_tile, shell_size2);
                                const int permuted_iimage_tile_end    = std::min(permuted_iimage_tile_start+iimage_tile, exp_nr_images);
                                const int ipsi_tile_end               = std::min(ipsi_tile_start+ipsi_tile,exp_nr_psi);
                                const int iover_rot_tile_end          = std::min(iover_rot_tile_start+iover_rot_tile,exp_nr_over_rot);
                                const int thread_exp_Mweight_sub_size = const_max_iimage_tile*exp_nr_trans*exp_nr_over_trans*exp_nr_over_rot;
                                const auto thread_exp_Mweight_sub_tid = thread_exp_Mweight_sub[tid].wptr(thread_exp_Mweight_sub_size);
                                //
                                for (int ipsi = ipsi_tile_start; ipsi < ipsi_tile_end; ipsi++)
                                {
                                    // TOOD...need to reset exp_Rot_significant and exp_Mcoarse_significant
                                    // because the significant point is updated
                                    if (!exp_Mcoarse_Rot_significant.isRotSignificant(iclass, idir, ipsi)) continue;
                                    auto exp_Mcoarse_significant_rptr = exp_Mcoarse_Rot_significant.isMcoarseSignificantRptrAll(iclass, idir, ipsi);
                                    // set thread_exp_Mweight_sub_tid
                                    // TODO... need put this function into the exp_Mweight for copy in and out part of it.
                                    for (int i = 0; i < thread_exp_Mweight_sub_size; i++)
                                        thread_exp_Mweight_sub_tid[i] = 0.;
                                    for (int permuted_iimage = permuted_iimage_tile_start; permuted_iimage < permuted_iimage_tile_end; permuted_iimage++)
                                    {
										typedef void iimage;	// because need to be explicit whether using iimagePermuteIndex or original_iimage

										const int unpermuted_iimage = imagePermutor.image(permuted_iimage);
#ifdef	EXP_MWEIGHT_NEW
                                        const auto  exp_Mweight_fine_sub_ptr = exp_Mweight_fine.mweightsForSomeSpinsAndSlidesOrNullptr(unpermuted_iimage, iclass, idir, ipsi);
                                        if (!exp_Mweight_fine_sub_ptr) continue;
                                        const auto& exp_Mweight_fine_sub = *exp_Mweight_fine_sub_ptr;
#else
                                        const auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(unpermuted_iimage, iclass, idir, ipsi);
                                        if (exp_Mweight_fine_sub.size() == 0) continue;
#endif
                                        auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid +
                                                                              (permuted_iimage - permuted_iimage_tile_start) * exp_nr_trans *
                                                                              exp_nr_over_trans * exp_nr_over_rot;

                                        // NOTE : the order of this for loop is same as getAllSquaredDifferencesFine()
                                        for (size_t i = 0; i < exp_Mweight_fine_sub.size(); i++)
                                        {
#ifdef	EXP_MWEIGHT_NEW
                                            double weight = exp_Mweight_fine_sub.value(i);
#else
                                            double weight = exp_Mweight_fine_sub[i].second;
#endif
                                            if (weight < exp_significant_weight.rptrAll()[unpermuted_iimage]) continue;
#ifdef	EXP_MWEIGHT_NEW
                                            int itrans,iover_rot,iover_trans;
                                            exp_Mweight_fine_sub.key(i, iover_rot, itrans, iover_trans);
                                            int ihidden_tid = itrans*exp_nr_over_trans + iover_trans;
                                            ihidden_tid = ihidden_tid*exp_nr_over_rot + iover_rot;
#else
                                            auto ihidden_tid = exp_Mweight_fine_sub[i].first;
#endif
                                            thread_exp_Mweight_sub_sub_tid[ihidden_tid] = weight/exp_sum_weight.rptrAll()[unpermuted_iimage];
                                            //
                                        }
                                    }
                                    //
                                    updateModel(thread_exp_Mweight_sub_tid,
                                                exp_Mcoarse_significant_rptr,
                                                thread_wsum_norm_correction_tid,
                                                tid, iclass, idir, ipsi,
                                                shell_n_start, shell_n_end,
                                                shell_norm_start, shell_norm_end,
                                                shell_scale_start, shell_scale_end,
                                                iover_rot_tile_start, iover_rot_tile_end,
												imagePermutor, permuted_iimage_tile_start, permuted_iimage_tile_end,
                                                iimage_sub_tile, itrans_sub_tile);
                                }// end loop ipsi
                            }// end loop iimage tile
                        }// end loop iover_rot_tile
                    }// end loop ipsi tile
                }// end of iclass
            }// end loop idir
        }// end loop n_tile
        //
    }//
    //
    void updateOtherParams()
    {
        thread_wsum_scale_correction_XA	.fill_with_first_touch(0.);
        thread_wsum_scale_correction_AA	.fill_with_first_touch(0.);
        thread_wsum_sigma2_noise		.fill_with_first_touch(0.);
        thread_wsum_norm_correction		.fill_with_first_touch(0.);

        if (/*tune*/false)
        {
#define VTUNE_ANA 1

#ifdef VTUNE_ANA
            //
            //
            TileCfg tilecfg;
            tilecfg.choiceTileForKNL(1);
            for (int i = 0; i < 5; i++)
            {
                double time_start = dtime();
                DoUpdateOtherParams::updateModelMain(tilecfg);

                double time_end = dtime();
                std::cout	<<"count : "<<i<<"time : "<<(time_end-time_start)<<std::endl;
            }
#endif
            EXIT_ABNORMALLY;
        }

        TileCfg tilecfg;
        DoUpdateOtherParams::updateModelMain(tilecfg);

        // reduce all the threads data
        exp_wsum_norm_correction		.fill(0.);
        exp_wsum_scale_correction_XA	.fill(0.);
        exp_wsum_scale_correction_AA	.fill(0.);

        for (int thread = 0; thread < maxthreads; thread++)
        {
            //
			{
				auto thread_wsum_norm_correction_aux = thread_wsum_norm_correction[thread].rptr(exp_nr_images);
				auto ptr = exp_wsum_norm_correction.wptrAll();
				for (int iimage = 0; iimage < exp_nr_images; iimage++)
					ptr[iimage] += thread_wsum_norm_correction_aux[iimage];
			}
            //
            for (int iimage = 0; iimage < exp_nr_images; iimage++)
            {
                int igroup = exp_metadata[iimage].GROUP_NO-1;
                auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
                auto thread_wsum_sigma2_noise_iimage = thread_wsum_sigma2_noise[thread][iimage].wptr(exp_current_Fsize2);
                auto thread_wsum_scale_correction_XA_iimage = thread_wsum_scale_correction_XA[thread][iimage].wptr(exp_current_Fsize2);
                auto thread_wsum_scale_correction_AA_iimage = thread_wsum_scale_correction_AA[thread][iimage].wptr(exp_current_Fsize2);
                //
                fourierShellTrans.transformBackFine(thread, thread_wsum_sigma2_noise_iimage, exp_current_Fsize2);
                fourierShellTrans.transformBackFine(thread, thread_wsum_scale_correction_XA_iimage, exp_current_Fsize2);
                fourierShellTrans.transformBackFine(thread, thread_wsum_scale_correction_AA_iimage, exp_current_Fsize2);
                for (int n = 0; n < exp_current_Fsize2; n++)
                {
                    // mapping the data to specturm
                    int ires = Mresol_fine.rptrAll()[n];
                    if (ires > -1){
                        wsum_sigma2_noise_igroup[ires] += thread_wsum_sigma2_noise_iimage[n];
                        if (do_scale_correction)
                        {
                            exp_wsum_scale_correction_XA[iimage].wptrAll()[ires] += thread_wsum_scale_correction_XA_iimage[n];
                            exp_wsum_scale_correction_AA[iimage].wptrAll()[ires] += thread_wsum_scale_correction_AA_iimage[n];
                        }
                    }
                }// end loop n
            }// end loop iimage
        }// end loop thread
        //

#ifdef DATA_STREAM
        global_data_stream.foutDouble(exp_wsum_norm_correction.wptrAll()[0], "updateOtherParams()_exp_wsum_norm_correction", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mlModel.wsum_sigma2_noise, mlModel.ori_Fsize, "updateOtherParams()_wsum_sigma2_noise1", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mlModel.wsum_sigma2_noise+(mlModel.nr_groups-1)*mlModel.ori_Fsize, mlModel.ori_Fsize, "updateOtherParams()_wsum_sigma2_noiseN", __FILE__, __LINE__);
        // todo,print different group
        if (do_scale_correction) {
            global_data_stream.foutDouble(exp_wsum_scale_correction_XA[0].wptrAll(), mlModel.ori_Fsize, "updateOtherParams()_exp_wsum_scale_correction_XA", __FILE__, __LINE__);
            global_data_stream.foutDouble(exp_wsum_scale_correction_AA[0].wptrAll(), mlModel.ori_Fsize, "updateOtherParams()_exp_wsum_scale_correction_AA", __FILE__, __LINE__);
        }
        global_data_stream.check();global_data_stream.flush();
#endif
    }
}

void updateOtherParams()
{
    //
    DoUpdateOtherParams::updateOtherParams();

    //
    thread_wsum_pdf_direction	.fill_with_first_touch(0.);
    thread_wsum_pdf_class		.fill_with_first_touch(0.);
    thread_wsum_sigma2_offset	.fill_with_first_touch(0.);
    thread_sumw_group			.fill_with_first_touch(0.);

#pragma omp parallel for collapse(3) schedule(dynamic)
    for (int iimage = 0; iimage < exp_nr_images; iimage++)
    {
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
            for (int idir = 0; idir < exp_nr_dir; idir++)
            {
                for (int ipsi = 0; ipsi < exp_nr_psi; ipsi++)
                {
					//
                    int tid = omp_get_thread_num();
                    double thread_wsum_pdf_iclass = 0.;
                    double thread_wsum_sigma2_offset_i = 0.;
                    auto thread_wsum_pdf_direction_tid = thread_wsum_pdf_direction[tid][iclass].wptr(exp_nr_dir);
                    auto thread_sumw_group_tid = thread_sumw_group[tid].wptr(exp_nr_images);
#ifdef EXP_MWEIGHT_NEW
					auto exp_Mweight_fine_sub_ptr = exp_Mweight_fine.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi);
					if (exp_Mweight_fine_sub_ptr) {
						auto& exp_Mweight_fine_sub = *exp_Mweight_fine_sub_ptr;
#else
					if (true) {
						auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(iimage, iclass, idir, ipsi);
#endif
						for (size_t i = 0; i < exp_Mweight_fine_sub.size(); i++)
						{
#ifdef EXP_MWEIGHT_NEW
							double weight = exp_Mweight_fine_sub.value(i);
#else
							double weight = exp_Mweight_fine_sub[i].second;
							int ihidden = exp_Mweight_fine_sub[i].first;
#endif
							// TODO,vectorize this???
							if (weight >= exp_significant_weight.rptrAll()[iimage])
							{
								// Normalise the weight (do this after the comparison with exp_significant_weight!)
								weight /= exp_sum_weight.rptrAll()[iimage];

								// Store sum of weights for this group
								thread_sumw_group_tid[iimage] += weight;

								// Store weights for this class and orientation
								thread_wsum_pdf_iclass += weight;

								thread_wsum_pdf_direction_tid[idir] += weight;

								int iover_rot, itrans, iover_trans;
#ifdef EXP_MWEIGHT_NEW
								exp_Mweight_fine_sub.key(i, iover_rot, itrans, iover_trans);
#else
								decodeRotOverTrans(ihidden, itrans, iover_trans, iover_rot);
#endif

								double offsetx = mlModel.prior_offsetx_class.rptrAll()[iclass] - exp_old_offsetx.rptrAll()[iimage] - samplingGrid.exp_over_trans_x[itrans][iover_trans];
								double offsety = mlModel.prior_offsety_class.rptrAll()[iclass] - exp_old_offsety.rptrAll()[iimage] - samplingGrid.exp_over_trans_y[itrans][iover_trans];

								//this may cause some false share!
								thread_wsum_sigma2_offset_i += weight * (offsetx*offsetx+offsety*offsety);
							}
						}
					}
                    // Fix the NUMA-Issua for padding size small than cache-line
                    auto thread_wsum_pdf_class_tid = thread_wsum_pdf_class[tid].wptr(nr_classes);
                    auto thread_wsum_sigma2_offset_tid = thread_wsum_sigma2_offset[tid].wptr(1);
                    thread_wsum_pdf_class_tid[iclass] += thread_wsum_pdf_iclass;
                    thread_wsum_sigma2_offset_tid[0] += thread_wsum_sigma2_offset_i;
                }// end loop ipsi
            } // end loop idir
        }// end of iclass
    } // end loop iimage

    for (int thread = 0; thread < maxthreads; thread++)
    {
        auto thread_wsum_sigma2_offset_aux = thread_wsum_sigma2_offset[thread].wptr(1);
        mlModel.wsum_sigma2_offset += thread_wsum_sigma2_offset_aux[0];
        auto thread_wsum_pdf_class_aux = thread_wsum_pdf_class[thread].wptr(nr_classes);
        for (int iclass = 0; iclass < nr_classes; iclass++)
        {
            mlModel.wsum_pdf_class.wptrAll()[iclass] += thread_wsum_pdf_class_aux[iclass];
            auto thread_wsum_pdf_direction_aux = thread_wsum_pdf_direction[thread][iclass].wptr(exp_nr_dir);
            for (int idir = 0; idir < exp_nr_dir; idir++){
                mlModel.wsum_pdf_direction[iclass][idir] += thread_wsum_pdf_direction_aux[idir];
            }
        }
        auto thread_sumw_group_tid = thread_sumw_group[thread].wptr(exp_nr_images);
        for (int iimage = 0; iimage < exp_nr_images; iimage++) {
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            mlModel.wsum_sumw_group.wptrAll()[igroup] += thread_sumw_group_tid[iimage];
        }
    }
#ifdef DATA_STREAM
    // global_data_stream.foutDouble(mlModel.wsum_sumw_group[0], "updateOtherParams()_wsum_sumw_group1", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_sumw_group[mlModel.nr_groups-1], "updateOtherParams()_wsum_sumw_groupN", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_sigma2_offset, "updateOtherParams()_wsum_sigma2_offset", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_pdf_class[0], "updateOtherParams()_wsum_pdf_class0", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_pdf_class[nr_classes-1], "updateOtherParams()_wsum_pdf_classN", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_pdf_direction, exp_nr_dir, "updateOtherParams()_wsum_pdf_direction0", __FILE__, __LINE__);
    // global_data_stream.foutDouble(mlModel.wsum_pdf_direction+(nr_classes-1)*mlModel.nr_directions, exp_nr_dir, "updateOtherParams()_wsum_pdf_directionN", __FILE__, __LINE__);
    // global_data_stream.check();global_data_stream.flush();
#endif
    /*** do on convertSquaredDifferencesToWeights()
	****do on convertSquaredDifferencesToWeights() ***///
#ifdef DATA_STREAM
    global_data_stream.foutDouble(exp_metadata[0].ROT, "exp_metadata[0].ROT", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].TILT, "exp_metadata[0].TILT", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].PSI, "exp_metadata[0].PSI", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].XOFF, "exp_metadata[0].XOFF", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].YOFF, "exp_metadata[0].YOFF", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].CLASS, "exp_metadata[0].CLASS", __FILE__, __LINE__);
    global_data_stream.foutDouble(exp_metadata[0].PMAX, "exp_metadata[0].PMAX", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}

// TODO : to fix the MAPModel.reduce() and large 3D map issue
namespace DoBackProjection
{
    int shell_size2;
    FDOUBLE COMPUTE_FLAG = (std::numeric_limits<FDOUBLE>::max)();

    inline void global_data_stream_inner(int iimage,int itrans,int iover_trans,double weight,
                                         FDOUBLE* Fimg_real,FDOUBLE* Fimg_imag,FDOUBLE* Fweight,int exp_current_Fsize2)
    {
#ifdef DATA_STREAM
        // global_data_stream.foutInt(iclass, "backProjection()_iclass", __FILE__, __LINE__);
        // global_data_stream.foutInt(idir, "backProjection()_idir", __FILE__, __LINE__);
        // global_data_stream.foutInt(ipsi, "backProjection()_ipsi", __FILE__, __LINE__);
        // global_data_stream.foutInt(iover_rot, "backProjection()_iover_rot", __FILE__, __LINE__);
        global_data_stream.foutInt(iimage, "backProjection()_iimage", __FILE__, __LINE__);
        global_data_stream.foutInt(itrans, "backProjection()_itrans", __FILE__, __LINE__);
        global_data_stream.foutInt(iover_trans, "backProjection()_iover_trans", __FILE__, __LINE__);
        global_data_stream.foutDouble(weight, "backProjection()_weight", __FILE__, __LINE__);
        global_data_stream.foutDouble(exp_significant_weight.wptrAll()[iimage], "backProjection()_exp_significant_weight[iimage]", __FILE__, __LINE__);
        global_data_stream.foutDouble(Fimg_real, exp_current_Fsize2, "backProjection()_backProjection()_Fimg_real", __FILE__, __LINE__);
        global_data_stream.foutDouble(Fimg_imag, exp_current_Fsize2, "backProjection()_backProjection()_Fimg_imag", __FILE__, __LINE__);
        global_data_stream.foutDouble(Fweight, exp_current_Fsize2, "backProjection()_backProjection()_Fweight", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
    }
    //
    void prepareNomaskImageAndMinvsigma2s()
    {
        // In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the exp_local_Minvsigma2s was omitted.
        // Set those back here
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            int igroup = exp_metadata[iimage].GROUP_NO-1;
            auto exp_local_Minvsigma2s_image = exp_local_Minvsigma2s[iimage].wptrAll();
            exp_local_Minvsigma2s_image[0] = 1. / (sigma2_fudge * mlModel.sigma2_noise[igroup][0]);
        }
#ifdef DATA_STREAM
        global_data_stream.foutDouble(exp_local_Minvsigma2s[0].wptrAll(),exp_nr_images*exp_current_size*(exp_current_size/2+1),
                                      "backProjection()_exp_local_Minvsigma2s is diff because change the array order...", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
#pragma omp parallel for
        for (int iimage = 0; iimage < exp_nr_images; iimage++)
        {
            auto exp_local_Fctfs_iimage = exp_local_Fctfs[iimage].wptrAll();
            auto exp_local_Minvsigma2s_iimage = exp_local_Minvsigma2s[iimage].wptrAll();

            if (!do_map)
            {
                for (int n = 0; n < shell_size2; n++)
                    exp_local_Minvsigma2s_iimage[n] = 1.0;
            }
            // Apply CTF to reference
            if (!do_ctf_correction)
            {
                for (int n = 0; n < shell_size2; n++)
                    exp_local_Fctfs_iimage[n] = 1.0;
            }
            if (do_scale_correction)
            {
                // For CTF-terms in BP
                int igroup = exp_metadata[iimage].GROUP_NO-1;
                double myscale = mlModel.scale_correction.rptrAll()[igroup];
                check_scale(myscale, igroup);
                for (int n = 0; n < shell_size2; n++)
                    exp_local_Fctfs_iimage[n] = exp_local_Fctfs_iimage[n] * myscale;
            }
            // NOTE : later code do some trick for Fimages and exp_local_Minvigma2s
            // Operation : exp_local_Minvsigma2s =  exp_local_Fctfs*exp_local_Minvsigma2s
            auto exp_local_Minvsigma2s_X_Fctfs_iimage = exp_local_Minvsigma2s[iimage].wptrAll();
            for (int n = 0; n < shell_size2; n++)
                exp_local_Minvsigma2s_X_Fctfs_iimage[n] = exp_local_Minvsigma2s_iimage[n]*exp_local_Fctfs_iimage[n];
            // NOTE : multiply nomask image by Minvsigma2s and Fctfs(exp_local_Minvsigma2s_X_Fctfs)
            // this will reduce some memory access in inner for loop
            auto Fimages_nomask_real = particleModel.Fimages_nomask_real[iimage].wptrAll();
            auto Fimages_nomask_imag = particleModel.Fimages_nomask_imag[iimage].wptrAll();
            for (int n = 0; n < shell_size2; n++) {
                Fimages_nomask_real[n] *= exp_local_Minvsigma2s_X_Fctfs_iimage[n];
                Fimages_nomask_imag[n] *= exp_local_Minvsigma2s_X_Fctfs_iimage[n];
            }
            // NOTE : multiply Fctf to exp_local_Minvsigma2s
            // then exp_local_Minvsigma2s will be equal to Minvsigma2s*Fctf^2
            auto exp_local_Minvsigma2s_X_Fctfs2_iimage = exp_local_Minvsigma2s[iimage].wptrAll();
            for (int n = 0; n < shell_size2; n++)
                exp_local_Minvsigma2s_X_Fctfs2_iimage[n] = exp_local_Minvsigma2s_X_Fctfs_iimage[n]*exp_local_Fctfs_iimage[n];
        }
    }
    //
    void backProjection(const double* thread_exp_Mweight_sub_tid,
                        const char* exp_Mcoarse_significant_rptr,
                        int tid,int iclass, int idir, int ipsi,
                        int shell_n_start, int shell_n_end,
                        int iover_rot_tile_start,int iover_rot_tile_end,
                        int iimage_tile_start,int iimage_tile_end,
                        int iimage_sub_tile, int itrans_sub_tile)
    {
        BackProjection_Kernel kernel(exp_current_size, shell_n_start, shell_n_end);
        FDOUBLE *exp_over_rot,*exp_over_tilt,*exp_over_psi;
        samplingGrid.setOrientation(sampler3d, exp_current_oversampling, idir, ipsi,
                                    exp_over_rot, exp_over_tilt, exp_over_psi);
        // flag whether need to re-compute projected class
        for (int iover_rot = iover_rot_tile_start; iover_rot < iover_rot_tile_end; iover_rot++)
        {
            auto Frefctf_real = thread_Frefctf_real[tid][iover_rot-iover_rot_tile_start].wptrAll();
            Frefctf_real[shell_n_start] = COMPUTE_FLAG;
        }
        //
        for (int iimage_sub_tile_start = iimage_tile_start; iimage_sub_tile_start < iimage_tile_end; iimage_sub_tile_start+=iimage_sub_tile)
        {
            for (int itrans_sub_tile_start = 0; itrans_sub_tile_start < exp_nr_trans; itrans_sub_tile_start+=itrans_sub_tile)
            {// 100(offset_range=10,step=2) or 400(offset_range=10,step=1)
                int iimage_sub_tile_end = std::min(iimage_sub_tile_start+iimage_sub_tile, iimage_tile_end);
                int itrans_sub_tile_end = std::min(itrans_sub_tile_start+itrans_sub_tile, exp_nr_trans);
                for (int iimage = iimage_sub_tile_start; iimage < iimage_sub_tile_end; iimage++){
                    for (int itrans = itrans_sub_tile_start; itrans < itrans_sub_tile_end; itrans++)
                    {
                        // TOOD.. TODO...TODO,need to reset because it use new significant point
                        if (!exp_Mcoarse_significant_rptr[iimage*exp_nr_trans+itrans]) continue;
                        auto exp_local_Minvsigma2s_X_Fctfs2_iimage = exp_local_Minvsigma2s[iimage].rptrAll();
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
                        auto Fimg_real 	=	thread_Fimg_real[tid][0].wptrAll();
                        auto Fimg_imag 	=	thread_Fimg_imag[tid][0].wptrAll();
#endif
#if defined(DOUBLE_TRANSLATION)
                        particleModel.getLargeShiftedNomaskImageOneTile(iimage, Fimg_real, Fimg_imag, shell_n_start, shell_n_end, itrans);
#elif defined(TRIPLE_TRANSLATION)
                        particleModel.getLargeShiftedNomaskImageDecompOneTile(iimage, Fimg_real, Fimg_imag, shell_n_start, shell_n_end, itrans, samplingGrid);
#endif
                        //
                        auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid +
                        (iimage - iimage_tile_start) * exp_nr_trans *
                        exp_nr_over_trans * exp_nr_over_rot;
                        kernel.acquireImage(
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
                                            Fimg_real, Fimg_imag,
#else
                                            particleModel.Fimages_nomask_real[iimage].rptrAll(),
                                            particleModel.Fimages_nomask_imag[iimage].rptrAll(),
#endif
                                            exp_local_Minvsigma2s_X_Fctfs2_iimage);
                        for(int iover_trans = 0;iover_trans < exp_nr_over_trans;iover_trans++)
                        {
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
							auto shifter = particleModel.smallShiftedABTable.acquireShifter(0,iover_trans, iover_trans);
#else
							auto shifter = particleModel.shiftImageAssistor.acquireShifter(itrans, iover_trans, itrans*exp_nr_over_trans+iover_trans);
#endif
							kernel.acquireTable(
                                shifter->aTable_rptr(),
                                shifter->bTable_rptr(),
                                shifter->tableSize());

                            for (int iover_rot = iover_rot_tile_start; iover_rot < iover_rot_tile_end; iover_rot++) {
                                //
                                double weight = thread_exp_Mweight_sub_sub_tid[(itrans*exp_nr_over_trans+iover_trans)*exp_nr_over_rot+iover_rot];
                                if (weight==0) continue;
                                //
                                auto Frefctf_real 	= thread_Frefctf_real[tid][iover_rot-iover_rot_tile_start].wptrAll();
                                auto Frefctf_imag 	= thread_Frefctf_imag[tid][iover_rot-iover_rot_tile_start].wptrAll();
                                auto Fweight   		= thread_Fweight[tid][iover_rot-iover_rot_tile_start]	  .wptrAll();
                                if (Frefctf_real[shell_n_start]==COMPUTE_FLAG)
                                {
                                    for (int n = shell_n_start; n < shell_n_end; n++)
                                        Frefctf_real[n] = Frefctf_imag[n] = Fweight[n] = 0.;
                                }
                                kernel.appendFref(Frefctf_real, Frefctf_imag, Fweight, weight);
                            }// end loop iover_rot
                            // Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
                            // Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
                            kernel.compute();
                            kernel.release();
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
							particleModel.smallShiftedABTable.releaseShifter(shifter);
#else
							particleModel.shiftImageAssistor.releaseShifter(shifter);
#endif
                        }// end loop iover_trans
                    }// end loop itrans
                }// end loop iimage
            }// end loop itrans_sub_tile
        }// end loop iimage_sub_tile
        // put all image back
        for (int iover_rot = iover_rot_tile_start; iover_rot < iover_rot_tile_end; iover_rot++)
        {
            auto Frefctf_real	= thread_Frefctf_real[tid][iover_rot-iover_rot_tile_start].rptrAll();
            auto Frefctf_imag	= thread_Frefctf_imag[tid][iover_rot-iover_rot_tile_start].rptrAll();
            auto Fweight		= thread_Fweight[tid][iover_rot-iover_rot_tile_start].rptrAll();
            if (Frefctf_real[shell_n_start] != COMPUTE_FLAG)
            {
                FDOUBLE A[3][3];
                Euler_angles2matrix(exp_over_rot[iover_rot], exp_over_tilt[iover_rot], exp_over_psi[iover_rot], A);
                mapModel.set2DFourierTransformOneTile(tid, Frefctf_real, Frefctf_imag, shell_n_start, shell_n_end,
                                                      exp_current_size, A, Fweight, fourierShellTrans.rptr_nIndexFine());
            }
        }
    }

    //
    void backProjectionMain(TileCfg& tileCfg)
    {
        shell_size2    = fourierShellTrans.fineFsize2();
        //
		maybeGoSerial("prepareNomaskImageAndMinvsigma2s");
        prepareNomaskImageAndMinvsigma2s();
        //
#ifdef DATA_STREAM
        // turn off data stream inside openmp for loop
        if(data_stream_node == node) global_data_stream.turnOff();
#endif
        // set tile
		maybeGoSerial("backProjectionMainLoop");
        int shell_N_tile 	= shell_size2;
        int N_tile,ipsi_tile,iover_rot_tile,iimage_tile,iimage_sub_tile,itrans_sub_tile;
        tileCfg.getTileForBackproject(shell_N_tile, N_tile, ipsi_tile, iover_rot_tile, iimage_tile, iimage_sub_tile, itrans_sub_tile);
        //
        for (int iclass = exp_iclass_min; iclass <= exp_iclass_max; iclass++)
        {
            for (int shell_n_start = 0; shell_n_start < shell_size2; shell_n_start+=shell_N_tile)
            {
#pragma omp parallel for collapse(4) schedule(dynamic)
                for (int idir = 0; idir < exp_nr_dir; idir++)
                {
                    for (int ipsi_tile_start = 0; ipsi_tile_start < exp_nr_psi; ipsi_tile_start += ipsi_tile)
                    {
                        for (int iover_rot_tile_start = 0; iover_rot_tile_start < exp_nr_over_rot; iover_rot_tile_start += iover_rot_tile)// 8
                        {
                            for (int iimage_tile_start = 0; iimage_tile_start < exp_nr_images; iimage_tile_start+=iimage_tile)
                            {
                                const int tid = omp_get_thread_num();
                                const int shell_n_end                 = std::min(shell_n_start+shell_N_tile, shell_size2);
                                const int ipsi_tile_end               = std::min(ipsi_tile_start+ipsi_tile, exp_nr_psi);
                                const int iimage_tile_end             = std::min(iimage_tile_start+iimage_tile, exp_nr_images);
                                const int iover_rot_tile_end          = std::min(iover_rot_tile_start+iover_rot_tile, exp_nr_over_rot);
                                const int thread_exp_Mweight_sub_size = const_max_iimage_tile*exp_nr_trans*exp_nr_over_trans*exp_nr_over_rot;
                                const auto thread_exp_Mweight_sub_tid = thread_exp_Mweight_sub[tid].wptr(thread_exp_Mweight_sub_size);

                                for (int ipsi = ipsi_tile_start; ipsi < ipsi_tile_end; ipsi++)
                                {
                                    // TOOD...need to reset exp_Rot_significant and exp_Mcoarse_significant
                                    // because the significant point is updated
                                    if (!exp_Mcoarse_Rot_significant.isRotSignificant(iclass, idir, ipsi)) continue;
                                    auto exp_Mcoarse_significant_rptr = exp_Mcoarse_Rot_significant.isMcoarseSignificantRptrAll(iclass, idir, ipsi);
                                    // set thread_exp_Mweight_sub_tid
                                    // TODO... need put this function into the exp_Mweight for copy in and out part of it.
                                    for (int i = 0; i < thread_exp_Mweight_sub_size; i++)
                                        thread_exp_Mweight_sub_tid[i] = 0.;
                                    for (int iimage = iimage_tile_start; iimage < iimage_tile_end; iimage++)
                                    {
#ifdef	EXP_MWEIGHT_NEW
                                        const auto  exp_Mweight_fine_sub_ptr = exp_Mweight_fine.mweightsForSomeSpinsAndSlidesOrNullptr(iimage, iclass, idir, ipsi);
                                        if (!exp_Mweight_fine_sub_ptr) continue;
                                        const auto& exp_Mweight_fine_sub = *exp_Mweight_fine_sub_ptr;
#else
                                        const auto& exp_Mweight_fine_sub = exp_Mweight_fine.wptr_sparse(iimage, iclass, idir, ipsi);
                                        if (exp_Mweight_fine_sub.size() == 0) continue;
#endif
                                        auto thread_exp_Mweight_sub_sub_tid = thread_exp_Mweight_sub_tid +
                                        (iimage - iimage_tile_start) * exp_nr_trans *
                                        exp_nr_over_trans * exp_nr_over_rot;

                                        // NOTE : the order of this for loop is same as getAllSquaredDifferencesFine()
                                        for (size_t i = 0; i < exp_Mweight_fine_sub.size(); i++)
                                        {
#ifdef	EXP_MWEIGHT_NEW
                                            double weight = exp_Mweight_fine_sub.value(i);
#else
                                            double weight = exp_Mweight_fine_sub[i].second;
#endif
                                            if (weight < exp_significant_weight.rptrAll()[iimage]) continue;
#ifdef	EXP_MWEIGHT_NEW
                                            int itrans,iover_rot,iover_trans;
                                            exp_Mweight_fine_sub.key(i, iover_rot, itrans, iover_trans);
                                            int ihidden_tid = itrans*exp_nr_over_trans + iover_trans;
                                            ihidden_tid = ihidden_tid*exp_nr_over_rot + iover_rot;
#else
                                            auto ihidden_tid = exp_Mweight_fine_sub[i].first;
#endif
                                            thread_exp_Mweight_sub_sub_tid[ihidden_tid] = weight/exp_sum_weight.rptrAll()[iimage];
                                            //
                                        }
                                    }
                                    //
                                    backProjection(thread_exp_Mweight_sub_tid, exp_Mcoarse_significant_rptr,
                                                   tid, iclass, idir, ipsi, shell_n_start, shell_n_end,
                                                   iover_rot_tile_start, iover_rot_tile_end,
                                                   iimage_tile_start, iimage_tile_end,
                                                   iimage_sub_tile, itrans_sub_tile);
                                }// end loop ipsi
                            } // end loop iimage_tile
                        }// end loop iover_rot_tile
                    }// enbd loop ipsi_tile
                } // end loop idir
            }// end loop n_tile
            // reduce iclass
            mapModel.reduceThreadMap(iclass);
        } // end loop iclass
        //
        // Reduce thread classes outside the backprojection()
#ifdef DATA_STREAM
        if(data_stream_node == node) global_data_stream.turnOn();
        int vol_size3 = mapModel.projector[0].pad_size;
        vol_size3 = vol_size3*vol_size3*(vol_size3/2+1);
        // global_data_stream.foutDouble(mapModel.backprojector[0].data_real.wptr(), vol_size3, "backProjection()_backRef_real_0", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[0].data_imag.wptr(), vol_size3, "backProjection()_backRef_imag_0", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[nr_classes-1].data_real.wptr(), vol_size3, "backProjection()_backRef_real_N-1", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[nr_classes-1].data_imag.wptr(), vol_size3, "backProjection()_backRef_imag_N-1", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[0].weight.wptr(), vol_size3, "backProjection()_backRef_weight_0", __FILE__, __LINE__);
        // global_data_stream.foutDouble(mapModel.backprojector[nr_classes-1].weight.wptr(), vol_size3, "backProjection()_backRef_weight_N-1", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.backprojector[0].weight.wptr(), vol_size3, "backProjection()_backRef_weight_0", __FILE__, __LINE__);
        global_data_stream.foutDouble(mapModel.backprojector[(nr_classes-1)].weight.wptr(), vol_size3, "backProjection()_backRef_weight_N-1", __FILE__, __LINE__);
        global_data_stream.check();global_data_stream.flush();
#endif
    }

}

void backProjection()
{
    if (/*tune*/false)
    {
#define VTUNE_ANA 1

#ifdef VTUNE_ANA
        //
        //
        TileCfg tilecfg;
        tilecfg.choiceTileForKNL(1);
        for (int i = 0; i < 5; i++)
        {
            double time_start = dtime();
            DoBackProjection::backProjectionMain(tilecfg);

            double time_end = dtime();
            std::cout	<<"count : "<<i<<"time : "<<(time_end-time_start)<<std::endl;
        }
#endif
        EXIT_ABNORMALLY;
    }

    //
    TileCfg tilecfg;
    // IS THIS NEEDED HERE?  tilecfg.choiceTileForKNL(1);
    DoBackProjection::backProjectionMain(tilecfg);
}


void maximization(bool update_tau2_with_fsc,bool is_whole_instead_of_half)
{
	NODE0ONLY std::cout << " Maximization ..." << std::endl;

	int first_local_class, last_local_class;
	first_local_class = 0; last_local_class=nr_classes-1;

#ifdef DO_RECONSTRUCT_EACH_NODE
	int nr_local_classes = divide_equally(nr_classes, nodes, node, first_local_class, last_local_class);
#endif

	auto beforeProcessing = [&](int iclass) {
#ifdef DATA_STREAM
		global_data_stream.foutInt(iclass, "maximization()_iclass", __FILE__, __LINE__);
		global_data_stream.foutDouble(mapModel.Irefs[iclass].wptr(), mapModel.Irefs[iclass].dimzyx, "maximization()_mymodel.Iref[iclass]", __FILE__, __LINE__);
		global_data_stream.foutInt(gridding_nr_iter, "maximization()_gridding_nr_iter", __FILE__, __LINE__);
		global_data_stream.foutInt(do_map, "maximization()_do_map", __FILE__, __LINE__);
		global_data_stream.foutDouble(tau2_fudge_factor, "maximization()_mymodel.tau2_fudge_factor", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.tau2_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.tau2_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.sigma2_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.sigma2_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.data_vs_prior_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.fsc_halves_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.fsc_halves_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.wsum_pdf_class.wptrAll()[iclass], "maximization()_wsum_model.pdf_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mapModel.minres_map, "maximization()_minres_map", __FILE__, __LINE__);
		global_data_stream.check(); global_data_stream.flush();
#endif
	};

	auto afterProcessing = [&](int iclass) {
#ifdef DATA_STREAM
		global_data_stream.foutInt(iclass, "maximization()_iclass", __FILE__, __LINE__);
		global_data_stream.foutDouble(mapModel.Irefs[iclass].wptr(), mapModel.Irefs[iclass].dimzyx, "maximization()_mymodel.Iref[iclass]", __FILE__, __LINE__);
		global_data_stream.foutInt(gridding_nr_iter, "maximization()_gridding_nr_iter", __FILE__, __LINE__);
		global_data_stream.foutInt(do_map, "maximization()_do_map", __FILE__, __LINE__);
		global_data_stream.foutDouble(tau2_fudge_factor, "maximization()_mymodel.tau2_fudge_factor", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.tau2_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.tau2_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.sigma2_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.sigma2_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.data_vs_prior_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.fsc_halves_class[iclass].wptr(ori_size/2+1), ori_size/2+1, "maximization()_mymodel.fsc_halves_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mlModel.wsum_pdf_class.wptrAll()[iclass], "maximization()_wsum_model.pdf_class[iclass]", __FILE__, __LINE__);
		global_data_stream.foutDouble(mapModel.minres_map, "maximization()_minres_map", __FILE__, __LINE__);
		global_data_stream.check(); global_data_stream.flush();
#endif
	};

	// Split this list up to get load balancing in the main loop
	//
	std::vector<int> classesToReconstruct;
	for (int iclass = first_local_class; iclass <= last_local_class; iclass++) {
		// TOOD : remember to remove this bug
		// if(mlModel.wsum_pdf_clsss[iclass] > 0.)
		if (mlModel.pdf_class.rptrAll()[iclass] > 0.) {
			classesToReconstruct.push_back(iclass);
		} else {
			beforeProcessing(iclass);
			mapModel.Irefs[iclass].fill(0);
			afterProcessing(iclass);
		}
	}

	// First reconstruct the images for each class
	int const nr_threads_for_reconstruction = classesToReconstruct.size() > 0 ? std::max<int>(1, omp_get_max_threads() / classesToReconstruct.size()) : 1;

	if (debug_flag) {
		NODE0ONLY std::cerr << __FILE__ << ":" << __LINE__
			<< " maximization() classesToReconstruct.size():" << classesToReconstruct.size()
			<< " nr_threads_for_reconstruction:" << nr_threads_for_reconstruction
			<< std::endl;
	}

#pragma omp parallel for
    for (int classesToReconstructIndex = 0; classesToReconstructIndex < classesToReconstruct.size(); classesToReconstructIndex++)
    {
		auto const iclass = classesToReconstruct[classesToReconstructIndex];

		beforeProcessing(iclass);

        // update
        mlModel.sigma2_class[iclass].zero();
        auto sigma2_iclass = mlModel.sigma2_class[iclass].wptr(ori_size/2+1);
        // update or not depend on update_tau2_with_fsc
        auto tau2_iclass = mlModel.tau2_class[iclass].wptr(ori_size/2+1);
        // update or not depend on update_tau2_with_fsc
        mlModel.data_vs_prior_class[iclass].zero();
        auto data_vs_prior_class_iclass = mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1);
        auto fsc_halves_class_iclass = mlModel.fsc_halves_class[iclass].rptr(ori_size/2+1);

        // void reconstruction(int iclass,int max_iter_preweight,bool do_map,double tau2_fudge,double* tau2,double* sigma2,double* data_vs_prior,
        //                     const double* fsc, /* only input*/, double normalise = 1., bool update_tau2_with_fsc = false,
        //                     bool is_whole_instead_of_half = false,int minres_map = -1)

        mapModel.reconstruction(iclass, gridding_nr_iter, do_map, tau2_fudge_factor, tau2_iclass, sigma2_iclass,
                                data_vs_prior_class_iclass, fsc_halves_class_iclass,
                                mlModel.wsum_pdf_class.rptrAll()[iclass], update_tau2_with_fsc, is_whole_instead_of_half, nr_threads_for_reconstruction);

		afterProcessing(iclass);
	}

#if defined(USEMPI) && defined(DO_RECONSTRUCT_EACH_NODE)
    mlModel.broadcastData();
    mapModel.broadcastData();
#endif

	// Then perform the update of all other model parameters
	maximizationOtherParameters();
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(mlModel.avg_norm_correction, "maximizationOtherParameters()_model_avg_norm_correction", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_class.wptrAll()[0], "maximizationOtherParameters()_model_pdf_class[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_class.wptrAll()[nr_classes-1], "maximizationOtherParameters()_model_pdf_class[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.scale_correction.wptr(mlModel.nr_groups), mlModel.nr_groups, "maximizationOtherParameters()_model_scale_correction", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[0].wptr(mlModel.nr_directions), mlModel.nr_directions, "maximizationOtherParameters()_model_pdf_direction[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.pdf_direction[nr_classes-1].wptr(mlModel.nr_directions), mlModel.nr_directions, "maximizationOtherParameters()_model_pdf_direction[nr_classes-1]", __FILE__, __LINE__);
    for (int igroup = 0; igroup < mlModel.nr_groups; igroup++) {
        global_data_stream.foutDouble(mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize), mlModel.ori_Fsize, "maximizationOtherParameters()_model_sigma2_noise_igroup", __FILE__, __LINE__);
    }
    global_data_stream.foutDouble(mlModel.sigma2_offset, "maximizationOtherParameters()_sigma2_offset", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.LL, "maximizationOtherParameters()_model_LL", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.ave_Pmax, "maximizationOtherParameters()_model_ave_Pmax", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.tau2_class[0].wptr(ori_size/2+1), ori_size/2+1, "maximizationOtherParameters()_model_tau2_class[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.tau2_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "maximizationOtherParameters()_model_tau2_class[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[0].wptr(ori_size/2+1), ori_size/2+1, "maximizationOtherParameters()_data_vs_prior_class[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mlModel.data_vs_prior_class[nr_classes-1].wptr(ori_size/2+1), ori_size/2+1, "maximizationOtherParameters()_data_vs_prior_class[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[0].wptr(), mapModel.Irefs[0].dimzyx, "maximizationOtherParameters()_classProjector.Irefs[0]", __FILE__, __LINE__);
    global_data_stream.foutDouble(mapModel.Irefs[nr_classes-1].wptr(), mapModel.Irefs[nr_classes-1].dimzyx, "maximizationOtherParameters()_classProjector.Irefs[nr_classes-1]", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
    
	// Keep track of changes in hidden variables
	hiddenVarMonitor.updateOverallChangesInHiddenVariables();
    
#ifdef DATA_STREAM
    global_data_stream.foutDouble(hiddenVarMonitor.current_changes_optimal_classes, "updateOverallChangesInHiddenVariables()_current_changes_optimal_classes", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.current_changes_optimal_orientations, "updateOverallChangesInHiddenVariables()_current_changes_optimal_orientations", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.current_changes_optimal_offsets, "updateOverallChangesInHiddenVariables()_current_changes_optimal_offsets", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.nr_iter_wo_large_hidden_variable_changes, "updateOverallChangesInHiddenVariables()_nr_iter_wo_large_hidden_variable_changes", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.smallest_changes_optimal_classes, "updateOverallChangesInHiddenVariables()_smallest_changes_optimal_classes", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.smallest_changes_optimal_offsets, "updateOverallChangesInHiddenVariables()_smallest_changes_optimal_offsets", __FILE__, __LINE__);
    global_data_stream.foutDouble(hiddenVarMonitor.smallest_changes_optimal_orientations, "updateOverallChangesInHiddenVariables()_smallest_changes_optimal_orientations", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}


void maximizationOtherParameters()
{
	// Calculate total sum of weights, and average CTF for each class (for SSNR estimation)
	FDOUBLE sum_weight = 0.;
	for (int iclass = 0; iclass < nr_classes; iclass++)
        sum_weight += mlModel.wsum_pdf_class.rptrAll()[iclass];
	// Update average norm_correction
	if (do_norm_correction)
	{
        mlModel.avg_norm_correction = mlModel.wsum_avg_norm_correction / sum_weight;
	}

    auto sum=[&](FDOUBLE* V,int N){
        FDOUBLE s = 0;
        for (int i = 0; i < N; i++) s += V[i];
        return s;
    };

    if ( do_scale_correction && !(iter==1 && do_firstiter_cc) )
    {
        std::vector<FDOUBLE> sorted(mlModel.nr_groups);
        //
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        {
            auto sumXA = sum(mlModel.wsum_signal_product_spectra[igroup].wptr(mlModel.ori_Fsize),mlModel.ori_Fsize);
            auto sumAA = sum(mlModel.wsum_reference_power_spectra[igroup].wptr(mlModel.ori_Fsize),mlModel.ori_Fsize);
            if (sumAA > 0.)
                mlModel.scale_correction.wptrAll()[igroup] = sumXA / sumAA;
            else
                mlModel.scale_correction.wptrAll()[igroup] = 1.;
            sorted[igroup] = mlModel.scale_correction.rptrAll()[igroup];
        }

        // TODO! Avoid extremities in scale estimates, because they lead to catastrophic events and instabilities in refinement
        // Let's exclude anything bigger than 5x the median or smaller than 1/5 of the median...
        // Use the median instead of the mean, because it is much more robust to outliers.
        std::sort(sorted.begin(), sorted.end());
        auto median = sorted[mlModel.nr_groups / 2];

        FDOUBLE avg_scale_correction = 0., nr_part = 0.;
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        {

            if (mlModel.scale_correction.rptrAll()[igroup] > 5. * median)
                mlModel.scale_correction.wptrAll()[igroup] = 5. * median;
            else if (mlModel.scale_correction.rptrAll()[igroup] < median / 5.)
                mlModel.scale_correction.wptrAll()[igroup] =  median / 5.;

            avg_scale_correction += (FDOUBLE)(mlModel.nr_particles_group.rptrAll()[igroup]) * mlModel.scale_correction.rptrAll()[igroup];
            nr_part += (FDOUBLE)(mlModel.nr_particles_group.rptrAll()[igroup]);
        }

        // Constrain average scale_correction to one.
        avg_scale_correction /= nr_part;
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        {
            mlModel.scale_correction.wptrAll()[igroup] /= avg_scale_correction;
//            #define DEBUG_UPDATE_SCALE
#ifdef DEBUG_UPDATE_SCALE
            std::cerr<< "Group "<<igroup+1<<": scale_correction= "<<mlModel.scale_correction[igroup]<<std::endl;
            for (int i = 0; i < mlModel.ori_Fsize; i++)
                if (mlModel.wsum_reference_power_spectra[igroup*mlModel.ori_Fsize+i]> 0.)
                    std::cerr << " i= " << i << " XA= " << mlModel.wsum_signal_product_spectra[igroup*mlModel.ori_Fsize+i]
                    << " A2= " << mlModel.wsum_reference_power_spectra[igroup*mlModel.ori_Fsize+i]
                    << " XA/A2= " << mlModel.wsum_signal_product_spectra[igroup*mlModel.ori_Fsize+i]/mlModel.wsum_reference_power_spectra[igroup*mlModel.ori_Fsize+i] << std::endl;
#endif
        }
    }

	// Update model.pdf_class vector (for each k)
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
        mlModel.pdf_class.wptrAll()[iclass] = mlModel.wsum_pdf_class.rptrAll()[iclass] / sum_weight;

		// for 2D also update priors of translations for each class!
        // None
        // TODO

        for (int idir = 0; idir < sampler3d.NrDir(); idir++)
        {
            mlModel.pdf_direction[iclass][idir] = mlModel.wsum_pdf_direction[iclass][idir] / sum_weight;
        }
	}

    mlModel.sigma2_offset = (mlModel.wsum_sigma2_offset) / (2. * sum_weight);

	// TODO: update estimates for sigma2_rot, sigma2_tilt and sigma2_psi!
    if(!(iter == 1 && do_firstiter_cc))
    {
        for (int igroup = 0; igroup < mlModel.nr_groups; igroup++)
        {
            auto sigma2_noise_igroup = mlModel.sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
            auto wsum_sigma2_noise_igroup = mlModel.wsum_sigma2_noise[igroup].wptr(mlModel.ori_Fsize);
            // Factor 2 because of the 2-dimensionality of the complex-plane
            for (int n = 0; n < mlModel.ori_Fsize; n++)
            {
                sigma2_noise_igroup[n] = wsum_sigma2_noise_igroup[n] / (2. * mlModel.wsum_sumw_group.rptrAll()[igroup] * Npix_per_shell.rptrAll()[n]);
            }
        }
    }


	// After the first iteration the references are always CTF-corrected
    if (do_ctf_correction)
    	refs_are_ctf_corrected = true;

	// Some statistics to output

    mlModel.LL = mlModel.wsum_LL;
    if ((iter==1 && do_firstiter_cc) || do_always_cc)
        mlModel.LL /= sum_weight; // this now stores the average ccf

    mlModel.ave_Pmax = mlModel.wsum_ave_Pmax / sum_weight;

    // After the first, special iteration, apply low-pass filter of -ini_high again
    if (iter == 1 && do_firstiter_cc)
    {
        mapModel.applyLowPassFilter();
        if (ini_high > 0.)
        {
            // Adjust the tau2_class and data_vs_prior_class, because they were calculated on the unfiltered maps
            // This is merely a matter of having correct output in the model.star file (these values are not used in the calculations)
            FDOUBLE radius = ori_size * pixel_size / ini_high;
            radius -= mapModel.width_fmask_edge / 2.;
            FDOUBLE radius_p = radius + mapModel.width_fmask_edge;

            for (int iclass = 0; iclass < nr_classes; iclass++)
            {
                auto tau2_aux = mlModel.tau2_class[iclass].wptr(ori_size/2+1);
                auto data_vs_prior_class_aux = mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1);
                for (int rr = 0; rr < ori_size/2+1; rr++)
                {
                    FDOUBLE r = (FDOUBLE)rr;
                    if (r < radius)
                        continue;
                    else if (r > radius_p)
                    {
                        tau2_aux[rr] = 0.;
                        data_vs_prior_class_aux[rr] = 0.;
                    }
                    else
                    {
                        FDOUBLE raisedcos = 0.5 - 0.5 * cos(rome_pi * (radius_p - r) / mapModel.width_fmask_edge);
                        tau2_aux[rr] *= raisedcos * raisedcos;
                        data_vs_prior_class_aux[rr] *= raisedcos * raisedcos;
                    }
                }
            }
        }

        if (true/*do_generate_seeds*/ && nr_classes > 1)
        {
            // In the first CC-iteration only a single reference was used
            // Now copy this one reference to all K references, for seed generation in the second iteration
            for (int iclass = 1; iclass < nr_classes; iclass++)
            {
                auto tau2_class_iclass = mlModel.tau2_class[iclass].wptr(ori_size/2+1);
                auto data_vs_prior_class_iclass = mlModel.data_vs_prior_class[iclass].wptr(ori_size/2+1);
                for (int i = 0; i < ori_size/2+1; i++) {
                    tau2_class_iclass[i] = mlModel.tau2_class[0][i];
                    data_vs_prior_class_iclass[i] = mlModel.data_vs_prior_class[0][i];
                }
                auto pdf_direction_iclass = mlModel.pdf_direction[iclass].wptr(mlModel.nr_directions);
                for (int i = 0; i < mlModel.nr_directions; i++) {
                    pdf_direction_iclass[i] = mlModel.pdf_direction[0][i];
                }

                auto& Iref_iclass = mapModel.Irefs[iclass];
                for (int i = 0; i < Iref_iclass.dimzyx; i++) {
                    Iref_iclass(0,0,i) = mapModel.Irefs[0](0,0,i);
                }
                mlModel.pdf_class.wptrAll()[iclass] = mlModel.pdf_class.rptrAll()[0] / nr_classes;
            }
            mlModel.pdf_class.wptrAll()[0] /= nr_classes;
        }
    }

}

static IntPerformanceCounter transformFourierArrayToShellIncreasedCoarse_performanceCounter("transformFourierArrayToShellIncreasedCoarse_performanceCounter");
static IntPerformanceCounter transformFourierArrayToShellIncreasedFine_performanceCounter  ("transformFourierArrayToShellIncreasedFine_performanceCounter");

void transformFourierArrayToShellIncreasedCoarse()
{
	transformFourierArrayToShellIncreasedCoarse_performanceCounter.count.v++;
	assert(coarse_size==exp_current_size);
	int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);
#pragma omp parallel for
    for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
		particleModel.shiftImageAssistor.transform(itrans, 0, itrans*exp_nr_over_trans + 0, exp_current_Fsize2, fourierShellTrans, true);
    }
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++) {
        int tid = omp_get_thread_num();
        fourierShellTrans.transformCoarse(tid, particleModel.Fimages_mask_coarse_real[iimage].wptrAll(), exp_current_Fsize2);
        fourierShellTrans.transformCoarse(tid, particleModel.Fimages_mask_coarse_imag[iimage].wptrAll(), exp_current_Fsize2);
        fourierShellTrans.transformCoarse(tid, exp_local_Fctfs[iimage].wptrAll(), exp_current_Fsize2);// TODO : If do not do ctf_correction
        fourierShellTrans.transformCoarse(tid, exp_local_Minvsigma2s[iimage].wptrAll(), exp_current_Fsize2);
    }
}

void transformFourierArrayToShellIncreasedFine()
{
	transformFourierArrayToShellIncreasedFine_performanceCounter.count.v++;
	assert(current_size==exp_current_size);
	const int exp_current_Fsize2 = exp_current_size*(exp_current_size/2+1);

#pragma omp parallel for
    for (int itrans = 0; itrans < exp_nr_trans; itrans++) {
        int tid = omp_get_thread_num();
#if defined(DOUBLE_TRANSLATION)
		particleModel.largeShiftedABTable.transform(itrans, 0, itrans, exp_current_Fsize2, fourierShellTrans, false);
#endif
        for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
			particleModel.shiftImageAssistor.transform(itrans, iover_trans, itrans*exp_nr_over_trans + iover_trans, exp_current_Fsize2, fourierShellTrans, false);
        }
    }

#if defined(TRIPLE_TRANSLATION)
    int exp_nr_single_trans = samplingGrid.exp_positive_shift_index.size();
#pragma omp parallel for
    for (int i_single_trans = 0; i_single_trans < exp_nr_single_trans; i_single_trans++) {
		particleModel.largeXShiftedABTable.transform(i_single_trans, 0, i_single_trans, exp_current_Fsize2, fourierShellTrans, false);
		particleModel.largeYShiftedABTable.transform(i_single_trans, 0, i_single_trans, exp_current_Fsize2, fourierShellTrans, false);
    }
#endif
    
#if defined(DOUBLE_TRANSLATION) || defined(TRIPLE_TRANSLATION)
#pragma omp parallel for
    for (int iover_trans = 0; iover_trans < exp_nr_over_trans; iover_trans++) {
		particleModel.smallShiftedABTable.transform(0, iover_trans, iover_trans, exp_current_Fsize2, fourierShellTrans, false);
    }
#endif
    
#pragma omp parallel for
    for (int iimage = 0; iimage < exp_nr_images; iimage++) {
        int tid = omp_get_thread_num();
        fourierShellTrans.transformFine(tid, particleModel.Fimages_mask_fine_real[iimage].wptrAll(), exp_current_Fsize2);
        fourierShellTrans.transformFine(tid, particleModel.Fimages_mask_fine_imag[iimage].wptrAll(), exp_current_Fsize2);
        fourierShellTrans.transformFine(tid, particleModel.Fimages_nomask_real   [iimage].wptrAll(), exp_current_Fsize2);
        fourierShellTrans.transformFine(tid, particleModel.Fimages_nomask_imag   [iimage].wptrAll(), exp_current_Fsize2);
        fourierShellTrans.transformFine(tid, exp_local_Fctfs                     [iimage].wptrAll(), exp_current_Fsize2);
        fourierShellTrans.transformFine(tid, exp_local_Minvsigma2s               [iimage].wptrAll(), exp_current_Fsize2);
    }
}

void readResult()
{
    if (continue_fn.find("_optimiser.star")!=std::string::npos) {
        std::string model_fn,sampling_fn;
        //
        std::tie(model_fn, sampling_fn) = readFromOptimizer(continue_fn,debug_flag);
        //
        mlModel.readFromModel(model_fn,mapModel,debug_flag);
        //
        sampler3d.readFromSampling(sampling_fn,debug_flag);
    }
    else if(continue_fn.find("_backup.back")) {
        std::string statusFn = continue_fn.substr(0,continue_fn.find("_backup"));
        statusTracer.recoveryStatus(statusFn);
    }
    else{
        ERROR_REPORT("Wrong continue file name "+continue_fn+",use *_optimiser.star or *_backup.back");
    }
}

void checkResult()
{
    // check result
    if(guidance_fn!="NULL" && guidance_fn!="BACK")
        statusTracer.checkStatus(guidance_fn+"_it"+num2str(iter));
    //
    if(guidance_fn=="BACK")
    {
        std::string fn 		= write_path+write_fn;
        std::string iterStr = num2str(iter);
        std::string statusFn= fn+"_it"+iterStr;
        statusTracer.backupStatus(statusFn);
    }
}
    
void writeResult(bool do_write_sampling, bool do_write_data, bool do_write_optimiser, bool do_write_model, bool do_write_mrc, int random_subset)
{
    // write result
    if(write_path != "NULL")
    {
        std::string fn = write_path+write_fn;
        std::string iterStr = num2str(iter);
        if (random_subset == 1 || random_subset == 2) {
            iterStr += "_half"+num2str(random_subset,1);
        }
        else{
            assert(random_subset==-1);
        }
        // Write Classes
        if (do_write_mrc) {
            std::string fn_class    = fn+"_it"+iterStr+"_class";
            mapModel.writeResult(fn_class);
            mlModel.writeOutBild(fn_class, mapModel, sampler3d);
        }
        //
        if (do_write_sampling) {
            std::string fn_sampling = fn+"_it"+iterStr+"_sampling";
            sampler3d.writeOutSampling(fn_sampling);
        }
        // Write Metadata
        if (do_write_data) {
            std::string fn_metadata = fn+"_it"+iterStr+"_data";
            metadata.writeToStar(fn_metadata);
        }
        //
        if (do_write_optimiser) {
            std::string fn_optimiser = fn+"_it"+iterStr+"_optimiser";
            writeOutOptimizer(fn_optimiser);
        }
        //
        if (do_write_model) {
            std::string fn_model = fn+"_it"+iterStr+"_model";
            mlModel.writeOutModel(fn_model, mapModel, metadata, random_subset);
        }
    }
}
    
void printMem(int set_nr_pool)
{
    std::string fn = write_path+write_fn;
    std::string iterStr	= num2str(iter);
    std::string fn_ram 	= fn+"_it"+iterStr+"_ram";
    std::ofstream ramFile;
    ramFile.open((fn_ram+".txt").c_str(), std::ios::out);

#define FOUT ramFile<<std::setw(40)<<std::left

    double perGb = 1/1024./1024./1024.;
    int ori_Fsize = (ori_size/2+1);
    int current_Fsize = (current_size/2+1);
    int current_Fsize2 = current_Fsize*current_size;
    //
    double images_data_size               = nr_local_images*sizeof(double)*ori_size*(ori_size/2+1)*perGb;
    double metadata_size                  = nr_global_images*sizeof(MetaDataElem)*perGb;
    double fix_size = images_data_size + metadata_size;
    FOUT<<"fixed size : "<<std::endl;
    FOUT<<"images_data size : "<<images_data_size<<" GB."<<std::endl;
    FOUT<<"metadata_size : "<<metadata_size<<" GB."<<std::endl;
    FOUT<<"Total : "<<fix_size<<" GB."<<std::endl;
    FOUT<<"------------------------------------------------"<<std::endl;
    mlModel.printSpaceInfo(ramFile);
    mapModel.printSpaceInfo(ramFile);
    //
    double exp_metadata_size             		=	set_nr_pool*sizeof(MetaDataElem)*perGb;
    double exp_imgs_size                 		=	set_nr_pool*ori_size*ori_size*sizeof(double)*perGb;
    double exp_power_imgs_size           		=	set_nr_pool*ori_Fsize*sizeof(double)*perGb;
    double exp_highres_Xi2_imgs_size     		=	set_nr_pool*sizeof(double)*perGb;
    double exp_min_diff2_size            		=	set_nr_pool*sizeof(double)*perGb;
    double exp_old_offsetx_size          		=	set_nr_pool*sizeof(double)*perGb;
    double exp_old_offsety_size					=	set_nr_pool*sizeof(double)*perGb;
    double exp_wsum_scale_correction_XA_size	=	set_nr_pool*ori_Fsize*sizeof(double)*perGb;
    double exp_wsum_scale_correction_AA_size	=	set_nr_pool*ori_Fsize*sizeof(double)*perGb;
    double exp_wsum_norm_correction_size		=	set_nr_pool*sizeof(double)*perGb;
    double exp_Fimgs_size                		= 	2*set_nr_pool*current_size*current_Fsize*sizeof(double)*perGb;
    double exp_Fimgs_nomask_size         		= 	2*set_nr_pool*current_size*current_Fsize*sizeof(double)*perGb;
    double exp_Fctfs_size                		= 	set_nr_pool*current_size*current_Fsize*sizeof(double)*perGb;
    double exp_significant_weight_size   		= 	set_nr_pool*sizeof(double)*perGb;
    double exp_sum_weight_size           		= 	set_nr_pool*sizeof(double)*perGb;
    double exp_Fimgs_shifted_nomask_size		= 	0;
    double exp_local_sqrtXi2_size				=	set_nr_pool*sizeof(double)*perGb;
    double exp_over_psi_size					=	sampler3d.NrPsi()*sampler3d.NrDir()*sampler3d.oversamplingFactorOrientations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_over_rot_size					=	sampler3d.NrPsi()*sampler3d.NrDir()*sampler3d.oversamplingFactorOrientations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_over_tilt_size					=	sampler3d.NrPsi()*sampler3d.NrDir()*sampler3d.oversamplingFactorOrientations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_over_trans_x						=	sampler3d.NrTrans()*sampler3d.oversamplingFactorTranslations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_over_trans_y						=	sampler3d.NrTrans()*sampler3d.oversamplingFactorTranslations(adaptive_oversampling)*sizeof(double)*perGb;
    double exp_Fimgs_shifted_size				=	2*set_nr_pool*sampler3d.NrTrans(adaptive_oversampling)*current_Fsize*sizeof(double)*perGb;
    double exp_local_Fctfs_size					=	set_nr_pool*current_Fsize2*sizeof(double)*perGb;
    double exp_local_Minvsigma2s_size			=	set_nr_pool*current_Fsize2*sizeof(double)*perGb;
    double exp_Rot_significant_size 			=	nr_classes*sampler3d.NrDir()*sampler3d.NrPsi()*sizeof(bool)*perGb;
    double exp_Mcoarse_significant_size			=	set_nr_pool*nr_classes*sampler3d.NrPoints(0)*sizeof(bool)*perGb;
    double exp_Mweight_coarse_size				=	set_nr_pool*nr_classes*sampler3d.NrPoints(0)*sizeof(double)*perGb;
    double exp_Mweight_fine_size				=	set_nr_pool*nr_classes*sampler3d.NrPoints(adaptive_oversampling)*sizeof(double)*perGb;
    //
    double unfixed_size = exp_metadata_size + exp_imgs_size + exp_power_imgs_size + exp_highres_Xi2_imgs_size + exp_min_diff2_size \
    + exp_old_offsetx_size + exp_old_offsety_size + exp_wsum_scale_correction_XA_size \
    + exp_wsum_scale_correction_AA_size + exp_wsum_norm_correction_size + exp_Fimgs_size + exp_Fimgs_nomask_size \
    + exp_Fctfs_size + exp_significant_weight_size + exp_sum_weight_size + exp_Fimgs_shifted_nomask_size \
    + exp_local_sqrtXi2_size + exp_over_psi_size + exp_over_rot_size + exp_over_tilt_size + exp_over_trans_x \
    + exp_over_trans_y + exp_Fimgs_shifted_size + exp_local_Fctfs_size + exp_local_Minvsigma2s_size + exp_Rot_significant_size \
    + exp_Mcoarse_significant_size + exp_Mweight_coarse_size + exp_Mweight_fine_size;
    int sugg_set_nr_pool = (6.5 - fix_size)/(unfixed_size/set_nr_pool);
    FOUT<<"------------------------------------------------------"<<std::endl;
    FOUT<<"suggestion nr_pool : "<<sugg_set_nr_pool<<std::endl;
    FOUT<<"real nr_pool : "<<set_nr_pool<<std::endl;
    FOUT<<"------------------------------------------------------"<<std::endl;
    // nr_pool = set_nr_pool;
    FOUT<<"unfixed size(Expectation step) : "<<std::endl;
    FOUT<<"exp_metadata_size : "<<exp_metadata_size<<" GB."<<std::endl;
    FOUT<<"exp_imgs_size : "<<exp_imgs_size<<" GB."<<std::endl;
    FOUT<<"exp_power_imgs_size : "<<exp_power_imgs_size<<" GB."<<std::endl;
    FOUT<<"exp_highres_Xi2_imgs_size : "<<exp_highres_Xi2_imgs_size<<" GB."<<std::endl;
    FOUT<<"exp_min_diff2_size : "<<exp_min_diff2_size<<" GB."<<std::endl;
    FOUT<<"exp_old_offsetx_size : "<<exp_old_offsetx_size<<" GB."<<std::endl;
    FOUT<<"exp_old_offsety_size : "<<exp_old_offsety_size<<" GB."<<std::endl;
    FOUT<<"exp_wsum_scale_correction_XA_size : "<<exp_wsum_scale_correction_XA_size<<" GB."<<std::endl;
    FOUT<<"exp_wsum_scale_correction_AA_size : "<<exp_wsum_scale_correction_AA_size<<" GB."<<std::endl;
    FOUT<<"exp_wsum_norm_correction_size  	: "<<exp_wsum_norm_correction_size<<" GB."<<std::endl;
    FOUT<<"exp_Fimgs_size : "<<exp_Fimgs_size<<" GB."<<std::endl;
    FOUT<<"exp_Fimgs_nomask_size : "<<exp_Fimgs_nomask_size<<" GB."<<std::endl;
    FOUT<<"exp_Fctfs_size : "<<exp_Fctfs_size<<" GB."<<std::endl;
    FOUT<<"exp_significant_weight_size : "<<exp_significant_weight_size<<" GB."<<std::endl;
    FOUT<<"exp_sum_weight_size : "<<exp_sum_weight_size<<" GB."<<std::endl;
    FOUT<<"exp_Fimgs_shifted_nomask_size : "<<exp_Fimgs_shifted_nomask_size<<" GB."<<std::endl;
    FOUT<<"exp_local_sqrtXi2_size : "<<exp_local_sqrtXi2_size<<" GB."<<std::endl;
    FOUT<<"exp_over_psi_size : "<<exp_over_psi_size<<" GB."<<std::endl;
    FOUT<<"exp_over_rot_size : "<<exp_over_rot_size<<" GB."<<std::endl;
    FOUT<<"exp_over_tilt_size : "<<exp_over_tilt_size<<" GB."<<std::endl;
    FOUT<<"exp_over_trans_x : "<<exp_over_trans_x<<" GB."<<std::endl;
    FOUT<<"exp_over_trans_y : "<<exp_over_trans_y<<" GB."<<std::endl;
    FOUT<<"exp_Fimgs_shifted_size : "<<exp_Fimgs_shifted_size<<" GB."<<std::endl;
    FOUT<<"exp_local_Fctfs_size : "<<exp_local_Fctfs_size<<" GB."<<std::endl;
    FOUT<<"exp_local_Minvsigma2s_size : "<<exp_local_Minvsigma2s_size<<" GB."<<std::endl;
    FOUT<<"exp_Rot_significant_size : "<<exp_Rot_significant_size<<" GB."<<std::endl;
    FOUT<<"exp_Mcoarse_significant_size : "<<exp_Mcoarse_significant_size<<" GB."<<std::endl;
    FOUT<<"exp_Mweight_coarse_size : "<<exp_Mweight_coarse_size<<" GB."<<std::endl;
    FOUT<<"exp_Mweight_fine_size : "<<exp_Mweight_fine_size<<" GB."<<std::endl;
    FOUT<<"Total : "<<unfixed_size<<" GB."<<std::endl;
    FOUT<<"------------------------------------------------------"<<std::endl;
    //
    double thread_Frefctf_size 					=	2*maxthreads*current_Fsize2*sizeof(double)*perGb;
    double thread_Fimg_nomask_size				=	2*maxthreads*ori_size*ori_size*sizeof(double)*perGb;
    double thread_Fimg_size						=	2*maxthreads*ori_size*ori_size*sizeof(double)*perGb;
    double thread_Fweight_size					=	maxthreads*ori_size*ori_size*sizeof(double)*perGb;
    double thread_max_weight_size				=	maxthreads*nr_pool*sizeof(double)*perGb;
    double thread_wsum_norm_correction_size		=	maxthreads*nr_pool*sizeof(double)*perGb;
    double thread_wsum_sigma2_noise_size		=	maxthreads*nr_pool*current_Fsize2*sizeof(double)*perGb;
    double thread_wsum_pdf_direction_size		=	maxthreads*nr_classes*exp_nr_dir*sizeof(double)*perGb;
    double thread_wsum_pdf_class_size			=	maxthreads*nr_classes*sizeof(double)*perGb;
    double thread_wsum_scale_correction_XA_size	=	maxthreads*nr_pool*current_Fsize2*sizeof(double)*perGb;
    double thread_wsum_scale_correction_AA_size	=	maxthreads*nr_pool*current_Fsize2*sizeof(double)*perGb;
    double thread_sumw_group_size				=	maxthreads*nr_pool*sizeof(double)*perGb;
    double thread_size = thread_Frefctf_size + thread_Fimg_nomask_size + thread_Fimg_size + thread_Fweight_size \
    + thread_max_weight_size + thread_wsum_norm_correction_size + thread_wsum_sigma2_noise_size \
    + thread_wsum_pdf_direction_size + thread_wsum_pdf_class_size + thread_wsum_scale_correction_XA_size \
    + thread_wsum_scale_correction_AA_size + thread_sumw_group_size;
    FOUT<<"thread size(Expectation step) : "<<std::endl;
    FOUT<<"thread_Frefctf_size : "<<thread_Frefctf_size<<" GB."<<std::endl;
    FOUT<<"thread_Fimg_nomask_size : "<<thread_Fimg_nomask_size<<" GB."<<std::endl;
    FOUT<<"thread_Fimg_size : "<<thread_Fimg_size<<" GB."<<std::endl;
    FOUT<<"thread_Fweight_size : "<<thread_Fweight_size<<" GB."<<std::endl;
    FOUT<<"thread_max_weight_size : "<<thread_max_weight_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_norm_correction_size : "<<thread_wsum_norm_correction_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_sigma2_noise_size : "<<thread_wsum_sigma2_noise_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_pdf_direction_size : "<<thread_wsum_pdf_direction_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_pdf_class_size : "<<thread_wsum_pdf_class_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_scale_correction_XA_size : "<<thread_wsum_scale_correction_XA_size<<" GB."<<std::endl;
    FOUT<<"thread_wsum_scale_correction_AA_size : "<<thread_wsum_scale_correction_AA_size<<" GB."<<std::endl;
    FOUT<<"thread_sumw_group_size : "<<thread_sumw_group_size<<" GB."<<std::endl;
    FOUT<<"Total : "<<thread_size<<" GB."<<std::endl;
    FOUT<<"------------------------------------------------------"<<std::endl;

    ramFile.close();
}


void debugStoreWeightedSums(){

    std::cout << " WARNING: norm_correction : "<< exp_metadata[0].NORM  << " for particle " << 0 <<std::endl;
    std::cout << " mymodel.current_size : " << current_size << " mymodel.ori_size= " << ori_size <<std::endl;
    std::cout << " coarse_size : " << coarse_size << std::endl;
    std::cout << " DIRECT_A2D_ELEM(exp_metadata2, my_image_no, exp_nr_imagas-1) : " <<exp_metadata[exp_nr_images-1].NORM << std::endl;
    // std::cout << " mymodel.avg_norm_correction : " << model_avg_norm_correction << std::endl;
    std::cout << " exp_wsum_norm_correction[ipart] : " << exp_wsum_norm_correction.rptrAll()[0] << std::endl;
    // std::cout << " old_norm_correction : " << old_norm_correction << std::endl;
    // std::cout << " wsum_model.avg_norm_correction : " << wsum_avg_norm_correction << std::endl;
    // std::cout << " group_id : " << group_id << " mymodel.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    // std::cout << " mymodel.sigma2_noise[group_id = 0] : " << model_sigma2_noise[0] << std::endl;
    // std::cout << " wsum_model.sigma2_noise[group_id = 0] : " << wsum_sigma2_noise[0] << std::endl;
    std::cout << " exp_power_imgs[my_image_no = 0] : " << exp_power_imgs[0][0] << std::endl;
    // std::cout << " exp_wsum_scale_correction_XA[ipart] : " << exp_wsum_scale_correction_XA[ipart] << " exp_wsum_scale_correction_AA[ipart] : " << exp_wsum_scale_correction_AA[ipart] << std::endl;
    // std::cout << " wsum_model.wsum_signal_product_spectra[group_id] : " << wsum_model.wsum_signal_product_spectra[group_id] << " wsum_model.wsum_reference_power_spectra[group_id] : " << wsum_model.wsum_reference_power_spectra[group_id] << std::endl;
    std::cout << " exp_min_diff2[ipart = 0] : " << exp_min_diff2.rptrAll()[0] << std::endl;
    // std::cout << " ml_model.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    std::cout << " exp_significant_weight[ipart = 0] : " << exp_significant_weight.rptrAll()[0] << std::endl;
    // std::cout << " exp_max_weight[ipart = 0] : " << exp_max_weight[0] << std::endl;

    std::cerr << " part_id : " << 0 << std::endl;
    std::cerr << " ipart : " << 0 << std::endl;
    std::cerr << " exp_min_diff2[ipart = 0] : " << exp_min_diff2.rptrAll()[0] << std::endl;
    // std::cerr << " logsigma2 : " << logsigma2 << std::endl;
    int group_id = 0;//mydata.getGroupId(part_id, 0);
    std::cerr << " group_id : " << group_id << std::endl;
    // std::cerr << " ml_model.scale_correction[group_id] : " << mymodel.scale_correction[group_id] << std::endl;
    std::cerr << " exp_significant_weight[ipart = 0] : " << exp_significant_weight.rptrAll()[0] << std::endl;
    // std::cerr << " exp_max_weight[ipart = 0]= " << exp_max_weight[0] << std::endl;
    // std::cerr << " ml_model.sigma2_noise[group_id = 0] : " << model_sigma2_noise[0] << std::endl;
}

} // end namespace MLoptimizer
