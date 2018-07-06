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

#define MAP2DOPTIMIZER_BASE_BEING_COMPILED
#include "../src/gtm_optimizer_kernel.h"
#include "../src/map2d_optimizer.h"
#include "../src/map2d_optimizer_kernel.h"
#include "../src/time.h"
#include "../src/tandemExecution.h"
#include "../src/option.h"

#if defined(INTEL_PGXSDK_USED)
#include "../intel_pgx_prealpha/release_win32-x86_64_cl_12.0/pgxsdk_1.0.0/include/pgxsdk1/pgx_sdk.h"
#endif


#include <cstdio>
#include <iostream>
int myrandom (int i) { return std::rand()%i;}

void MLProgram() {
    
    Map2dOptimizer::setupMap2dOptimizer();
    Map2dOptimizer::readImages();
	TUNING_FLUSH
    Map2dOptimizer::prepare();
	TUNING_FLUSH
    Map2dOptimizer::iterate();
	TUNING_FLUSH
    Map2dOptimizer::writeClassesAndMetadata();
	TUNING_FLUSH
    Map2dOptimizer::destroyMap2dOptimizer();
    
}


int main(int argc, char * argv[]) {

    Option option;
    // general option
    option.addOption("-i",                  "Input metadata file with images needed to align"                               					);
    option.addOption("-o",                  "Output metadata"                                                               					);
    option.addOption("-K",                  "Number of classes needed to classify!!!!! please use Uppercase —K instead of lowercase -k !!!!!" 	);
    option.addOption("-iter",               "Maximum number of iterations to perform"                                       					);
    option.addOption("-angpix",             "Pixel size (in Angstroms)!!!!! please use -angpix instead of -pixel !!!!!"        					);
    option.addOption("-particle_diameter",  "Particle_diameter",													   					"-1"	);
    option.addOption("-tau2_fudge",         "Regularisation parameter (values higher than 1 give more weight to the data,2 for 2D)",	"2"		);
    option.addOption("-oversampling",       "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)","1"		);
    // advanced option
    option.addOption("-pool",               "Number of images to be processed together for each EM step",              					"50"	);
    option.addOption("-random_seed",        "Number of the random seed generator",                                     					"33"	);
    option.addOption("-offset_step",        "Sampling rate (before oversampling) for origin offsets (in pixels)",      					"2" 	);
    option.addOption("-offset_range",       "Search range for origin offsets (in pixels)",                             					"10"	);
    option.addOption("-psi_step",           "Sampling rate (before oversampling) for the in-plane angle",             				 	"10"	);
    //
    option.addOption("-ctf",                "Perform CTF correction?",                                                  		        "0" 	);
    option.addOption("-only_flip_phases",	"Only perform CTF phase-flipping? (default is full amplitude-correction)",					"0"		);
    option.addOption("-ctf_phase_flipped",	"Have the data been CTF phase-flipped?",													"0"		);
    option.addOption("-ctf_corrected_ref", 	"Have the input references been CTF-amplitude corrected?",									"0"		);
    option.addOption("-norm",               "Do norm correction for image data",                                        				"0" 	);
    //
    option.addOption("-zero_mask", 			"Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)","0");
    option.addOption("-flatten_solvent", 	"Perform masking on the references as well",												"0"		);
    // testing option
    option.addOption("-datastream_in", 		"data stream for comparing",																"NULL"	);
    option.addOption("-datastream_out", 	"write data stream to file",																"NULL"  );
    option.addOption("-continue",           "continue work from specified iter.(remember to change -i input star)",    					"0"  	);
    option.addOption("-testguidance",       "keep each iteration same as guidefile",                                 					"NULL"	);
    option.addOption("-unittest",           "ml2d unit testing (test must be compiled in) : \n \
               Map2dOptimizer_Kernel_correctness,Map2dOptimizer_Kernel_performance \n \
               GTM_Optimizer_Kernel_correctness,GTM_Optimizer_Kernel_performance  ",          			  		  						"no_test");
    // performance tuning options
	option.addOption("-kernel3x3",          "use the kernel 3x3 optimization 0,1 ",							       	   					"0"		);
#if defined(TUNING)
	option.addOption("-checker_ftrace",		"writes a ftrace file showing execution times (must be compiled in) : \n",					""   	);
#endif

    if (argc < 3) {
        option.printHelp();
        std::cerr << "example:" << std::endl;
        std::cerr << "rome_ml2d  -i ../../smalldataset/data8_160/data8_160 -o c:/temp/rome_ml2d_classes -n 100 -K 20 -angpix 4. -iter 10 -pool  30 >  ../../smalldataset/data8_160/emit_test_output_stab.txt" << std::endl;
        std::cerr << std::endl;
        EXIT_ABNORMALLY;
    }
    
    option.readCommandLine(argc, argv);
    
#if defined(TUNING)
	{
		std::string ftrace = option.getStrOption("-checker_ftrace");
		if (ftrace != "") {
			Tuning::setFTraceFilename(ftrace.c_str());
			TUNING_SCOPE_STEP_BEGIN(setFTraceFilename_set)
			TUNING_SCOPE_STEP_END
			TUNING_SCOPE_STEP_BEGIN(setFTraceFilename_sleep)
			sleepInSecs(0.5);
			TUNING_SCOPE_STEP_END
			if (false) {
				std::cerr << "exit(0) after TUNING_SCOPE_END" << std::endl;
				exit(0);
			}
		}
	}
#endif

	{
		int kernel3x3 = option.getIntOption("-kernel3x3");
		Map2dOptimizer_Kernel::setKernel3x3(0 != kernel3x3);
	}

	std::string unittest = option.getStrOption("-unittest");
    
	if (unittest != "no_test") {
		if (unittest == "Map2dOptimizer_Kernel_correctness") {
			Map2dOptimizer_Kernel::unitTestCorrectness();
			return 0;
		}
		if (unittest == "Map2dOptimizer_Kernel_performance") {
			Map2dOptimizer_Kernel::unitTestPerformance();
			return 0;
		}
		if (unittest == "GTM_Optimizer_Kernel_correctness") {
			GTM_Optimizer_Kernel::unitTestCorrectness();
			return 0;
		}
		if (unittest == "GTM_Optimizer_Kernel_performance") {
			GTM_Optimizer_Kernel::unitTestPerformance();
			return 0;
		}
		if (unittest == "shiftImageInFourierTransformUnitTestCorrectness") {
			shiftImageInFourierTransformUnitTestCorrectness();
			return 0;
		}
		if (unittest == "shiftImageInFourierTransformUnitTestPerformance") {
			shiftImageInFourierTransformUnitTestPerformance();
			return 0;
		}
	}
    
#ifdef USEMPI
    MPI::Init();
    if(MPI::COMM_WORLD.Get_rank()  == 0)
        option.printValue();
#else
    option.printValue();
#endif
    
#if defined(INTEL_PGXSDK_USED)
	{
		pgxsdk1::Stats stats;
		pgxsdk1::getStats(stats);
		std::cerr << "stats.currentTime:" << stats.currentTime << std::endl;
		static bool pinAttached(false);
		std::string pinfolder = "D:\\bevin\\pgx\\_install\\release_win32-x86_64_cl_12.0\\bin\\";
		pgxsdk1::attachPin2(pinAttached, pinfolder);
		std::cerr << "pgxsdk1::attachPin pinAttached:" << pinAttached << std::endl;
	}
#endif

    double t1 = dtime();
    //
    std::string star_fn						=	pathRemoveSuffix(option.getStrOption("-i"))+".star";
    MapOptimizer_base_new::star_fn			=	star_fn;
    MapOptimizer_base_new::mrcs_dir			=	pathGetDir(star_fn);
    std::string output_fn					=   pathRemoveSuffix(option.getStrOption("-o"));
    MapOptimizer_base_new::write_fn			=   pathGetFilename(output_fn);
    MapOptimizer_base_new::write_path			=	pathGetDir(output_fn);
    MapOptimizer_base_new::ref_fn				=   "NULL";
    //
    MapOptimizer_base_new::nr_classes			=	option.getIntOption("-K");
    MapOptimizer_base_new::nr_pool			=	option.getIntOption("-pool");
    MapOptimizer_base_new::pixel_size			=	option.getFloatOption("-angpix");
    MapOptimizer_base_new::nr_iter			=	option.getIntOption("-iter");
    MapOptimizer_base_new::random_seed		=	option.getIntOption("-random_seed");
    MapOptimizer_base_new::offset_step 		=	option.getIntOption("-offset_step");
    MapOptimizer_base_new::offset_range		=	option.getIntOption("-offset_range");
    MapOptimizer_base_new::psi_step			=	option.getIntOption("-psi_step");
    MapOptimizer_base_new::tau2_fudge_factor  =   option.getFloatOption("-tau2_fudge");
    MapOptimizer_base_new::particle_diameter  =   option.getFloatOption("-particle_diameter");
    MapOptimizer_base_new::adaptive_oversampling=  option.getIntOption("-oversampling");
    //
    MapOptimizer_base_new::do_ctf_correction 	=   option.getBoolOption("-ctf");
    MapOptimizer_base_new::only_flip_phases	=	option.getBoolOption("-only_flip_phases");
    MapOptimizer_base_new::ctf_phase_flipped	=	option.getBoolOption("-ctf_phase_flipped");
    MapOptimizer_base_new::do_norm_correction	=   option.getBoolOption("-norm");
    MapOptimizer_base_new::do_zero_mask	  	=	option.getBoolOption("-zero_mask");
    MapOptimizer_base_new::do_solvent			=	option.getBoolOption("-flatten_solvent");
    //
    MapOptimizer_base_new::data_stream_in_fn	=	option.getStrOption("-datastream_in");
    MapOptimizer_base_new::data_stream_out_fn	=	option.getStrOption("-datastream_out");
    Map2dOptimizer::set_iter(option.getIntOption("-continue"));
    MapOptimizer_base_new::guidance_fn 	= option.getStrOption("-testguidance");
    
    MLProgram();
    
    double t2 = dtime();
    
    if (MapOptimizer_base_new::node  == 0)
        std::cout<<"ML2D costs : "<<(t2-t1)<<std::endl;
    
#ifdef USEMPI
    MPI::Finalize();
#endif

    return 0;
}
