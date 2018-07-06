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

#include "../src/util.h"		// used for building precompiled headers on Windows

#include "../src/checker.h"
#include "../src/exp_model.h"
#include "../src/findBestPoints.h"
#include "../src/map_model.h"
#include "../src/map3d_optimizer.h"
#include "../src/map3d_optimizer_new.h"
#include "../src/map3d_optimizer_kernel.h"
#include "../src/option.h"

#include "../src/map3d_optimizer_new.h"
#define Map3dOptimizer Map3dOptimizer_new

//#include "../src/map3d_autorefinement.h"
//#define Map3dOptimizer Map3dAutoRefinement

int main(int argc, char * argv[]) {
    Option option;
	// 
	option.addOption("-@",                  "Read the options, one per line, from the specified file",                                                        "NULL");
	// general option
    option.addOption("-i",                  "Input metadata file with images needed to align",                                                                "NULL");
    option.addOption("-continue",        	"continue work from specified file.( from *_optimiser.star or *_backup.back file)",    						   	  "NULL");
    option.addOption("-o",                  "Output metadata"                                                                                                       );
    option.addOption("-ref",                "3d reference file name(*.mrc)"                                                                                         );
    option.addOption("-particle_diameter",  "Particle_diameter"                                                                                                     );
    option.addOption("-K",                  "Number of classes needed to classify!!!!! please use Uppercase —K instead of lowercase -k !!!!!" 						);
    option.addOption("-iter",               "Maximum number of iterations to perform"                                                                               );
    option.addOption("-angpix",             "Pixel size (in Angstroms)!!!!! please use -angpix instead of -pixel !!!!!"									            );
    option.addOption("-ini_high",           "Resolution (in Angstroms) to which to limit refinement in the first iteration"                                         );
    option.addOption("-oversampling",       "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)"                             );
    option.addOption("-healpix_order",      "Healpix order for the angular sampling (before oversampling) on the (3D) sphere: hp2=15deg, hp3=7.5deg, etc"           );
    option.addOption("-tau2_fudge",         "Regularisation parameter (values higher than 1 give more weight to the data)"                                          );
    option.addOption("-sym",                "Symmetry group"						                                                                                );
    // advanced option
    option.addOption("-pool",               "Particles processing together each EM expectation.",                                                              "50" );
    option.addOption("-random_seed",        "Seed of randomly shuffle the image data and noise generation",                                                    "33" );
    option.addOption("-offset_step",        "The offset step of image shift searching",                                                                        "2"  );
    option.addOption("-offset_range",       "The offset range of image shift searching,(-10~10)",                                                              "10" );
    //
    option.addOption("-ctf",                "Perform CTF correction?",                                                  		                               "0" 	);
    option.addOption("-only_flip_phases",	"Only perform CTF phase-flipping? (default is full amplitude-correction)",										   "0"	);
    option.addOption("-ctf_phase_flipped",	"Have the data been CTF phase-flipped?",																		   "0"	);
    option.addOption("-ctf_corrected_ref", 	"Have the input references been CTF-amplitude corrected?",														   "0"	);
    option.addOption("-intact_ctf_first_peak", "Ignore CTFs until their first peak?",																		   "0"	);
    option.addOption("-firstiter_cc",       "Perform CC-calculation in the first iteration (use this if references are not on the absolute intensity scale)",  "0"  );
    option.addOption("-always_cc",          "Perform CC-calculation in all iterations (useful for faster denovo model generation?)",                           "0"  );
    option.addOption("-scale",              "Perform intensity-scale corrections on image groups?",                                                            "0"  );
    option.addOption("-norm",               "Perform normalisation-error correction?",                                                                         "0"  );
    option.addOption("-zero_mask",          "Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)",      "0"  );
    option.addOption("-flatten_solvent",    "Perform masking on the references as well?",                                                                      "0"  );
    option.addOption("-solvent_mask", 		"User-provided mask for the references (default is to use spherical mask with particle_diameter)",				  "NULL");
    option.addOption("-no_map",             "Do not use map estimitor",                                                                                        "0"  );
    option.addOption("-sigma_ang", 			"Stddev on all three Euler angles for local angular searches (of +/- 3 stddev)",								   "0"	);
    
    option.addOption("-datastream_in", 		"data stream for comparing",																					  "NULL");
    option.addOption("-datastream_out", 	"write data stream to file",																					  "NULL");
    option.addOption("-testguidance",       "keep each iteration same as guidefile(use 'BACK' to backup the data)",     									  "NULL");
    option.addOption("-debug",     			"print out debug info",                                                                                       		""  );
    
	option.addOption("-checker_goSerial",	"When to force serial execution - used for debugging nondeterminism", "");
	option.addOption("-checker_goParallel", "When to return to parallel execution - used for debugging nondeterminism", "");
	option.addOption("-checker_i",          "File to compare against",                                                                                         ""   );
    option.addOption("-checker_o",          "File to write for future comparisons",                                                                            ""   );
    option.addOption("-checker_test",       "The string describing the test being measured",                                                 "default_checker_test" );
    option.addOption("-checker_ftrace",     "File to write an ftrace to",                                                                                       ""  );
    

    
    if (argc < 3) {
        option.printHelp();
        std::cerr << "example:" << std::endl;
        std::cerr << "./bin/rome_map3d -i dataset/class19.star -o dataset/class19_ml2d -K 10 -iter 30 -angpix 1.74 > dataset/ml2d_output.txt" << std::endl;
        std::cerr << std::endl;
		EXIT_ABNORMALLY;
    }
    
    option.readCommandLine(argc, argv);

	std::string optionInclude = option.getStrOption("-@");
	if (optionInclude != "NULL") {
		option.readIncludeFile(optionInclude);
	}

	std::string checker_goSerial   = option.getStrOption("-checker_goSerial");
	std::string checker_goParallel = option.getStrOption("-checker_goParallel");
	if (checker_goSerial  .size() > 0) setGoSerial  (checker_goSerial);
	if (checker_goParallel.size() > 0) setGoParallel(checker_goParallel);
	maybeGoSerial("immediately");

    std::string checker_i_fn		=   option.getStrOption("-checker_i");
    std::string checker_o_fn		=   option.getStrOption("-checker_o");
    std::string checker_test		=   option.getStrOption("-checker_test");
	std::string checker_ftrace		=   option.getStrOption("-checker_ftrace");
	Checker::Benchmark* benchmark = nullptr;
	if (checker_i_fn.size() > 0 || checker_o_fn.size() > 0) {
		benchmark = sNew(Checker::Benchmark);
		benchmark->setFiles(checker_i_fn, checker_o_fn);
		benchmark->setDefaultColumnHeader(checker_test);
	}
	if (checker_ftrace.size() > 0) {
		Checker::setFtraceFile(checker_ftrace);
	}

#ifdef USEMPI
    MPI::Init();
    if(MPI::COMM_WORLD.Get_rank()  == 0)
        option.printValue();
#else
    option.printValue();
#endif
    
    double t1 = dtime();
    
    std::string set_star_fn         =   pathRemoveSuffix(option.getStrOption("-i"))+".star";
    std::string output_fn           =   pathRemoveSuffix(option.getStrOption("-o"));
    
    // set parameters
    Map3dOptimizer::star_fn					=   set_star_fn;
    Map3dOptimizer::write_path				=   pathGetDir(output_fn);
    Map3dOptimizer::write_fn				=   pathGetFilename(output_fn);
    Map3dOptimizer::ref_fn					=   option.getStrOption("-ref");
    Map3dOptimizer::nr_classes				=   option.getIntOption("-K");
    Map3dOptimizer::nr_pool					=   option.getIntOption("-pool");
    Map3dOptimizer::pixel_size				=   option.getFloatOption("-angpix");
    Map3dOptimizer::nr_iter					=   option.getIntOption("-iter");
    Map3dOptimizer::random_seed				=   option.getIntOption("-random_seed");
    Map3dOptimizer::offset_step				=   option.getIntOption("-offset_step");
    Map3dOptimizer::offset_range			=   option.getIntOption("-offset_range");
    
    Map3dOptimizer::do_ctf_correction		=   option.getBoolOption("-ctf");
    Map3dOptimizer::only_flip_phases		=	option.getBoolOption("-only_flip_phases");
    Map3dOptimizer::ctf_phase_flipped		=	option.getBoolOption("-ctf_phase_flipped");
    Map3dOptimizer::refs_are_ctf_corrected	=	option.getBoolOption("-ctf_corrected_ref");
    Map3dOptimizer::intact_ctf_first_peak	=	option.getBoolOption("-intact_ctf_first_peak");
    
    Map3dOptimizer::ini_high				=   option.getFloatOption("-ini_high");
    Map3dOptimizer::adaptive_oversampling	=   option.getIntOption("-oversampling");
    Map3dOptimizer::sampler3d_healpix_order	=   option.getIntOption("-healpix_order");
    Map3dOptimizer::tau2_fudge_factor		=   option.getFloatOption("-tau2_fudge");
    Map3dOptimizer::sampler3d_fn_sym		=   option.getStrOption("-sym");
    Map3dOptimizer::particle_diameter		=   option.getFloatOption("-particle_diameter");
    Map3dOptimizer::do_firstiter_cc			=   option.getBoolOption("-firstiter_cc");
    Map3dOptimizer::do_always_cc			=   option.getBoolOption("-always_cc");
    Map3dOptimizer::do_scale_correction		=   option.getBoolOption("-scale");
    Map3dOptimizer::do_norm_correction		=   option.getBoolOption("-norm");
    Map3dOptimizer::do_zero_mask			=   option.getBoolOption("-zero_mask");
    Map3dOptimizer::do_solvent             	=   option.getBoolOption("-flatten_solvent");
    Map3dOptimizer::do_map                 	=   !option.getBoolOption("-no_map");
    Map3dOptimizer::sigma2_angle			=	option.getFloatOption("-sigma_ang");
    Map3dOptimizer::mask_fn					=	option.getStrOption("-solvent_mask");
    Map3dOptimizer::do_shifts_onthefly		=	true;
    
    Map3dOptimizer::data_stream_in_fn		=	option.getStrOption("-datastream_in");
    Map3dOptimizer::data_stream_out_fn		=	option.getStrOption("-datastream_out");
    Map3dOptimizer::continue_fn				=	option.getStrOption("-continue");
    Map3dOptimizer::guidance_fn       		=   option.getStrOption("-testguidance");
    Map3dOptimizer::debug_flag				=	option.getBoolOption("-debug");
    
    Map3dOptimizer::setupMLoptimizer();
    
    Map3dOptimizer::prepare();
    
    Map3dOptimizer::iterate();
    
    Map3dOptimizer::destroyMLoptimizer();

    
    double t2 = dtime();
    
    if(Map3dOptimizer::node  == 0)
        std::cout<<"ML2D costs : "<<(t2-t1)
#ifdef DOUBLE_TRANSLATION
			<< " using DOUBLE_TRANSLATION"
#endif
#ifdef TRIPLE_TRANSLATION
			<< " using TRIPLE_TRANSLATION"
#endif
			<< std::endl;
    
	if(benchmark) sDelete(benchmark);

#ifdef USEMPI
    MPI::Finalize();
#endif
    
    return 0;
}
