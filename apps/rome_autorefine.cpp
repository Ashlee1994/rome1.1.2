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
#include "../src/option.h"
#include "../src/map3d_autorefinement.h"

int main(int argc, char * argv[]) {
    Option option;
	// 
	option.addOption("-@",                  "Read the options, one per line, from the specified file",                                                        "NULL");
	// general option
    option.addOption("-i",                  "Input metadata file with images needed to align"                                                                       );
    option.addOption("-o",                  "Output metadata"                                                                                                       );
    option.addOption("-ref",                "3d reference file name(*.mrc)"                                                                                         );
    option.addOption("-particle_diameter",  "Particle_diameter"                                                                                                     );
    option.addOption("-low_resol_join_halves", "Resolution (in Angstrom) up to which the two random half-reconstructions will not be independent to prevent diverging orientations");
    // option.addOption("-K",                  "Number of classes needed to classify!!!!! please use Uppercase —K instead of lowercase -k !!!!!" 						);
    // option.addOption("-iter",               "Maximum number of iterations to perform"                                                                               );
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
    option.addOption("-split_random_halves", "Refine two random halves of the data completely separately",													   "0"	);
    option.addOption("-auto_local_healpix_order", "Minimum healpix order (before oversampling) from which autosampling procedure will use local searches",	   "4"	);
    //
    option.addOption("-datastream_in", 		"data stream for comparing",																					  "NULL");
    option.addOption("-datastream_out", 	"write data stream to file",																					  "NULL");
    option.addOption("-continue",        	"continue work from specified file.( from *_optimiser.star or *_backup.back file)",    						   	  "NULL");
    option.addOption("-testguidance",       "keep each iteration same as guidefile(use 'BACK' to backup the data)",     									  "NULL");
    
	option.addOption("-checker_goSerial",	"When to force serial execution - used for debugging nondeterminism", "");
	option.addOption("-checker_goParallel", "When to return to parallel execution - used for debugging nondeterminism", "");
	option.addOption("-checker_i",          "File to compare against",                                                                                         ""   );
    option.addOption("-checker_o",          "File to write for future comparisons",                                                                            ""   );
    option.addOption("-checker_test",       "The string describing the test being measured",                                                 "default_checker_test" );
    option.addOption("-checker_ftrace",     "File to write an ftrace to",                                                                                       ""  );
    
    if (argc < 3) {
        option.printHelp();
        std::cerr << "example:" << std::endl;
        std::cerr << "./bin/rome_autorefine ....." << std::endl;
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
    
    Map3dAutoRefinement::do_auto_refine			=	true;
    
    // set parameters
    Map3dAutoRefinement::star_fn				=   set_star_fn;
    Map3dAutoRefinement::write_path				=   pathGetDir(output_fn);
    Map3dAutoRefinement::write_fn				=   pathGetFilename(output_fn);
    Map3dAutoRefinement::ref_fn					=   option.getStrOption("-ref");
    Map3dAutoRefinement::nr_classes				=   1;//option.getIntOption("-K");
    Map3dAutoRefinement::nr_pool				=   option.getIntOption("-pool");
    Map3dAutoRefinement::pixel_size				=   option.getFloatOption("-angpix");
    Map3dAutoRefinement::nr_iter				=   50;//option.getIntOption("-iter");
    Map3dAutoRefinement::random_seed			=   option.getIntOption("-random_seed");
    Map3dAutoRefinement::offset_step			=   option.getIntOption("-offset_step");
    Map3dAutoRefinement::offset_range			=   option.getIntOption("-offset_range");
    
    Map3dAutoRefinement::do_ctf_correction		=   option.getBoolOption("-ctf");
    Map3dAutoRefinement::only_flip_phases		=	option.getBoolOption("-only_flip_phases");
    Map3dAutoRefinement::ctf_phase_flipped		=	option.getBoolOption("-ctf_phase_flipped");
    Map3dAutoRefinement::refs_are_ctf_corrected	=	option.getBoolOption("-ctf_corrected_ref");
    Map3dAutoRefinement::intact_ctf_first_peak	=	option.getBoolOption("-intact_ctf_first_peak");
    
    Map3dAutoRefinement::ini_high				=   option.getFloatOption("-ini_high");
    Map3dAutoRefinement::adaptive_oversampling	=   option.getIntOption("-oversampling");
    Map3dAutoRefinement::sampler3d_healpix_order=   option.getIntOption("-healpix_order");
    Map3dAutoRefinement::tau2_fudge_factor		=   option.getFloatOption("-tau2_fudge");
    Map3dAutoRefinement::sampler3d_fn_sym		=   option.getStrOption("-sym");
    Map3dAutoRefinement::particle_diameter		=   option.getFloatOption("-particle_diameter");
    Map3dAutoRefinement::do_firstiter_cc		=   option.getBoolOption("-firstiter_cc");
    Map3dAutoRefinement::do_always_cc			=   option.getBoolOption("-always_cc");
    Map3dAutoRefinement::do_scale_correction	=   option.getBoolOption("-scale");
    Map3dAutoRefinement::do_norm_correction		=   option.getBoolOption("-norm");
    Map3dAutoRefinement::do_zero_mask			=   option.getBoolOption("-zero_mask");
    Map3dAutoRefinement::do_solvent            	=   option.getBoolOption("-flatten_solvent");
    Map3dAutoRefinement::do_map                	=   !option.getBoolOption("-no_map");
    Map3dAutoRefinement::mask_fn				=	option.getStrOption("-solvent_mask");
    Map3dAutoRefinement::do_shifts_onthefly		=	true;
    
    Map3dAutoRefinement::do_split_random_halves	=	option.getBoolOption("-split_random_halves");
    Map3dAutoRefinement::low_resol_join_halves	=	option.getFloatOption("-low_resol_join_halves");
    Map3dAutoRefinement::autosampling_hporder_local_searches = option.getFloatOption("-auto_local_healpix_order");
    
    Map3dAutoRefinement::data_stream_in_fn		=	option.getStrOption("-datastream_in");
    Map3dAutoRefinement::data_stream_out_fn		=	option.getStrOption("-datastream_out");
    Map3dAutoRefinement::continue_fn			=	option.getStrOption("-continue");
    Map3dAutoRefinement::guidance_fn       		=   option.getStrOption("-testguidance");
    
    Map3dAutoRefinement::setupMLoptimizer();
    
    Map3dAutoRefinement::prepare();
    
    Map3dAutoRefinement::iterate();
    
    Map3dAutoRefinement::destroyMLoptimizer();

    
    double t2 = dtime();
    
    if(Map3dAutoRefinement::node  == 0)
        std::cout<<"ML2D costs : "<<(t2-t1)
#ifdef DOUBLE_TRANSLATION
			<< " using DOUBLE_TRANSLATION"
#endif
#ifdef TRIPLE_TRANSLATION
			<< " using TRIPLE_TRANSLATION"
#endif
			<< std::endl;
    
	sDelete(benchmark);

#ifdef USEMPI
    MPI::Finalize();
#endif
    
    return 0;
}
