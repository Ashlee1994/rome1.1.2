/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Authors: "Bevin Brett"
 * Intel Corporation
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
#include "../src/map_optimizer_base_new.h"
#include <algorithm>

#ifdef USE_TBB
#include "tbb/tbbmalloc_proxy.h"
    // preliminary laptop measurements show this slower than the default allocator
#endif
#include <list>

int rs = MapOptimizer_base_new::random_seed;

void TestOutputMweight::print(std::ostream & os) {
	std::sort(&rds[0], &rds[size]);
	os << "TestOutputMweight{" << std::endl;
	for (size_t i = 0; i < size; i++) {
		os << rds[i].rot_trans_over << " => " << rds[i].diff2 << std::endl;
	}
	os << "}";
}

void __declspec(noinline) tandemCompareFailed() {
    return;
}


namespace MapOptimizer_base_new {

	bool all_versions_agree_emit_test_output_possible  = true;
	bool all_versions_agree_emit_test_output_continues = true;

#define SEP
#define ELT(TYPE, NAME, INIT) TYPE NAME INIT;
    MAP2DOPTIMIZER_BASE_VARS
#undef ELT
#undef SEP

const int maxthreads =
    omp_get_max_threads()
    + 2     // make sure we can cope when omp doesn't use all the threads for a parallel region
    ;


void readOneImage(const size_t bufLen, float* buffer, MetaDataElem const & metaDataElem) {
    TUNING_SCOPE_STEP(readOneImage)

     assert(bufLen == ori_size*ori_size);

#ifndef RANDDATA

    size_t image_id = metaDataElem.IMAGE.INDEX;  // For some reason, IMAGE_NAME is a double

    FILE* filehandle = fopen((mrcs_dir+metaDataElem.IMAGE.NAME).c_str(),"rb");
    if (!filehandle) {
        std::cerr<<"fopen failed opening "<<mrcs_dir+metaDataElem.IMAGE.NAME<<std::endl;
		EXIT_ABNORMALLY;
    }

    long offset  = (256+long(image_id-1)*ori_size*ori_size)*sizeof(float);
    if (fseek(filehandle,offset,SEEK_SET)) {
        std::cerr<<"read file failed."<<std::endl;
		EXIT_ABNORMALLY;
    }

    if (fread((char*)buffer, ori_size*ori_size*sizeof(float), 1, filehandle) == NULL) {
        std::cerr<<"read file failed."<<std::endl;
		EXIT_ABNORMALLY;
    }

    fclose(filehandle);

    if (false) normalizeData(buffer, ori_size*ori_size);

#else
    for (int i = 0; i < ori_size*ori_size; i++) {
        buffer[i] = (float)rand()/RAND_MAX;
    }
    normalizeData(buffer,ori_size*ori_size);
#endif

}

void setupMap2dOptimizer() {

	// Initialize all the MAP2DOPTIMIZER_BASE_VARS in the order they appear in that definition
	//
    // random_seed = set_random_seed;
    // Also randomize random-number-generator for perturbations on the angles
    
#ifdef USEMPI
    // MPI::Init();
    nodes = MPI::COMM_WORLD.Get_size();
    node = MPI::COMM_WORLD.Get_rank();
    // MPI::Get_processor_name(nodeName,nodeNameLen);
#else
    nodes = 1;
    node  = 0;
#endif
    
    NODE0ONLY std::cout<<"maxthreads = "<<maxthreads<<std::endl;
    
	MetaDataTable metadata;
    metadata.readFromStar(star_fn);
    nr_global_images = metadata.numberOfParticles();
    
    //int mrcsHead[256];
    Mrcs::MrcsHead mrcsHead;
    Mrcs::readMrcsHead(metadata[0].IMAGE.NAME, mrcsHead);
    NODE0ONLY std::cout<<"x = "<<mrcsHead.NC<<",y = "<<mrcsHead.NC<<",z = "<<mrcsHead.NS<<",mode = "<<mrcsHead.MODE<<std::endl;
    
    ori_size = mrcsHead.NC;

#ifndef SAME_IMAGE_DIVISION
    //divided like not mpi version,just to easily debug
    nr_local_images = divide_equally_pool(nr_global_images,nr_pool,nodes,node, first_local_image, last_local_image);
#else
    nr_local_images = divide_equally(nr_global_images,nodes,node, first_local_image, last_local_image);
#endif
    
    sigma2_offset = 3*3;
    
    // Initialise the sampling object (sets prior mode and fills translations and rotations inside sampling object)
    sampling2d.initialize(offset_step,offset_range,psi_step);
	
    dontShare_Random_generator.init(random_seed);
    
#ifdef DATA_STREAM
    global_data_stream.init(data_stream_out_fn, data_stream_in_fn);
    global_data_stream.foutInt(&nr_iter, 1, "nr_iter", __FILE__, __LINE__);
    global_data_stream.foutInt(&nr_classes, 1, "nr_classes", __FILE__, __LINE__);
    global_data_stream.foutDouble(&pixel_size, 1, "pixel_size", __FILE__, __LINE__);
    global_data_stream.foutInt(&random_seed, 1, "random_seed", __FILE__, __LINE__);
    double ini_high = 1;
    global_data_stream.foutDouble(&ini_high, 1, "ini_high", __FILE__, __LINE__);
    global_data_stream.foutDouble(&tau2_fudge_factor, 1, "tau2_fudge_factor", __FILE__, __LINE__);
    global_data_stream.foutDouble(&particle_diameter, 1, "particle_diameter", __FILE__, __LINE__);
    global_data_stream.foutInt(&adaptive_oversampling, 1, "adaptive_oversampling", __FILE__, __LINE__);
    global_data_stream.foutInt(&sampling2d.healpix_order, 1, "sampling.healpix_order", __FILE__, __LINE__);
    global_data_stream.foutDouble(&sampling2d.psi_step, 1, "sampling.psi_step", __FILE__, __LINE__);
    global_data_stream.foutDouble(-91, "sampling.limit_tilt", __FILE__, __LINE__);
    global_data_stream.foutDouble(&offset_range, 1, "sampling.offset_range", __FILE__, __LINE__);
    global_data_stream.foutDouble(&offset_step, 1, "sampling.offset_step", __FILE__, __LINE__);
    global_data_stream.foutDouble(0.5, "sampling.perturbation_factor", __FILE__, __LINE__);
    global_data_stream.check();global_data_stream.flush();
#endif
}

void destroyMap2dOptimizer() {
   //
}

void readImages() {
}

void prepare() {
}

void iterate() {
}

void __declspec(noinline) update_metadata(MetaDataTable& metadata, int my_first_image, MetaDataElem const * exp_metadata, int exp_nr_images) {
    // memcpy(metadata+my_first_image, exp_metadata, sizeof(MyMetaData)*exp_nr_images);
    for (int iimage = 0; iimage < exp_nr_images; iimage++) {
        metadata[my_first_image+iimage] = exp_metadata[iimage];
    }
}

void set_sigma2_offset(double to) {
	sigma2_offset = to;
}

void resetRandomlyPerturbedSampling(int iter) {
    dontShare_Random_generator.init(random_seed + iter);
    // Reset the random perturbation for this sampling
    sampling2d.resetRandomlyPerturbedSampling();
}


void writeClassesAndMetadata(std::string fn_class,std::string fn_metadata, MetaDataTable const & metadata, double* model_Irefs, std::ostream& testos) {
	ListOfImagesFromDoubleVector listOfImages(model_Irefs, nr_classes, ori_size, ori_size*ori_size);
	writeClassesAndMetadata(fn_class, fn_metadata, metadata, listOfImages, testos);
}

void writeClassesAndMetadata(std::string fn_class, std::string fn_metadata, MetaDataTable const & metadata, ListOfImages & model_Irefs, std::ostream& testos) {
    
    NODE0ONLY
    {
        // A. Write Classes
        // int mrcsHead[256];
        // mrcsHead[0] = ori_size;
		// mrcsHead[1] = ori_size;
		// mrcsHead[2] = nr_classes;
		// mrcsHead[3] = 2;
        // writeMrcsData(fn_class, mrcsHead, model_Irefs);
		class FloatImages : public ::FloatImages {
		public:
			FloatImages(ListOfImages & v) : v(v), fv(v.nr_images()) {}
			virtual int nr_images() { return v.nr_images(); }
			virtual int imageSide() { return v.imageSide(); }
			virtual float* image_ptr(size_t i) {
				auto& fi = fv[i];
				if (fi.size() == 0) {
					auto n = imageSide()*imageSide();
					auto d = v.rImage(i);
					fi.resize(n);
					auto f = fi.data();
					for (int i = 0; i < n; i++) f[i] = float(d[i]);
				}
				return fi.data();
			}
		private:
			ListOfImages & v;
			std::vector<std::vector<float>> fv;
		} floatImages(model_Irefs);
        Mrcs::writeMrcsData(fn_class, floatImages);
        
		if (emit_test_output()) {
	       for (int img_id = 0; img_id < nr_global_images; img_id++) {
	    		testos << "~~~TEST OUTPUT: metadata.writeToStar image " << img_id << std::endl;
	    		testos << "~~~ " << metadata[img_id].CTF_VOLTAGE 		<< std::endl;
	    		testos << "~~~ " << metadata[img_id].CTF_DEFOCUS_U		<< std::endl;
	    		testos << "~~~ " << metadata[img_id].CTF_DEFOCUS_V		<< std::endl;
	    		testos << "~~~ " << metadata[img_id].CTF_DEFOCUS_ANGLE	<< std::endl;
	    		testos << "~~~ " << metadata[img_id].CTF_CS				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].CTF_Q0				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].IMAGE              << std::endl;
	    		testos << "~~~ " << metadata[img_id].NORM				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].GROUP_NO			<< std::endl;
	    		testos << "~~~ " << metadata[img_id].XOFF				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].YOFF				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].ROT				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].TILT				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].PSI				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].CLASS				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].DLL				<< std::endl;
	    		testos << "~~~ " << metadata[img_id].NR_SIGN			<< std::endl;
	    		testos << "~~~ " << metadata[img_id].PMAX				<< std::endl;
	    	}
		}
        
		metadata.writeToStar(fn_metadata);

    }
}

}	// namespace MapOptimizer_base_new


