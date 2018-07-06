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
#pragma once

#if defined(ROME_MAP3D)
#error "map2d used in non ROME2D project"
#endif

#include "./initialize.h"	// ListOfImages
#include "./util.h"
#include "./memory.h"

// There are several partially orthogonal choices
//
// There are three implementations of the optimizer
//		Map2dOptimizer_old
//		Map2dOptimizer_original
//		Map2dOptimizer_new
//	They are all run via
//		map2d_optimizer.{h,cpp}
//	In theory, any subset can be run at the same time, but only the following tested
//		old 
//		original
//		new
//		original + new		// tandem comparison support has been built into these
// The following selects which of these three to run
//
namespace MapOptimizer_base_new {
	static const bool version_do_old      = false;
    static const bool version_do_new      = true;
	static const bool comparing_versions  = (int(version_do_old) + int(version_do_new) > 1);
}


namespace MapOptimizer_base_new {
	extern bool all_versions_agree_emit_test_output_possible ;
	extern bool all_versions_agree_emit_test_output_continues;

	inline bool emit_test_output() {
		return 
			true 
			&& all_versions_agree_emit_test_output_possible;
	}

	inline bool emit_test_output_prolific() {
		return 
			true											// change this as you desire
			&& all_versions_agree_emit_test_output_possible 
			&& all_versions_agree_emit_test_output_continues;
	}
}




//#define MY_MACHINE_HAS_AVX512
	//
	// This affects ml_optimizer_kernel and maybe other places where there is a better algorithm
	// if this hardware is available

#include "./sampling.h"
#include "./mrcs.h"
#include "./image.h"
#include "./fft_fftw3.h"
#include "./ctf.h"
#include "./time.h"
#include "./mpi.h"
#include "./metadata.h"
#include "./initialize.h"
#include "./string.h"
#include "./math.h"
#include "./ml_model.h"
#include "./map_model.h"
#include "./exp_model.h"

#include "./tandemExecution.h"

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif


namespace MapOptimizer_base_new {

	// Other Utilities
	//

	//#define LOADBALANCEANALYZER

//==

//--

//#define CHECK_FOR_NOT_WRITING_WHAT_SHOULD_BE_CONST
#if defined(MAP2DOPTIMIZER_BASE_BEING_COMPILED) || !defined(CHECK_FOR_NOT_WRITING_WHAT_SHOULD_BE_CONST)
#define ROME_ML2D_CONST 
#else
#define ROME_ML2D_CONST const
#endif

#if defined(MAP2DOPTIMIZER_ORIGINAL_BEING_COMPILED) || !defined(CHECK_FOR_NOT_WRITING_WHAT_SHOULD_BE_CONST)
#define MAP2DOPTIMIZER_CONST 
#else
#define MAP2DOPTIMIZER_CONST const
#endif

extern const int maxthreads;

// These variables are in the base, and should only be read by the specific optimizers
// TBD - enforce this
//
#define MAP2DOPTIMIZER_BASE_VARS \
	ELT(int			  , random_seed				, )				SEP \
	ELT(int			  , node					, )				SEP \
	ELT(int			  , nodes					, )				SEP \
    ELT(int			  , nr_global_images		, )				SEP \
    ELT(std::string   , mrcs_dir  				, )				SEP \
    ELT(std::string   , star_fn					, )				SEP \
    ELT(std::string   , write_path				, )				SEP \
    ELT(std::string   , write_fn				, )				SEP \
    ELT(std::string   , ref_fn					, )				SEP \
    ELT(std::string   , guidance_fn				, )             SEP \
    ELT(std::string   , data_stream_in_fn		, )             SEP \
    ELT(std::string   , data_stream_out_fn		, )             SEP \
    ELT(DataStream    , global_data_stream		, )				SEP \
	ELT(int			  , ori_size				, )				SEP \
    ELT(int			  , nr_pool					, )				SEP \
    ELT(int			  , nr_local_images			, )				SEP \
    ELT(int			  , first_local_image		, )				SEP \
	ELT(int			  , last_local_image		, )				SEP \
    ELT(double		  , sigma2_offset			, )				SEP \
    ELT(int			  , nr_classes				, )				SEP \
	ELT(int			  , adaptive_oversampling	, )				SEP \
    ELT(double        , offset_step				, )				SEP \
    ELT(double		  , offset_range			, )				SEP \
    ELT(double		  , psi_step				, )				SEP \
	ELT(HealpixSampler, sampling2d				, )				SEP \
    ELT(double		  , pixel_size				, )				SEP \
    ELT(double		  , particle_diameter		, )	/*in Ang*/	SEP \
    /*WARNING : some defalut change,if you wonder why get different result*/								\
    ELT(double		  , tau2_fudge_factor		, )	/*1.5 => 2*/SEP \
    ELT(bool		  , do_ctf_correction		, )	/*T => F */	SEP \
    ELT(bool		  , only_flip_phases		, )				SEP \
    ELT(bool		  , ctf_phase_flipped		, )				SEP \
    ELT(bool		  , do_norm_correction		, )	/*T => F */	SEP \
    ELT(bool		  , do_solvent				, )	/*T => F */	SEP \
    ELT(bool		  , do_zero_mask			, )	/*T => F */	SEP \
    ELT(int			  , nr_iter					, = 20)				// end of macro

#define SEP 
#define ELT(TYPE, NAME, INIT) extern TYPE ROME_ML2D_CONST NAME;
	MAP2DOPTIMIZER_BASE_VARS
#undef ELT
#undef SEP

    inline int   getPixelFromResolution(double resol) { return (int)(resol * pixel_size * ori_size); }
    inline double getResolution(int ipix)	          { return ipix/(pixel_size * ori_size);		 }
    inline double getResolutionAngstrom(int ipix)	  { return (ipix==0) ? 999. : (pixel_size * ori_size)/(double)ipix; }

    static const int  width_mask_edge		= 5;
    static const int  width_fmask_edge		= 2;
    static const double adaptive_fraction   = 0.999f;
    
    static const bool do_map                = true;
    
    // 2D or 3D interpolation parameter
    static const int gridding_nr_iter       = 10;
    
	// These are done from within the algorithm, even though they affect the above data
	// so when doing tandem execution, only one of them should execute
	//
	void update_metadata(MetaDataTable& metadata, int my_first_image, MetaDataElem const * exp_metadata, int exp_nr_images);
	void set_sigma2_offset(double to);
	void resetRandomlyPerturbedSampling(int iter);

	// top level flow
	// currently in rom_ml2d.cpp
	//
	void setupMap2dOptimizer();
	void readImages();
	void prepare();
	void iterate();
	void writeClassesAndMetadata(std::string fn_class,std::string fn_metadata, MetaDataTable const & metadata, double* model_Irefs,        std::ostream& testos);
	void writeClassesAndMetadata(std::string fn_class,std::string fn_metadata, MetaDataTable const & metadata, ListOfImages & model_Irefs, std::ostream& testos);
		// strangely this is implemented in rome2d.cpp
	void destroyMap2dOptimizer();

	void readOneImage(const size_t bufLen, float* buffer, MetaDataElem const & metaDataElem);

#define MAP2DOPTIMIZER_STATIC_SCALARS_EXPLICITLY_DEFINED \
	ELT(int,	iter)										SEP \
	ELT(double, ini_high)									SEP \
	ELT(double, current_resolution)							SEP \
	ELT(int,    current_size)								SEP \
    ELT(int,    coarse_size)								 // end of macro

#define MAP2DOPTIMIZER_STATIC_SCALARS_EXPECTATION_DEFINED \
    ELT(int,    exp_first_image)							SEP \
    ELT(int,    exp_last_image)                             SEP \
    ELT(int,    exp_nr_images)                              SEP \
    ELT(int,    exp_iclass_min)                             SEP \
    ELT(int,    exp_iclass_max)                             SEP \
    ELT(int,    exp_ipass)                                  SEP \
    ELT(int,    exp_current_size)                           SEP \
    ELT(int,    exp_current_oversampling)                   SEP \
    ELT(int,    exp_nr_over_psi)							SEP \
	ELT(int,    exp_nr_over_trans)							SEP \
	ELT(int,    exp_nr_trans)								SEP \
	ELT(int,	exp_nr_trans_adaptive_oversampling)			SEP \
	ELT(int,	exp_nr_psi_adaptive_oversampling)	            // end of macro

#define MAP2DOPTIMIZER_STATIC_SCALARS_IMPLICITLY_DEFINED \
    ELT(bool    , refs_are_ctf_corrected					/*, false	 */ )   SEP \
    ELT(bool    , has_converged								/*, false	 */ )   SEP \
    ELT(double  , smallest_changes_optimal_offsets          /*, 999.     */ )   SEP \
    ELT(double  , smallest_changes_optimal_orientations     /*, 999.     */ )   SEP \
    ELT(int     , smallest_changes_optimal_classes          /*, 9999999  */ )   SEP \
    ELT(double  , current_changes_optimal_orientations      /*, 999.     */ )   SEP \
    ELT(double  , current_changes_optimal_offsets           /*, 999.     */ )   SEP \
    ELT(double  , current_changes_optimal_classes           /*, 9999999. */ )   SEP \
    ELT(double  , sum_changes_optimal_orientations          /*, 0.       */ )   SEP \
    ELT(double  , sum_changes_optimal_offsets               /*, 0.       */ )   SEP \
    ELT(double  , sum_changes_optimal_classes               /*, 0.       */ )   SEP \
    ELT(double  , sum_changes_count                         /*, 0.       */ )   SEP \
    ELT(int     , nr_iter_wo_large_hidden_variable_changes  /*, 0        */ )   // end of macro

#define MAP2DOPTIMIZER_STATIC_SCALARS \
	MAP2DOPTIMIZER_STATIC_SCALARS_EXPLICITLY_DEFINED SEP		\
	MAP2DOPTIMIZER_STATIC_SCALARS_EXPECTATION_DEFINED SEP		\
	MAP2DOPTIMIZER_STATIC_SCALARS_IMPLICITLY_DEFINED			// end of macro

#define MAP2DOPTIMIZER_WSUM_SCALARS \
	ELT(double, LL)											SEP \
	ELT(double, ave_Pmax)									SEP \
	ELT(double, avg_norm_correction)						SEP \
	ELT(double, sumw_group)									SEP \
	ELT(double, sigma2_offset)								SEP \
    ELT(int,    data_size)										// end of macro
    
#define MAP2DOPTIMIZER_MODEL_SCALARS \
    ELT(double, LL)											SEP \
    ELT(double, avg_norm_correction)                        SEP \
    ELT(double, ave_Pmax)                                   SEP \
    ELT(int, Frefs_pad_size)                                    // end of macro

#define MAP2DOPTIMIZER_EXPITER_SCALARS \
	ELT(int,				Frefs_Rot_len_per_class)		SEP \
	ELT(Exp_Mweight,		Mweight)						SEP \
    ELT(Exp_Mcoarse_significant, Mcoarse_significant)	    // end of macro

#define MAP2DOPTIMIZER_EXPITER_ARRAYS \
	ELTONE( PerPsiDouble,   over_rot_psi)					SEP \
    ELTONE( PerTransDouble, over_trans_x)					SEP \
    ELTONE( PerTransDouble, over_trans_y)					SEP \
    ELTVEC( double,			highres_Xi2_imgs)				SEP \
    ELTVEC( double,			min_diff2)						SEP \
    ELTVEC( double,			old_offsetx)					SEP \
    ELTVEC( double,			old_offsety)					SEP \
    ELTVEC( double,			sum_weight)						SEP \
    ELTVEC( double,			significant_weight)				SEP \
    ELTVEC( double,			wsum_norm_correction)			SEP \
    ELTVoV( ImageDouble,	imgs)							SEP \
    ELTVoV( ModelData,		power_imgs)						SEP \
    ELTVoV( FimgsData,		Fimgs_real)						SEP \
    ELTVoV( FimgsData,		Fimgs_imag)						SEP \
    ELTVoV( FimgsData,		Fimgs_nomask_real)				SEP \
    ELTVoV( FimgsData,		Fimgs_nomask_imag)					// end of macro

#define MAP2DOPTIMIZER_EXPITER_ARRAYS_BOOL \
    ELTVoV( NrPsiChar,	    Rot_significant)					// end of macro


#define MAP2DOPTIMIZER_STATIC_ARRAYS_MyMetaData \
    ELTVEC( MetaDataElem,	exp_metadata)						// end of macro

#define MAP2DOPTIMIZER_STATIC_ARRAYS_INT \
    ELTVEC( int,		Npix_per_shell)						SEP \
    ELTVEC( int,		Mresol_coarse)						SEP \
    ELTVEC( int,		Mresol_fine)							// end of macro

#define MAP2DOPTIMIZER_STATIC_ARRAYS_DOUBLE \
    ELTVoV( FimgsData,	exp_Fctfs_writable)					SEP \
    ELTVoV( FimgsData,	exp_local_Fctfs_writable)			SEP \
    ELTVoV( FimgsData,	exp_local_Minvsigma2s)					// end of macro

#define MAP2DOPTIMIZER_STATIC_ARRAYS	\
    MAP2DOPTIMIZER_STATIC_ARRAYS_MyMetaData					SEP \
	MAP2DOPTIMIZER_STATIC_ARRAYS_INT						SEP \
	MAP2DOPTIMIZER_STATIC_ARRAYS_DOUBLE						// end of macro

#define MAP2DOPTIMIZER_MODEL_ARRAYS \
    ELTVEC( ModelData,      data_vs_prior_class)			SEP \
																\
    ELTVEC( ImageDouble,	Irefs)							SEP \
    ELTVoV( FrefPadData,    Frefs_pad_real)					SEP \
    ELTVoV( FrefPadData,    Frefs_pad_imag)					SEP \
    ELTVoV( ModelData,		tau2_class)						SEP \
    ELTONE( ModelData,		sigma2_noise)					SEP \
    ELTVoV( ModelData,		sigma2_class)					SEP \
    /* NOT USED ELTVEC( ModelData,		fsc_halves_class) */ SEP \
    ELTONE( PerClassDouble, pdf_class)						SEP \
    ELTONE( PerClassDouble, pdf_direction)					SEP \
    ELTONE( PerClassDouble, prior_offsetx_class)			SEP \
    ELTONE( PerClassDouble, prior_offsety_class)			// end of macro
															
#define MAP2DOPTIMIZER_WSUM_ARRAYS \
    ELTONE( ModelData,      sigma2_noise)					SEP \
    ELTONE( PerClassDouble, pdf_class)						SEP \
    ELTONE( PerClassDouble, prior_offsetx_class)			SEP \
    ELTONE( PerClassDouble, prior_offsety_class)			SEP \
    ELTONE( PerClassDouble, pdf_direction)					SEP \
	ELTVEC( WSUMModel::Data, data_real)					    SEP \
	ELTVEC( WSUMModel::Data, data_imag)					    SEP \
    ELTVEC( WSUMModel::Data, weight)						    // end of macro

	class LoadBalanceAnalyzer {
	public:
		LoadBalanceAnalyzer(const char* func, const char* file, int line, size_t numberOfTasks) 
		  : 
		  func(func), file(file), line(line), numberOfTasks(numberOfTasks),
			timeStarted(dtime()), 
			ompMaxThreads(omp_get_max_threads()), 
				// the real value, since measuring balancing across the threads
			perThread(ompMaxThreads)
		{
		}
		~LoadBalanceAnalyzer() {
			std::cout << "LoadBalanceAnalyzer " << func << " " << file << ":" << line << " numberOfTasks:" << numberOfTasks << std::endl;
			std::cerr << "LoadBalanceAnalyzer " << func << " " << file << ":" << line << " numberOfTasks:" << numberOfTasks << std::endl;
			double minStart = +std::numeric_limits<double>::infinity();
			double maxStart = -std::numeric_limits<double>::infinity();
			double minEnd   = +std::numeric_limits<double>::infinity();
			double maxEnd   = -std::numeric_limits<double>::infinity();
			for (int tid = 0; tid < ompMaxThreads; tid++) {
				auto & pt = perThread[tid];
				auto s = pt.started;
				auto e = pt.ended;
				minStart = std::min(minStart,s);
				maxStart = std::max(maxStart,s);
				minEnd   = std::min(minEnd,  e);
				maxEnd   = std::max(maxEnd,  e);
			}
			auto pct = [&](double t) { return int(100*(t-minStart)/(maxEnd-minStart)); };
			for (int tid = 0; tid < ompMaxThreads; tid++) {
				auto & pt = perThread[tid];
				if (pt.started > pt.ended) continue;
				auto f = [&](std::ostream&os) {
					os << " " << tid << " " << pct(pt.started) << "% .. " << pct(pt.ended);
					if (pt.iterations) os << " iterations:" << pt.iterations << " itMin:" << pt.iterationMinDuration << " itMax:" << pt.iterationMaxDuration;
					os << std::endl;
				};
				f(std::cout); f(std::cerr);
			}
		}
		void note() {
			int tid = omp_get_thread_num();
			double now = dtime();
			note(tid,now);
		}
		void iterationBegin() {
			int tid = omp_get_thread_num();
			double now = dtime();
			note(tid,now);
			auto & pt = perThread[tid];
			pt.iterationStarted = now;
			pt.iterations++;
		}
		void iterationEnd() {
			int tid = omp_get_thread_num();
			double now = dtime();
			note(tid,now);
			auto & pt = perThread[tid];
			auto duration = now - pt.iterationStarted;
			pt.iterationMaxDuration = std::max(pt.iterationMaxDuration, duration);
			pt.iterationMinDuration = std::min(pt.iterationMinDuration, duration);
		}
	private:
		void note(int tid, double now) {
			auto & pt = perThread[tid];
			pt.started = std::min(now,pt.started);
			pt.ended   = std::max(now,pt.ended);
		}
		const char* const func; const char* const file; int const line; size_t const numberOfTasks;
		double const timeStarted;
		const int ompMaxThreads;
		struct PerThread {
			double started;
			double ended;
			double iterationMinDuration;
			double iterationMaxDuration;
			double iterationStarted;
			size_t iterations;
			PerThread() { 
				started = iterationMinDuration = +std::numeric_limits<double>::infinity();
				ended   = iterationMaxDuration = -std::numeric_limits<double>::infinity();
				iterations = 0;
			}
		};
		std::vector<PerThread> perThread;
	};

	class Exp_Mbase {
	public:
		Exp_Mbase() : _size(0), _capacityPerClass(0), _sizePerClass(-1), _nr_images(0), _nr_classes(0) {}
		Exp_Mbase(Exp_Mbase const & rhs) : _nr_images(rhs._nr_images), _nr_classes(rhs._nr_classes), _capacityPerClass(rhs._capacityPerClass), _size(rhs._size) {}

		int nr_images	    () const { return _nr_images;    }
		int nr_classes      () const { return _nr_classes;   }
		int capacityPerClass() const { return _capacityPerClass; }
		int size            () const { return _size;         }

		int sizePerClass() const {
			assert(0 < _sizePerClass);
			return _sizePerClass;
		}

		void setSizePerClass(int to) {
			assert(0 < to && to <= capacityPerClass());
			_sizePerClass = to;
		}

	protected:
		void initBase(int nr_images, int nr_classes, int capacityPerClass) {
			_nr_images = nr_images;
			_nr_classes = nr_classes;
			_capacityPerClass = capacityPerClass;
			_size = nr_images*nr_classes*capacityPerClass;
		}

		void finiBase() {
			_size = 0;
			_capacityPerClass = 0;
			_nr_classes = 0;
			_nr_images = 0;
		}

		int  _nr_images;
		int  _nr_classes;
		int  _capacityPerClass;
		int  _sizePerClass;
		int  _size;

		int index(int iimage, int iclass, int withinClass) const {
			assert(0 <= iimage && iimage < _nr_images);
			assert(0 <= iclass && iclass < _nr_classes);
			assert(0 <= withinClass   && withinClass < _capacityPerClass);
			int index = iimage;
			index     = index*_nr_classes   + iclass;
			index     = index*_capacityPerClass + withinClass;
			assert(index < _size);
			return index;
		}
	};

	static const double Exp_Mweight_unsetWeight = -666.666;
	class Exp_Mweight : public Exp_Mbase {
		// metadata of the weight,structure is iclass->rot->iover_rot->trans->iover_trans

	public:
        // Double precision is needed because use involves summing and storing the square of the difference of image elements
		//
        typedef double Elt;

	private:
		Elt*  _ptr;
		char* _validImageClass;
		std::vector<bool> _minimasComputed;
		void resetMinimasComputed() {
			for (int i = 0; i < _nr_images; i++) _minimasComputed[i] = false;
		}
		int validImageClassIndex(int iimage, int iclass) const {
			assert(0 <= iimage && iimage < _nr_images);
			assert(0 <= iclass && iclass < _nr_classes);
			return iimage * _nr_classes + iclass;
		}
	public:
		Exp_Mweight() : Exp_Mbase(), _ptr(NULL), _minimasComputed(_nr_images) {
			resetMinimasComputed();
		}

		Exp_Mweight const & operator=(Exp_Mweight const & rhs) {
			if (_nr_images != rhs._nr_images || _nr_classes != rhs._nr_classes || _capacityPerClass != rhs._capacityPerClass) {
				if (_ptr) fini();
				init(rhs._nr_images, rhs._nr_classes, rhs._capacityPerClass);
			}
            for (int i = 0; i < _size;i++) _ptr[i] = rhs._ptr[i];
			return *this;
		}

		void computeMinimas(int iimage, bool to = true) {
			_minimasComputed[iimage] = to;
		}

		bool operator==(Exp_Mweight const & rhs) const {
			if (_size != rhs._size) return false;
			for (int i = 0; i < _size; i++) {
				auto lv = _ptr[i];
				auto rv = rhs._ptr[i];
				if (std::abs(lv-Exp_Mweight_unsetWeight) < 0.001) continue;		// these are used as "uninit" values in the new code
				if (std::abs(rv-Exp_Mweight_unsetWeight) < 0.001) continue;
                if (!nearEnoughTemplate(lv,rv)) {
                    std::cerr.precision(10);
                    std::cerr<<"diff in Exp_Mweight : "<<lv<<" "<<rv<<std::endl;
                    return false;
                }
			}
			return true;
		}

		void print(std::ostream & os) const {
			os << "Exp_Mweight { size: " << _size;
			for (int i = 0; i < std::min(_size, 5); i++) os << " " << i << ":" << _ptr[i];
			if (_size > 5) os << "...";
			os << "}";
		}

		void init(int nr_images, int nr_classes, int capacityPerClass) {
			assert(!_ptr);
			Exp_Mbase::initBase(nr_images, nr_classes, capacityPerClass);
			_ptr = mallocCacheAligned(Elt,_size);
			_minimasComputed.resize(nr_images);
			// Clearing them all is unnecessary and expensive.
			// Instead keep track of which are valid, and make sure they are written before being used.
            //		for (int i = 0; i < _size; i++) _ptr[i] = 0.;
			_validImageClass = vNew(char,nr_images*nr_classes);
			invalidateAll();
		}

		void invalidateAll() {
			assert(!!_ptr);
			resetMinimasComputed();
			for (int i = 0; i < _nr_images*_nr_classes; i++) _validImageClass[i] = 0;
		}

		void changeCapacityPerClass(int capacityPerClass) {
			assert(_ptr);
			aFree(_ptr);
			init(_nr_images, _nr_classes, capacityPerClass);
		}

		void fini() {
			vDelete(_validImageClass);
			aFree(_ptr); 
			Exp_Mbase::finiBase();
		}

		Elt* ptrForWrite(int iimage, int iclass, int withinClass)	{
			assert(!!_ptr);
			assert(!_minimasComputed[iimage]);
			_validImageClass[validImageClassIndex(iimage,iclass)] = 1;
			return _ptr+index(iimage, iclass, withinClass);  
		}
		Elt* ptrForModify(int iimage, int iclass, int withinClass)	{
			assert(!!_ptr);
			assert(!_minimasComputed[iimage]);
			assert(_validImageClass[validImageClassIndex(iimage,iclass)]);
			return _ptr+index(iimage, iclass, withinClass);  
		}
		Elt const * ptr(int iimage, int iclass, int withinClass) const {
			assert(!!_ptr);
			return _ptr+index(iimage, iclass, withinClass);  
		}
		bool isSet(int iimage, int iclass) const { return _validImageClass[validImageClassIndex(iimage,iclass)]; }
		Elt get(int iimage, int iclass, int withinClass, bool allowNegative=false) const {
			assert(!!_ptr); 
			//  if (!allowNegative && _ptr[index] < 0) {
			//  	std::cerr << "Exp_Mweight get uninited value " << std::endl;
			//  	exit(0);
			//  }
			assert(_validImageClass[validImageClassIndex(iimage,iclass)]);
			return _ptr[index(iimage, iclass, withinClass)]; 
		}
		void set(int iimage, int iclass, int withinClass, Elt to) {
			if (to < 0 || isnan(to)) {
				std::cerr << "Exp_Mweight set to " << to << std::endl;
				std::cout << "Exp_Mweight set to " << to << std::endl;
				EXIT_ABNORMALLY;
			}
			assert(!!_ptr); 
			assert(!_minimasComputed[iimage]);
			assert(_validImageClass[validImageClassIndex(iimage,iclass)]);
			_ptr[index(iimage, iclass, withinClass)] = to;   
		}
		void unsetAll(int iimage, int iclass, Elt to) {
			//  if (to >= 0 || isnan(to)) {
			//  	std::cerr << "Exp_Mweight unset to " << to << std::endl;
			//  	exit(0);
			//  }
			assert(!!_ptr);
			resetMinimasComputed();
			auto p = ptrForWrite(iimage, iclass, 0);
		    for (int i = 0; i < _capacityPerClass; i++) p[i] = to;
		}

		// Used for testing
		Exp_Mweight(Exp_Mweight const & rhs) : Exp_Mbase(rhs) {
			_ptr = mallocCacheAligned(Elt,_size);
			for (int i = 0; i < _size; i++) _ptr[i] = rhs._ptr[i];
		}
		void swapPtrs(Exp_Mweight & with) {
			if (_size != with._size) {
				std::cerr << "Wrong Exp_Mweight._size" << std::endl;
				std::cout << "Wrong Exp_Mweight._size" << std::endl;
				EXIT_ABNORMALLY;
			}
			std::swap(_ptr, with._ptr);
			std::swap(_validImageClass, with._validImageClass);
		}
		void check(const char* why, Exp_Mweight const & with) {
			assert(_size == with._size);
			// Get an average to help decide whether an error is important
			double sum(0);
			for (int i = 0; i < _size; i++) {
				sum += std::abs(_ptr[i]);
			}
			const Elt average = Elt(sum / _size);
			
			static size_t errors(0);
			static size_t warnings(0);
			for (int i = 0; i < _size; i++) {
				auto cv =      _ptr[i];
				auto ov = with._ptr[i];
				auto diff = std::abs(cv - ov);
				if (diff == 0) continue;
				auto max = std::max(std::abs(cv), std::abs(ov));
				if (diff*100.0 < max) continue;
				bool warningOnly = (diff*100.0 < average);
				if (warnings++ < 5) {
					std::cerr << "Exp_Mweight::check: Diff " << why << " " 
						<< " [i:" << i << "] is " << cv << " instead of " << ov 
						<< (warningOnly ? "  <<<<<<<<<<<<<<<<<<<<<<<<<<<<" 
										: "  ****************************" )
						<< " average copy[i] is " << average
						<< " so diff/avg=" << diff/average
						<< std::endl;
				}
				if (warningOnly) continue;
				if (++errors > 5) {
					std::cout << "Exp_Mweight::check more than 5 errors" << std::endl;
					EXIT_ABNORMALLY;
				}
			}
		}
	};

	inline std::ostream& operator<<(std::ostream& os, Exp_Mweight const & rhs) {
		rhs.print(os);
		return os;
	}

	class Exp_Mcoarse_significant : public Exp_Mbase {
	public:
		Exp_Mcoarse_significant() : Exp_Mbase(), _ptr(NULL) {}

		Exp_Mcoarse_significant const & operator=(Exp_Mcoarse_significant const & rhs) {
			if (_nr_images != rhs._nr_images || _nr_classes != rhs._nr_classes || _capacityPerClass != rhs._capacityPerClass) {
				if (_ptr) fini();
				init(rhs._nr_images, rhs._nr_classes, rhs._nr_psi, rhs._nr_trans);
			}
            for (int i = 0; i < _size; i++) _ptr[i] = rhs._ptr[i];
			return *this;
		}

		bool operator==(Exp_Mcoarse_significant const & rhs) const {
			if (_size != rhs._size) return false;
			size_t diffs = 0;
			for (int i = 0; i < _size; i++) {
				if (_ptr[i] != rhs._ptr[i]) diffs++; 
			}
			return diffs <= _size/100 || diffs <= 2;	// allow some error,  0:0 should return true
		}

		void print(std::ostream & os) const {
			os << "Exp_Mcoarse_significant { size: " << _size;
			for (int i = 0; i < std::min(_size, 5); i++) os << " " << i << ":" << int(_ptr[i]);
			if (_size > 5) os << "...";
			os << "}";
		}

		void init(int nr_images, int nr_classes, int nr_psi, int nr_trans) {
			Exp_Mbase::initBase(nr_images, nr_classes, nr_psi*nr_trans);
			_nr_psi = nr_psi; _nr_trans = nr_trans;
			setSizePerClass(_nr_psi*_nr_trans);
			_ptr = mallocCacheAligned(char,_size);
            for (int i = 0; i < _size; i++) _ptr[i] = 0;
		}
		void fini() {
			aFree(_ptr);
			Exp_Mbase::finiBase();
		}
		bool get(int iimage, int iclass, int irot, int itrans) {
			return _ptr[index(iimage, iclass, irot, itrans)];
		}

		void set(Exp_Mweight const & mweight, int iimage, double significant_weight) {
			assert(0 <= iimage);
			assert(iimage < _nr_images);

			if (_nr_classes != mweight.nr_classes()) {
				std::cerr << "Exp_Mcoarse_significant::set _nr_classes mismatch" << std::endl;
				std::cout << "Exp_Mcoarse_significant::set _nr_classes mismatch" << std::endl;
				EXIT_ABNORMALLY;
			}
			for (int iclass = 0; iclass < _nr_classes; iclass++) {
				char* const significantPtr = _ptr + index(iimage, iclass, 0,0);
				auto weightPtr = mweight.ptr(iimage,iclass,0);

				// This seems weird, but the current code initialized the mweight
				// to its maximum size which happens when adaptive_oversampling selects it
				// but then only uses the first portion of the per class data
				//
				if (_sizePerClass > _capacityPerClass) {
					std::cerr << "Exp_Mcoarse_significant::set sizePerClass _capacityPerClass mismatch" 
						<< "sizePerClass:" << _sizePerClass
						<< "E_capacityPerClass:" << _capacityPerClass
						<< std::endl;
				}
				if (_sizePerClass > mweight.capacityPerClass()) {
					std::cerr << "Exp_Mcoarse_significant::set _capacityPerClass mismatch" 
						<< "Exp_Mcoarse_significant _capacityPerClass:" << _capacityPerClass
						<< "Exp_Mweight _capacityPerClass:" << mweight.capacityPerClass()
						<< std::endl;
					std::cout << "Exp_Mcoarse_significant::set _capacityPerClass mismatch"
						<< "Exp_Mcoarse_significant _capacityPerClass:" << _capacityPerClass
						<< "Exp_Mweight _capacityPerClass:" << mweight.capacityPerClass()
						<< std::endl;
					EXIT_ABNORMALLY;
				}
				for (int i = 0; i < _sizePerClass; i++) {
                    // NOTE : use > instead of >= will cause error
					significantPtr[i] = (weightPtr[i] >= significant_weight);
				}
            }
		}

		bool isAnyTransSet(int iimage, int iclass, int irot) const {
			auto p = _ptr + index(iimage, iclass, irot, 0);
			for (int itrans = 0; itrans < _nr_trans; itrans++) {
				if (p[itrans]) return true;
			}
			return false;
		}

	private:
		int index(int iimage, int iclass, int irot, int itrans) const {
			assert(0 <= irot && irot < _nr_psi);
			assert(0 <= itrans && itrans < _nr_trans);
			return Exp_Mbase::index(iimage, iclass, irot*_nr_trans+itrans);
		}
		int _nr_psi;
		int _nr_trans;
	    char *_ptr;
	};

	inline std::ostream& operator<<(std::ostream& os, Exp_Mcoarse_significant const & rhs) {
		rhs.print(os);
		return os;
	}
}

class TestOutputMweight {
public:
	TestOutputMweight() : size(0) {}
	static bool shouldDo(int iclass, int iimage) {
		return (iclass ^ iimage) % 16 == 0;
	}
	void note(int rot_trans_over, double diff2) {
		if (rot_trans_over % 13 != 0) return;
		if (size >= capacity) return;
		rds[size].rot_trans_over = rot_trans_over;
		rds[size].diff2			 = diff2;
		size++;
	}
	void print(std::ostream & os);
private:
	static const size_t capacity = 256;
	size_t size;
	struct RD {int rot_trans_over; double diff2; bool operator<(RD const & rhs) const { return rot_trans_over < rhs.rot_trans_over; } };
	RD rds[capacity];
};

void tandemCompareFailed();

template <typename T>
bool tandemCompareSingleton(T & tentative, T const & known, int & errorCount, const char* name) {
	if (nearEnoughTemplate(tentative,known)) return true;
	errorCount++;
	if (errorCount < 4) {
        std::cerr.precision(20);
		std::cerr << "tandemCompareSingleton "
			<< name << " tentative:" << tentative << " != known:" << known 
			<< "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
			<< std::endl;
		nearEnoughTemplate(tentative,known,true);
		tandemCompareFailed();
	}
	return false;
}

static bool tandemCompareSingleton(int & tentative, bool const & knownBool, int & errorCount, const char* N) {
	// Needed to deal with VectorOfInt matching bool*
	int knownInt = int(knownBool);
	return tandemCompareSingleton(tentative, knownInt, errorCount, N);
}

namespace Map2dOptimizer {
    
#ifdef USEMPI
    extern MAP2DOPTIMIZER_CONST int nodeNameLen;
    extern MAP2DOPTIMIZER_CONST char nodeName[MPI_MAX_PROCESSOR_NAME];
#endif
    
}

#undef MAP2DOPTIMIZER_CONST
