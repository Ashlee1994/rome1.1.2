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

#include "./map2d_optimizer.h"
#include "./map2d_optimizer_old.h"
#include "./map2d_optimizer_new.h"


namespace Map2dOptimizer {

	static const bool use_old			 = MapOptimizer_base_new::version_do_old;
	static const bool use_new			 = MapOptimizer_base_new::version_do_new;
	static const bool comparing_versions = MapOptimizer_base_new::comparing_versions;

	class TandemPhaseAlgorithms : public TandemExecution::Algorithms {
	public:
		TandemPhaseAlgorithms() {}
		virtual bool compare(unsigned int tag) {
			std::cout << "TandemPhaseAlgorithms::compare" << std::endl;
			if (!comparing_versions) return true;
			return Map2dOptimizer_new::waypointCompare();
		}
	};

	// major phases
	//
	void set_iter(int to) {
//		if (use_old)		Map2dOptimizer_old     ::set_iter(to);
		if (use_new)	    Map2dOptimizer_new     ::set_iter(to);
	}

	void setupMap2dOptimizer() {
		TUNING_SCOPE_STEP(setupMap2dOptimizer)
							MapOptimizer_base_new  ::setupMap2dOptimizer();
		if (use_old)		Map2dOptimizer_old     ::setupMap2dOptimizer();
		if (use_new)	    Map2dOptimizer_new     ::setupMap2dOptimizer();
		// will be destroyed in reverse order
	}

	void destroyMap2dOptimizer() {
		TUNING_SCOPE_STEP(destroyMap2dOptimizer)
		// reverse order to construction
							MapOptimizer_base_new  ::destroyMap2dOptimizer();
		if (use_old)		Map2dOptimizer_old	   ::destroyMap2dOptimizer();
		if (use_new)		Map2dOptimizer_new     ::destroyMap2dOptimizer();
	}

	void readImages() {
		TUNING_SCOPE_STEP(readImages)
							MapOptimizer_base_new  ::readImages();
		if (use_old)		Map2dOptimizer_old     ::readImages();
		if (use_new)		Map2dOptimizer_new     ::readImages();
	}

	void prepare() {
		TUNING_SCOPE_STEP(prepare)
								MapOptimizer_base_new  ::prepare();
		if (!comparing_versions) {
			if (use_old)		Map2dOptimizer_old     ::prepare();
			if (use_new)		Map2dOptimizer_new     ::prepare();
		} else {
			assert(!use_old);
			std::cout << "prepare started <<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
			TandemPhaseAlgorithms a;
			if (use_new)	   Map2dOptimizer_new     ::prepareAlgorithm(&a);
			a.run();
			std::cout << "prepare finished <<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
		}
	}

	void iterate() {
								MapOptimizer_base_new  ::iterate();
		if (!comparing_versions) {
			if (use_old)		Map2dOptimizer_old     ::iterate();
			if (use_new)		Map2dOptimizer_new     ::iterate();
		} else {
			assert(!use_old);
			std::cout << "iterate started <<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
			TandemPhaseAlgorithms a;
			if (use_new)	   Map2dOptimizer_new     ::iterateAlgorithm(&a);
			a.run();
			std::cout << "iterate finished <<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
		}
	}

	void writeClassesAndMetadata() {
		TUNING_SCOPE_STEP(writeClassesAndMetadata)
		// (use_old)       Map2dOptimizer_old     ::writeClassesAndMetadata();	// iterate does this
		if (use_new)	   Map2dOptimizer_new     ::writeClassesAndMetadata();
	}

}