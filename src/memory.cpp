/***************************************************************************
 *
 * Authors: "Yongbei(Glow) Ma,Jiayi (Timmy) Wu, Youdong (Jack) Mao"
 * Dana-Farber Cancer Institute, Harvard Medical School and Peking University
 * Authors: "Brett, Bevin"
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

#include "memory.h"

static size_t populationCount;

//#define DEBUG_LEAKS
#if defined(DEBUG_LEAKS)
static size_t populationShowTrigger = 10000;
struct BirthPlace {
    const char* file; size_t line;
    BirthPlace(const char* file, size_t line) : file(file), line(line) {}
    bool operator<(BirthPlace const & rhs) const {
        auto fd = strcmp(file, rhs.file);
        return (fd != 0) ? (fd < 0) : (line < rhs.line);
    }
};

struct CensusCountable {
    const char* file; size_t line;  size_t tick;
    CensusCountable(const char* file, size_t line, size_t tick) : file(file), line(line), tick(tick) {}
};
std::map<void*, CensusCountable> population;
#endif

size_t census() {
#if defined(DEBUG_LEAKS)
    return population.size();
#else
    return populationCount;
#endif
}



void showPopulation() {
#if defined(DEBUG_LEAKS)
    std::map<BirthPlace,size_t> perBP;
    for (auto i = population.begin(); i != population.end(); i++) {
        BirthPlace bp(i->second.file, i->second.line);
        auto c = perBP.find(bp);
        if (c == perBP.end()) { perBP.insert(std::make_pair(bp,0)); c = perBP.find(bp); }
        c->second++;
    }
    std::cerr << "Census taken! populationCount:" << populationCount << std::endl;
    for (auto i = perBP.begin(); i != perBP.end(); i++) {
        std::cerr << "    " << i->first.file << ":" << i->first.line << " " << i->second << std::endl;
    }
#endif
}

void censusWkr(void* p, const char* file, size_t line, bool died) {
    if (died && !p) return;
    
    if (died)
#pragma omp atomic
        populationCount--;
    else
#pragma omp atomic
        populationCount++;
    
#if defined(DEBUG_LEAKS)
#pragma omp critical
    {
        // Note: Ptrs can be freed on different threads than the creator
        //
        static size_t watchTick = 0;
        static void*  watchP    = nullptr;
        static size_t tick;
        
        auto i = population.find(p);
        if (i == population.end()) {
            if (died) {
                std::cerr << "Undocumented alien died" << std::endl;
            } else {
                population.insert(i, std::make_pair(p,CensusCountable(file,line,++tick)));
                if (watchTick == tick) {
                    std::cerr << "Identity theft victim p:" << (void*)p << " being made at tick:" << tick << std::endl;
                    watchP = p;
                }
            }
        } else {
            if (died) {
                population.erase(i);
                if (watchP == p) {
                    std::cerr << "watchP:" << (void*)p << " died" << std::endl;
                    watchP = nullptr;
                }
            } else {
                std::cerr << "Identify theft p:" << (void*)p << " previously made at tick:" << i->second.tick << std::endl;
                exit(1);
            }
        }
        
        if (!died && populationCount > populationShowTrigger) {
            if (populationShowTrigger < 1000000) populationShowTrigger *= 10; else populationShowTrigger += 1000000;
            showPopulation();
        }
    }
#endif
}


namespace Heap {
    
    static double* malloc2DDoubles(const char* name, size_t sizeOfItem, size_t s) {
        assert(s % sizeOfItem == 0);
        auto p = mallocDoubles(s);
        return p;
    }
    static double* mallocZeroed2DDoubles(const char* name, size_t sizeOfItem, size_t s) {
        assert(s % sizeOfItem == 0);
        auto p = mallocZeroedDoubles(s);
        return p;
    }
    
    double* allocDoublesForImages(const char* name, size_t nr_images, size_t s) { return malloc2DDoubles(name, nr_images, s); }
    double* allocDoublesForClasses(const char* name, size_t nr_classes, size_t s) { return malloc2DDoubles(name, nr_classes, s); }
    double* allocZeroedDoublesForClasses(const char* name, size_t nr_classes, size_t s) { return mallocZeroed2DDoubles(name, nr_classes, s); }
};
