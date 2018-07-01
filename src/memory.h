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

#ifndef MEMORY_H_
#define MEMORY_H_

#include <cstddef> /* NULL size_t*/
#include <cassert>
#include <cstring> /* memset */
// Census support
//
size_t census();
void showPopulation();

void censusWkr(void* p, const char* file, size_t line, bool died);

template<typename T>
T* bornTemplate(T* p, const char* file, size_t line) {
    censusWkr(p, file, line, false);
    return p;
}
#define born(p) bornTemplate(p,__FILE__,__LINE__)

template<typename T>
T* diedTemplate(T* p, const char* file, size_t line) {
    censusWkr(p, file, line, true);
    return p;
}
#define died(p) diedTemplate(p, __FILE__,__LINE__)


// Heap support
//
template<typename T>
T* mallocCacheAlignedTemplate(size_t len, const char* file, size_t line) {
    return bornTemplate((T*)_mm_malloc(sizeof(T)*len, 64), file, line);
}
#define mallocCacheAligned(T,L) mallocCacheAlignedTemplate<T>(L, __FILE__, __LINE__)


template<typename T>
void freeCacheAligned(T*&p) {
    _mm_free(died(p));
    p = NULL;
}

namespace Heap {
    template <class ScalarType>
    ScalarType* allocScalars(size_t s, const char* file, size_t line) {
        return mallocCacheAlignedTemplate<ScalarType>(s, file, line);
    }
    template <class ScalarType>
    ScalarType* allocInitedScalars(size_t s, ScalarType v, const char* file, size_t line) {
        auto p = allocScalars<ScalarType>(s, file, line);
        for (size_t i = 0; i < s; i++) p[i] = v;
        return p;
    }
    template <class ScalarType>
    ScalarType* allocZeroedScalars(size_t s, const char* file, size_t line) {
        return allocInitedScalars<ScalarType>(s, ScalarType(0), file, line);
    }
    template <class ScalarType>
    void freeScalars(ScalarType* & p) {
        freeCacheAligned(p); p = NULL;
    }
    static auto allocChars = allocScalars < char > ;
    static auto allocFloats = allocScalars < float > ;
    static auto allocDoubles = allocScalars < double > ;
    static auto allocZeroedDoubles = allocZeroedScalars < double > ;
    static auto freeDoubles = freeScalars < double > ;
    static auto freeChars = freeScalars < char > ;
    static auto freeFloats = freeScalars < float > ;
    static auto freeInts = freeScalars < int > ;
    double* allocDoublesForImages(const char* name, size_t nr_images, size_t s);
    double* allocDoublesForClasses(const char* name, size_t nr_classes, size_t s);
    double* allocZeroedDoublesForClasses(const char* name, size_t nr_classes, size_t s);
    void traceDoubles(double* ptr, const char* name, size_t sizeOfItem, size_t s);
    void traceFloats(float * ptr, const char* name, size_t sizeOfItem, size_t s);
    void untrace(void  * ptr);
};

#define mallocDoubles(len)			Heap::allocDoubles      (len,		__FILE__, __LINE__)
#define mallocFloats(len)			Heap::allocFloats       (len,		__FILE__, __LINE__)
#define mallocValDoubles(len,val)   Heap::allocInitedDoubles(len, val,	__FILE__, __LINE__)
#define mallocValFloats(len,val)    Heap::allocInitedFloats (len, val,	__FILE__, __LINE__)
#define mallocZeroedDoubles(len)	Heap::allocZeroedDoubles(len,		__FILE__, __LINE__)
#define mallocZeroedFloats(len)   	Heap::allocZeroedFloats (len,		__FILE__, __LINE__)


#ifdef USEMCDRAM
// high bandwidth memory
#include  <hbwmalloc.h>   // hbwmalloc interface
#endif
//
static void *
allocMemory(size_t __size, size_t __align,bool __inHBM = false)
{
    void *__mallocedMemory;
#ifdef USEMCDRAM
    if (__inHBM) int ret = hbw_posix_memalign((void**) &__mallocedMemory, 64, __size);
    else __mallocedMemory = (void*)_mm_malloc(__size, __align);
#else
    __mallocedMemory = (void*)_mm_malloc(__size, __align);
#endif
    return __mallocedMemory;
}
//
static void
freeMemory(void *__p,bool __inHBM = false)
{
#ifdef USEMCDRAM
    if (__inHBM) hbw_free(__p);
    else _mm_free(__p);
#else
    _mm_free(__p);
#endif
}

#endif /* defined(MEMORY_H_) */
