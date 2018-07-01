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

#ifndef UTIL_H_
#define UTIL_H_

#include <assert.h>

#include <cmath>
#include <algorithm>
#include <map>
#include <omp.h>
#include <vector>

#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>



#if 0
namespace std {
#define STD_TO_STRING // TOOD
#ifndef STD_TO_STRING
	std::string to_string(double d);
#endif
}
#endif

template <class T>
T square(T v) {
	return v*v;
}


static int divRoundUp(int lhs, int rhs) { return (lhs + rhs - 1) / rhs; }

// convenient vector operations
//
template<typename T>
static T sumvec(T *vec, int length)
{
	T sum = 0;
	for (int i = 0; i < length; i++) {
		sum += fabs(vec[i]);
	}
	return sum;
}

template<typename T>
static T sumvec(T *real_vec, T *imag_vec, int length)
{
    T sum = 0;
    for (int i = 0; i < length; i++) {
        sum += fabs(real_vec[i]);
        sum += fabs(imag_vec[i]);
    }
    return sum;
}

template<typename T>
static void printvec(std::ostream& out, T* vec, int length)
{
	for (int i = 0; i < std::min(100, length); i++)
		out << vec[i] << " ";
	out << std::endl;
};

template<typename T>
static bool checkvec(std::ostream& out, T* vec1, T* vec2, int length){
	size_t count = 0;
	for (int i = 0; i < length; i++){
		if ((vec1[i] - vec2[i]) > 1e-5)
		{
			/*MASTERNODE*/ if (count == 0) out << "#### large error(index,value1,value2) : (" << i << "," << vec1[i] << "," << vec2[i] << ") ";
			count++;
		}
	}
	/*MASTERNODE*/ if (count > 0) out << " (error_count,data_length,percentage)(" << count << "," << length << "," << (double)count / length*100. << "%) ";
	if (count > 0) return false;
	else return true;
};

template<typename T1,typename T2>
static bool checkDiffTypeVec(std::ostream& out, T1* vec1, T2* vec2, int length){
    size_t count = 0;
    for (int i = 0; i < length; i++){
        if ((vec1[i] - vec2[i]) > 1e-5)
        {
            /*MASTERNODE*/ if (count == 0) out << "#### large error(index,value1,value2) : (" << i << "," << vec1[i] << "," << vec2[i] << ") ";
            count++;
        }
    }
    /*MASTERNODE*/ if (count > 0) out << " (error_count,data_length,percentage)(" << count << "," << length << "," << (double)count / length*100. << "%) ";
    if (count > 0) return false;
    else return true;
};


class DoSomePerIter {
public:
	int count, limit, currIter;
	DoSomePerIter() : currIter(-1) {}
	template <typename Callable>
	void note(const bool doIt, int iter, Callable callable)
	{
		if (doIt && omp_get_thread_num() == 0)
#pragma omp critical
		{
			if (currIter != iter) { currIter = iter; count = 0; limit = 1; }
			if (count++ >= limit) {
				if (limit < 4) limit++; else limit *= 2;
				callable(count);
			}
		}
	}
};

#define EXIT_ABNORMALLY exitAbnormally(__FILE__,__LINE__)
void exitAbnormally(const char* file, int line);


struct DoublePair {
	double lo;
	double hi;
	DoublePair() : lo(0), hi(0) {}
	DoublePair(double lo, double hi) : lo(lo), hi(hi) {}
	DoublePair(std::string);
	std::string string() const;
	bool operator==(DoublePair const & rhs) const { return lo == rhs.lo && hi == rhs.hi; }
	bool operator!=(DoublePair const & rhs) const { return !(*this == rhs); }
	double mid() const { return (lo + hi) / 2.0; }
	static DoublePair midFrac(double mid, double frac) { auto diff = std::abs(frac*mid); return DoublePair(mid - diff, mid + diff); }
};


template <typename T> void destroy(T*&p) { delete p; p = NULL; }


// By using this as a member or as a base class, the code will get an error message
// if an attempt is made to copy the object.  Useful when copying doesn't make sense.
//
class NoCopy {
	NoCopy(NoCopy const & rhs) = delete;				// Private and unimplemented so can't be used
	NoCopy& operator=(NoCopy const & rhs) = delete;
public:
    NoCopy()=default;
};

// Much faster than a vector<char> or string
//
struct FastCharBuf {
	FastCharBuf(size_t initialCapacity)
		: capacity(std::max(size_t(64), initialCapacity)),
		ptr(new char[capacity]),
		size(0)
	{}
	~FastCharBuf() { delete[] ptr; }

	size_t capacity;
	char*  ptr;
	size_t size;

	void push_back(char c) {
		if (size == capacity) grow();
		ptr[size++] = c;
	}

	void grow() {
		capacity += capacity / 2;
		char* p = new char[capacity];
		memcpy(p, ptr, size);
		delete[] ptr;
		ptr = p;
	}
};


// Floating point comparisons for equality are too risky in many cases
// so this provides a slightly less code-generation-sensitive test
//
template <typename T>
bool nearEnoughTemplate(T const & tentative, T const & known, bool info = false) {
	if (info) { std::cerr << "nearEnoughTemplate<T> " << tentative << " != " << known << std::endl; }
	return tentative == known;
}
template <> bool nearEnoughTemplate(double const & tentative, double const & known, bool info);
template <> bool nearEnoughTemplate(float const & tentative, float const & known, bool info);

// Some convenience io support
//
class ifstreamCheckingExistence : public std::ifstream {
public:
	ifstreamCheckingExistence(const char* fileName);
};

class ofstreamCheckingCreated : public std::ofstream {
public:
	ofstreamCheckingCreated(const char* fileName);
};

#endif
