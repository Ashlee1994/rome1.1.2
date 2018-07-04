#include "util.h"

#include <algorithm>
#include <ostream>
#include <sstream>
//#include "primitives_for_each_os.h"

#if 0
namespace std {
#ifndef STD_TO_STRING
	std::string to_string(double d) {
		std::ostringstream os;
		os << d;
		return os.str();
	}
#endif
}
#endif

void exitAbnormally(const char* file, int line) {
	std::cerr << "exitAbnormally called" << std::endl; 	exit(1);
}

DoublePair::DoublePair(std::string s) {
	std::istringstream is(s);
	is >> lo;
	char colon; is >> colon;
	is >> hi;
}

std::string DoublePair::string() const {
	return std::to_string((long double)lo) + ":" + std::to_string((long double)hi);
}

ifstreamCheckingExistence::ifstreamCheckingExistence(const char* fileName) : std::ifstream(fileName) {
	if (!fail()) return;
	std::cerr << "Failed to open " << fileName << std::endl;
	//std::cerr << "Current directory: " << currentDirectory() << std::endl;
	EXIT_ABNORMALLY;
}

ofstreamCheckingCreated::ofstreamCheckingCreated(const char* fileName) : std::ofstream(fileName) {
	if (!fail()) return;
	std::cerr << "Failed to create " << fileName << std::endl;
	//std::cerr << "Current directory: " << currentDirectory() << std::endl;
	EXIT_ABNORMALLY;
}

template <> bool nearEnoughTemplate(double const & tentative, double const & known, bool info) {
	auto max = std::max(std::abs(tentative), std::abs(known));
	if (info) {
		std::cerr
			<< "nearEnoughTemplate<float> " << tentative << " ?= " << known
			<< ": Rel err (" << std::abs(tentative - known) / max << ") too large ( > 1.0e-1)"
			<< std::endl;
		std::cerr << std::endl;
	}
	return (max <= 1.0e-2) || std::abs(tentative - known) / max < 1.0e-1;
}

template <> bool nearEnoughTemplate(float const & tentative, float const & known, bool info) {
	double t = tentative;
	double k = known;
	return nearEnoughTemplate(t, k, info);
}
