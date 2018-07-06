#ifdef _WIN32

#pragma once

#include "./util.h"

class BmpWriter : public NoCopy {
public:
	typedef size_t V;
	const V _width;
	const V _height;
	BmpWriter(std::string const  bmpFilename, V width, V height);
	~BmpWriter();
	void draw(FloatRaster& floatRaster);	// writes every pixel

private:
	class Impl;
	Impl* const _impl;
};

void rasterToBmp(std::string const filename, FloatRaster& raster);

void rasterToBmpTest();

#endif