#ifdef _WIN32

#include "./util.h"

#include "./bmpWriter.h"
#include "./mrcs.h"
#include "./jpegEtc.h"

void rasterToBmp(std::string const filename, FloatRaster& raster) {
	BmpWriter bmpWriter(filename, raster._width, raster._height);
	bmpWriter.draw(raster);
}


static Lock     lock;
static JPEGEtc* jpegEtc;


class Raster_DrawPhoto : public JPEGEtcGrey {
public:
	Raster_DrawPhoto(
		size_t width,
		size_t height,
		FloatRaster& floatRaster)
	  : JPEGEtcGrey(SRCLOC),
		width(width), height(height), floatRaster(floatRaster) {
		maxIntensity = -std::numeric_limits<float>::max();
		minIntensity = std::numeric_limits<float>::max();
		for (FloatRaster::V h = 0; h < height; h++) {
			for (FloatRaster::V w = 0; w < width; w++) {
				auto f = floatRaster(w, h);
				maxIntensity = max(maxIntensity, f);
				minIntensity = min(minIntensity, f);
			}
		}
		if (minIntensity == maxIntensity) maxIntensity += 1.0f;
	}
	virtual float intensity(int x, int y) const {
		auto w = toIntegral<FloatRaster::V>(x);
		auto h = toIntegral<FloatRaster::V>(y);
		auto f = floatRaster(w, h);
		return (f - minIntensity) / (maxIntensity - minIntensity);
	}
private:
	const FloatRaster::V width;
	const FloatRaster::V height;
	FloatRaster&         floatRaster;
	float                maxIntensity;
	float                minIntensity;
};



class BmpWriter::Impl : public Citizen {
public:
	Impl(SRCLOCARG) : Citizen(srcLoc), acquire(lock, srcLoc) {
		if (!jpegEtc) jpegEtc = JPEGEtc::alloc(SRCLOC);
	}
	~Impl() {}
	Lock::ForScope	acquire;
	Filename		bmpFilename;
	JPEGEtc::Photo*	photo;
};

BmpWriter::BmpWriter(SRCLOCARG, Filename const & bmpFilename, V width, V height)
  : Citizen(srcLoc), _width(width), _height(height),
	_impl(snew<Impl>(srcLoc))
{
	_impl->bmpFilename = bmpFilename;
	_impl->photo = jpegEtc->makePhoto(_width, _height);
}

BmpWriter::~BmpWriter() {
	_impl->photo->saveAsBMP(_impl->bmpFilename);
	jpegEtc->freePhoto(_impl->photo);
	sdeleteC(_impl);
}

void BmpWriter::draw(FloatRaster& floatRaster) {
	Raster_DrawPhoto brush(_width, _height, floatRaster);
	_impl->photo->setAllPixels(brush);
}

void BmpWriter::draw(FloatDotter& floatDotter) {
	float maxCoeff = -1.0;
	for (size_t i = 0; i < floatDotter._numberOfDots; i++) {
		FloatDotter::V w,h;
		float coeff;
		floatDotter.get(coeff,w,h, i);
		if (maxCoeff < coeff) {
			maxCoeff = coeff;
		}
		_impl->photo->setOneGreyPixel(coeff,toIntegral<int>(w), toIntegral<int>(h));
	}
}


void rasterToBmpTest() {
	
	auto xfp = mrcsTestRead();
	auto & xf = *xfp;

	ListOfSizeTRange subset;
	subset.push_back({0,xf.numberOfImages()});

	Images images(SRCLOC);
	xf.getImages(images, subset);

	for (auto ii = images.begin(); ii != images.end(); ii++) {
		const size_t 
			w = toIntegral<size_t>(images.xSize().v()), 
			h = toIntegral<size_t>(images.ySize().v());
		auto& image = images[ii];
		Filename bmpFilename = "C:/TEMP/3DRecon_rasterToBmpTest_";
		bmpFilename += std::to_string(ii.v());
		BmpWriter bmpWriter(SRCLOC, bmpFilename, w, h);
		bmpWriter.draw(image.rawPixels().raster());
	}

	sdelete(xfp);
}

#endif