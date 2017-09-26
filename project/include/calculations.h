#pragma once

#include "io.h"

#define MIN_SHIFT -15
#define MAX_SHIFT 15
#define BORDER_PERCENT 0.1 // % of picture to omit in finding border
#define BRIGHTNESS_BOUND 256

typedef unsigned long long ull;

enum metrics { STANDARD_DEVIATION, CROSS_CORRELATION };

ull calculateStandardDeviation(const Image &imgA, const Image &imgB, const short horShift, const short verShift);
ull calculateCrossCorrelation(const Image &imgA, const Image &imgB, const short horShift, const short verShift);
std::tuple<short, short> calculateOptimalShift(const Image &imgA, const Image &imgB, metrics metricType);
Image combineImages(const Image &imgA, const Image &imgB, const Image &imgC,
					const std::tuple<short, short> shiftAB, const std::tuple<short, short> shiftAC);
uint find_median_value_of(uint * const hist, const uint stop_sum);
Image add_border(Image src_image, const uint border_radius);
Image cut_border(Image src_image, const uint border_radius);
Image mirror_border(Image src_image, const uint border_radius);
template<typename T, typename R> void shoveInBounds(T& a, T& b, T& c, const R Min, const R Max) {
	T min = static_cast<T>(Min);
	T max = static_cast<T>(Max);
	a = (a > max) ? max : (a < min) ? min : a;
	b = (b > max) ? max : (b < min) ? min : b;
	c = (c > max) ? max : (c < min) ? min : c;
}

class ConvolutionFilter
{
public:
	ConvolutionFilter(const Matrix<double> &kernel);
	std::tuple<uint, uint, uint> operator()(const Image &filterZone) const;
	const uint radius;
private:
	Matrix<double> filterKernel;
};

class MedianFilter
{
public:
	MedianFilter(const int &Radius);
	std::tuple<uint, uint, uint> operator()(const Image &filterZone) const;
	const uint radius;
};
