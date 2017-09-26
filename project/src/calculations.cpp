#include "calculations.h"

using std::get;

ull calculateStandardDeviation(const Image &imgA, const Image &imgB, const short horShift, const short verShift) {
	ull metric = 0;
	const uint rows = imgA.n_rows - abs(verShift),
			   cols = imgA.n_cols - abs(horShift),
			   deltaAhor = (horShift < 0) ? 0 : horShift,
			   deltaAver = (verShift < 0) ? 0 : verShift,
			   deltaBhor = (horShift < 0) ? -horShift : 0,
			   deltaBver = (verShift < 0) ? -verShift : 0;

	for (uint i = 0; i < rows; ++i) {
		for (uint j = 0; j < cols; ++j) {
			metric += pow(static_cast<int>(get<0>(imgA(i + deltaAver, j + deltaAhor))) -
						  static_cast<int>(get<0>(imgB(i + deltaBver, j + deltaBhor))), 2);
		}
	}
	return metric / (imgA.n_rows * imgA.n_cols);
}

ull calculateCrossCorrelation(const Image &imgA, const Image &imgB, const short horShift, const short verShift) {
	ull metric = 0;
	const uint rows = imgA.n_rows - abs(verShift),
			   cols = imgA.n_cols - abs(horShift),
			   deltaAhor = (horShift < 0) ? 0 : horShift,
		   	   deltaAver = (verShift < 0) ? 0 : verShift,
			   deltaBhor = (horShift < 0) ? -horShift : 0,
			   deltaBver = (verShift < 0) ? -verShift : 0;

	for (uint i = 0; i < rows; ++i) {
		for (uint j = 0; j < cols; ++j) {
			metric += get<0>(imgA(i + deltaAver, j + deltaAhor)) *
					  get<0>(imgB(i + deltaBver, j + deltaBhor));
		}
	}
	return metric;
}

std::tuple<short, short> calculateOptimalShift(const Image &imgA, const Image &imgB, metrics metricType) {
	ull optMetrics, metric;
	std::tuple<short, short> optShift;
	short optHorShift, optVerShift;

	if (metricType == STANDARD_DEVIATION) {
		optMetrics = std::numeric_limits<ull>::max();
		for (short horShift = MIN_SHIFT; horShift <= MAX_SHIFT; horShift++)
		{
			for (short verShift = MIN_SHIFT; verShift <= MAX_SHIFT; verShift++)
			{
				metric = calculateStandardDeviation(imgA, imgB, horShift, verShift);
				if (metric < optMetrics)
				{
					optMetrics = metric;
					optHorShift = horShift;
					optVerShift = verShift;
				}
			}
		}
	}
	if (metricType == CROSS_CORRELATION) {
		optMetrics = 0;
		for (short horShift = MIN_SHIFT; horShift <= MAX_SHIFT; horShift++)
		{
			for (short verShift = MIN_SHIFT; verShift <= MAX_SHIFT; verShift++)
			{
				metric = calculateCrossCorrelation(imgA, imgB, horShift, verShift);
				if (metric > optMetrics)
				{
					optMetrics = metric;
					optHorShift = horShift;
					optVerShift = verShift;
				}
			}
		}
	}
	optShift = std::make_tuple(optHorShift, optVerShift);
	return optShift;
}

Image combineImages(const Image &imgA, const Image &imgB, const Image &imgC,
				    const std::tuple<short, short> shiftAB, const std::tuple<short, short> shiftAC) {
	Image resImage(imgA.n_rows, imgA.n_cols);
	uint R, G, B;

	for (int i = 0; i < static_cast<int>(imgA.n_rows); ++i) {
		for (int j = 0; j < static_cast<int>(imgA.n_cols); ++j) {
			B = (i - get<1>(shiftAB) < 0 || i - get<1>(shiftAB) >= static_cast<int>(imgA.n_rows) ||
				 j - get<0>(shiftAB) < 0 || j - get<0>(shiftAB) >= static_cast<int>(imgA.n_cols)) ?
			   	 get<0>(imgA(i, j)) : get<0>(imgB(i - get<1>(shiftAB), j - get<0>(shiftAB)));
			G = get<0>(imgA(i, j));
			R = (i - get<1>(shiftAC) < 0 || i - get<1>(shiftAC) >= static_cast<int>(imgA.n_rows) ||
				 j - get<0>(shiftAC) < 0 || j - get<0>(shiftAC) >= static_cast<int>(imgA.n_cols)) ?
				 get<0>(imgA(i, j)) : get<0>(imgC(i - get<1>(shiftAC), j - get<0>(shiftAC)));
			resImage(i, j) = std::make_tuple(R, G, B);
		}
	}
	return resImage;
}

ConvolutionFilter::ConvolutionFilter(const Matrix<double> &kernel) : radius(kernel.n_rows), filterKernel(kernel) {}

std::tuple<uint, uint, uint> ConvolutionFilter::operator()(const Image &filterZone) const
{
	double R, G, B, SR, SG, SB;

	SR = SG = SB = 0;
	for (uint i = 0; i < radius; ++i) {
		for (uint j = 0; j < radius; ++j) {
			std::tie(R, G, B) = filterZone(i, j);
			SR += R * filterKernel(i, j);
			SG += G * filterKernel(i, j);
			SB += B * filterKernel(i, j);
		}
	}

	shoveInBounds(SR, SG, SB, 0, 255);

	return std::make_tuple(SR, SG, SB);
}

MedianFilter::MedianFilter(const int &Radius) : radius(Radius) {}

std::tuple<uint, uint, uint> MedianFilter::operator()(const Image &filterZone) const
{
	uint R, G, B, median = (radius * radius) / 2;
	std::vector<uint> pixel[3];

	for (uint i = 0; i < radius; ++i) {
		for (uint j = 0; j < radius; ++j) {
			std::tie(R, G, B) = filterZone(i, j);
			pixel[0].push_back(R);
			pixel[1].push_back(G);
			pixel[2].push_back(B);
		}
	}
	for (uint i = 0; i < 3; ++i)
		std::sort(pixel[i].begin(), pixel[i].end());

	return std::make_tuple(pixel[0][median], pixel[1][median], pixel[2][median]);
}

uint find_median_value_of(uint* const hist, const uint stop_sum) {
	uint sum = 0;
	for (uint i = 0; i < BRIGHTNESS_BOUND - 1; ++i) {
		sum += hist[i];
		if (sum >= stop_sum)
			return i;
	}
	return 0;
}

Image add_border(Image src_image, const uint border_radius) {
	Image new_image(src_image.n_rows + 2 * border_radius, src_image.n_cols + 2 * border_radius);
	for (uint i = 0; i < src_image.n_rows; ++i) {
		for (uint j = 0; j < src_image.n_cols; ++j) {
			new_image(i + border_radius, j + border_radius) = src_image(i, j);
		}
	}
	return new_image;
}

Image cut_border(Image src_image, const uint border_radius) {
	if (2 * border_radius + 1 > src_image.n_rows || 2 * border_radius + 1 > src_image.n_rows)
		throw std::string("radius is too big for this picture");
	Image new_image(src_image.n_rows - 2 * border_radius, src_image.n_cols - 2 * border_radius);
	for (uint i = 0; i < new_image.n_rows; ++i) {
		for (uint j = 0; j < new_image.n_cols; ++j) {
			new_image(i, j) = src_image(i + border_radius, j + border_radius);
		}
	}
	return new_image;
}

Image mirror_border(Image src_image, const uint border_radius) {
 // Mirrowing edge squares of image
	for (uint i = 0; i < border_radius; ++i) {
		for (uint j = 0; j < border_radius; ++j) {
		 // Left upper
			src_image(i, j) = src_image(2 * border_radius - 1 - j, 2 * border_radius - 1 - i);
		 // Right upper
			src_image(i, src_image.n_cols - 1 - j) = src_image(src_image.n_rows - 2 * border_radius + j, src_image.n_cols - 2 * border_radius + i);
		 // Left down
			src_image(src_image.n_rows - 1 - i, j) = src_image(src_image.n_rows - 2 * border_radius + j, 2 * border_radius - 1 - i);
		 // Right down
			src_image(src_image.n_rows - 1 - i, src_image.n_cols - 1 - j) = src_image(src_image.n_rows - 2 * border_radius + j, src_image.n_cols - 2 * border_radius + i);
		}
	}
 // Mirrowing side borders of image
	for (uint i = 0; i < src_image.n_rows - 2 * border_radius; ++i) {
		for (uint j = 0; j < border_radius; ++j) {
		 // Left
			src_image(i + border_radius, j) = src_image(i + border_radius, 2 * border_radius - 1 - j);
		 // Right
			src_image(i + border_radius, src_image.n_cols - 1 - j) = src_image(i + border_radius, src_image.n_cols - 2 * border_radius + j);
		}
	}
	for (uint i = 0; i < src_image.n_cols - 2 * border_radius; ++i) {
		for (uint j = 0; j < border_radius; ++j) {
		 // Up
			src_image(j, i + border_radius) = src_image(2 * border_radius - 1 - j, i + border_radius);
		 // Down
			src_image(src_image.n_rows - 1 - j, i + border_radius) = src_image(src_image.n_rows - 2 * border_radius + j, i + border_radius);
		}
	}
	return src_image;
}
