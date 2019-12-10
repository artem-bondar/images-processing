#include <string>
#include "align.h"
#include "calculations.h"

Image align(Image srcImage, bool isPostprocessing, string postprocessingType, double fraction, bool isMirror,
            bool isInterp, bool isSubpixel, double subScale)
{
	const uint res_rows = srcImage.n_rows / 3, res_cols = srcImage.n_cols,
			   border_rows = res_rows * BORDER_PERCENT, border_cols = res_cols * BORDER_PERCENT;
	std::tuple<short, short> optimalShiftGB, optimalShiftGR;
	metrics metricForUse = STANDARD_DEVIATION; // As a result, standard deviation works better on average,
											   // only on 2 of 11 pictures cross-corelation is quite better,
											   // while on another 9 results are totally worse
	Image channel_R, channel_G, channel_B, resImage(res_rows, res_cols);

 // 1-2 rows from bottom can be lost
	channel_B = srcImage.submatrix(0, 0, res_rows, res_cols);
	channel_G = srcImage.submatrix(res_rows, 0, res_rows, res_cols);
	channel_R = srcImage.submatrix(2 * res_rows, 0, res_rows, res_cols);

	optimalShiftGB = calculateOptimalShift(channel_G.submatrix(border_rows, border_cols, res_rows - border_rows, res_cols - border_cols),
		channel_B.submatrix(border_rows, border_cols, res_rows - border_rows, res_cols - border_cols), metricForUse);
	optimalShiftGR = calculateOptimalShift(channel_G.submatrix(border_rows, border_cols, res_rows - border_rows, res_cols - border_cols),
		channel_R.submatrix(border_rows, border_cols, res_rows - border_rows, res_cols - border_cols), metricForUse);

	resImage = combineImages(channel_G, channel_B, channel_R, optimalShiftGB, optimalShiftGR);

	if (isPostprocessing)
	{
		if (postprocessingType == "--gray-world")
			resImage = gray_world(resImage);
		else if (postprocessingType == "--unsharp")
		{
			if (isMirror)
			{
				resImage = add_border(resImage, 3);
				resImage = mirror_border(resImage, 3);
			}
			resImage = unsharp(resImage);
			if (isMirror)
			{
				resImage = cut_border(resImage, 3);
			}
		}
		else if (postprocessingType == "--autocontrast")
			resImage = autocontrast(resImage, fraction);
		else if (postprocessingType == "--white-balance")
			resImage = white_balance(resImage);
	}

    return resImage;
}

Image sobel_x(Image src_image) {
    Matrix<double> kernel = {{-1, 0, 1},
                             {-2, 0, 2},
                             {-1, 0, 1}};
    return custom(src_image, kernel);
}

Image sobel_y(Image src_image) {
    Matrix<double> kernel = {{ 1,  2,  1},
                             { 0,  0,  0},
                             {-1, -2, -1}};
    return custom(src_image, kernel);
}

Image unsharp(Image src_image) {
	Matrix<double> kernel = { { -1. / 6, -2. / 3, -1. / 6 },
							  { -2. / 3, 13. / 3, -2. / 3 },
							  { -1. / 6, -2. / 3, -1. / 6 } };
    return custom(src_image, kernel);
}

Image gray_world(Image src_image) {
	double S, SR, SB, SG, R, G, B;

	S = SR = SB = SG = 0;

	for (uint i = 0; i < src_image.n_rows; ++i) {
		for (uint j = 0; j < src_image.n_cols; ++j) {
			SR += get<0>(src_image(i, j));
			SG += get<1>(src_image(i, j));
			SB += get<2>(src_image(i, j));
		}
	}

	SR /= (src_image.n_rows * src_image.n_cols);
	SG /= (src_image.n_rows * src_image.n_cols);
	SB /= (src_image.n_rows * src_image.n_cols);
	S = (SR + SG + SB) / 3;

	for (uint i = 0; i < src_image.n_rows; ++i) {
		for (uint j = 0; j < src_image.n_cols; ++j) {
			R = (get<0>(src_image(i, j)) * S / SR);
			G = (get<1>(src_image(i, j)) * S / SG);
			B = (get<2>(src_image(i, j)) * S / SB);
			shoveInBounds(R, G, B, 0, 255);
			src_image(i, j) = std::make_tuple(R, G, B);
		}
	}

    return src_image;
}

Image resize(Image src_image, double scale) {
    return src_image;
}

Image custom(Image src_image, Matrix<double> kernel) {
    return src_image.unary_map(ConvolutionFilter(kernel));
}

Image autocontrast(Image src_image, double fraction) {
	// fraction is given from [0.0, 0.4],
	// meanwhile it was checked in main(),
	// it wasn't checked in parse_args(),
	// so I added that check there.

	uint R, G, B, yMin = 0, yMax = BRIGHTNESS_BOUND - 1, brightnessHistogramm[BRIGHTNESS_BOUND] = { 0 };
	int omitMinBrightnessAmount, omitMaxBrightnessAmount; 
	double rCoeff = 0.2125, gCoeff = 0.7154, bCoeff = 0.0721,
		   multiplier;

	for (uint i = 0; i < src_image.n_rows; ++i) {
		for (uint j = 0; j < src_image.n_cols; ++j) {
			std::tie(R, G, B) = src_image(i, j); // quite safe cast, never gets > 255 + epsilon
			brightnessHistogramm[static_cast<int>(rCoeff * R + gCoeff * G + bCoeff * B)]++;
		}
	}
	omitMinBrightnessAmount = omitMaxBrightnessAmount = (src_image.n_rows * src_image.n_cols) * fraction;
	for (uint i = 0; i < BRIGHTNESS_BOUND; i++)
	{
		omitMinBrightnessAmount -= brightnessHistogramm[i];
		if (omitMinBrightnessAmount < 0)
		{// special bound magic
			yMin = i ? i - 1 : i;
			break;
		}
	}
	for (int i = BRIGHTNESS_BOUND - 1; i >= 0; i--)
	{
		omitMaxBrightnessAmount -= brightnessHistogramm[i];
		if (omitMaxBrightnessAmount < 0)
		{// special bound magic
			yMax = (i == BRIGHTNESS_BOUND - 1) ? BRIGHTNESS_BOUND - 1 : i + 1;
			break;
		}
	}
	
	if (yMin != yMax)
 // else means that all pixels are on same brightness, nothing to do there
	{									   // cast requered for accuracy or will receive rounded int
		multiplier = (BRIGHTNESS_BOUND - 1) / static_cast<double>(yMax - yMin);
		for (uint i = 0; i < src_image.n_rows; ++i) {
			for (uint j = 0; j < src_image.n_cols; ++j) {
				std::tie(R, G, B) = src_image(i, j);
			 // uint - uint < 0 not to become white, when black needed
				R = ((R > yMin) ? (R - yMin) : 0) * multiplier;
				G = ((G > yMin) ? (G - yMin) : 0) * multiplier;
				B = ((B > yMin) ? (B - yMin) : 0) * multiplier;
				shoveInBounds(R, G, B, 0, 255);
				src_image(i, j) = std::make_tuple(R, G, B);
			}
		}
	}
	return src_image;
}

Image white_balance(Image src_image) {
 // As said in the article, there will be 3 stages:
 // 1. White Object Purification
 // 2. White Point Detection
 // 3. White Balance Adjustment
 //
 // Some copypaste from article will be added as commentary

	bool Rflag = false, Gflag = false, Bflag = false;
	int probableWhitePixels = 0,
		referenceWhitePixels = 0,
		colorcast = 0,
		R, G, B,
		Y, Cr, Cb,
		Rmin = 0, Gmin = 0, Bmin = 0,
		Rmax = BRIGHTNESS_BOUND - 1, Gmax = BRIGHTNESS_BOUND - 1, Bmax = BRIGHTNESS_BOUND - 1,
		Rhist[BRIGHTNESS_BOUND] = { 0 }, Ghist[BRIGHTNESS_BOUND] = { 0 }, Bhist[BRIGHTNESS_BOUND] = { 0 },
		omitMinRpixels, omitMaxRpixels, omitMinGpixels, omitMaxGpixels, omitMinBpixels, omitMaxBpixels,
		Rsum = 0, Gsum = 0, Bsum = 0,
		YSum = 0, CrSum = 0, CbSum = 0;
	double
	 // fraction of histogram
		fraction = 0.05,
	 // Coefficients for RGB/YCrCb conversion according to Rec. 601
		Kr = 0.299, Kg = 0.587, Kb = 0.114,
		Rcoeff, Gcoeff, Bcoeff,
		Rscale, Gscale, Bscale,
		Rfactor = 1, Gfactor = 1, Bfactor = 1,
		Rgwa, Ggwa, Bgwa,
		Ravg, Gavg, Bavg,
		Rw, Gw, Bw,
		YAvg, CrAvg, CbAvg,
	 // 0, 4, 4 - values out of range of that the probable white pixel takes
		YBright = 0, CrBright = 4, CbBright = 4,
		Yl, Crl, Cbl, // lower bound
		Yu, Cru, Cbu, // upper bound
		Yw;
	Image Ihist(src_image.n_rows, src_image.n_cols);
	Matrix<std::tuple<int, int, int>> IyccHist(src_image.n_rows, src_image.n_cols);

 // 1. White Object Purification
 // Apply the histogram equalization to each RGB channel separately
 //
 // Building histogram for each channel
	for (uint i = 0; i < src_image.n_rows; ++i) {
		for (uint j = 0; j < src_image.n_cols; ++j) {
			std::tie(R, G, B) = src_image(i, j);
			Rhist[R]++;
			Ghist[G]++;
			Bhist[B]++;
		}
	}
 // Finding RGB min and max values
	omitMinRpixels = omitMaxRpixels = omitMinGpixels = omitMaxGpixels = omitMinBpixels = omitMaxBpixels = (src_image.n_rows * src_image.n_cols) * fraction;
	for (uint i = 0; i < BRIGHTNESS_BOUND; i++)
	{
		omitMinRpixels -= Rhist[i];
		omitMinGpixels -= Ghist[i];
		omitMinBpixels -= Bhist[i];
		if (!Rflag && omitMinRpixels < 0)
		{
			Rmin = i;
			Rflag = true;
		}
		if (!Gflag && omitMinGpixels < 0)
		{
			Gmin = i;
			Gflag = true;
		}
		if (!Bflag && omitMinBpixels < 0)
		{
			Bmin = i;
			Bflag = true;
		}
		if (Rflag && Gflag && Bflag)
			break;
	}
	Rflag = Gflag = Bflag = false;
	for (int i = BRIGHTNESS_BOUND - 1; i >= 0; i--)
	{
		omitMaxRpixels -= Rhist[i];
		omitMaxGpixels -= Ghist[i];
		omitMaxBpixels -= Bhist[i];
		if (!Rflag && omitMaxRpixels < 0)
		{
			Rmax = i;
			Rflag = true;
		}
		if (!Gflag && omitMaxGpixels < 0)
		{
			Gmax = i;
			Gflag = true;
		}
		if (!Bflag && omitMaxBpixels < 0)
		{
			Bmax = i;
			Bflag = true;
		}
		if (Rflag && Gflag && Bflag)
			break;
	}
 // Calculating coefficients
	Rcoeff = static_cast<double>(BRIGHTNESS_BOUND - 1) / (Rmax - Rmin);
	Gcoeff = static_cast<double>(BRIGHTNESS_BOUND - 1) / (Gmax - Gmin);
	Bcoeff = static_cast<double>(BRIGHTNESS_BOUND - 1) / (Bmax - Bmin);
 // Histogram equalization
	for (uint i = 0; i < src_image.n_rows; ++i) {
		for (uint j = 0; j < src_image.n_cols; ++j) {
			std::tie(R, G, B) = src_image(i, j);
			R = (R - Rmin) * Rcoeff;
			G = (G - Gmin) * Gcoeff;
			B = (B - Bmin) * Bcoeff;
			shoveInBounds(R, G, B, 0, 255);
			Ihist(i, j) = std::make_tuple(R, G, B);
		}
	}
 // Cover image space form RGB to YCrCb, I Hist(Y Hist, Cr Hist, Cb Hist)
	for (uint i = 0; i < Ihist.n_rows; ++i) {
		for (uint j = 0; j < Ihist.n_cols; ++j) {
			std::tie(R, G, B) = Ihist(i, j);
			Y = Kr * R + Kg * G + Kb * B;
			Cr = R - Y;
			Cb = B - Y;
			IyccHist(i, j) = std::make_tuple(Y, Cr, Cb);
		}
	}
 // 2. White Point Detection
 // Detection of the probable white pixels
	for (uint i = 0; i < IyccHist.n_rows; ++i) {
		for (uint j = 0; j < IyccHist.n_cols; ++j) {
			std::tie(Y, Cr, Cb) = IyccHist(i, j);
			if (Y >= 210 &&
				(Cr >= -3 && Cr <= 3) &&
				(Cb >= -3 && Cb <= 3))
			{
				probableWhitePixels++;
				YSum += Y;
				CrSum += Cr;
				CbSum += Cb;
			 // Among the probably white pixels detect the brightest pixel
		   	 // with highest Y Hist value, and Cr Hist, Cb Hist values should be closest to zero
				if (Y > YBright || (Y == static_cast<int>(YBright) && abs(Cr) + abs(Cb) < abs(CrBright) + abs(CbBright)))
				{
					YBright = Y;
					CrBright = Cr;
					CbBright = Cb;
				}
			}
		}
	}
 // If no probable white pixel is detected, we stop for the white balancing process
	if (!probableWhitePixels)
		return src_image;
 // Calculate the average probable white pixels
	YAvg = static_cast<double>(YSum) / static_cast<double>(probableWhitePixels);
	CrAvg = static_cast<double>(CrSum) / static_cast<double>(probableWhitePixels);
	CbAvg = static_cast<double>(CbSum) / static_cast<double>(probableWhitePixels);

 // Select the minimum and maximum between YBright, CrBright, CbBright,
 //	and YAvg, CrAvg, CbAvg respectively
	Yl = std::min(YBright, YAvg);
	Yu = std::max(YBright, YAvg);
	Crl = std::min(CrBright, CrAvg);
	Cru = std::max(CrBright, CrAvg);
	Cbl = std::min(CbBright, CbAvg);
	Cbu = std::max(CbBright, CbAvg);
 // Calculate the average of reference white pixels as W(Rw , Gw , Bw)
	for (uint i = 0; i < IyccHist.n_rows; ++i) {
		for (uint j = 0; j < IyccHist.n_cols; ++j) {
			std::tie(Y, Cr, Cb) = IyccHist(i, j);
			if ((Y >= Yl && Y <= Yu) &&
				(Cr >= Crl && Cr <= Cru) &&
				(Cb >= Cbl && Cb <= Cbu))
			{
				std::tie(R, G, B) = src_image(i, j);
				referenceWhitePixels++;
				Rsum += R;
				Gsum += G;
				Bsum += B;
			}
		}
	}

	Rw = static_cast<double>(Rsum) / static_cast<double>(referenceWhitePixels);
	Gw = static_cast<double>(Gsum) / static_cast<double>(referenceWhitePixels);
	Bw = static_cast<double>(Bsum) / static_cast<double>(referenceWhitePixels);

 // 3. White Balance Adjustment
 // Calculate the scale factors
	Yw = Kr * Rw + Kg * Gw + Kb * Bw;

	Rscale = Yw / Rw;
	Gscale = Yw / Gw;
	Bscale = Yw / Bw;

 // Inverse conversion: YCrCb to RGB
	Ravg = YAvg + CrAvg;
	Gavg = YAvg - (Kb / Kg) * CbAvg - (Kr / Kg) * CrAvg;
	Bavg = YAvg + CbAvg;

	Rgwa = YAvg / Ravg;
	Ggwa = YAvg / Gavg;
	Bgwa = YAvg / Bavg;

 // Calculate the color cast type 'colorcast' as 0, 1, 2, or 3 (default value = 0)
 // Calculate scale factors
	if (Bavg + 3 >= Gavg && Bavg > Ravg) // bluish cast
	{
		colorcast = 1;
		Rfactor = Rscale;
		Gfactor = Gscale;
		Bfactor = Bgwa;
	}
	if (Gavg + 3 > Ravg && Ravg > Bavg) // greenish cast
	{
		colorcast = 2;
		Rfactor = Rscale;
		Gfactor = Ggwa;
		Bfactor = Bscale;
	}
	if (Ravg > Gavg && Gavg > Bavg) // reddish cast
	{
		colorcast = 3;
		Rfactor = Rgwa;
		Gfactor = Gscale;
		Bfactor = Bscale;
	}

	if (!colorcast) // no color cast needed
		return src_image;

 // Apply the scale factors
	for (uint i = 0; i < src_image.n_rows; ++i) {
		for (uint j = 0; j < src_image.n_cols; ++j) {
			std::tie(R, G, B) = src_image(i, j);
			R *= Rfactor;
			G *= Gfactor;
			B *= Bfactor;
			shoveInBounds(R, G, B, 0, 255);
			src_image(i, j) = std::make_tuple(R, G, B);
		}
	}
	return src_image;
}

Image gaussian(Image src_image, double sigma, int radius)  {
    return src_image;
}

Image gaussian_separable(Image src_image, double sigma, int radius) {
    return src_image;
}

Image median(Image src_image, int radius) {
	return src_image.unary_map(MedianFilter(radius));
}

Image median_linear(Image src_image, int radius) {
	bool moveRight = true;
	uint R, G, B,
		Rmed, Gmed, Bmed,
		Rhist[BRIGHTNESS_BOUND] = { 0 }, Ghist[BRIGHTNESS_BOUND] = { 0 }, Bhist[BRIGHTNESS_BOUND] = { 0 };
	const uint stop_sum = 2 * radius * (radius + 1) + 1;
	int deltaJadd, deltaJdel;
	Image dst_image = src_image.deep_copy();
	if (2 * static_cast<uint>(radius) + 1 > src_image.n_rows)
		throw std::string("radius is too big for this picture");

 // Initialization of histogram
	for (int i = 0; i < 2 * radius + 1; ++i) {
		for (int j = 0; j < 2 * radius + 1; ++j) {
			std::tie(R, G, B) = src_image(i, j);
			Rhist[R]++;
			Ghist[G]++;
			Bhist[B]++;
		}
	}
 // Map:
 //                 Moving direction
 // 000000...00000        ->
 // 000XCC...CCY00        <-
 // 000YCC...CCZ00        ->
 // 000ZCC...CCY00        <-
 // 
 // Legend:
 // 0 - not counted
 // C - counted in J-circle
 // Y - counted in I-circle before move down
 // Z - counted in I-circle after move down
 // X - counted here:
	Rmed = find_median_value_of(Rhist, stop_sum);
	Gmed = find_median_value_of(Ghist, stop_sum);
	Bmed = find_median_value_of(Bhist, stop_sum);
	dst_image(radius, radius) = std::make_tuple(Rmed, Gmed, Bmed);

	uint j = 0;
	for (uint i = radius; i < src_image.n_rows - radius; ++i, moveRight = !moveRight) {
		deltaJadd = moveRight ? radius : -radius;
		deltaJdel = moveRight ? -radius - 1 : radius + 1;
		for (uint J = radius + 1; J < src_image.n_cols - radius - 1; ++J) {
		 // Direction
			j = moveRight ? J : src_image.n_cols - 1 - J;
		 // Move histogram right or left
		 // C - counted here:
			for (int k = -radius; k <= radius; ++k) {
				std::tie(R, G, B) = src_image(i + k, j + deltaJadd);
				Rhist[R]++;
				Ghist[G]++;
				Bhist[B]++;
				std::tie(R, G, B) = src_image(i + k, j + deltaJdel);
				Rhist[R]--;
				Ghist[G]--;
				Bhist[B]--;
			}
			Rmed = find_median_value_of(Rhist, stop_sum);
			Gmed = find_median_value_of(Ghist, stop_sum);
			Bmed = find_median_value_of(Bhist, stop_sum);
			dst_image(i, j) = std::make_tuple(Rmed, Gmed, Bmed);
		}
	 // Y - counted here:
		j = moveRight ? src_image.n_cols - radius - 1: radius;
		for (int k = -radius; k <= radius; ++k) {
			std::tie(R, G, B) = src_image(i + k, j + deltaJadd);
			Rhist[R]++;
			Ghist[G]++;
			Bhist[B]++;
			std::tie(R, G, B) = src_image(i + k, j + deltaJdel);
			Rhist[R]--;
			Ghist[G]--;
			Bhist[B]--;
		}
		Rmed = find_median_value_of(Rhist, stop_sum);
		Gmed = find_median_value_of(Ghist, stop_sum);
		Bmed = find_median_value_of(Bhist, stop_sum);
		dst_image(i, j) = std::make_tuple(Rmed, Gmed, Bmed);
	 // Move histogram down
	 // Z - counted here:
		if (i == src_image.n_rows - radius - 1)
			break;
		for (int k = -radius; k <= radius; ++k) {
			std::tie(R, G, B) = src_image(i + radius + 1, j + k);
			Rhist[R]++;
			Ghist[G]++;
			Bhist[B]++;
			std::tie(R, G, B) = src_image(i - radius, j + k);
			Rhist[R]--;
			Ghist[G]--;
			Bhist[B]--;
		}
		Rmed = find_median_value_of(Rhist, stop_sum);
		Gmed = find_median_value_of(Ghist, stop_sum);
		Bmed = find_median_value_of(Bhist, stop_sum);
		dst_image(i, j) = std::make_tuple(Rmed, Gmed, Bmed);
	}
    return dst_image;
}

Image median_const(Image src_image, int radius) {
	bool moveRight = false;
	uint R, G, B,
		Rmed, Gmed, Bmed,
		size = 2 * radius + 1,
		Rhist[BRIGHTNESS_BOUND] = { 0 }, Ghist[BRIGHTNESS_BOUND] = { 0 }, Bhist[BRIGHTNESS_BOUND] = { 0 },
		columnsRhist[src_image.n_cols][BRIGHTNESS_BOUND] = { { 0 } },
		columnsGhist[src_image.n_cols][BRIGHTNESS_BOUND] = { { 0 } },
		columnsBhist[src_image.n_cols][BRIGHTNESS_BOUND] = { { 0 } };
	const uint stop_sum = 2 * radius * (radius + 1) + 1;
	int deltaJadd, deltaJdel;
	Image dst_image = src_image.deep_copy();
	if (2 * static_cast<uint>(radius) + 1 > src_image.n_rows)
		throw std::string("radius is too big for this picture");

 // Initialization of columns histograms
	for (uint i = 0; i < size; ++i) {
		for (uint j = 0; j < src_image.n_cols; ++j) {
			std::tie(R, G, B) = src_image(i, j);
			columnsRhist[j][R]++;
			columnsGhist[j][G]++;
			columnsBhist[j][B]++;
		}
	}
 // Initialization of kernel histogram
	for (uint i = 0; i < size; ++i) {
		for (uint j = 0; j < BRIGHTNESS_BOUND - 1; ++j) {
			Rhist[j] += columnsRhist[i][j];
			Ghist[j] += columnsGhist[i][j];
			Bhist[j] += columnsBhist[i][j];
		}
	}
	// Map:
	//                 Moving direction
	// 000000...00000        ->
	// 000XFF...FFF00        <-
	// 000CCC...CCY00        ->
	// 000YCC...CCC00        <-
	// 
	// Legend:
	// 0 - not counted
	// F - counted in first row circle
	// C - counted in another row circle
	// Y - counted in I-circle before move down
	// X - counted here:
	Rmed = find_median_value_of(Rhist, stop_sum);
	Gmed = find_median_value_of(Ghist, stop_sum);
	Bmed = find_median_value_of(Bhist, stop_sum);
	dst_image(radius, radius) = std::make_tuple(Rmed, Gmed, Bmed);

 // Because of the same columns histograms position to kernel histogram
 // First row is calculated in another method, without first step
	deltaJadd = radius;
	deltaJdel = -radius - 1;
	for (int i = radius + 1; i < static_cast<int>(src_image.n_cols) - radius; ++i) {
	 // Move histogram right
	 // F - counted here:
		for (int j = 0; j < BRIGHTNESS_BOUND - 1; ++j) {
			Rhist[j] += columnsRhist[i + deltaJadd][j];
			Ghist[j] += columnsGhist[i + deltaJadd][j];
			Bhist[j] += columnsBhist[i + deltaJadd][j];
			Rhist[j] -= columnsRhist[i + deltaJdel][j];
			Ghist[j] -= columnsGhist[i + deltaJdel][j];
			Bhist[j] -= columnsBhist[i + deltaJdel][j];
		}
		Rmed = find_median_value_of(Rhist, stop_sum);
		Gmed = find_median_value_of(Ghist, stop_sum);
		Bmed = find_median_value_of(Bhist, stop_sum);
		dst_image(radius, i) = std::make_tuple(Rmed, Gmed, Bmed);
	}
 // First Y - counted here:
 // Move histograms down, both kernel and columns
	if (radius == static_cast<int>(src_image.n_rows) - radius - 1)
		return dst_image;
	for (int k = -radius; k <= radius; ++k) {
		std::tie(R, G, B) = src_image(radius + 1, src_image.n_cols - radius - 1 + k);
		columnsRhist[src_image.n_cols - radius - 1 + k][R]++;
		columnsGhist[src_image.n_cols - radius - 1 + k][G]++;
		columnsBhist[src_image.n_cols - radius - 1 + k][B]++;
		Rhist[R]++;
		Ghist[G]++;
		Bhist[B]++;
		std::tie(R, G, B) = src_image(0, src_image.n_cols - radius - 1 + k);
		columnsRhist[src_image.n_cols - radius - 1 + k][R]--;
		columnsGhist[src_image.n_cols - radius - 1 + k][G]--;
		columnsBhist[src_image.n_cols - radius - 1 + k][B]--;
		Rhist[R]--;
		Ghist[G]--;
		Bhist[B]--;
	}
	Rmed = find_median_value_of(Rhist, stop_sum);
	Gmed = find_median_value_of(Ghist, stop_sum);
	Bmed = find_median_value_of(Bhist, stop_sum);
	dst_image(radius + 1, src_image.n_cols - radius - 1) = std::make_tuple(Rmed, Gmed, Bmed);

 // Default cycle and algorithm in two steps starts from here:
	uint j = 0;
	for (uint i = radius + 1; i < src_image.n_rows - radius; ++i, moveRight = !moveRight) {
		deltaJadd = moveRight ? radius : -radius;
		deltaJdel = moveRight ? -radius - 1 : radius + 1;
		for (uint J = radius + 1; J < src_image.n_cols - radius; ++J) {
		 // Direction
			j = moveRight ? J : src_image.n_cols - 1 - J;
		 // C - counted here:
		 // Step 1 - update proper column histogram
			std::tie(R, G, B) = src_image(i + radius, j + deltaJadd);
			columnsRhist[j + deltaJadd][R]++;
			columnsGhist[j + deltaJadd][G]++;
			columnsBhist[j + deltaJadd][B]++;
			std::tie(R, G, B) = src_image(i - radius - 1, j + deltaJadd);
			columnsRhist[j + deltaJadd][R]--;
			columnsGhist[j + deltaJadd][G]--;
			columnsBhist[j + deltaJadd][B]--;
		 // Step 2 - move kernel histogram right or left
			for (int k = 0; k < BRIGHTNESS_BOUND - 1; ++k) {
				Rhist[k] += columnsRhist[j + deltaJadd][k];
				Ghist[k] += columnsGhist[j + deltaJadd][k];
				Bhist[k] += columnsBhist[j + deltaJadd][k];
				Rhist[k] -= columnsRhist[j + deltaJdel][k];
				Ghist[k] -= columnsGhist[j + deltaJdel][k];
				Bhist[k] -= columnsBhist[j + deltaJdel][k];
			}
			Rmed = find_median_value_of(Rhist, stop_sum);
			Gmed = find_median_value_of(Ghist, stop_sum);
			Bmed = find_median_value_of(Bhist, stop_sum);
			dst_image(i, j) = std::make_tuple(Rmed, Gmed, Bmed);
		}
	 // Y - counted here:
	 // Move histograms down, both kernel and columns
		if (i == src_image.n_rows - radius - 1)
			break;
		j = moveRight ? src_image.n_cols - radius - 1 : radius;
		for (int k = -radius; k <= radius; ++k) {
			std::tie(R, G, B) = src_image(i + radius + 1, j + k);
			columnsRhist[j + k][R]++;
			columnsGhist[j + k][G]++;
			columnsBhist[j + k][B]++;
			Rhist[R]++;
			Ghist[G]++;
			Bhist[B]++;
			std::tie(R, G, B) = src_image(i - radius, j + k);
			columnsRhist[j + k][R]--;
			columnsGhist[j + k][G]--;
			columnsBhist[j + k][B]--;
			Rhist[R]--;
			Ghist[G]--;
			Bhist[B]--;
		}
		Rmed = find_median_value_of(Rhist, stop_sum);
		Gmed = find_median_value_of(Ghist, stop_sum);
		Bmed = find_median_value_of(Bhist, stop_sum);
		dst_image(radius, i) = std::make_tuple(Rmed, Gmed, Bmed);
	}
    return dst_image;
}

Image canny(Image src_image, int threshold1, int threshold2) {
    return src_image;
}
