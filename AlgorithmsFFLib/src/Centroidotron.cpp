//
// Created by anichols on 3/12/24.
//

#include "Centroidotron.h"

#include "EigenKernelUtils.h"
#include "EigenSparseUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MathUtils.h"
#include "MsUtils.h"
#include "ParallelUtils.h"
#include "PeakIntegratomatic.h"

#include <boost/range/iterator_range.hpp>

#include <sort_together.hpp>

Err Centroidotron::init(
        double peakWidth,
        int hashingPrecision,
        int filterLength
        ) {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(
            peakWidth,
            0.0,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;
    m_peakWidth = peakWidth;

    e = ErrorUtils::isAboveThreshold(
            filterLength,
            0,
            ErrorUtilsParam::ExcludeThreshold
    ); ree;
    m_filterLength = filterLength;

    e = ErrorUtils::isAboveThreshold(
            hashingPrecision,
            0,
            ErrorUtilsParam::ExcludeThreshold); ree;
    m_hashingPrecision = hashingPrecision;

    ERR_RETURN
}

ScanPoints Centroidotron::runningAverage(const ScanPoints &scanPoints) {
    return ScanPoints();
}

Err Centroidotron::smoothIntensities(
        const QVector<double> &intzVals,
        QVector<double> *intzValsSmoothed
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(intzVals); ree;

    Eigen::VectorX<double> vec = EigenUtils::convertQVectorToEigenVector(intzVals);

    const Eigen::VectorXd mexicanHatKernel = EigenKernelUtils::buildMexicanHatFilter1D<double>(9, m_peakWidth);

    const Eigen::VectorXd vecSmoothed = EigenKernelUtils::convolveVectorWithKernel(vec, mexicanHatKernel);

    *intzValsSmoothed = EigenUtils::convertEigenVectorToQVector(vecSmoothed);

    ERR_RETURN
}

Err Centroidotron::smoothIntensities(
        const ScanPoints &scanPoints,
        ScanPoints *scanPointsSmoothed
        ) {

    ERR_INIT

    Eigen::SparseVector<float> vec;
    e = EigenSparseUtils::pointsToSparseVector(
            scanPoints,
            1200.0f,
            m_hashingPrecision,
            &vec
            ); ree

    const Eigen::VectorX<float> mexicanHatKernel = EigenKernelUtils::buildMexicanHatFilter1D<float>(
            m_filterLength,
            m_peakWidth
            );

    int filterLen = 3;
    Eigen::VectorX<float> movingAverage(filterLen);
    movingAverage.setOnes();
    movingAverage /= static_cast<float>(filterLen - 1);

    const Eigen::VectorX<float> gaussKernel = EigenKernelUtils::buildGaussianFilter1D<float>(
            m_filterLength,
            m_peakWidth
    );

    Eigen::SparseVector<float> vecSmoothed = vec;
    for (int i = 0; i < 1; i++) {
        vecSmoothed = EigenKernelUtils::convolveVectorWithKernel(
                vecSmoothed,
                gaussKernel
        );
    }

    e = EigenSparseUtils::sparseVectorToPoints(
            vecSmoothed,
            m_hashingPrecision,
            scanPointsSmoothed
            ); ree;

    ERR_RETURN
}

namespace {

    double minSnr_;
    int fixedPeaksKeep_;
    double mzTol_;
    bool centroid_;
    int nScales = 10;
    std::vector<double> scalings;

    // Function for determining the scales we want to sample for the CWT calculation
    //
    // This function assumes the profile data is approximately (locally) equally spaced. If the spacing
    // is much different then issues may arise. If exceptions are being thrown this would be the place to
    // start looking. Check the Xspacing. For instance, if you pass a peak list to this function you'll
    // get some unpredictable behavior.
    void getScales(
            const std::vector <double> &mzData,
            const std::vector <double> &intensityData,
            std::vector<std::vector<std::vector<int>>> &nPoints,
            std::vector<double> &widths
            ) {

        int mzLength = static_cast<int>(mzData.size());
        std::vector <double> Xspacing(mzLength,0.0);

        double lastXspacing = 0;

        // get the sampling rate as a function of m/z
        // this assumes that zero-intensity points flank the "islands" of data
        for (int i = 1; i < mzLength - 1; i++) {

            if ( intensityData[i] != 0.0 || intensityData[i-1] != 0.0 ) {
                Xspacing[i] = mzData[i] - mzData[i-1];

                if (Xspacing[i] <= 0) {
                    throw std::runtime_error("[CwtPeakDetector::getScales] m/z profile data are unsorted or contain duplicates");
                }
                if (Xspacing[i] > 10) {
                    qDebug() << "[CwtPeakDetector::getScales] m/z profile data seems to lack flanking zeros between peak profiles";
                    throw std::runtime_error("[CwtPeakDetector::getScales] m/z profile data seems to lack flanking zeros between peak profiles");
                }

                lastXspacing = Xspacing[i];
            }
            else {
                Xspacing[i] = lastXspacing;
            }
        }

        Xspacing[0] = Xspacing[1];
        Xspacing[mzLength-1] = Xspacing[mzLength-2];

        // I'm going to take a small average around each point b/c I worry that
        // there may be cases where the Xspacing isn't as consistent as the files I've tested.
        // This might also be dangerous if a spectrum is extremely sparse and neighboring points
        // are separated by several hundred m/z, for example. However, the sampling rate is different
        // by about a factor of 3-4 from low to high m/z. I think we should be okay here.
        int window_size = 10;
        if (window_size > mzLength) {
            window_size = mzLength / 2;
        }

        int hf_window = window_size / 2;
        window_size = 2 * hf_window; // ensures consistency since original window_size could be odd

        for (int i = 1; i < mzLength - 1; i++) {

            int scalesToInclude = nScales;
            if ( intensityData[i] < 0.75*intensityData[i-1] || intensityData[i] < 0.75*intensityData[i+1] ) scalesToInclude = 1;

            int windowLow = i - hf_window; // inclusive
            if (windowLow<0) {
                windowLow = 0;
            }
            int windowHigh = i + hf_window; // not inclusive
            if (windowHigh>mzLength-1){
                windowHigh = mzLength - 1;
            }

            // need to recalculate the window size in case we are at one of the edges
            int nTot = windowHigh - windowLow; // recall that windowHigh is not inclusive, so leave off the +1

            // get average
            double sum = accumulate( Xspacing.begin() + windowLow, Xspacing.begin() + windowHigh, 0.0);
            widths[i] = sum / double(nTot);

            // figure out the number of wavelet points you'll need to sample for each m/z point
            for (int j=0; j < scalesToInclude; ++j) {

                double maxMZwindow = widths[i] * scalings[j] * 3.0; // this returns the max possible m/z away from the current point where a wavelet may still contribute to the correlation

                int nPointsLeft = 0;
                int counter = i;
                while (--counter >= 0){
                    if ( mzData[i] - mzData[counter] > maxMZwindow ) {
                        break;
                    }
                    nPointsLeft++;
                }

                int nPointsRight = 0;
                counter = i;
                while (++counter < mzLength) {
                    if ( mzData[counter] - mzData[i] > maxMZwindow ) break;
                    nPointsRight++;
                }

               std::cout << "LSDkfjsdlfjdslF " << j << " " << i << " " << nPoints[0][j].size() << " " << mzLength << std::endl;
                nPoints[0][j][i] = nPointsLeft;
//                nPoints[1][j][i] = nPointsRight;
            }
        }
    }

    void ricker2d(
            const std::vector<double> &Pad_mz,
            const int col,
            const int rickerPointsLeft,
            const int rickerPointsRight,
            const double A,
            const double wsq,
            const double centralMZ,
            std::vector <double> &total
            ){

        if (rickerPointsRight + rickerPointsLeft >= (int) total.size()) {
            throw std::runtime_error("[CwtPeakDetector::ricker2d] invalid input parameters");
        }

        for (int i = col-rickerPointsLeft, cnt=0, end = col+rickerPointsRight; i <= end; i++, cnt++) {
            double vec = Pad_mz[i] - centralMZ;
            double tsq = vec * vec;
            double mod = 1.0 - tsq / wsq;
            double gauss = exp( -1.0 * tsq / (2.0 * wsq) );
            total[cnt] = A * mod * gauss;
        }

    }

    void calcCorrelation(
            const std::vector<double> &mz,
            const std::vector<double> &intensity,
            const std::vector<std::vector<std::vector<int>>> & waveletPoints,
            const std::vector<double> &widths,
            std::vector<std::vector<double>> &matrix){

        int mzLength = static_cast<int>(mz.size());

        // setup padded arrays for facilitating the wavelet computation
        int paddingPoints = 500; // should be more than sufficient
        std::vector<double> PadMz(mzLength+2*paddingPoints,0.0); // x array with padding on the front and back
        std::vector<double> PadIntensity(mzLength+2*paddingPoints,0.0); // y array with padding on the front and back
        for (int j = paddingPoints; j < mzLength + paddingPoints; j++) {
            PadMz[j] = mz[j-paddingPoints];
            PadIntensity[j] = intensity[j-paddingPoints];
        }
        std::vector<double> waveletData(paddingPoints,0.0);

        // calculate correlation between wavelet and spectrum data, populate correlation matrix
        for (int i = 0; i<nScales ; i++) {

            double currentScaling = scalings[i];
            int matrixIndex = 2; // This depends on the starting index in the loop immediately below!!

            for (int j = 1; j < mzLength-1; j++) {

                // exclude points that are unlikely to have high correlation
                if ( i > 0 ) { // calculate first row no matter what, as this is important for the noise calculation

                    if ( intensity[j] < 0.75*intensity[j-1] || intensity[j] < 0.75*intensity[j+1] ) {
                        matrixIndex += 2;
                        continue;
                    }
                }

                int nPointsLeft = waveletPoints[0][i][j];
                int nPointsRight = waveletPoints[1][i][j];
                int paddedCol = j + paddingPoints;
                double width = widths[j]*currentScaling;
                double param1 = 2.0 / ( sqrt(3.0 * width) * (sqrt( sqrt(3.141519) ) ) ); // ricker wavelet parameter
                double param2 = width * width; // ricker wavelet parameter

                ricker2d(PadMz, paddedCol, nPointsLeft, nPointsRight, param1, param2, PadMz[paddedCol], waveletData);

                int startPoint = paddedCol - nPointsLeft;
                for (int k = 0; k < nPointsLeft+nPointsRight+1 ; k++)
                    matrix[i][matrixIndex] += waveletData[k] * PadIntensity[startPoint+k];

                matrixIndex++;

                if (j == mzLength-1) {
                    break; // jump out before final iteration
                }

                double moverzShift = ( PadMz[ paddedCol ] + PadMz[ paddedCol + 1 ] ) / 2.0;

                // calculate the correlation at the midpoint between two m/z points, as well. This is why
                // the length of the correlation matrix is (almost) twice that of the number of m/z points.
                ricker2d(PadMz, paddedCol, nPointsLeft, nPointsRight, param1, param2, moverzShift, waveletData);

                for (int k = 0; k < nPointsLeft+nPointsRight+1 ; k++){
                    matrix[i][matrixIndex] += waveletData[k] * PadIntensity[startPoint+k];
                }

                matrixIndex++;

            } // end for over mzPoints

        } // end for over 10 scales

    }

    double scoreAtPercentile(
            const double perc,
            const std::vector<double> &dataSorted,
            const int nTot
            ){

        //using nTot - 1, which is what's done in scipy.scoreatpercentile
        double nBelow = (double)(nTot-1) * perc / 100.0;

        if (ceil(nBelow) == nBelow) {// whole number
            return dataSorted[(int)nBelow]; // consistent with scipy.scoreatpercentile
        }
        else {
            // fraction method used in scipy.scoreatpercentile
            // linear interpolation between two points
            double loInd = floor(nBelow);
            double hiInd = ceil(nBelow);
            double fraction = nBelow - loInd;
            double lo = dataSorted[(int)loInd];
            double hi = dataSorted[(int)hiInd];
            return lo + (hi - lo) * fraction;
        }
    }

    double convertColToMZ(
            const std::vector<double> &mzs,
            const int Col
            ){

        int mapIndex = Col / 2;
        if (Col % 2 == 1) {
            return (mzs[ mapIndex ]+mzs[mapIndex+1]) / 2.0;
        }

        return mzs[ mapIndex ];
    }

    int getColHighBound(
            const std::vector<double> &mzs,
            const double target
            ) {

        const std::vector<double>::const_iterator it = upper_bound(mzs.begin(), mzs.end(), target);
        const std::vector<double>::const_iterator decrementedIt = it - 1; // decrement by one to get the value that is less than target

        return decrementedIt - mzs.begin();
    }

    int getColLowBound(
            const std::vector<double> &mzs,
            const double target
            ){

        const auto it = lower_bound(mzs.begin(),mzs.end(),target);
        return it - mzs.begin();
    }


    void getPeakLines(
            const std::vector<std::vector<double>> &corrMatrix,
            const std::vector<double> &x,
            std::vector<ridgeLine> &allLines,
            std::vector <double> & snrs
            ) {

        int corrMatrixLength = corrMatrix[0].size();

        std::vector<int> colMaxes(corrMatrixLength,0);
        for (int i=0; i<corrMatrixLength; ++i) {

            double corrMax = 0.0;
            for (int j=0; j<nScales; ++j) {

                if ( corrMatrix[j][i] > corrMax ) {
                    corrMax = corrMatrix[j][i];
                    colMaxes[i] = j;
                }
            }
        }

        // step 2, setup bins of 300 points and calculate noise threshold within each bin
        double noise_per = 95.0;
        int window_size = 300;
        if (window_size > corrMatrixLength) window_size = corrMatrixLength / 2;
        int hf_window = window_size / 2;
        window_size = 2 * hf_window; // ensures consistency since original window_size could be odd

        int nNoiseBins = corrMatrixLength / window_size + 1;
        std::vector<double> noises(nNoiseBins,0.0);
        for (int i=0; i < nNoiseBins; ++i) {

            int windowLow = i * window_size; // inclusive
            int windowHigh = windowLow + window_size; // not inclusive
            if ( i == nNoiseBins - 1 ) windowHigh = corrMatrixLength;
            int nTot = windowHigh - windowLow; // don't need +1 because windowHigh is not inclusive

            std::vector<double> sortedData(nTot,0.0);
            for (int j = 0; j < nTot; j++) {
                sortedData[j] = corrMatrix[0][windowLow + j]; // first row of correlation matrix
            }

            // sort correlation data on first row within window using STL sort function
            sort(sortedData.begin(),sortedData.end());

            double noise = scoreAtPercentile( noise_per, sortedData, nTot );
            if ( noise < 1.0 ) noise = 1.0;

            noises[i] = noise;

        }

        std::vector<double> interpolatedXpoints(corrMatrixLength,0.0);
        for (int i=0; i<corrMatrixLength; ++i) {
            interpolatedXpoints[i] = convertColToMZ( x, i );
        }

        // step 3, find local maxima that are separated by at least mzTol_
        for (int i=2; i<corrMatrixLength-2; ++i) {

            double correlationVal = corrMatrix[colMaxes[i]][i];

            if ( correlationVal < corrMatrix[colMaxes[i-1]][i-1] ||
                 correlationVal < corrMatrix[colMaxes[i-2]][i-2] ||
                 correlationVal < corrMatrix[colMaxes[i+1]][i+1] ||
                 correlationVal < corrMatrix[colMaxes[i+2]][i+2]
                 ) {
                continue;
            }

            double mzCol = convertColToMZ( x, i );
            double lowTol = mzCol - mzTol_;
            double highTol = mzCol + mzTol_;

            // get the indices for the lower and upper bounds
            int lowBound = getColLowBound(interpolatedXpoints,lowTol);
            int highBound = getColHighBound(interpolatedXpoints,highTol);

            double maxCorr = 0.0;
            int maxCol = 0;
            for (int j=lowBound; j <= highBound; ++j) {
                int row = colMaxes[j];
                if ( corrMatrix[row][j] > maxCorr ) {
                    maxCorr = corrMatrix[row][j];
                    maxCol = j;
                }
            }

            int noiseBin = maxCol / window_size;
            if ( noiseBin > nNoiseBins - 1 ) noiseBin = nNoiseBins - 1;
            double snr = maxCorr / noises[noiseBin];

            if ( snr < minSnr_ ) continue;

            int nLines = static_cast<int>(allLines.size());
            if ( nLines > 0 )
            {
                double mzNewLine = convertColToMZ(x,maxCol);
                double mzPrevLine = convertColToMZ(x,allLines[nLines-1].Col);
                double mzDiff = mzNewLine - mzPrevLine;
                double corrPrev = corrMatrix[allLines[nLines-1].Row][allLines[nLines-1].Col];
                if ( mzDiff > mzTol_ ) {
                    ridgeLine newLine;
                    newLine.Col = maxCol;
                    newLine.Row = colMaxes[maxCol];
                    allLines.push_back(newLine);
                    snrs.push_back(snr);
                }
                else if ( maxCorr > corrPrev ) {
                    // remove last ridge line
                    allLines.pop_back();
                    snrs.pop_back();
                    // add new ridge line
                    ridgeLine newLine;
                    newLine.Col = maxCol;
                    newLine.Row = colMaxes[maxCol];
                    allLines.push_back(newLine);
                    snrs.push_back(snr);
                }
            }
            else {
                ridgeLine newLine;
                newLine.Col = maxCol;
                newLine.Row = colMaxes[maxCol];
                allLines.push_back(newLine);
                snrs.push_back(snr);
            }
        }
    }

    void refinePeaks(
            const std::vector<double> &noisyX,
            const std::vector<double> &noisyY,
            const std::vector<ridgeLine> &lines,
            const std::vector<double> &widths,
            std::vector<double> &smoothX,
            std::vector <double> &smoothY,
            std::vector<double> &snrs
            ){

        if (lines.size() == 0) {
            return;
        }

        for (int i=0, iend=lines.size(); i < iend; ++i) {

            double mzCol = convertColToMZ(noisyX,lines[i].Col);

            int row = lines[i].Row;
            double currentScaling = scalings[row];
            double offset = currentScaling * widths[lines[i].Col / 2];

            // get the indices for the lower and upper bounds that encapsulate the peak
            int startFittingPoint = getColLowBound(noisyX,mzCol-offset);
            int endFittingPoint = getColHighBound(noisyX,mzCol+offset);

            double maxIntensity = 0.0;
            double intensityAccumulator = 0.0;
            double mzCentroid = 0.0;
            double bestMZ = 0.0;

            // take weighted average of points in the peak to get centroid m/z
            if (centroid_) {
                for (int j = startFittingPoint; j <= endFittingPoint; ++j)
                {
                    intensityAccumulator += noisyY[j];
                    mzCentroid += noisyY[j]*noisyX[j];
                    if (noisyY[j] >= maxIntensity)
                        maxIntensity = noisyY[j];
                }
                bestMZ = mzCentroid / intensityAccumulator;
            }
                // sum up the intensity and find the highest point. If there are multiple
                // points with the maxIntensity value, take the one with the highest m/z.
            else {
                for (int j = startFittingPoint; j <= endFittingPoint; ++j) {
                    intensityAccumulator += noisyY[j];
                    if (noisyY[j] >= maxIntensity)
                    {
                        maxIntensity = noisyY[j];
                        bestMZ = noisyX[j];
                    }
                }
            }

            if (smoothX.empty() || bestMZ != smoothX.back()) {
                smoothX.push_back(bestMZ);
                smoothY.push_back(maxIntensity);
            }

        }

        // This is an intensity threshold filter for removing peaks with intensity of 1
        // and the best matched wavelet at the lowest scale (indicitive of noise)
        for (int k = smoothX.size()-1; k > 0; --k) {
            if ( smoothY[k] < 2.0 && lines[k].Row < 1 ) {
                smoothX.erase(smoothX.begin()+k);
                smoothY.erase(smoothY.begin()+k);
                snrs.erase(snrs.begin()+k);
            }
        }

        // peaks are occasionally out of order
        pwiz::util::sort_together(
                smoothX,
                std::vector<boost::iterator_range<std::vector<double>::iterator>> { smoothY, snrs }
                );

        // possible to list the same peak if two lines are drawn on the same peak
        // and fall back to the same max intensity value (or a very similar max intensity value)
        for (int k = smoothX.size()-1; k > 0; --k)
        {
            if ( smoothX[k] - smoothX[k-1] < mzTol_ )
            {
                int removePeakIndex = smoothY[k] > smoothY[k-1] ? k-1 : k;
                smoothX.erase(smoothX.begin()+removePeakIndex);
                smoothY.erase(smoothY.begin()+removePeakIndex);
                snrs.erase(snrs.begin()+removePeakIndex);
            }
        }


        // If required, trim number of peaks to requested size.
        // This is used by Turbocharger.
        // it's important to do this after peak refinement, since replicated peaks may be removed at that point.
        if ( fixedPeaksKeep_ > 0 )
        {
            int sizeSNR = snrs.size();
            if ( sizeSNR > fixedPeaksKeep_ )
            {
                double updatedPercentLinesToFilter = 100.0*( 1.0 - ( (double)fixedPeaksKeep_ / (double)snrs.size() ) );
                std::vector<double> sortedSnrs = snrs;
                sort( sortedSnrs.begin(), sortedSnrs.end() );
                double cutoff = scoreAtPercentile(updatedPercentLinesToFilter,sortedSnrs,sortedSnrs.size());

                for (std::vector<double>::iterator it = snrs.begin(); it != snrs.end(); ) {
                    if ( *it < cutoff ) {
                        int index1 = it - snrs.begin();
                        smoothX.erase(smoothX.begin()+index1);
                        smoothY.erase(smoothY.begin()+index1);
                        snrs.erase(it);
                    }
                    else {
                        ++it;
                    }
                }
            }
        }
    }

}//namespace
void Centroidotron::proteoWizDetect(
        const std::vector<double> &x_,
        const std::vector<double> &y_,
        std::vector<double> &xPeakValues,
        std::vector<double> &yPeakValues,
        std::vector<Peak> *peaks
        ) {

    if (x_.size() != y_.size()){
        throw std::runtime_error("[CwtPeakDetector::detect()] x and y arrays must be the same size");
    }

    //cout << "peakPicking on " << x.size() << " points" << endl;

    if (x_.size() <= 2) {
        return;
    }

    // local array copies
    std::vector<double> x(x_);
    std::vector<double> y(y_);

//    // ensure data is m/z sorted
    pwiz::util::sort_together(x, y);

    std::vector<double> binnedX;
    binnedX.reserve(x.size());
    std::vector<double> binnedY;
    binnedY.reserve(y.size());

    // bin identical m/z values
    if (x.size() > 1)
    {
        binnedX.push_back(x[0]);
        binnedY.push_back(y[0]);
        for (size_t i = 1; i < x.size(); ++i){
            if (fabs(binnedX.back() - x[i]) < 1e-6) {
                for (; i < x.size() && fabs(binnedX.back() - x[i]) < 1e-6; ++i){
                    binnedY.back() += y[i];
                }
                --i;
            }
            else {
                binnedX.push_back(x[i]);
                binnedY.push_back(y[i]);
            }
        }

        swap(x, binnedX);
        swap(y, binnedY);
    }

    int mzLength = x.size(); // number of data points in spectrum
    if (mzLength <= 2) return;
    int corrMatrixLength = 2 * mzLength - 1; // number of data points in a row of the correlation matrix

    // Data arrays
    std::vector<std::vector<double>> corrMatrix(nScales, std::vector<double>(corrMatrixLength, 0.0)); // correlation matrix
    std::vector<std::vector<bool>> locMaxs(nScales, std::vector<bool>(corrMatrixLength, false)); // local maxima matrix
    std::vector <double> widths(mzLength,0.0);
    std::vector<std::vector<std::vector<int>>> waveletPoints(2, std::vector<std::vector<int>>(nScales, std::vector<int>(x_.size(), 0)));
    getScales(x, y, waveletPoints, widths);

//    calcCorrelation(x, y, waveletPoints, widths, corrMatrix); // calculate the correlation matrix

//    // step 1: find maxima in each column
//    // step 2: apply sliding window with fixed width to generate list of (row,col) maxima (i.e., "lines")
//    // step 3: filter the list of maxima with SNR
//    std::vector<ridgeLine> allLines; // identified/filtered lines
//    std::vector <double> snrs; // snr associated with each line
//    getPeakLines(corrMatrix, x, allLines, snrs);
//
//    // refine the peak positions and remove peaks using fixedPeaksKeep_, if needed.
//    xPeakValues.reserve(allLines.size()), yPeakValues.reserve(allLines.size());
//    refinePeaks( x, y, allLines, widths, xPeakValues, yPeakValues, snrs);

}

namespace {

    Err buildFirstDerivativeMzVec(
            const ScanPoints &scanPoints,
            QVector<double> *firstDirMzVec
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(scanPoints); ree;

        for (int i = 0; i < scanPoints.size() - 1; i++) {
            firstDirMzVec->push_back(scanPoints.at(i+1).x() - scanPoints.at(i).x());
        }

        ERR_RETURN
    }

    Err addZeroPoints(
            const ScanPoints &scanPoints,
            ScanPoints *scanPointsWithZerosOnPeakEdges
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scanPoints); ree;

        QVector<double> mzValsFirstDeriv;
        e = buildFirstDerivativeMzVec(
                scanPoints,
                &mzValsFirstDeriv
                ); ree;

        const double medianDiffFirst20 = MathUtils::median(mzValsFirstDeriv.mid(0, 20));

        const double rejectFactor = 1.5;
        double currentDelta = medianDiffFirst20;

        for (int i = 0; i < mzValsFirstDeriv.size(); i++) {

            if (mzValsFirstDeriv.at(i) < currentDelta * rejectFactor) {
                scanPointsWithZerosOnPeakEdges->push_back(scanPoints.at(i));
                currentDelta = mzValsFirstDeriv.at(i);
                continue;
            }

            scanPointsWithZerosOnPeakEdges->push_back(scanPoints.at(i));

            scanPointsWithZerosOnPeakEdges->push_back({static_cast<float>(scanPoints.at(i).x() + currentDelta), 0.0f});

            if (i + 1 < scanPoints.size()) {
                scanPointsWithZerosOnPeakEdges->push_back({static_cast<float>(scanPoints.at(i+1).x() - currentDelta), 0.0f});
            }

        }

        scanPointsWithZerosOnPeakEdges->push_back({static_cast<float>(scanPoints.back().x() + currentDelta), 0.0f});

        ERR_RETURN
    }

    double mzWeightedAverageFromScanPoints(const ScanPoints &scanPoints) {

        double runningSum = 0.0;
        double numerator = 0.0;
        for (const ScanPoint &sp : scanPoints) {
            runningSum += sp.y();
            numerator += sp.x() * sp.y();
        }

        return numerator / runningSum;
    }

}//namespace
Err Centroidotron::centroidScan(
        const ScanPoints &_scanPoints,
        ScanPoints *centroidedScanPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(_scanPoints); ree;
    centroidedScanPoints->clear();

    ScanPoints scanPoints;
    e = addZeroPoints(_scanPoints, &scanPoints); ree;

    QPair<QVector<double>, QVector<double>> mzVsIntensityVecs = ParallelUtils::unZip(scanPoints);
    const QVector<double> &mzVals = mzVsIntensityVecs.first;
    const QVector<double> &intensityVals = mzVsIntensityVecs.second;
    const QVector<float> intensityValsFloat(intensityVals.begin(), intensityVals.end());

    PeakIntegratomaticParameters params;
    params.filterLength = m_filterLength;
    params.stopThresholdFraction = 0.5; //TODO make this settable
    params.signalToNoiseRatio = 1.0; //TODO make this settable
    params.smoothCount = 0; //NOTE: set to zero because we're smoothing above.
    params.sigma = 1.0;

    PeakIntegratomatic peakIntegratomatic;
    e = peakIntegratomatic.init(params);

    QVector<QPair<PeakIntegrationIndexes, float>> peakIntegrationIndexesVsIntensity;
    const int intMaxVal = std::numeric_limits<int>::max();
    e = peakIntegratomatic.simpleIntegrator(
            intensityValsFloat,
            intMaxVal,
            intMaxVal,
            &peakIntegrationIndexesVsIntensity
            ); ree;

    ScanPoints scanPointsLimits;
    ScanPoints apexes;
    for (const QPair<PeakIntegrationIndexes, float> &res : peakIntegrationIndexesVsIntensity) {
        const PeakIntegrationIndexes &pii = res.first;
        scanPointsLimits.push_back({static_cast<float>(mzVals.at(pii.first)), static_cast<float>(mzVals.at(pii.second))});

        const ScanPoints scanPointsLimitss = scanPoints.mid(pii.first, (pii.second - pii.first) + 1);
        const auto centroid = static_cast<float>(mzWeightedAverageFromScanPoints(scanPointsLimitss));
        apexes.push_back({centroid, centroid});

    }

    ERR_RETURN

}

namespace {

    double mexicanHatWavelet(double u) {
        double coeff = 2 / std::sqrt(3 * std::sqrt(std::acos(-1)));
        double u2 = u * u;
        return coeff * (1 - u2) * std::exp(-u2 / 2);
    }

}//namespace
Err Centroidotron::performCWT(
        const ScanPoints &scanPoints,
        int minScale,
        int maxScale,
        ScanPoints *processedScanPoints) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    e = ErrorUtils::isTrue(minScale < maxScale); ree;
    e = ErrorUtils::isAboveThreshold(
            minScale,
            0,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;

    processedScanPoints->clear();

    QPair<QVector<double>, QVector<double>> mzVsIntensity = ParallelUtils::unZip(scanPoints);
    const QVector<double> &mzVals = mzVsIntensity.first;
    const QVector<double> &intensityVals = mzVsIntensity.second;
    Eigen::VectorX<double> intensityValsVec = EigenUtils::convertQVectorToEigenVector(intensityVals);

    int length = intensityVals.size();
    Eigen::MatrixX<double> cwtResult(maxScale - minScale + 1, length);

    for (int scale = minScale; scale <= maxScale; ++scale) {
        const double delta = scale / (double)length;
        Eigen::VectorX<double> wavelet(length);

        for (int i = 0; i < length; ++i) {
            double u = (i - length / 2.0) * delta;
            wavelet(i) = mexicanHatWavelet(u) / std::sqrt(delta);
        }

        Eigen::VectorX<double> convolved(length);
        convolved.setZero();

        for (int i = 0; i < length; ++i) {
            for (int j = 0; j < length; ++j)
            {
                int k = (i - j + length) % length;
                convolved(i) += wavelet(k) * intensityValsVec.coeff(j);
            }
        }

        cwtResult.row(scale - minScale) = convolved;
    }

    for (int row = 0; row < cwtResult.rows(); row++) {
        const Eigen::VectorX<double> v = cwtResult.row(row);
        const QVector<double> qV = EigenUtils::convertEigenVectorToQVector(v);

        QVector<QPointF> rowScanPoints;
        e = ParallelUtils::zip(
                mzVals,
                qV,
                &rowScanPoints
                ); ree;

        const QString ogScanFileName = QString("/home/anichols/Downloads/ogScan_%1.csv").arg(row);
        e = MsUtils::writePointsToCSV(rowScanPoints, ogScanFileName); ree;

    }


    ERR_RETURN
}














