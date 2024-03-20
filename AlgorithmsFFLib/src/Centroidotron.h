//
// Created by anichols on 3/12/24.
//

#ifndef PYTHIADIACPP_CENTROIDOTRON_H
#define PYTHIADIACPP_CENTROIDOTRON_H

#include "AlgorithmsFFLib_Exports.h"

#include "GlobalSettings.h"
#include "Error.h"

#include <complex>

using namespace Error;

struct Peak
        {
        double x; /// x value of a signal peak (or centroid)
        double y; /// y value of a signal peak (or centroid), aka intensity/abundance/amplitude
        double start; // x value where the peak's profile starts
        double stop; // x value where the peak's profile stops
        double area; // area under the profile between start and stop
        };

typedef struct {
    int Col;
    int Row;
} ridgeLine;


class ALGORITHMSFFLIB_EXPORTS Centroidotron {

    friend class CentroidotronTests;

public:

    Centroidotron() = default;
    ~Centroidotron() = default;

    Err init(
            double peakWidth,
            int hashingPrecision,
            int filterLength
            );

    Err centroidScan(
            const ScanPoints &scanPoints,
            ScanPoints *centroidedScanPoints
            );

private:

    ScanPoints runningAverage(const ScanPoints &scanPoints);

    Err smoothIntensities(
            const QVector<double> &intzVals,
            QVector<double> *intzValsSmoothed
            );

    Err smoothIntensities(
            const ScanPoints &scanPoints,
            ScanPoints *scanPointsSmoothed
            );

    Err performCWT(
            const ScanPoints &scanPoints,
            int minScale,
            int maxScale,
            ScanPoints *processedScanPoints
            );

    static void proteoWizDetect(
            const std::vector<double>& x_,
            const std::vector<double>& y_,
            std::vector<double>& xPeakValues,
            std::vector<double>& yPeakValues,
            std::vector<Peak> *peaks
            );

private:

    double m_peakWidth;
    int m_hashingPrecision;
    int m_filterLength;


};


#endif //PYTHIADIACPP_CENTROIDOTRON_H
