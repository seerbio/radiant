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
