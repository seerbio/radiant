//
// Created by Drucifer on 4/27/2022.
//

#include "MsUtils.h"

#include "BiophysicalCalcs.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "IsotopicDistributionBuilder.h"
#include "ParallelUtils.h"

namespace {
    const auto sortAscMz = [](const QPointF &l, const QPointF &r){return l.x() < r.x();};
}
ExtractPoints MsUtils::extractPointsFromPoints(
        const QVector<QPointF> &_points,
        const QVector<QPointF> &_pointsToExtract,
        double extractionPPM
        ) {

    ExtractPoints extractPointsOutput;
    extractPointsOutput.mzFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});
    extractPointsOutput.intensityFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});

    QVector<QPointF> pointsToExtract = _pointsToExtract;
    std::sort(pointsToExtract.begin(), pointsToExtract.end(), sortAscMz);

    QVector<QPointF> points = _points;
    std::sort(points.begin(), points.end(), sortAscMz);

    int currentExtractionIndex = 0;
    double extractionPointX = pointsToExtract.at(currentExtractionIndex).x();
    double extractionPointY = pointsToExtract.at(currentExtractionIndex).y();

    extractPointsOutput.mzFoundVsSearched[currentExtractionIndex].ry() = extractionPointX;
    extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].ry() = extractionPointY;

    double xPPM = MathUtils::calculatePPM(extractionPointX, extractionPPM);
    double xLo = extractionPointX - xPPM;
    double xHi = extractionPointX + xPPM;

    for (const QPointF &xPoint : points) {

        const double xVal = xPoint.x();

        while (xVal > xHi) {

            currentExtractionIndex++;

            if (currentExtractionIndex > _pointsToExtract.size() - 1) {
                break;
            }

            extractionPointX = pointsToExtract.at(currentExtractionIndex).x();
            extractionPointY = pointsToExtract.at(currentExtractionIndex).y();

            xPPM = MathUtils::calculatePPM(extractionPointX, extractionPPM);
            xLo = extractionPointX - xPPM;
            xHi = extractionPointX + xPPM;

            extractPointsOutput.mzFoundVsSearched[currentExtractionIndex].ry() = extractionPointX;
            extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].ry() = extractionPointY;
        }

        if (xLo <= xVal && xVal <= xHi) {

            extractPointsOutput.mzFoundVsSearched[currentExtractionIndex]
                    =  xPoint.y() > extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].x()
                    ? QPointF(xVal , extractionPointX)
                    : extractPointsOutput.mzFoundVsSearched[currentExtractionIndex];

            extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].rx() =  std::max(
                    xPoint.y(),
                    extractPointsOutput.intensityFoundVsSearched[currentExtractionIndex].x()
                    );
        }
    }

    while (extractPointsOutput.intensityFoundVsSearched.back().y() < 0) {
        extractPointsOutput.mzFoundVsSearched.pop_back();
        extractPointsOutput.intensityFoundVsSearched.pop_back();
    }

    return extractPointsOutput;
}

ScanPoints MsUtils::extractPointsFromPoints(
        const QVector <QPointF> &points,
        const QVector<double> &extractionPoints,
        double extractionPPM,
        bool removeZeroPoints
) {

    QVector<QPointF> extractQPoints;
    std::transform(
            extractionPoints.begin(),
            extractionPoints.end(),
            std::back_inserter(extractQPoints),
            [](double mz){return QPointF(mz, 1.0);}
    );

    const ExtractPoints ep = extractPointsFromPoints(
            points,
            extractQPoints,
            extractionPPM
            );

    QVector<MS2Ion> extractedPoints;
    for (int i = 0; i < ep.mzFoundVsSearched.size(); i++) {
        const double mzFound = ep.mzFoundVsSearched.at(i).x();
        const double intensityFound = ep.intensityFoundVsSearched.at(i).x();

        if (mzFound < 0 && removeZeroPoints) {
            continue;
        }

        extractedPoints.push_back({mzFound, intensityFound});
    }

    return extractedPoints;
}

namespace {


     Err normalizeTheoPointsToActual(
            const ExtractPoints &extractPoints,
            double mzHiPassCutoff,
            QVector<QPointF> *normalizedPoints
            ) {


         ERR_INIT

         const QVector<QPointF> &mzVals = extractPoints.mzFoundVsSearched;
         const QVector<QPointF> &intzVals = extractPoints.intensityFoundVsSearched;

         const QPair<QVector<double>, QVector<double>> splitPointsMz = ParallelUtils::unZip(mzVals);
         const QPair<QVector<double>, QVector<double>> splitPointsIntz = ParallelUtils::unZip(intzVals);

         QVector<QPointF> founds;
         e = ParallelUtils::zip(
                 splitPointsMz.first,
                 splitPointsIntz.first,
                 &founds
         ); ree;

         QVector<QPointF> theos;
         e = ParallelUtils::zip(
                 splitPointsMz.second,
                 splitPointsIntz.second,
                 &theos
         ); ree;

         e = ErrorUtils::isEqual(
                 founds.size(),
                 theos.size()
         ); ree;

         QVector<double> sumFounds;
         QVector<double> sumTheos;
         for (int i = 0; i < founds.size(); ++i) {

             const QPointF &found = founds.at(i);
             const QPointF &theo = theos.at(i);

             if (found.x() > mzHiPassCutoff) {
                 sumFounds.push_back(found.y());
             }

             if (theo.x() > mzHiPassCutoff) {
                 sumTheos.push_back(theo.y());
             }
         }

         const double meanFounds = MathUtils::mean(sumFounds);
         const double meanTheos = MathUtils::mean(sumTheos);

         const double scaleFactor = MathUtils::tZero(meanTheos) ? 0.0 : meanFounds / meanTheos;

         QVector<QPointF> normalizedFoundPoints;
         for (int i = 0; i < splitPointsMz.first.size(); i++) {

             const double mz = splitPointsMz.first.at(i);
             const double intz = splitPointsIntz.second.at(i) * scaleFactor;
             normalizedFoundPoints.push_back({mz, intz});
         }

         *normalizedPoints = normalizedFoundPoints;

         ERR_RETURN
    }


}//namespace
Err MsUtils::buildDeletionPoints(
        const ExtractPoints &ep,
        double mzHiPassCutoff,
        QVector<QPointF> *delPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ep.mzFoundVsSearched); ree;
    e = ErrorUtils::isEqual(
            ep.intensityFoundVsSearched.size(),
            ep.mzFoundVsSearched.size()
            ); ree;

    QVector<QPointF> founds;
    e = normalizeTheoPointsToActual(
            ep,
            mzHiPassCutoff,
            &founds
            );

    ERR_RETURN
}

namespace {

    QVector<QPointF> calculateMzPointsFromMzStart(
        double mzStart,
        int numberOfIonsToGenerate,
        int charge,
        bool generateInLeftDirection
        ) {

        const int direction = generateInLeftDirection ? -1 : 1;
        const double chargeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / charge;

        QVector<QPointF> points;
        for (int i = 1; i <= numberOfIonsToGenerate; i++) {

            const double mz = mzStart + (i * chargeDistance * direction);
            points.push_back({mz, 1});
        }

        return points;
    }

    int calcNumberOfIonsToGenerateLeft(int charge) {

        int numberOfIonsToGenerateLeft = 3;

        switch (charge) {

            case 1:
                numberOfIonsToGenerateLeft = 1;
                break;
            case 2:
                numberOfIonsToGenerateLeft = 2;
                break;
            default:
                numberOfIonsToGenerateLeft = 3;
                break;
        }

        return numberOfIonsToGenerateLeft;
    }

    void extractPoints(
            const QVector<QPointF> &scanPoints,
            int charge,
            double mzCenterPointVal,
            double ppmTol,
            ExtractPoints *extractedPointsLeft,
            ExtractPoints *extractedPointsRight
            ) {

        const int numberOfIonsToGenerateLeft = calcNumberOfIonsToGenerateLeft(charge);

        const int numberOfIonsToGenerateRight = 10;

        const QVector<QPointF> leftPoints = calculateMzPointsFromMzStart(
                mzCenterPointVal,
                numberOfIonsToGenerateLeft,
                charge,
                true
        );

        const QVector<QPointF> rightPoints = calculateMzPointsFromMzStart(
                mzCenterPointVal,
                numberOfIonsToGenerateRight,
                charge,
                false
        );

        *extractedPointsLeft
                = MsUtils::extractPointsFromPoints(scanPoints, leftPoints, ppmTol);

        *extractedPointsRight
                = MsUtils::extractPointsFromPoints(scanPoints, rightPoints, ppmTol);
    }

    double iterateFoundPoints(
            const QPointF &mzCenterPoint,
            const ExtractPoints &extractedPoints
            ) {

        double intensitySum = 0.0;

        double currentIntensity = mzCenterPoint.y();
        for (int i = 0; i < extractedPoints.mzFoundVsSearched.size(); i++) {

            const QPointF pointMz = extractedPoints.mzFoundVsSearched.at(i);
            const QPointF pointsIntensity = extractedPoints.intensityFoundVsSearched.at(i);

            const double mzFound = pointMz.x();
            const double intensityFound = pointsIntensity.x();

            if (mzFound < 0 || intensityFound > currentIntensity) {
                break;
            }

            intensitySum += intensityFound;
            currentIntensity = intensityFound;
        }

        return intensitySum;
    }

    Err chargeScorer(
            const QPointF &mzCenterPoint,
            const QVector<QPointF> &scanPoints,
            int charge,
            double ppmTol,
            double *intensitySum
            ) {

        ERR_INIT

        *intensitySum = 0.0;

        ExtractPoints extractedPointsLeft;
        ExtractPoints extractedPointsRight;
        extractPoints(
                scanPoints,
                charge,
                mzCenterPoint.x(),
                ppmTol,
                &extractedPointsLeft,
                &extractedPointsRight
                );

        if (extractedPointsRight.mzFoundVsSearched.isEmpty()) {
            *intensitySum = -1.0;
            ERR_RETURN
        }

        *intensitySum += iterateFoundPoints(
                mzCenterPoint,
                extractedPointsLeft
                );

        *intensitySum += iterateFoundPoints(
                mzCenterPoint,
                extractedPointsRight
        );


        ERR_RETURN
    }


}//namespace
Err MsUtils::chargeDeterminator(
        const QPointF &mzCenterPoint,
        const QVector<QPointF> &scanPoints,
        double ppmTol,
        int chargeMin,
        int chargeMax,
        int *charge
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;

    *charge = -1;

    double currentBestScore = 0.0;

    for (int chrg = chargeMin; chrg <= chargeMax; chrg++) {

        double intensitySum;

        e = chargeScorer(
                mzCenterPoint,
                scanPoints,
                chrg,
                ppmTol,
                &intensitySum
                ); ree;

        if (intensitySum > currentBestScore) {
            currentBestScore = intensitySum;
            *charge = chrg;
        }

    }

    ERR_RETURN
}

namespace {

    int getCenterPointIndex(
            const QVector<QPointF> &points,
            const QPointF &mzCenterPoint
            ) {

        QVector<double> mzVals;
        std::transform(
                points.begin(),
                points.end(),
                std::back_inserter(mzVals),
                [](const QPointF &p){return p.x();}
        );

        return MathUtils::closest(mzVals, mzCenterPoint.x());
    }

    Err buildExtractedPoints(
            const QPointF &mzCenterPoint,
            const QVector<QPointF> &scanPoints,
            double ppmTol,
            int charge,
            int *startCenterPointIdxOG,
            QVector<QPointF> *extractedPointsFilled
            ) {

        ERR_INIT

        const int leftPointCount = 5;
        const int rightPointCount = 11;

        const QVector<double> mzTheoPoints
                = BiophysicalCalcs::calculateIsotopesFromMz(mzCenterPoint.x(), charge, leftPointCount, rightPointCount);

        const bool removeZeroPoints = false;
        QVector<QPointF> extractedScanPoints = MsUtils::extractPointsFromPoints(
                scanPoints,
                mzTheoPoints,
                ppmTol,
                removeZeroPoints
        );

        QVector<QPointF> extractedPointsFilledTemp(mzTheoPoints.size(), {0,0});
        *extractedPointsFilled = extractedPointsFilledTemp;

        std::copy(extractedScanPoints.begin(), extractedScanPoints.end(), extractedPointsFilled->begin());

        *startCenterPointIdxOG = getCenterPointIndex(*extractedPointsFilled, mzCenterPoint);

        e = ErrorUtils::isTrue(
                MathUtils::tZero(extractedPointsFilled->at(*startCenterPointIdxOG).x() - mzCenterPoint.x())
        ); ree;

        ERR_RETURN
    }

    Eigen::VectorX<double> convertQPointsYToVec(const QVector<QPointF> &vec) {

        QVector<double> yVals;
        yVals.reserve(vec.size());

        std::transform(
                vec.begin(),
                vec.end(),
                std::back_inserter(yVals),
                [](const QPointF &p){return p.y();}
                );

        return EigenUtils::convertQVectorToEigenVector(yVals);
    }

QVector<QPointF> buildSubractionPoints(
        const QVector<QPointF> &extractedPointsFilled,
        int startCenterPointIdx,
        int monoOffset
        ) {

        QVector<QPointF> subtractPoints;

        const int return4Points = 4;
        for (int i = startCenterPointIdx - monoOffset; i < extractedPointsFilled.size(); i++) {

            const QPointF &p = extractedPointsFilled.at(i);
            if (p.x() < 0) {
                continue;
            }

            subtractPoints.push_back(p);

            if (subtractPoints.size() >= return4Points) {
                break;
            }

        }


        return subtractPoints;
    }

}//namespace
Err MsUtils::monoIsotopeDeterminator(
        const QPointF &mzCenterPoint,
        const QVector<QPointF> &scanPoints,
        double ppmTol,
        int charge,
        int *monoIsoOffset,
        QVector<QPointF> *subtractionPoints,
        double *bestCosineSim
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;

    *bestCosineSim = 0;
    *monoIsoOffset = 1000;

    int startCenterPointIdxOG;
    QVector<QPointF> extractedPointsFilled;
    e = buildExtractedPoints(
            mzCenterPoint,
            scanPoints,
            ppmTol,
            charge,
            &startCenterPointIdxOG,
            &extractedPointsFilled
            ); ree;

    const Eigen::VectorX<double> extractedPointsFilledVec = convertQPointsYToVec(extractedPointsFilled);

    const double roughMass = mzCenterPoint.x() * charge;
    QVector<double> isoDistribution
            = IsotopicDistributionBuilder::getIsotopicDistribution(roughMass);

    const int maxIsotopeToUse = 7;
    isoDistribution.resize(std::min(maxIsotopeToUse, isoDistribution.size()));

    int startCenterPointIdx = startCenterPointIdxOG;

    int cycle = 0;
    while (startCenterPointIdx >= 0) {

        QVector<double> isoDistributionPointsFilled(extractedPointsFilled.size(), 0);
        std::copy(isoDistribution.begin(), isoDistribution.end(), isoDistributionPointsFilled.begin() + startCenterPointIdx);

        const Eigen::VectorX<double> isoDistributionPointsFilledVec
                = EigenUtils::convertQVectorToEigenVector(isoDistributionPointsFilled);

        const double cosineSim = EigenUtils::cosineSimilarity(
                extractedPointsFilledVec,
                isoDistributionPointsFilledVec
                );

        if (cosineSim > *bestCosineSim) {
            *bestCosineSim = cosineSim;
            *monoIsoOffset = cycle;
        }

        cycle++;
        startCenterPointIdx--;
    }

    *subtractionPoints = buildSubractionPoints(
            extractedPointsFilled,
            startCenterPointIdxOG,
            *monoIsoOffset
            );

    ERR_RETURN
}
