//
// Created by Drucifer on 4/27/2022.
//

#include "MsUtils.h"

#include "ErrorUtils.h"
#include "GlobalSettings.h"
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

    return extractPointsOutput;
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

    QPointF getMaxPoint(const QVector<QPointF> &scanPoints) {

        const auto maxPointLogic = [](const QPointF &l, const QPointF &r){
            return l.y() < r.y();
        };

        const QPointF maxPoint = *std::max_element(
                scanPoints.begin(),
                scanPoints.end(),
                maxPointLogic
                );

        return maxPoint;
    }

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

    Err chargeScorer(
            const QPointF &mzCenterPoint,
            const QVector<QPointF> &scanPoints,
            int charge,
            double ppmTol,
            double *intensitySum
            ) {

        ERR_INIT

        *intensitySum = 0.0;

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

        const int numberOfIonsToGenerateRight = 10;

        const QVector<QPointF> leftPoints = calculateMzPointsFromMzStart(
                mzCenterPoint.x(),
                numberOfIonsToGenerateLeft,
                charge,
                true
                );

        const QVector<QPointF> rightPoints = calculateMzPointsFromMzStart(
                mzCenterPoint.x(),
                numberOfIonsToGenerateRight,
                charge,
                false
        );

        const ExtractPoints extractedPointsLeft
                = MsUtils::extractPointsFromPoints(scanPoints, leftPoints, ppmTol);

        const ExtractPoints extractedPointsRight
                = MsUtils::extractPointsFromPoints(scanPoints, rightPoints, ppmTol);

        if (extractedPointsRight.mzFoundVsSearched.isEmpty()) {
            *intensitySum = -1.0;
            ERR_RETURN
        }


        double currentIntensity = mzCenterPoint.y();
        for (int i = 0; i < extractedPointsRight.mzFoundVsSearched.size(); i++) {

            const QPointF pointMz = extractedPointsRight.mzFoundVsSearched.at(i);
            const QPointF pointsIntensity = extractedPointsRight.intensityFoundVsSearched.at(i);

            const double mzFound = pointMz.x();
            const double intensityFound = pointsIntensity.x();

            if (mzFound < 0 || intensityFound > currentIntensity) {
                break;
            }

            *intensitySum += intensityFound;
            currentIntensity = intensityFound;
        }

        currentIntensity = mzCenterPoint.y();
        for (int i = 0; i < extractedPointsLeft.mzFoundVsSearched.size(); i++) {

            const QPointF pointMz = extractedPointsLeft.mzFoundVsSearched.at(i);
            const QPointF pointsIntensity = extractedPointsLeft.intensityFoundVsSearched.at(i);

            const double mzFound = pointMz.x();
            const double intensityFound = pointsIntensity.x();

            if (mzFound < 0 || intensityFound > currentIntensity) {
                break;
            }

            *intensitySum += intensityFound;
            currentIntensity = intensityFound;
        }

        ERR_RETURN
    }


}
Err MsUtils::chargeDeterminator(
        const QPointF &mzCenterPoint,
        const QVector<QPointF> &scanPoints,
        double ppmTol,
        int *charge
        ) {

    ERR_INIT

    *charge = -1;

    double currentBestScore = 0.0;

    const int chargeMin = 1;
    const int chargeMax = 5;

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

