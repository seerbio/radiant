//
// Created by Drucifer on 4/27/2022.
//

#include "MsUtils.h"

#include "ErrorUtils.h"
#include "IsotopicDistributionBuilder.h"
#include "ParallelUtils.h"


ExtractPoints MsUtils::extractPointsFromPoints(
        const QVector<QPointF> &_points,
        const QVector<QPointF> &_pointsToExtract,
        double extractionPPM
        ) {

    ExtractPoints extractPointsOutput;
    extractPointsOutput.mzFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});
    extractPointsOutput.intensityFoundVsSearched = QVector<QPointF>(_pointsToExtract.size(), {-1.0,-1.0});

    QVector<QPointF> pointsToExtract = _pointsToExtract;
    std::sort(pointsToExtract.begin(), pointsToExtract.end(), [](const QPointF &l, const QPointF &r){return l.x() < r.x();});

    QVector<QPointF> points = _points;
    std::sort(points.begin(), points.end(), [](const QPointF &l, const QPointF &r){return l.x() < r.x();});

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

QVector<QPointF> MsUtils::extractPointsFromPoints(
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

    QVector<QPointF> extractedPoints;
    for (int i = 0; i < ep.mzFoundVsSearched.size(); i++) {
        const double mzSearched = ep.mzFoundVsSearched.at(i).y();
        const double intensityFound = ep.intensityFoundVsSearched.at(i).x();

        if (removeZeroPoints && (intensityFound < 0 || MathUtils::tZero(intensityFound))) {
            continue;
        }

        extractedPoints.push_back({mzSearched, std::max(intensityFound, 0.0)});
    }

    return extractedPoints;
}
