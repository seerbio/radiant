//
// Created by Drucifer on 4/27/2022.
//

#ifndef PYTHIACPP_MSUTILS_H
#define PYTHIACPP_MSUTILS_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"

#include <QPointF>

using namespace Error;

struct ALGORITHMSFFLIB_EXPORTS ExtractPoints {
    QVector<QPointF> mzFoundVsSearched;
    QVector<QPointF> intensityFoundVsSearched;
};


class ALGORITHMSFFLIB_EXPORTS MsUtils {


public:

    static ExtractPoints extractPointsFromPoints(
            const QVector<QPointF> &points,
            const QVector<QPointF> &extractionPoints,
            double extractionPPM
    );

    static QVector<QPointF> extractPointsFromPoints(
            const QVector<QPointF> &points,
            const QVector<double> &extractionPoints,
            double extractionPPM,
            bool removeZeroPoints = true
    );

    static Err buildDeletionPoints(
            const ExtractPoints &ep,
            double mzHiPassCutoff,
            QVector<QPointF> *delPoints
            );

    static Err chargeDeterminator(
            const QPointF &mzCenterPoint,
            const QVector<QPointF> &scanPoints,
            double ppmTol,
            int chargeMin,
            int chargeMax,
            int *charge
            );

    static int getCenterPointIndex(
            const QVector<QPointF> &points,
            const QPointF &mzCenterPoint
    );

    static Err monoIsotopeDeterminator(
            const QPointF &mzCenterPoint,
            const QVector<QPointF> &scanPoints,
            double ppmTol,
            int charge,
            int *monoIsoOffset,
            QVector<QPointF> *subtractionPoints,
            double *bestCosineSim
    );

    static Err writePointsToCSV(
            const QVector<QPointF> &points,
            const QString &destFilePath
    );

};


#endif //PYTHIACPP_MSUTILS_H
