//
// Created by Drucifer on 4/27/2022.
//

#ifndef PYTHIACPP_MSUTILS_H
#define PYTHIACPP_MSUTILS_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"

#include <QPointF>

using namespace Error;

struct ALGORITHMSLIB_EXPORTS ExtractPoints {
    QVector<QPointF> mzFoundVsSearched;
    QVector<QPointF> intensityFoundVsSearched;
};


namespace MsUtilsNamespace{

    const extern std::function<void(void)> ALGORITHMSLIB_EXPORTS sortAscMz;
    const extern std::function<void(void)> ALGORITHMSLIB_EXPORTS sortAscIntensity;
}


class ALGORITHMSLIB_EXPORTS MsUtils {


public:

    static ExtractPoints extractPointsFromPoints(
            const QVector<QPointF> &points,
            const QVector<QPointF> &extractionPoints,
            double extractionPPM
    );

    static ScanPoints extractPointsFromPoints(
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

    static Err calcXCorr(
            const QVector<double> &foundMassVals,
            const QVector<double> &theoMassVals,
            const QVector<double> &foundIntensityVals,
            const QVector<double> &theoIntensityVals,
            double ppmTol,
            double *xCorr
    );


};


#endif //PYTHIACPP_MSUTILS_H
