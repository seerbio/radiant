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

    static Err monoIsotopeDeterminator(
            const QPointF &mzCenterPoint,
            const QVector<QPointF> &scanPoints,
            double ppmTol,
            int charge,
            int *monoIsoOffset,
            QVector<QPointF> *subtractionPoints,
            double *bestCosineSim
    );



};


#endif //PYTHIACPP_MSUTILS_H
