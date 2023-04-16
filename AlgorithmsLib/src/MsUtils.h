//
// Created by Drucifer on 4/27/2022.
//

#ifndef PYTHIACPP_MSUTILS_H
#define PYTHIACPP_MSUTILS_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
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
            const QVector<QPointF> &vecToExtract,
            const QVector<QPointF> &extractionPoints,
            double extractionPPM
    );

    static Err buildDeletionPoints(
            const ExtractPoints &ep,
            double mzHiPassCutoff,
            QVector<QPointF> *delPoints
            );


};


#endif //PYTHIACPP_MSUTILS_H
