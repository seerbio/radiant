//
// Created by anichols on 6/5/24.
//

#ifndef DECONVOLVOTRON_H
#define DECONVOLVOTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS Deconvolvotron {

public:

    Deconvolvotron();
    ~Deconvolvotron() = default;

    Err init(int precision);

    Err deconvolve(
        const QMap<IdStr, QVector<QPointF>> &aMatrixPoints,
        const QVector<QPointF> &bVecPoints,
        QVector<QPair<IdStr, Score>> *idStrVsScore
        );

private:

    Err buildXValsSet(
        const QMap<IdStr, QVector<QPointF>>& aMatrixPoints,
        const QVector<QPointF>& bVecPoints,
        QVector<int> *xVals
        ) const;

private:

    int m_precision;

};



#endif //DECONVOLVOTRON_H
