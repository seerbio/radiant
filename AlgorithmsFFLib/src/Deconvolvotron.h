//
// Created by anichols on 6/5/24.
//

#ifndef DECONVOLVOTRON_H
#define DECONVOLVOTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"

using namespace Error;

struct DeconvolvotronResult {
    DiscScore discScore = -1.0;
    double pVal = -1.0;
    double pValFrameFtest = -1.0;
};

class ALGORITHMSFFLIB_EXPORTS Deconvolvotron {

public:

    Deconvolvotron();
    ~Deconvolvotron() = default;

    Err init(int precision);

    Err deconvolve(
        const QVector<QPair<IdStr, QVector<QPointF>>> &aMatrixPoints,
        const QVector<QPointF> &bVecPoints,
        QVector<QPair<IdStr, DeconvolvotronResult>> *idStrVsScore
        ) const;

private:

    Err buildXValsSet(
        const QVector<QPair<IdStr, QVector<QPointF>>> &aMatrixPoints,
        const QVector<QPointF>& bVecPoints,
        QVector<double> *xValsReturn
        ) const;

private:

    int m_precision;

};

#endif //DECONVOLVOTRON_H
