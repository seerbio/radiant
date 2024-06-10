//
// Created by anichols on 6/5/24.
//

#ifndef DECONVOLVOTRON_H
#define DECONVOLVOTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"

class CandidateScores;

using namespace Error;

struct DeconvolvotronResult {
    DiscScore discScore = -1.0;
    double pVal = -1.0;
    double pValFrameFtest = -1.0;
};

class ALGORITHMSFFLIB_EXPORTS Deconvolvotron {

    friend class PlayGroundTests;

public:

    Deconvolvotron();
    ~Deconvolvotron() = default;

    Err init(
        int mzGroupingPrecision,
        double ppmExtractionTolerance
        );

    Err deconvolve(
        const QVector<QPair<CandidateScores*, QVector<QPointF>>> &aMatrixPoints,
        const QVector<QPointF> &bVecPoints,
        QVector<QPair<CandidateScores*, DeconvolvotronResult>> *candScoresVsScore
        ) const;

private:

    Err buildXValsSet(
        const QVector<QPair<CandidateScores*, QVector<QPointF>>> &aMatrixPoints,
        QMap<int, QVector<double>> *hashedXValsVsMzValsGrouped
        ) const;

    static Err buildAMatrixAndBVecTestAccess(
        const QVector<QPair<CandidateScores*, QVector<QPointF>>> &aMatrixPoints,
        const QVector<QPointF>& bVecPoints,
        const QMap<int, QVector<double>> &hashedXValsVsMzValsGrouped,
        double ppmExtractTol,
        QVector<QVector<double>> *aMat,
        QVector<double> *bVec
        );

private:

    int m_mzGroupingPrecision;
    double m_ppmExtractionTolerance;

};

#endif //DECONVOLVOTRON_H
