//
// Created by anichols on 4/21/23.
//

#ifndef PYTHIADIACPP_TANDEMSPECTRADECONVOLVOTRON_H
#define PYTHIADIACPP_TANDEMSPECTRADECONVOLVOTRON_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"
#include "PythiaParameterReader.h"

using namespace Error;

struct TandemDeconvolverResult {
    DiscScore discScore = -1.0;
    double tTestVal = -1.0;
    double pVal = -1.0;
};


class ALGORITHMSLIB_EXPORTS TandemSpectraDeconvolvotron {

public:

    TandemSpectraDeconvolvotron();
    ~TandemSpectraDeconvolvotron() = default;

    Err init(
            int precision,
            double mzMax,
            int iterationsMax,
            double stopTolerance
            );

    Err deconvolveTandemSpectra(
            const ScanPoints &scanPoints,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPredictions,
            QMap<PeptideStringWithMods, TandemDeconvolverResult> *pepSeqVsWeight,
            double *fStat,
            double *pValFTest
            ) const;

private:



private:

    int m_iterationsMax;
    int m_precision;
    double m_mzMax;
    double m_stopTolerance;

    bool m_isInit;

};

#endif //PYTHIADIACPP_TANDEMSPECTRADECONVOLVOTRON_H
