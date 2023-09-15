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
    double frameFStat = -1.0;
    double pValFrameFtest = -1.0;
    double frameError = -1.0;
    int scanNumberCandidateCount = -1;
};


class ALGORITHMSLIB_EXPORTS TandemSpectraDeconvolvotron {

public:

    TandemSpectraDeconvolvotron();
    ~TandemSpectraDeconvolvotron() = default;

    Err init(
            int precision,
            double mzMax,
            double ppmTol,
            int iterationsMax,
            double stopTolerance,
            double pValThreshold
            );

    Err deconvolveTandemSpectra(
            const ScanPoints &scanPoints,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPredictions,
            QMap<PeptideStringWithMods, TandemDeconvolverResult> *pepSeqVsWeight
            );

private:

    int m_iterationsMax;
    int m_precision;
    double m_mzMax;
    double m_stopTolerance;
    double m_pValThreshold;
    double m_ppmTol;

    bool m_isInit;

};

#endif //PYTHIADIACPP_TANDEMSPECTRADECONVOLVOTRON_H
