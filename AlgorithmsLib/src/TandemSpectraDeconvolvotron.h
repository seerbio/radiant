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

class MS2Ion;

class ALGORITHMSLIB_EXPORTS TandemSpectraDeconvolvotron {

public:

    TandemSpectraDeconvolvotron();
    ~TandemSpectraDeconvolvotron() = default;

    Err setParameters(
            int precision,
            double mzMax,
            int iterationsMax,
            double stopTolerance
            );

    Err deconvolveTandemSpectra(
            const ScanPoints &scanPoints,
            const QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> &tandemPredictions,
            QMap<PeptideSequenceChargeKey, double> *pepSeqVsWeight
            );

private:

    int m_iterationsMax;
    int m_precision;
    double m_mzMax;
    double m_stopTolerance;

};

#endif //PYTHIADIACPP_TANDEMSPECTRADECONVOLVOTRON_H
