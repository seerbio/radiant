//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H

#include "AlgorithmsFFLib_Exports.h"

#include "CandidateScores.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;


class MsCalibratomatic;
class TargetDecoyPairParallelInput;


class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePairScoretron {

public:

    TargetDecoyCandidatePairScoretron();
    ~TargetDecoyCandidatePairScoretron() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanTimeMS1,
            MsReaderPointerAcc *msReaderPointerAcc,
            QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames
            );

    Err scoreTargetDecoyPairs(
            int topNMS2Ions,
            const QPair<double, double> &scanTimeMinMax,
            const MsCalibratomatic &msCalibratomatic,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
            QVector<CandidateScores> *candidateScoresVec
            );

    bool isInit();

    Err setPythiaParameters(const PythiaParameters &pythiaParameters);

private:

    Err buildParallelInput(
            int topNMS2Ions,
            const QPair<double, double> &scanTimeMinMax,
            const MsCalibratomatic &msCalibratomatic,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
            QVector<TargetDecoyPairParallelInput> *input
            );


private:

    PythiaParameters m_pythiaParameters;
    MsReaderPointerAcc *m_msReaderPointerAcc;
    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *m_diaTargetFrames;
    QMap<ScanNumber, ScanPoints> m_ms1Frame;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
