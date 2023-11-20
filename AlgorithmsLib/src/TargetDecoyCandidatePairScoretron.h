//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;

class TargetDecoyPairParallelInput;
class MsCalibratomatic;


class ALGORITHMSLIB_EXPORTS TargetDecoyCandidatePairScoretron {

public:

    TargetDecoyCandidatePairScoretron() = default;
    ~TargetDecoyCandidatePairScoretron() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanTimeMS1,
            MsReaderPointerAcc *msReaderPointerAcc,
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
            TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
            );

    Err scoreTargetDecoyPairs(
            int topNMS2Ions,
            const MsCalibratomatic &msCalibratomatic,
            QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers
            );

    bool isInit();

    Err setPythiaParameters(const PythiaParameters &pythiaParameters);

private:

    Err buildParallelInput(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
            int topNMS2Ions,
            const MsCalibratomatic &msCalibratomatic,
            QVector<TargetDecoyPairParallelInput> *input
            );


private:

    PythiaParameters m_pythiaParameters;
    MsReaderPointerAcc *m_msReaderPointerAcc;
    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints*>> *m_diaTargetFrames;
    QMap<ScanNumber, ScanPoints*> m_ms1Frame;
    TargetDecoyCandidatePairManager *m_targetDecoyCandidatePairManager;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
