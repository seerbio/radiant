//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;

class TargetDecoyPairParallelInput;


class ALGORITHMSLIB_EXPORTS TargetDecoyCandidatePairScoretron {

public:

    TargetDecoyCandidatePairScoretron() = default;
    ~TargetDecoyCandidatePairScoretron() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QVector<MsScanInfo> &msScanInfos,
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
            TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
            );

    Err scoreTargetDecoyPairs();

private:

    Err buildParallelInput(
            double randomSelectionFraction,
            QVector<TargetDecoyPairParallelInput> *input
            );



private:

    PythiaParameters m_pythiaParameters;
    QVector<MsScanInfo> m_msScanInfos;
    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *m_diaTargetFrames;
    TargetDecoyCandidatePairManager *m_targetDecoyCandidatePairManager;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRSCORETRON_H
