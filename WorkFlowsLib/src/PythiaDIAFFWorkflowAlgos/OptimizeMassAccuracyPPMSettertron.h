//
// Created by andrewnichols on 9/28/24.
//

#ifndef OPTIMIZEMASSACCURACYPPMSETTERTRON_H
#define OPTIMIZEMASSACCURACYPPMSETTERTRON_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretron.h"

using namespace Error;

class WORKFLOWSLIB_EXPORTS OptimizeMassAccuracyPPMSettertron {

public:

    OptimizeMassAccuracyPPMSettertron();
    ~OptimizeMassAccuracyPPMSettertron() = default;

    Err initExec(
        MsReaderPointerAcc *msReaderPointerAcc,
        MsCalibratomatic *msCalibratomatic,
        PythiaParameters *pythiaParameters,
        TargetDecoyCandidatePairScoretron2 *targetDecoyCandidatePairScoretron,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPairPntrs
        );

    [[nodiscard]] QVector<float> weights() const;

private:

    int calculateNumberOfTranches() const;

    Err optimizePPM();

private:

    MsReaderPointerAcc *m_msReaderPointerAcc;
    MsCalibratomatic *m_msCalibratomatic;
    PythiaParameters *m_pythiaParameters;
    TargetDecoyCandidatePairScoretron2 *m_targetDecoyCandidatePairScoretron;

    QVector<TargetDecoyCandidatePair*> *m_targetDecoyPairPntrs;
    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> m_candidateScorePairs;

    QVector<float> m_weights;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints>> m_diaTargetFrames;
    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> m_diaTargetFramesPntrs;

};



#endif //OPTIMIZEMASSACCURACYPPMSETTERTRON_H
