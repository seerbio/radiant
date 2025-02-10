//
// Created by andrewnichols on 9/27/24.
//

#ifndef MSCALIBRATOMATICSETTERTRON_H
#define MSCALIBRATOMATICSETTERTRON_H

#include "WorkFlowsLib_Exports.h"

#include "CandidateScores.h"
#include "Error.h"
#include "FeatureFinderHillBuilder.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsReaderPointerAcc.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretron.h"


using namespace Error;

class WORKFLOWSLIB_EXPORTS MsCalibratomaticSettertron {

public:

    MsCalibratomaticSettertron();
    ~MsCalibratomaticSettertron() = default;

    Err init(
        PythiaParameters *pythiaParameters,
        MsReaderPointerAcc *msReaderPointerAcc,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager,
        TargetDecoyCandidatePairScoretron2 *targetDecoyCandidatePairScoretron,
        bool excludeDecoys
        );

    Err buildCalibration(MsCalibratomatic *msCalibratomatic);

private:

    int calculateNumberOfTranches() const;

    Err honeIRTAndMassCalibration(
        QVector<CandidateScores*> *candidateScoresVecScoredPntrs,
        int topNCandidates,
        int topCandidatesMass
    );

    Err setMsCalibratomaticMetrics();

    Err recalibrateMzVals(
        const MSLevelEnum &msLevel,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
        ) const;

private:

    MsReaderPointerAcc *m_msReaderPointerAcc;
    TargetDecoyCandidatePairManager *m_targetDecoyCandidatePairManager;
    TargetDecoyCandidatePairScoretron2 *m_targetDecoyCandidatePairScoretron;
    PythiaParameters *m_pythiaParameters;

    MsCalibratomatic m_msCalibratomatic;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints>> m_diaTargetFrames;
    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> m_diaTargetFramesPntrs;

    QVector<TargetDecoyCandidatePair*> m_targetDecoyPairPntrs;
    QVector<TargetDecoyCandidatePair*> m_targetDecoyCandidatePairsTopScores;
    QHash<TargetDecoyCandidatePair*, bool> m_entered;

    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> m_candidateScorePairs;

    QVector<float> m_scanTimeStDevs;
    QVector<float> m_ionMobilityStDevs;
    QVector<float> m_ms2PPMStDevs;
    QVector<float> m_weights;

    bool m_excludeDecoys;

};



#endif //MSCALIBRATOMATICSETTERTRON_H
