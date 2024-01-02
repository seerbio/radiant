//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretron.h"

using namespace Error;

class CandidateScores;
class MsScanInfo;
class MsReaderPointerAcc;
class TargetDecoyCandidatePair;

enum class QValueScoreType {
    DiscriminantScore,
    NNClassifierScore
};

struct KarnnNNTarget {
    QString seq;
    float nnScore = 0.0;
    bool isDecoy = false;
    QVector<float> scoreVec;
    int index = -1;
};

class WORKFLOWSLIB_EXPORTS PythiaDIAFFWorkflow {

public:

    friend class PythiaDIAFFWorkflowTests;

    PythiaDIAFFWorkflow();
    ~PythiaDIAFFWorkflow() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &fastaUri
            );

    Err processFile(const QString &msDataFilePath);


private:

    Err buildCalibration(
            MsReaderPointerAcc *msReaderPointerAcc,
            QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
            QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1,
            QVector<CandidateScores*> *candidateScoresForTrainings
            );

    Err buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            const QVector<MsScanInfo> &msScanInfos,
            double selectionFraction,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
            );

    Err setDiscriminantScoreForCandidates(
        bool useExtendedScores,
        bool useNeuralNetworkScores
        );

    Err recalibrateMzVals(
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1,
        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
        MsReaderPointerAcc *msReaderPointerAcc
        );

    Err optimizeParameters(
            const QVector<CandidateScores*> &candidateScoresTrainings,
            MsReaderPointerAcc *msReaderPointerAcc
            );

    Err mainAnalysis(
            MsReaderPointerAcc *msReaderPointerAcc,
            int *targetCountBelowFDRThresholdOnePercent
            );

    Err buildCandidateScoresPtrs(QVector<CandidateScores*> *candidateScoresPntrs);

    Err removeInterferingCandidates(
            int ionsSharedToReject,
            QVector<CandidateScores*> *candidates
            );

    Err applyNeuralNetClassifier(
            const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            int seed,
            QVector<CandidateScores*> *candidateScoreClassifier
            );

    Err updateProteinGroupAnnotation(
        const QString &fastaFilePath,
        QVector<CandidateScores*> *candidateScores
        );


    static Err setQValueForCandidates(
            const QValueScoreType &qValueScoreType,
            QVector<CandidateScores*> *candidateScores
            );


private:

    TargetDecoyCandidatePairManager m_targetDecoyCandidatePairManager;
    TargetDecoyCandidatePairScoretron m_targetDecoyCandidatePairScoretron;
    MsCalibratomatic m_msCalibratomatic;

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_fastaUri;

    int m_minTopNMs2Ions;

    QVector<CandidateScores> m_candidateScores;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
