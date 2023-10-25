//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsFrame.h"
#include "MsReaderPointerAcc.h"
#include "ProteinDigestomatic.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"


using namespace Error;

class MsReaderParquet;
class CandidatePeptide;
class ScoredCandidate;
class TargetDecoyCandidatePairScoretron;


class WORKFLOWSLIB_EXPORTS PythiaDIAWorkflow {

public:

    PythiaDIAWorkflow();
    ~PythiaDIAWorkflow() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &fastaUri
            );

    Err processFile(const QString &msDataFilePath);

private:

    Err deisotopeScans(MsReaderPointerAcc *msReaderPointerAcc);

    Err buildCalibration(TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron);

    Err setTargetDecoyCandidateScores(
            TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
            int topNMS2Ions,
            double calibrationTrainingFraction,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers,
            QMap<QString, int> *fdrVsCount
            );

    Err setDiscriminantScoreForCandidates(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            int theoMzIonsSize
            );

    static Err setQValueForCandidates(const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs);

    Err optimizeParameters(TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron);

    Err mainAnalysis(
            MsReaderPointerAcc *msReaderPointerAcc,
            QVector<ScoredCandidate> *scoredCandidatesAll
            );

    Err removeInterferingCandidates(
            MsReaderPointerAcc *msReaderPointerAcc,
            const QVector<ScoredCandidate> &scoredCandidatesTargetsFDRThresholded,
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            QVector<ScoredCandidate> *scoredCandidatesAllUpdated
            );

    Err applyNeuralNetClassifier(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            const QVector<ScoredCandidate> &scoredCandidatesCulled,
            const QPair<double, double> &scanTimeMinMax,
            QVector<ScoredCandidate> *scoredCandidatesClassifier
            );

    Err updateProteinGroupAnnotation(
            const QString &fastaFilePath,
            QVector<ScoredCandidate> *scoredCandidates
    );

private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_fastaUri;

    TargetDecoyCandidatePairManager m_targetDecoyCandidatePairManager;
    QVector<MsScanInfo> m_msScanInfos;
    QPair<double, double> m_scanTimeMinMax;

    MsCalibratomatic m_msCalibratomatic;

    const int m_minTopNMs2Ions;
    const bool m_byIonsOnly;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
