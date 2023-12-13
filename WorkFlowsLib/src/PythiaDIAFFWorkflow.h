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


class MsScanInfo;
class MsReaderPointerAcc;
class TargetDecoyCandidatePair;

enum class QValueScoreType {
    DiscriminantScore,
    NNClassifierScore
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

    Err buildCalibration(MsReaderPointerAcc *msReaderPointerAcc);

    Err buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            const QVector<MsScanInfo> &msScanInfos,
            double selectionFraction,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
            );

    Err setDiscriminantScoreForCandidates(
        const QPair<double, double> &scanTimeMinMax,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int theoMzIonsSize,
        QVector<CandidateScores> *candidateScores
        );

    static Err setQValueForCandidates(
            const QValueScoreType &qValueScoreType,
            QVector<CandidateScores> *candidateScores
            );


private:

    TargetDecoyCandidatePairManager m_targetDecoyCandidatePairManager;
    TargetDecoyCandidatePairScoretron m_targetDecoyCandidatePairScoretron;
    MsCalibratomatic m_msCalibratomatic;

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_fastaUri;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
