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

    Err buildCalibration(
            double calibrationTrainingFraction,
            bool useExtendedScores,
            TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron
            );

    Err recalibrateMzVals(
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrame,
            QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1,
            TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
            MsReaderPointerAcc *msReaderPointerAcc
            );

    Err setTargetDecoyCandidateScores(
        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
        int topNMS2Ions,
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

    static Err setQValueForCandidates(QVector<CandidateScores> *candidateScores);

    Err optimizeParameters(TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron);

    Err mainAnalysis(
            TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
            QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers
            );

    Err removeInterferingCandidates(
            MsReaderPointerAcc *msReaderPointerAcc,
            const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointers,
            QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointersUpdatedTargets
            );

    Err applyNeuralNetClassifier(
            const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointers,
            const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointersFDRFiltered,
            const QPair<double, double> &scanTimeMinMax,
            bool reportDecoys,
            QVector<CandidateScores> *candidateScoreClassifier
            );

    Err updateProteinGroupAnnotation(
            const QString &fastaFilePath,
            QVector<CandidateScores> *candidateScores
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
