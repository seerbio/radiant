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


using namespace Error;

class MsReaderParquet;
class ScoredCandidate;


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

    Err buildCalibration(MsReaderPointerAcc *msReaderPointerAcc);

    Err buildCandidates(
            int topNMs2Ions,
            double selectionListFraction,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide
    );

    Err buildCandidates(
            const QVector<PeptideStringWithMods> &inclusionList,
            int topNMs2Ions,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide
    );

    Err optimizeParameters(MsReaderPointerAcc *msReaderPointerAcc);

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

    FragLibReader m_fragLibReader;
    QVector<MsScanInfo> m_msScanInfos;

    MsCalibratomatic m_msCalibratomatic;

    const int m_minTopNMs2Ions;
    const bool m_byIonsOnly;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
