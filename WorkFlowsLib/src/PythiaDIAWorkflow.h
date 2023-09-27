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
#include "MsReaderPointerFactory.h"
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
            const QString &fragLibUri
            );

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &iRTReCalFilePath
    );

    Err processFile(const QString &msDataFilePath);


private:

    Err deisotopeScans(MsReaderParquet *msReaderParquet);

    Err buildCalibration(MsReaderParquet *msReaderParquet);

    Err extractionLoopLogic(
            const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptide,
            double fdrThreshold,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            int topNMs2IonsMainAnalysis,
            MsReaderParquet *msReaderParquet,
            QVector<ScoredCandidate> *scoredCandidatesAll,
            QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded
            );

    Err extractTargetDecoyData(
            const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptideCalibration,
            const PythiaParameters &pythiaParameters,
            MsReaderParquet *msReaderParquet,
            QVector<ScoredCandidate> *combinedResults
    );

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

    Err optimizeParameters(MsReaderParquet *msReaderParquet);

    Err mainAnalysis(
            MsReaderParquet *msReaderParquet,
            QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded,
            QVector<ScoredCandidate> *scoredCandidatesAll
            );

    Err removeInterferingCandidates(
            MsReaderParquet *msReaderParquet,
            const QVector<ScoredCandidate> &scoredCandidatesTargetsFDRThresholded,
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            QVector<ScoredCandidate> *scoredCandidatesAllUpdated
            );

    Err applyNeuralNetClassifier(
            const QVector<ScoredCandidate> &scoredCandidatesCulled,
            MsReaderParquet *msReaderParquet,
            QVector<ScoredCandidate> *scoredCandidatesClassifier
            );

    Err returnAllCandidatesScoredFullFragIons(
            MsReaderParquet *msReaderParquet,
            QVector<ScoredCandidate> *scoredCandidatesAllTemp
            );

    Err updateProteinGroupAnnotation(
            const QString &fastaFilePath,
            QVector<ScoredCandidate> *scoredCandidates
    );

    static Err buildUniqueMsInfoScanKeyVsScanPoints(
            MsReaderParquet *msReaderParquet,
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
            QMap<ScanNumber, ScanTime> *scanNumberVsScanTime,
            QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
            );

private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_iRTReCalFilePath;

    FragLibReader m_fragLibReader;
    QVector<MsScanInfo> m_msScanInfos;

    MsCalibratomatic m_msCalibratomatic;

    const int m_minTopNMs2Ions;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
