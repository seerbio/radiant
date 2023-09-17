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

namespace NeuralNetDataNamespace {
    const QString SCORES = QStringLiteral("scores");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString IS_TEST = QStringLiteral("isTest");
    const QString PEP_STR = QStringLiteral("peptideStringWithMods");
    const QString CHARGE = QStringLiteral("charge");
    const QString TARGET_KEY = QStringLiteral("targetKey");
}

struct NeuralNetData : public ParquetReaderInputBase {

public:

    QVector<double> scores;
    bool isDecoy = false;
    bool isTest = false;
    Charge charge = -1;
    QString targetKey;
    PeptideStringWithMods peptideStringWithMods;

    QMap<QString, QVariant> map() override {

        using namespace NeuralNetDataNamespace;

        return {
                {SCORES, QVariant(qVectorToQByteArray(scores))},
                {IS_DECOY, QVariant(isDecoy)},
                {IS_TEST, QVariant(isTest)},
                {CHARGE, QVariant(charge)},
                {TARGET_KEY, QVariant(targetKey)},
                {PEP_STR, QVariant(peptideStringWithMods)}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        ERR_INIT

        using namespace NeuralNetDataNamespace;

        const QStringList keysToCheck = {
                SCORES,
                IS_DECOY,
                IS_TEST,
                CHARGE,
                TARGET_KEY,
                PEP_STR
        };

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        scores = bytesArrayToQVector<double>(dataMap.value(SCORES).toByteArray());
        isDecoy = dataMap.value(IS_DECOY).toBool();
        isTest = dataMap.value(IS_TEST).toBool();
        peptideStringWithMods = dataMap.value(PEP_STR).toString();
        charge = dataMap.value(CHARGE).toInt();
        targetKey = dataMap.value(TARGET_KEY).toString();

        ERR_RETURN
    }
};

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
            MsReaderParquet *msReaderParquet
            );

    Err returnAllCandidatesScoredFullFragIons(
            MsReaderParquet *msReaderParquet,
            QVector<ScoredCandidate> *scoredCandidatesAllTemp
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
