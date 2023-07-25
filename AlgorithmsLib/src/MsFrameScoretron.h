//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRON_H
#define PYTHIADIACPP_MSFRAMESCORETRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FeatureFinderHillBuilder.h"
#include "FragLibIonRTree.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MS2ChargeDeconvolvotron.h"
#include "MsFrame.h"
#include "ParquetReader.h"
#include "PeakIntegratomatic.h"
#include "PSMsReader.h"
#include "PythiaParameterReader.h"


using namespace Error;

class TandemDeconvolverResult;
class TurboXIC;

namespace ScoredCandidateNamespace {

    const QString FREQ_PCT_SUM = QStringLiteral("frequencyPercentSum");
    const QString FREQ_PCT_SUM_POSSIBLE = QStringLiteral("frequencyPercentSumBestPossible");
    const QString KL_DIV = QStringLiteral("klDivergence");
    const QString COS_SIM = QStringLiteral("cosineSim");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString PEP_STR_W_MODS = QStringLiteral("peptideStringWithMods");
    const QString CHARGE = QStringLiteral("charge");
    const QString MASS = QStringLiteral("mass");
    const QString SCAN_NUM = QStringLiteral("scanNumber");
    const QString SCAN_TIME = QStringLiteral("scanTime");

    const QStringList keysToCheck = {
            FREQ_PCT_SUM,
            FREQ_PCT_SUM_POSSIBLE,
            KL_DIV,
            COS_SIM,
            IS_DECOY,
            PEP_STR_W_MODS,
            CHARGE,
            MASS,
            SCAN_NUM,
            SCAN_TIME
    };
}

struct FILEREADERSLIB_EXPORTS ScoredCandidate : public ParquetReaderInputBase {

    double frequencyPercentSum = -1.0;
    double frequencyPercentSumBestPossible = -1.0;
    double klDivergence = -1.0;
    double cosineSim = -1.0;
    bool isDecoy = false;
    PeptideStringWithMods peptideStringWithMods;
    Charge charge = -1;
    double mass = -1.0;
    ScanNumber scanNumber = -1;
    ScanTime scanTime = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace ScoredCandidateNamespace;

        return {
                {FREQ_PCT_SUM, QVariant(frequencyPercentSum)},
                {FREQ_PCT_SUM_POSSIBLE, QVariant(frequencyPercentSumBestPossible)},
                {KL_DIV, QVariant(klDivergence)},
                {COS_SIM, QVariant(cosineSim)},
                {IS_DECOY, QVariant(isDecoy)},
                {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
                {CHARGE, QVariant(charge)},
                {MASS, QVariant(mass)},
                {SCAN_NUM, QVariant(scanNumber)},
                {SCAN_TIME, QVariant(scanTime)}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace ScoredCandidateNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        frequencyPercentSum = dataMap.value(FREQ_PCT_SUM).toDouble();
        frequencyPercentSumBestPossible = dataMap.value(FREQ_PCT_SUM_POSSIBLE).toDouble();
        klDivergence = dataMap.value(KL_DIV).toDouble();
        cosineSim = dataMap.value(COS_SIM).toDouble();
        isDecoy = dataMap.value(IS_DECOY).toBool();
        peptideStringWithMods = dataMap.value(PEP_STR_W_MODS).toString();
        charge = dataMap.value(CHARGE).toInt();
        mass = dataMap.value(MASS).toDouble();
        scanNumber = dataMap.value(SCAN_NUM).toInt();
        scanTime = dataMap.value(SCAN_TIME).toDouble();

        ERR_RETURN
    }
};


class ALGORITHMSLIB_EXPORTS MsFrameScoretron {

public:

    friend class MissingPeptideManualTroubleshooter;
    friend class MsFrameScoretronProcessormaticTests;
    friend class MsFrameScoretronTests;

    MsFrameScoretron() = default;
    ~MsFrameScoretron() = default;

    Err init(
            const PythiaParameters &params,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const QString &fragLibBackgroundFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
            );

    Err init(
            const PythiaParameters &params,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const QString &fragLibBackgroundFilePath,
            const QString &iRTRecalibrationFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
    );

    Err scoreFrameCandidates(QMap<FrameIndex , QVector<ScoredCandidate>> *frameIndexVsScoredCandidates);


private:

    Err initMsFrame(
            const QString &msDataFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
            );

    Err buildMsFrameApexTurboXIC(TurboXIC *turboXic);

    Err iterateApexScanPoints(
            const QMap<FrameIndex, ScanPoints> &apexScanPoints,
            QMap<FrameIndex , QVector<ScoredCandidate>> *frameIndexVsScoredCandidates
            );

    Err extractFragLibIonsForScanPoints(
            const ScanPoints &scanPoints,
            ScanTime scanTime,
            FragLibIonRTree *fragLibIonRTree,
            QMap<PeptideId, QVector<FragLibIon>> *peptideIdVsFragLibIonsForFrameIndexOutput
            ) const;

    Err extractScores(
            const QMap<PeptideId, QVector<FragLibIon>> &peptideIdVsFragLibIonsForFrameIndex,
            FrameIndex frameIndex,
            QVector<ScoredCandidate> *scoredCandidates
    );


private:

    QMap<PeptideStringWithMods, MS2IonsSeparated> m_fragPreds;
    QMap<PeptideStringWithMods, bool> m_fragPredsIsDecoy;
    QMap<PeptideStringWithMods, double> m_fragPredsMass;
    QMap<PeptideStringWithMods, IRT> m_fragPredsIRT;

    QMap<PeptideStringWithMods, MS2IonsSeparated> m_fragPredsBackground;
    QMap<PeptideStringWithMods, bool> m_fragPredsBackgroundIsDecoy;

    FragLibIonRTree m_fragLibIonRTree;
    FragLibIonRTree m_fragLibIonRTreeBackground;

    MS2ChargeDeconvolvotron m_ms2ChargeDeconvolvotron;

    MsFrame m_msFrame;
    QVector<FeatureFinderHill> m_featureFinderHills;
    FeatureFinderHillBuilder m_featureFinderHillBuilder;

    PythiaParameters m_params;
    QString m_msDataFilePath;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    QPair<double, double> m_mzTargetStartStop;
    double m_mzStartStopMean;
};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
