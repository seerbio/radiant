//
// Created by anichols on 9/5/23.
//

#ifndef PYTHIADIACPP_CANDIDATEPROCESSERTRON_H
#define PYTHIADIACPP_CANDIDATEPROCESSERTRON_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "PeakIntegratomatic.h"
#include "PythiaParameterReader.h"
#include "TurboXIC.h"

using namespace Error;

namespace ScoredCandidateNamespace {

    const QString COS_SIM_SUM = QStringLiteral("cosineSimSum");
    const QString IS_DECOY = QStringLiteral("isDecoy");
    const QString PEP_STR_W_MODS = QStringLiteral("peptideStringWithMods");
    const QString CHARGE = QStringLiteral("charge");
    const QString MASS = QStringLiteral("mass");
    const QString SCAN_NUM = QStringLiteral("scanNumber");
    const QString SCAN_TIME = QStringLiteral("scanTime");
    const QString SCAN_ION_CNT = QStringLiteral("scanIonCount");
    const QString X_CORR = QStringLiteral("xCorr");
    const QString SCAN_TIME_PRED = QStringLiteral("scanTimePredicted");
    const QString IRT_PRED = QStringLiteral("iRTPredicted");

    const QString MZ_SRCH_V = QStringLiteral("mzSearchedVec");
    const QString THEO_INTZ_V = QStringLiteral("theoIntensityVec");
    const QString MZ_FND_MEAN_V = QStringLiteral("mzFoundMeanVec");
    const QString MZ_FND_STDEV_V = QStringLiteral("mzFoundStDevVec");
    const QString INTZ_FND_MAX_V = QStringLiteral("intensityFoundMaxVec");
    const QString FRAME_IND_MAX_DIV_ANCH_V = QStringLiteral("frameIndexMaxDiffFromAnchorVec");
    const QString COS_SIM_SUM_ANCH_V = QStringLiteral("cosineSimToAnchorVec");
    const QString PK_PNT_CNT_FND_V = QStringLiteral("peakPointCountFoundVec");
    const QString FRAG_FRQ_V = QStringLiteral("fragmenFrequencyVec");
    const QString TARGET_KEY = QStringLiteral("targetKey");

    const QString KL_DIV_SUM = QStringLiteral("klDivSum");
    const QString KL_DIV_TO_ANCHOR = QStringLiteral("klDivToAnchorVec");
    const QString KL_DIV_SPECTRUM = QStringLiteral("klDivSpectrum");
    const QString COSINE_SIM_SPECTRUM = QStringLiteral("cosineSimSpectrum");
    const QString COSINE_SIM_MS1 = QString("cosineSimMS1");
    const QString THEO_FRAG_CNT = QString("theoFragmentCount");
    const QString DISC_SCORE = QString("discriminateScore");
    const QString Q_VAL = QString("qValue");
    const QString DECOY_RATIO = QString("decoyRatio");

    const QStringList keysToCheck = {
            COS_SIM_SUM,
            IS_DECOY,
            PEP_STR_W_MODS,
            CHARGE,
            MASS,
            SCAN_NUM,
            SCAN_TIME,
            SCAN_ION_CNT,
            MZ_SRCH_V,
            THEO_INTZ_V,
            MZ_FND_MEAN_V,
            MZ_FND_STDEV_V,
            INTZ_FND_MAX_V,
            FRAME_IND_MAX_DIV_ANCH_V,
            COS_SIM_SUM_ANCH_V,
            PK_PNT_CNT_FND_V,
            FRAG_FRQ_V,
            SCAN_TIME_PRED,
            IRT_PRED,
            TARGET_KEY,
            X_CORR,
            KL_DIV_SUM,
            KL_DIV_TO_ANCHOR,
            KL_DIV_SPECTRUM,
            COSINE_SIM_SPECTRUM,
            COSINE_SIM_MS1,
            THEO_FRAG_CNT,
            DISC_SCORE,
            Q_VAL,
            DECOY_RATIO
    };
}

struct FILEREADERSLIB_EXPORTS ScoredCandidate : public ParquetReaderInputBase {

    PeptideStringWithMods peptideStringWithMods;
    double cosineSimSum = -1.0;
    bool isDecoy = false;
    Charge charge = -1;
    double mass = -1.0;
    ScanNumber scanNumber = -1;
    ScanTime scanTime = -1.0;
    int scanIonCount = -1;
    double xCorr = -1.0;
    double klDivSum = -1.0;
    double klDivSpectrum = -1.0;
    double cosineSimSpectrum = -1.0;
    double cosineSimMS1 = -1.0;
    ScanTime scanTimePredicted = -1.0;
    double iRTPredicted = -1.0;
    QVector<double> mzSearchedVec;
    QVector<double> theoIntensityVec;
    QVector<double> mzFoundMeanVec;
    QVector<double> mzFoundStDevVec;
    QVector<double> intensityFoundMaxVec;
    QVector<int> frameIndexMaxDiffFromAnchorVec;
    QVector<double> cosineSimToAnchorVec;
    QVector<double> klDivToAnchorVec;
    QVector<int> peakPointCountFoundVec;
    QVector<double> fragmentFrequencyVec;
    QString targetKey;
    int theoFragmentCount = -1;
    double discriminateScore = -1.0;
    double qValue = -1.0;
    double decoyRatio = -1.0;

//    double discScore = -1.0;
//    double pVal = -1.0;
//    double frameError = -1.0;
//    double tTestVal = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace ScoredCandidateNamespace;

        return {
                {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
                {COS_SIM_SUM, QVariant(cosineSimSum)},
                {IS_DECOY, QVariant(isDecoy)},
                {CHARGE, QVariant(charge)},
                {MASS, QVariant(mass)},
                {SCAN_NUM, QVariant(scanNumber)},
                {SCAN_TIME, QVariant(scanTime)},
                {SCAN_ION_CNT, QVariant(scanIonCount)},
                {MZ_SRCH_V, QVariant(qVectorToQByteArray(mzSearchedVec))},
                {THEO_INTZ_V, QVariant(qVectorToQByteArray(theoIntensityVec))},
                {MZ_FND_MEAN_V, QVariant(qVectorToQByteArray(mzFoundMeanVec))},
                {MZ_FND_STDEV_V, QVariant(qVectorToQByteArray(mzFoundStDevVec))},
                {INTZ_FND_MAX_V, QVariant(qVectorToQByteArray(intensityFoundMaxVec))},
                {FRAME_IND_MAX_DIV_ANCH_V, QVariant(qVectorToQByteArray(frameIndexMaxDiffFromAnchorVec))},
                {COS_SIM_SUM_ANCH_V, QVariant(qVectorToQByteArray(cosineSimToAnchorVec))},
                {PK_PNT_CNT_FND_V, QVariant(qVectorToQByteArray(peakPointCountFoundVec))},
                {FRAG_FRQ_V, QVariant(qVectorToQByteArray(fragmentFrequencyVec))},
                {SCAN_TIME_PRED, QVariant(scanTimePredicted)},
                {IRT_PRED , QVariant(iRTPredicted)},
                {TARGET_KEY, QVariant(targetKey)},
                {X_CORR, QVariant(xCorr)},
                {KL_DIV_SUM, QVariant(klDivSum)},
                {KL_DIV_TO_ANCHOR, QVariant(qVectorToQByteArray(klDivToAnchorVec))},
                {KL_DIV_SPECTRUM, QVariant(klDivSpectrum)},
                {COSINE_SIM_SPECTRUM, QVariant(cosineSimSpectrum)},
                {COSINE_SIM_MS1, QVariant(cosineSimMS1)},
                {THEO_FRAG_CNT, QVariant(theoFragmentCount)},
                {DISC_SCORE, QVariant(discriminateScore)},
                {Q_VAL, QVariant(qValue)},
                {DECOY_RATIO, QVariant(decoyRatio)}
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

        cosineSimSum = dataMap.value(COS_SIM_SUM).toDouble();
        isDecoy = dataMap.value(IS_DECOY).toBool();
        peptideStringWithMods = dataMap.value(PEP_STR_W_MODS).toString();
        charge = dataMap.value(CHARGE).toInt();
        mass = dataMap.value(MASS).toDouble();
        scanNumber = dataMap.value(SCAN_NUM).toInt();
        scanTime = dataMap.value(SCAN_TIME).toDouble();
        scanIonCount = dataMap.value(SCAN_ION_CNT).toInt();
        scanTimePredicted = dataMap.value(SCAN_TIME_PRED).toDouble();
        iRTPredicted = dataMap.value(IRT_PRED).toDouble();
        mzSearchedVec = bytesArrayToQVector<double>(dataMap.value(MZ_SRCH_V).toByteArray());
        theoIntensityVec = bytesArrayToQVector<double>(dataMap.value(THEO_INTZ_V).toByteArray());
        mzFoundMeanVec = bytesArrayToQVector<double>(dataMap.value(MZ_FND_MEAN_V).toByteArray());
        mzFoundStDevVec = bytesArrayToQVector<double>(dataMap.value(MZ_FND_STDEV_V).toByteArray());
        intensityFoundMaxVec = bytesArrayToQVector<double>(dataMap.value(INTZ_FND_MAX_V).toByteArray());
        frameIndexMaxDiffFromAnchorVec = bytesArrayToQVector<int>(dataMap.value(FRAME_IND_MAX_DIV_ANCH_V).toByteArray());
        cosineSimToAnchorVec = bytesArrayToQVector<double>(dataMap.value(COS_SIM_SUM_ANCH_V).toByteArray());
        peakPointCountFoundVec = bytesArrayToQVector<int>(dataMap.value(PK_PNT_CNT_FND_V).toByteArray());
        fragmentFrequencyVec = bytesArrayToQVector<double>(dataMap.value(FRAG_FRQ_V).toByteArray());
        targetKey = dataMap.value(TARGET_KEY).toString();
        xCorr = dataMap.value(X_CORR).toDouble();
        klDivSum = dataMap.value(KL_DIV_SUM).toDouble();
        klDivSpectrum = dataMap.value(KL_DIV_SPECTRUM).toDouble();
        cosineSimMS1 = dataMap.value(COSINE_SIM_MS1).toDouble();
        cosineSimSpectrum = dataMap.value(COSINE_SIM_SPECTRUM).toDouble();
        klDivToAnchorVec = bytesArrayToQVector<double>(dataMap.value(KL_DIV_TO_ANCHOR).toByteArray());
        theoFragmentCount = dataMap.value(THEO_FRAG_CNT).toInt();
        discriminateScore = dataMap.value(DISC_SCORE).toDouble();
        qValue = dataMap.value(Q_VAL).toDouble();
        decoyRatio = dataMap.value(DECOY_RATIO).toDouble();

        ERR_RETURN
    }
};


class ALGORITHMSLIB_EXPORTS CandidateProcessertron {

public:

    CandidateProcessertron();
    ~CandidateProcessertron() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
            const MsFrame &msFrame,
            const MsFrame &msFrameMS1,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const QMap<MzHashed, FrequencyPercent> &fragmentFrequencies,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey
            );

    Err init(
            const PythiaParameters &pythiaParameters,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
            const MsFrame &msFrame,
            const MsFrame &msFrameMS1,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const QMap<MzHashed, FrequencyPercent> &fragmentFrequencies,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QMap<PeptideStringWithMods, ScanTime> &fragPredsPredictedScanTime
    );

    Err processCandidateTarget(
            const CandidatePeptide &candidatePeptideTarget,
            ScoredCandidate *scoredCandidate
            );

    Err processCandidateDecoy(
            const CandidatePeptide &candidatePeptideDecoy,
            ScanTime scanTime,
            ScoredCandidate *scoredCandidateDecoy
            );

private:

    Err filterSummedVecPeakIntegrationsByPredictedScanTime(
            ScanTime predictedScanTime,
            double windowWidthMinutes,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
    );

    Err findCandidateIntegrations(
            const QVector<double> &summedMatToVec,
            int minFoundMzPeaks,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
            );

    Err findCandidateIntegrations(
            const QVector<double> &summedMatToVec,
            int minFoundMzPeaks,
            double scanTime,
            double scanTimeTolerance,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
    );

    Err buildScores(
            const CandidatePeptide &candidatePeptide,
            const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
            const QVector<double> &summedMatVecToVec,
            ScoredCandidate *scoredCandidate
            );

    Err calculateEmpiricalMzStats(
            const CandidatePeptide &candidatePeptide,
            PeakIntegrationIndexes &bestPeakIntegrationIndexes,
            QVector<double> *mzMeanValsFound,
            QVector<double> *stdMeanValsFound,
            QVector<double> *mzValsSearched,
            QVector<double> *ppmDifference,
            QVector<double> *theoApexIntensity
            );

private:

    PythiaParameters m_pythiaParameters;
    QMap<MzHashed, XICPoints> m_mzHashedVsXICPoints;
    QMap<MzHashed, QVector<double>> m_mzHashedVsIonPresence;

    PeakIntegratomatic m_peakIntegratomatic;
    MsFrame m_msFrame;
    MsFrame m_msFrameMS1;
    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;
    TurboXIC m_turboXICMS1;
    QMap<MzHashed, FrequencyPercent> m_fragmentFrequencies;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;

    QMap<PeptideStringWithMods, ScanTime> m_fragPredsPredictedScanTime;

    const int m_peakWidthMin;

};


#endif //PYTHIADIACPP_CANDIDATEPROCESSERTRON_H
