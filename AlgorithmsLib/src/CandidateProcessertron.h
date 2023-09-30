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

    const QString COS_SIM_SUM_100 = QStringLiteral("cosineSimSum100");
    const QString COS_SIM_SUM_45 = QStringLiteral("cosineSimSum45");
    const QString COS_SIM_SUM_20 = QStringLiteral("cosineSimSum20");
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
    const QString TARGET_KEY = QStringLiteral("targetKey");

    const QString KL_DIV_SUM = QStringLiteral("klDivSum");
    const QString KL_DIV_TO_ANCHOR = QStringLiteral("klDivToAnchorVec");
    const QString KL_DIV_SPECTRUM = QStringLiteral("klDivSpectrum");
    const QString COSINE_SIM_SPECTRUM = QStringLiteral("cosineSimSpectrum");
    const QString COSINE_SIM_100_MS1 = QString("cosineSim100MS1");
    const QString COSINE_SIM_45_MS1 = QString("cosineSim45MS1");
    const QString COSINE_SIM_20_MS1 = QString("cosineSim20MS1");
    const QString THEO_FRAG_CNT = QString("theoFragmentCount");
    const QString DISC_SCORE = QString("discriminateScore");
    const QString Q_VAL = QString("qValue");
    const QString DECOY_RATIO = QString("decoyRatio");
    const QString MATRIX_WEIGHT = QString("matrixWeight");
    const QString MATRIX_PVAL = QString("matrixPVal");
    const QString MATRIX_ERROR = QString("matrixError");
    const QString SCAN_NUM_CAND_CNT = QString("scanNumberCandidateCount");
    const QString CLASSIFIER_SCORE = QString("classifierScore");
    const QString PROTEIN_GRP = QString("proteinGroup");

    const QStringList keysToCheck = {
            COS_SIM_SUM_100,
            COS_SIM_SUM_45,
            COS_SIM_SUM_20,
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
            SCAN_TIME_PRED,
            IRT_PRED,
            TARGET_KEY,
            X_CORR,
            KL_DIV_SUM,
            KL_DIV_TO_ANCHOR,
            KL_DIV_SPECTRUM,
            COSINE_SIM_SPECTRUM,
            COSINE_SIM_100_MS1,
            COSINE_SIM_45_MS1,
            COSINE_SIM_20_MS1,
            THEO_FRAG_CNT,
            DISC_SCORE,
            Q_VAL,
            DECOY_RATIO,
            MATRIX_WEIGHT,
            MATRIX_PVAL,
            MATRIX_ERROR,
            SCAN_NUM_CAND_CNT,
            CLASSIFIER_SCORE,
            PROTEIN_GRP
    };
}

struct FILEREADERSLIB_EXPORTS ScoredCandidate : public ParquetReaderInputBase {

    PeptideStringWithMods peptideStringWithMods;
    double cosineSimSum100 = -1.0;
    double cosineSimSum45 = -1.0;
    double cosineSimSum20 = -1.0;
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
    double cosineSim100MS1 = -1.0;
    double cosineSim45MS1 = -1.0;
    double cosineSim20MS1 = -1.0;
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
    QString targetKey;
    int theoFragmentCount = -1;
    double discriminateScore = -1.0;
    double qValue = 1.0;
    double decoyRatio = -1.0;
    double matrixWeight = -1.0;
    double matrixPValue = 1.0;
    double matrixError = 1.0;
    int scanNumberCandidateCount = -1;
    double classifierScore = -1.0;
    QString proteinGroup;

    QMap<QString, QVariant> map() override {

        using namespace ScoredCandidateNamespace;

        return {
                {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
                {COS_SIM_SUM_100, QVariant(cosineSimSum100)},
                {COS_SIM_SUM_45, QVariant(cosineSimSum45)},
                {COS_SIM_SUM_20, QVariant(cosineSimSum20)},
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
                {SCAN_TIME_PRED, QVariant(scanTimePredicted)},
                {IRT_PRED , QVariant(iRTPredicted)},
                {TARGET_KEY, QVariant(targetKey)},
                {X_CORR, QVariant(xCorr)},
                {KL_DIV_SUM, QVariant(klDivSum)},
                {KL_DIV_TO_ANCHOR, QVariant(qVectorToQByteArray(klDivToAnchorVec))},
                {KL_DIV_SPECTRUM, QVariant(klDivSpectrum)},
                {COSINE_SIM_SPECTRUM, QVariant(cosineSimSpectrum)},
                {COSINE_SIM_100_MS1, QVariant(cosineSim100MS1)},
                {COSINE_SIM_45_MS1, QVariant(cosineSim45MS1)},
                {COSINE_SIM_20_MS1, QVariant(cosineSim20MS1)},
                {THEO_FRAG_CNT, QVariant(theoFragmentCount)},
                {DISC_SCORE, QVariant(discriminateScore)},
                {Q_VAL, QVariant(qValue)},
                {DECOY_RATIO, QVariant(decoyRatio)},
                {MATRIX_WEIGHT, QVariant(matrixWeight)},
                {MATRIX_PVAL, QVariant(matrixPValue)},
                {MATRIX_ERROR, QVariant(matrixError)},
                {SCAN_NUM_CAND_CNT, QVariant(scanNumberCandidateCount)},
                {PROTEIN_GRP, QVariant(proteinGroup)},
                {CLASSIFIER_SCORE, QVariant(classifierScore)}
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

        cosineSimSum100 = dataMap.value(COS_SIM_SUM_100).toDouble();
        cosineSimSum45 = dataMap.value(COS_SIM_SUM_45).toDouble();
        cosineSimSum20 = dataMap.value(COS_SIM_SUM_20).toDouble();
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
        targetKey = dataMap.value(TARGET_KEY).toString();
        xCorr = dataMap.value(X_CORR).toDouble();
        klDivSum = dataMap.value(KL_DIV_SUM).toDouble();
        klDivSpectrum = dataMap.value(KL_DIV_SPECTRUM).toDouble();
        cosineSim100MS1 = dataMap.value(COSINE_SIM_100_MS1).toDouble();
        cosineSim45MS1 = dataMap.value(COSINE_SIM_45_MS1).toDouble();
        cosineSim20MS1 = dataMap.value(COSINE_SIM_20_MS1).toDouble();
        cosineSimSpectrum = dataMap.value(COSINE_SIM_SPECTRUM).toDouble();
        klDivToAnchorVec = bytesArrayToQVector<double>(dataMap.value(KL_DIV_TO_ANCHOR).toByteArray());
        theoFragmentCount = dataMap.value(THEO_FRAG_CNT).toInt();
        discriminateScore = dataMap.value(DISC_SCORE).toDouble();
        qValue = dataMap.value(Q_VAL).toDouble();
        decoyRatio = dataMap.value(DECOY_RATIO).toDouble();
        matrixWeight = dataMap.value(MATRIX_WEIGHT).toDouble();
        matrixPValue = dataMap.value(MATRIX_PVAL).toDouble();
        matrixError = dataMap.value(MATRIX_ERROR).toDouble();
        scanNumberCandidateCount = dataMap.value(SCAN_NUM_CAND_CNT).toInt();
        classifierScore = dataMap.value(CLASSIFIER_SCORE).toDouble();
        proteinGroup = dataMap.value(PROTEIN_GRP).toString();

        ERR_RETURN
    }
};


class ALGORITHMSLIB_EXPORTS CandidateProcessertron {

public:

    CandidateProcessertron();
    ~CandidateProcessertron() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            int topNMS2Ions,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints100,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints45,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints20,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
            const MsFrame &msFrame,
            const MsFrame &msFrameMS1,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey
            );

    Err init(
            const PythiaParameters &pythiaParameters,
            int topNMS2Ions,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints100,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints45,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints20,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
            const MsFrame &msFrame,
            const MsFrame &msFrameMS1,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
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
    QMap<MzHashed, XICPoints> m_mzHashedVsXICPoints100;
    QMap<MzHashed, XICPoints> m_mzHashedVsXICPoints45;
    QMap<MzHashed, XICPoints> m_mzHashedVsXICPoints20;
    QMap<MzHashed, QVector<double>> m_mzHashedVsIonPresence;

    PeakIntegratomatic m_peakIntegratomatic;
    MsFrame m_msFrame;
    MsFrame m_msFrameMS1;
    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;
    TurboXIC m_turboXICMS1;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;

    QMap<PeptideStringWithMods, ScanTime> m_fragPredsPredictedScanTime;

    int m_topNMS2Ions;

    const int m_peakWidthMin;

};


#endif //PYTHIADIACPP_CANDIDATEPROCESSERTRON_H
