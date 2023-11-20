//
// Created by anichols on 10/20/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORES_H
#define PYTHIADIACPP_CANDIDATESCORES_H


#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"


namespace CandidateScoresNamespace {

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
    const QString SCAN_TIME_PRED = QStringLiteral("scanTimePredicted");
    const QString IRT_PRED = QStringLiteral("iRTPredicted");
    const QString MZ_SRCH_V = QStringLiteral("mzSearchedVec");
    const QString THEO_INTZ_V = QStringLiteral("theoIntensityVec");
    const QString MZ_FND_MEAN_V = QStringLiteral("mzFoundMeanVec");
    const QString MZ_FND_STDEV_V = QStringLiteral("mzFoundStDevVec");
    const QString INTZ_FND_MAX_V = QStringLiteral("intensityFoundMaxVec");
    const QString COS_SIM_SUM_ANCH_V = QStringLiteral("cosineSimToAnchorVec");
    const QString TARGET_KEY = QStringLiteral("targetKey");
    const QString KL_DIV_SUM = QStringLiteral("klDivSum");
    const QString KL_DIV_SPECTRUM = QStringLiteral("klDivSpectrum");
    const QString COSINE_SIM_SPECTRUM = QStringLiteral("cosineSimSpectrum");
    const QString COSINE_SIM_100_MS1 = QStringLiteral("cosineSim100MS1");
    const QString COSINE_SIM_100_MS1_ISO1 = QStringLiteral("cosineSim100MS1Iso1");
    const QString COSINE_SIM_100_MS1_ISO2 = QStringLiteral("cosineSim100MS1Iso2");
    const QString COSINE_SIM_45_MS1 = QStringLiteral("cosineSim45MS1");
    const QString COSINE_SIM_20_MS1 = QStringLiteral("cosineSim20MS1");
    const QString THEO_FRAG_CNT = QStringLiteral("theoFragmentCount");
    const QString DISC_SCORE = QStringLiteral("discriminateScore");
    const QString Q_VAL = QStringLiteral("qValue");
    const QString DECOY_RATIO = QStringLiteral("decoyRatio");
    const QString MATRIX_WEIGHT = QStringLiteral("matrixWeight");
    const QString MATRIX_PVAL = QStringLiteral("matrixPVal");
    const QString MATRIX_ERROR = QStringLiteral("matrixError");
    const QString SCAN_NUM_CAND_CNT = QStringLiteral("scanNumberCandidateCount");
    const QString CLASSIFIER_SCORE = QStringLiteral("classifierScore");
    const QString PROTEIN_GRP = QStringLiteral("proteinGroup");
    const QString PEAK_SHAPE_RATIO_1 = QStringLiteral("peakShapeRatio1");
    const QString PEAK_SHAPE_RATIO_2 = QStringLiteral("peakShapeRatio2");
    const QString PEAK_SHAPE_RATIO_3 = QStringLiteral("peakShapeRatio3");
    const QString B2_CORR = QStringLiteral("b2Corr");
    const QString B3_CORR = QStringLiteral("b3Corr");
    const QString B2B3_COSINE_SIM_SUM = QStringLiteral("b2b3CosineSimSum");
    const QString SHADOWS_COSINE_SIM_SUM = QStringLiteral("shadowsCosineSimSum");
    const QString COS_SIM_SUM_ANCH_SHADOW_V = QStringLiteral("cosineSimShadowsToAnchorVec");
    const QString SHADOW_INTZ_RATIO_VEC = QStringLiteral("shadowIntensityRatioVec");
    const QString ALL_MAX_IND_CNT = QStringLiteral("allignedMaxIndexesCount");


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
            COS_SIM_SUM_ANCH_V,
            SCAN_TIME_PRED,
            IRT_PRED,
            TARGET_KEY,
            KL_DIV_SUM,
            KL_DIV_SPECTRUM,
            COSINE_SIM_SPECTRUM,
            COSINE_SIM_100_MS1,
            COSINE_SIM_100_MS1_ISO1,
            COSINE_SIM_100_MS1_ISO2,
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
            PROTEIN_GRP,
            PEAK_SHAPE_RATIO_1,
            PEAK_SHAPE_RATIO_2,
            PEAK_SHAPE_RATIO_3,
            B2_CORR,
            B3_CORR,
            B2B3_COSINE_SIM_SUM,
            SHADOWS_COSINE_SIM_SUM,
            COS_SIM_SUM_ANCH_SHADOW_V,
            SHADOW_INTZ_RATIO_VEC,
            ALL_MAX_IND_CNT
    };
}

class ALGORITHMSLIB_EXPORTS CandidateScores : public ParquetReaderInputBase {

public:

    CandidateScores() = default;
    ~CandidateScores() = default;

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
    double klDivSum = -1.0;
    double klDivSpectrum = 10.0;
    double cosineSimSpectrum = -1.0;
    double cosineSim100MS1 = -1.0;
    double cosineSim45MS1 = -1.0;
    double cosineSim20MS1 = -1.0;
    double cosineSim100MS1Iso1 = -1.0;
    double cosineSim100MS1Iso2 = -1.0;
    ScanTime scanTimePredicted = -1.0;
    double iRTPredicted = -1.0;
    QVector<double> mzSearchedVec;
    QVector<double> theoIntensityVec;
    QVector<double> mzFoundMeanVec;
    QVector<double> mzFoundStDevVec;
    QVector<double> intensityFoundMaxVec;
    QVector<double> cosineSimToAnchorVec;
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
    double peakShapeRatio1 = -1.0;
    double peakShapeRatio2 = -1.0;
    double peakShapeRatio3 = -1.0;
    double b2Corr = -1.0;
    double b3Corr = -1.0;
    double b2b3CosineSimSum = -1.0;
    QVector<double> cosineSimShadowsToAnchorVec;
    QVector<double> shadowsIntensityRatioVec;
    double shadowsCosineSimSum = -1.0;
    int allignedMaxIndexesCount = -1;

    void clear() {
        peptideStringWithMods= "";
        mzSearchedVec.clear();
        theoIntensityVec.clear();
        mzFoundMeanVec.clear();
        mzFoundStDevVec.clear();
        intensityFoundMaxVec.clear();
        cosineSimToAnchorVec.clear();
        cosineSimShadowsToAnchorVec.clear();
        shadowsIntensityRatioVec.clear();
        targetKey = "";
        proteinGroup = "";
        cosineSimSum100 = -1.0;
        cosineSimSum45 = -1.0;
        cosineSimSum20 = -1.0;
        isDecoy = false;
        charge = -1;
        mass = -1.0;
        scanNumber = -1;
        scanTime = -1.0;
        scanIonCount = -1;
        klDivSum = -1.0;
        klDivSpectrum = -1.0;
        cosineSimSpectrum = -1.0;
        cosineSim100MS1 = -1.0;
        cosineSim45MS1 = -1.0;
        cosineSim20MS1 = -1.0;
        cosineSim100MS1Iso1 = -1.0;
        cosineSim100MS1Iso2 = -1.0;
        scanTimePredicted = -1.0;
        iRTPredicted = -1.0;
        theoFragmentCount = -1;
        discriminateScore = -1.0;
        qValue = 1.0;
        decoyRatio = -1.0;
        matrixWeight = -1.0;
        matrixPValue = 1.0;
        matrixError = 1.0;
        scanNumberCandidateCount = -1;
        classifierScore = -1.0;
        peakShapeRatio1 = -1.0;
        peakShapeRatio2 = -1.0;
        peakShapeRatio3 = -1.0;
        b2Corr = -1.0;
        b3Corr = -1.0;
        b2b3CosineSimSum = -1.0;
        shadowsCosineSimSum = -1.0;
        allignedMaxIndexesCount = -1;
    }

    QMap<QString, QVariant> map() override {

        using namespace CandidateScoresNamespace;

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
                {COS_SIM_SUM_ANCH_V, QVariant(qVectorToQByteArray(cosineSimToAnchorVec))},
                {SCAN_TIME_PRED, QVariant(scanTimePredicted)},
                {IRT_PRED , QVariant(iRTPredicted)},
                {TARGET_KEY, QVariant(targetKey)},
                {KL_DIV_SUM, QVariant(klDivSum)},
                {KL_DIV_SPECTRUM, QVariant(klDivSpectrum)},
                {COSINE_SIM_SPECTRUM, QVariant(cosineSimSpectrum)},
                {COSINE_SIM_100_MS1, QVariant(cosineSim100MS1)},
                {COSINE_SIM_100_MS1_ISO1, QVariant(cosineSim100MS1Iso1)},
                {COSINE_SIM_100_MS1_ISO2, QVariant(cosineSim100MS1Iso2)},
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
                {PEAK_SHAPE_RATIO_1, QVariant(peakShapeRatio1)},
                {PEAK_SHAPE_RATIO_2, QVariant(peakShapeRatio2)},
                {PEAK_SHAPE_RATIO_3, QVariant(peakShapeRatio3)},
                {CLASSIFIER_SCORE, QVariant(classifierScore)},
                {B2_CORR, QVariant(b2Corr)},
                {B3_CORR, QVariant(b3Corr)},
                {B2B3_COSINE_SIM_SUM, QVariant(b2b3CosineSimSum)},
                {SHADOWS_COSINE_SIM_SUM, QVariant(shadowsCosineSimSum)},
                {COS_SIM_SUM_ANCH_SHADOW_V, QVariant(qVectorToQByteArray(cosineSimShadowsToAnchorVec))},
                {SHADOW_INTZ_RATIO_VEC, QVariant(qVectorToQByteArray(shadowsIntensityRatioVec))},
                {ALL_MAX_IND_CNT, QVariant(allignedMaxIndexesCount)}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace CandidateScoresNamespace;

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
        cosineSimToAnchorVec = bytesArrayToQVector<double>(dataMap.value(COS_SIM_SUM_ANCH_V).toByteArray());
        targetKey = dataMap.value(TARGET_KEY).toString();
        klDivSum = dataMap.value(KL_DIV_SUM).toDouble();
        klDivSpectrum = dataMap.value(KL_DIV_SPECTRUM).toDouble();
        cosineSim100MS1 = dataMap.value(COSINE_SIM_100_MS1).toDouble();
        cosineSim100MS1Iso1 = dataMap.value(COSINE_SIM_100_MS1_ISO1).toDouble();
        cosineSim100MS1Iso2 = dataMap.value(COSINE_SIM_100_MS1_ISO2).toDouble();
        cosineSim45MS1 = dataMap.value(COSINE_SIM_45_MS1).toDouble();
        cosineSim20MS1 = dataMap.value(COSINE_SIM_20_MS1).toDouble();
        cosineSimSpectrum = dataMap.value(COSINE_SIM_SPECTRUM).toDouble();
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
        peakShapeRatio1 = dataMap.value(PEAK_SHAPE_RATIO_1).toDouble();
        peakShapeRatio2 = dataMap.value(PEAK_SHAPE_RATIO_2).toDouble();
        peakShapeRatio3 = dataMap.value(PEAK_SHAPE_RATIO_3).toDouble();

        b2Corr = dataMap.value(B2_CORR).toDouble();
        b3Corr = dataMap.value(B3_CORR).toDouble();
        b2b3CosineSimSum = dataMap.value(B2B3_COSINE_SIM_SUM).toDouble();

        shadowsCosineSimSum = dataMap.value(SHADOWS_COSINE_SIM_SUM).toDouble();
        cosineSimShadowsToAnchorVec = bytesArrayToQVector<double>(dataMap.value(COS_SIM_SUM_ANCH_SHADOW_V).toByteArray());
        shadowsIntensityRatioVec = bytesArrayToQVector<double>(dataMap.value(SHADOW_INTZ_RATIO_VEC).toByteArray());
        allignedMaxIndexesCount = dataMap.value(ALL_MAX_IND_CNT).toInt();

        ERR_RETURN
    }

};


#endif //PYTHIADIACPP_CANDIDATESCORES_H
