//
// Created by anichols on 10/20/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORES_H
#define PYTHIADIACPP_CANDIDATESCORES_H


#include "AlgorithmsFFLib_Exports.h"

#include "CSVReader.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "TargetDecoyCandidatePair.h"

class ALGORITHMSFFLIB_EXPORTS CandidateScores {

public:

    enum Features {
        CosineSimSum100 = 0,
        AllignedMaxIndexesCount,
        CosineSim100MS1,
        CosineSimSpectrumCubed,
        KlDivSpectrumCubeRoot,
        CosineSimSum45,
        CosineSimSum20,
        CosineSimSumTop6,
        CosineSimSumBottom6,
        TopBottomRatio,
        TopBottomRatioNorm,//10
        ChargeNorm,
        Mass,
        ScanTimeDelta,
        ScanTimeRange,
        ScanTimePd,
        ScanIonCount,
        MzNorm,
        KlDivSpectrum,
        CosineSimSpectrum,
        CosineSim45MS1,//20
        CosineSim20MS1,
        CosineSim100MS1PreMono,
        CosineSim100MS1Iso1,
        CosineSim100MS1Iso2,
        PeptideLengthNorm,
        ScanTimePredicted,
        TheoFragmentCount,
        TotalIntensityLog,
        PeakShapeRatio1,
        PeakShapeRatio2,//30
        PeakShapeRatio3,
        ShadowsCosineSimSum,
        IRTPredicted,
        CosineSimToAnchor1,
        CosineSimToAnchor2,
        CosineSimToAnchor3,
        CosineSimToAnchor4,
        CosineSimToAnchor5,
        CosineSimToAnchor6,
        CosineSimToAnchor7,//40
        CosineSimToAnchor8,
        CosineSimToAnchor9,
        CosineSimToAnchor10,
        CosineSimToAnchor11,
        CosineSimToAnchor12,
        CosineSimShadowsToAnchor1,
        CosineSimShadowsToAnchor2,
        CosineSimShadowsToAnchor3,
        CosineSimShadowsToAnchor4,
        CosineSimShadowsToAnchor5,//50
        CosineSimShadowsToAnchor6,
        CosineSimShadowsToAnchor7,
        CosineSimShadowsToAnchor8,
        CosineSimShadowsToAnchor9,
        CosineSimShadowsToAnchor10,
        CosineSimShadowsToAnchor11,
        CosineSimShadowsToAnchor12,
        ShadowsIntensityRatio1,
        ShadowsIntensityRatio2,
        ShadowsIntensityRatio3,//60
        ShadowsIntensityRatio4,
        ShadowsIntensityRatio5,
        ShadowsIntensityRatio6,
        ShadowsIntensityRatio7,
        ShadowsIntensityRatio8,
        ShadowsIntensityRatio9,
        ShadowsIntensityRatio10,
        ShadowsIntensityRatio11,
        ShadowsIntensityRatio12,
        MzSearched1,//70
        MzSearched2,
        MzSearched3,
        MzSearched4,
        MzSearched5,
        MzSearched6,
        MzSearched7,
        MzSearched8,
        MzSearched9,
        MzSearched10,
        MzSearched11,//80
        MzSearched12,
        TheoIntensity1,
        TheoIntensity2,
        TheoIntensity3,
        TheoIntensity4,
        TheoIntensity5,
        TheoIntensity6,
        TheoIntensity7,
        TheoIntensity8,
        TheoIntensity9,//90
        TheoIntensity10,
        TheoIntensity11,
        TheoIntensity12,
        MzFoundMean1,
        MzFoundMean2,
        MzFoundMean3,
        MzFoundMean4,
        MzFoundMean5,
        MzFoundMean6,
        MzFoundMean7,//100
        MzFoundMean8,
        MzFoundMean9,
        MzFoundMean10,
        MzFoundMean11,
        MzFoundMean12,
        IntensityFoundMax1,
        IntensityFoundMax2,
        IntensityFoundMax3,
        IntensityFoundMax4,
        IntensityFoundMax5,//110
        IntensityFoundMax6,
        IntensityFoundMax7,
        IntensityFoundMax8,
        IntensityFoundMax9,
        IntensityFoundMax10,
        IntensityFoundMax11,
        IntensityFoundMax12,
        MzPeakLengthsNorm1,
        MzPeakLengthsNorm2,
        MzPeakLengthsNorm3,//120
        MzPeakLengthsNorm4,
        MzPeakLengthsNorm5,
        MzPeakLengthsNorm6,
        MzPeakLengthsNorm7,
        MzPeakLengthsNorm8,
        MzPeakLengthsNorm9,
        MzPeakLengthsNorm10,
        MzPeakLengthsNorm11,
        MzPeakLengthsNorm12,
        ColumnApexIndexRatiosToAnchor1,//130
        ColumnApexIndexRatiosToAnchor2,
        ColumnApexIndexRatiosToAnchor3,
        ColumnApexIndexRatiosToAnchor4,
        ColumnApexIndexRatiosToAnchor5,
        ColumnApexIndexRatiosToAnchor6,
        ColumnApexIndexRatiosToAnchor7,
        ColumnApexIndexRatiosToAnchor8,
        ColumnApexIndexRatiosToAnchor9,
        ColumnApexIndexRatiosToAnchor10,
        ColumnApexIndexRatiosToAnchor11,//140
        ColumnApexIndexRatiosToAnchor12,
        AminoAcidCountA,
        AminoAcidCountC,
        AminoAcidCountD,
        AminoAcidCountE,
        AminoAcidCountF,
        AminoAcidCountG,
        AminoAcidCountH,
        AminoAcidCountI,
        AminoAcidCountK,//150
        AminoAcidCountL,
        AminoAcidCountM,
        AminoAcidCountN,
        AminoAcidCountP,
        AminoAcidCountQ,
        AminoAcidCountR,
        AminoAcidCountS,
        AminoAcidCountT,
        AminoAcidCountV,
        AminoAcidCountW,//160
        AminoAcidCountY,
        MzFoundStDev1,
        MzFoundStDev2,
        MzFoundStDev3,
        MzFoundStDev4,
        MzFoundStDev5,
        MzFoundStDev6,
        MzFoundStDev7,
        MzFoundStDev8,
        MzFoundStDev9,//170
        MzFoundStDev10,
        MzFoundStDev11,
        MzFoundStDev12,
        MzAccuracy1,
        MzAccuracy2,
        MzAccuracy3,
        MzAccuracy4,
        MzAccuracy5,
        MzAccuracy6,
        MzAccuracy7,//180
        MzAccuracy8,
        MzAccuracy9,
        MzAccuracy10,
        MzAccuracy11,
        MzAccuracy12,
        FeaturesSize
    };

    CandidateScores() = default;
    ~CandidateScores() = default;

    TargetDecoyCandidatePair *targetDecoyCandidatePair = nullptr;
    QString targetKey;
    QString proteinGroup;
    bool isDecoy = false;
    ScanNumber scanNumber = -1;
    ScanTime scanTime = -1.0;
    ScanTime scanTimePredicted = -1.0;
    double classifierScore = -1.0;
    double discriminantScore = -1.0;
    double qValue = 1.0;
    double decoyRatio = -1.0;

    QVector<float> featuresArray;

    [[nodiscard]] QVector<float>* featuresArrayRef();

    void initFeaturesArray();

};


namespace CandidateScoresReaderRowNamespace {

    const QString COS_SIM_SUM_100 = QStringLiteral("CosineSimSum100");
    const QString ALL_MAX_IND_CNT = QStringLiteral("AllignedMaxIndexesCount");
    const QString COS_SIM_SUM_MS1_100 = QStringLiteral("CosineSim100MS1");
    const QString COS_SIM_SPEC_CUBED = QStringLiteral("CosineSimSpectrumCubed");
    const QString KL_DIV_SPEC_CUBE_RT = QStringLiteral("KlDivSpectrumCubeRoot");
    const QString COS_SIM_SUM_45 = QStringLiteral("CosineSimSum45");
    const QString COS_SIM_SUM_20 = QStringLiteral("CosineSimSum20");
    const QString COS_SIM_SUM_TOP_6 = QStringLiteral("CosineSimSumTop6");
    const QString COS_SIM_SUM_BOTTOM_6 = QStringLiteral("CosineSimSumBottom6");
    const QString TOP_BOTTOM_RATIO = QStringLiteral("TopBottomRatio");
    const QString TOP_BOTTOME_RATIO_NORM = QStringLiteral("TopBottomRatioNorm");
    const QString CHARGE_NORM = QStringLiteral("ChargeNorm");
    const QString MASS = QStringLiteral("Mass");
    const QString SCAN_TIME_DELTA = QStringLiteral("ScanTimeDelta");
    const QString SCAN_TIME_RANGE = QStringLiteral("ScanTimeRange");
    const QString SCAN_TIME_PD = QStringLiteral("ScanTimePd");
    const QString SCAN_ION_CNT = QStringLiteral("ScanIonCount");
    const QString MZ_NORM = QStringLiteral("MzNorm");
    const QString KL_DIV_SPEC = QStringLiteral("KlDivSpectrum");
    const QString COSINE_SIM_SPEC = QStringLiteral("CosineSimSpectrum");
    const QString COSINE_SIM_SUM_MS1_45 = QStringLiteral("CosineSim45MS1");
    const QString COSINE_SIM_SUM_MS1_20 = QStringLiteral("CosineSim20MS1");
    const QString COSINE_SIM_SUM_MS1_PRE_MONO = QStringLiteral("CosineSim100MS1PreMono");
    const QString COSINE_SIM_SUM_MS1_ISO_1 = QStringLiteral("CosineSim100MS1Iso1");
    const QString COSINE_SIM_SUM_MS1_ISO_2 = QStringLiteral("CosineSim100MS1Iso2");
    const QString PEP_LEN_NORM = QStringLiteral("PeptideLengthNorm");
    const QString SCAN_TIME_PRED = QStringLiteral("ScanTimePredicted");
    const QString THEO_FRAG_CNT = QStringLiteral("TheoFragmentCount");
    const QString TOT_INT_LOG = QStringLiteral("TotalIntensityLog");
    const QString PEAK_RATIO_1 = QStringLiteral("PeakShapeRatio1");
    const QString PEAK_RATIO_2 = QStringLiteral("PeakShapeRatio2");
    const QString PEAK_RATIO_3 = QStringLiteral("PeakShapeRatio3");
    const QString SHADOW_COSINE_SIM_SUM = QStringLiteral("ShadowsCosineSimSum");
    const QString IRT_PRED = QStringLiteral("IRTPredicted");
    const QString COS_SIM_ANCH_1 = QStringLiteral("CosineSimToAnchor1");
    const QString COS_SIM_ANCH_2 = QStringLiteral("CosineSimToAnchor2");
    const QString COS_SIM_ANCH_3 = QStringLiteral("CosineSimToAnchor3");
    const QString COS_SIM_ANCH_4 = QStringLiteral("CosineSimToAnchor4");
    const QString COS_SIM_ANCH_5 = QStringLiteral("CosineSimToAnchor5");
    const QString COS_SIM_ANCH_6 = QStringLiteral("CosineSimToAnchor6");
    const QString COS_SIM_ANCH_7 = QStringLiteral("CosineSimToAnchor7");
    const QString COS_SIM_ANCH_8 = QStringLiteral("CosineSimToAnchor8");
    const QString COS_SIM_ANCH_9 = QStringLiteral("CosineSimToAnchor9");
    const QString COS_SIM_ANCH_10 = QStringLiteral("CosineSimToAnchor10");
    const QString COS_SIM_ANCH_11 = QStringLiteral("CosineSimToAnchor11");
    const QString COS_SIM_ANCH_12 = QStringLiteral("CosineSimToAnchor12");
    const QString COS_SIM_ANCH_SHADOW_1 = QStringLiteral("CosineSimShadowsToAnchor1");
    const QString COS_SIM_ANCH_SHADOW_2 = QStringLiteral("CosineSimShadowsToAnchor2");
    const QString COS_SIM_ANCH_SHADOW_3 = QStringLiteral("CosineSimShadowsToAnchor3");
    const QString COS_SIM_ANCH_SHADOW_4 = QStringLiteral("CosineSimShadowsToAnchor4");
    const QString COS_SIM_ANCH_SHADOW_5 = QStringLiteral("CosineSimShadowsToAnchor5");
    const QString COS_SIM_ANCH_SHADOW_6 = QStringLiteral("CosineSimShadowsToAnchor6");
    const QString COS_SIM_ANCH_SHADOW_7 = QStringLiteral("CosineSimShadowsToAnchor7");
    const QString COS_SIM_ANCH_SHADOW_8 = QStringLiteral("CosineSimShadowsToAnchor8");
    const QString COS_SIM_ANCH_SHADOW_9 = QStringLiteral("CosineSimShadowsToAnchor9");
    const QString COS_SIM_ANCH_SHADOW_10 = QStringLiteral("CosineSimShadowsToAnchor10");
    const QString COS_SIM_ANCH_SHADOW_11 = QStringLiteral("CosineSimShadowsToAnchor11");
    const QString COS_SIM_ANCH_SHADOW_12 = QStringLiteral("CosineSimShadowsToAnchor12");
    const QString SHAD_INTS_RATIO_1 = QStringLiteral("ShadowsIntensityRatio1");
    const QString SHAD_INTS_RATIO_2 = QStringLiteral("ShadowsIntensityRatio2");
    const QString SHAD_INTS_RATIO_3 = QStringLiteral("ShadowsIntensityRatio3");
    const QString SHAD_INTS_RATIO_4 = QStringLiteral("ShadowsIntensityRatio4");
    const QString SHAD_INTS_RATIO_5 = QStringLiteral("ShadowsIntensityRatio5");
    const QString SHAD_INTS_RATIO_6 = QStringLiteral("ShadowsIntensityRatio6");
    const QString SHAD_INTS_RATIO_7 = QStringLiteral("ShadowsIntensityRatio7");
    const QString SHAD_INTS_RATIO_8 = QStringLiteral("ShadowsIntensityRatio8");
    const QString SHAD_INTS_RATIO_9 = QStringLiteral("ShadowsIntensityRatio9");
    const QString SHAD_INTS_RATIO_10 = QStringLiteral("ShadowsIntensityRatio10");
    const QString SHAD_INTS_RATIO_11 = QStringLiteral("ShadowsIntensityRatio11");
    const QString SHAD_INTS_RATIO_12 = QStringLiteral("ShadowsIntensityRatio12");
    const QString MZ_SEARCHED_1 = QStringLiteral("MzSearched1");
    const QString MZ_SEARCHED_2 = QStringLiteral("MzSearched2");
    const QString MZ_SEARCHED_3 = QStringLiteral("MzSearched3");
    const QString MZ_SEARCHED_4 = QStringLiteral("MzSearched4");
    const QString MZ_SEARCHED_5 = QStringLiteral("MzSearched5");
    const QString MZ_SEARCHED_6 = QStringLiteral("MzSearched6");
    const QString MZ_SEARCHED_7 = QStringLiteral("MzSearched7");
    const QString MZ_SEARCHED_8 = QStringLiteral("MzSearched8");
    const QString MZ_SEARCHED_9 = QStringLiteral("MzSearched9");
    const QString MZ_SEARCHED_10 = QStringLiteral("MzSearched10");
    const QString MZ_SEARCHED_11 = QStringLiteral("MzSearched11");
    const QString MZ_SEARCHED_12 = QStringLiteral("MzSearched12");
    const QString THEO_INTS_1 = QStringLiteral("TheoIntensity1");
    const QString THEO_INTS_2 = QStringLiteral("TheoIntensity2");
    const QString THEO_INTS_3 = QStringLiteral("TheoIntensity3");
    const QString THEO_INTS_4 = QStringLiteral("TheoIntensity4");
    const QString THEO_INTS_5 = QStringLiteral("TheoIntensity5");
    const QString THEO_INTS_6 = QStringLiteral("TheoIntensity6");
    const QString THEO_INTS_7 = QStringLiteral("TheoIntensity7");
    const QString THEO_INTS_8 = QStringLiteral("TheoIntensity8");
    const QString THEO_INTS_9 = QStringLiteral("TheoIntensity9");
    const QString THEO_INTS_10 = QStringLiteral("TheoIntensity10");
    const QString THEO_INTS_11 = QStringLiteral("TheoIntensity11");
    const QString THEO_INTS_12 = QStringLiteral("TheoIntensity12");
    const QString MZ_FND_MEAN_1 = QStringLiteral("MzFoundMean1");
    const QString MZ_FND_MEAN_2 = QStringLiteral("MzFoundMean2");
    const QString MZ_FND_MEAN_3 = QStringLiteral("MzFoundMean3");
    const QString MZ_FND_MEAN_4 = QStringLiteral("MzFoundMean4");
    const QString MZ_FND_MEAN_5 = QStringLiteral("MzFoundMean5");
    const QString MZ_FND_MEAN_6 = QStringLiteral("MzFoundMean6");
    const QString MZ_FND_MEAN_7 = QStringLiteral("MzFoundMean7");
    const QString MZ_FND_MEAN_8 = QStringLiteral("MzFoundMean8");
    const QString MZ_FND_MEAN_9 = QStringLiteral("MzFoundMean9");
    const QString MZ_FND_MEAN_10 = QStringLiteral("MzFoundMean10");
    const QString MZ_FND_MEAN_11 = QStringLiteral("MzFoundMean11");
    const QString MZ_FND_MEAN_12 = QStringLiteral("MzFoundMean12");
    const QString INTS_FND_MAX_1 = QStringLiteral("IntensityFoundMax1");
    const QString INTS_FND_MAX_2 = QStringLiteral("IntensityFoundMax2");
    const QString INTS_FND_MAX_3 = QStringLiteral("IntensityFoundMax3");
    const QString INTS_FND_MAX_4 = QStringLiteral("IntensityFoundMax4");
    const QString INTS_FND_MAX_5 = QStringLiteral("IntensityFoundMax5");
    const QString INTS_FND_MAX_6 = QStringLiteral("IntensityFoundMax6");
    const QString INTS_FND_MAX_7 = QStringLiteral("IntensityFoundMax7");
    const QString INTS_FND_MAX_8 = QStringLiteral("IntensityFoundMax8");
    const QString INTS_FND_MAX_9 = QStringLiteral("IntensityFoundMax9");
    const QString INTS_FND_MAX_10 = QStringLiteral("IntensityFoundMax10");
    const QString INTS_FND_MAX_11 = QStringLiteral("IntensityFoundMax11");
    const QString INTS_FND_MAX_12 = QStringLiteral("IntensityFoundMax12");
    const QString MZ_PK_LEN_NORM_1 = QStringLiteral("MzPeakLengthsNorm1");
    const QString MZ_PK_LEN_NORM_2 = QStringLiteral("MzPeakLengthsNorm2");
    const QString MZ_PK_LEN_NORM_3 = QStringLiteral("MzPeakLengthsNorm3");
    const QString MZ_PK_LEN_NORM_4 = QStringLiteral("MzPeakLengthsNorm4");
    const QString MZ_PK_LEN_NORM_5 = QStringLiteral("MzPeakLengthsNorm5");
    const QString MZ_PK_LEN_NORM_6 = QStringLiteral("MzPeakLengthsNorm6");
    const QString MZ_PK_LEN_NORM_7 = QStringLiteral("MzPeakLengthsNorm7");
    const QString MZ_PK_LEN_NORM_8 = QStringLiteral("MzPeakLengthsNorm8");
    const QString MZ_PK_LEN_NORM_9 = QStringLiteral("MzPeakLengthsNorm9");
    const QString MZ_PK_LEN_NORM_10 = QStringLiteral("MzPeakLengthsNorm10");
    const QString MZ_PK_LEN_NORM_11 = QStringLiteral("MzPeakLengthsNorm11");
    const QString MZ_PK_LEN_NORM_12 = QStringLiteral("MzPeakLengthsNorm12");
    const QString COL_APX_IND_RATIO_TO_ANCH_1 = QStringLiteral("ColumnApexIndexRatiosToAnchor1");
    const QString COL_APX_IND_RATIO_TO_ANCH_2 = QStringLiteral("ColumnApexIndexRatiosToAnchor2");
    const QString COL_APX_IND_RATIO_TO_ANCH_3 = QStringLiteral("ColumnApexIndexRatiosToAnchor3");
    const QString COL_APX_IND_RATIO_TO_ANCH_4 = QStringLiteral("ColumnApexIndexRatiosToAnchor4");
    const QString COL_APX_IND_RATIO_TO_ANCH_5 = QStringLiteral("ColumnApexIndexRatiosToAnchor5");
    const QString COL_APX_IND_RATIO_TO_ANCH_6 = QStringLiteral("ColumnApexIndexRatiosToAnchor6");
    const QString COL_APX_IND_RATIO_TO_ANCH_7 = QStringLiteral("ColumnApexIndexRatiosToAnchor7");
    const QString COL_APX_IND_RATIO_TO_ANCH_8 = QStringLiteral("ColumnApexIndexRatiosToAnchor8");
    const QString COL_APX_IND_RATIO_TO_ANCH_9 = QStringLiteral("ColumnApexIndexRatiosToAnchor9");
    const QString COL_APX_IND_RATIO_TO_ANCH_10 = QStringLiteral("ColumnApexIndexRatiosToAnchor10");
    const QString COL_APX_IND_RATIO_TO_ANCH_11 = QStringLiteral("ColumnApexIndexRatiosToAnchor11");
    const QString COL_APX_IND_RATIO_TO_ANCH_12 = QStringLiteral("ColumnApexIndexRatiosToAnchor12");
    const QString AA_A = QStringLiteral("AminoAcidCountA");
    const QString AA_C = QStringLiteral("AminoAcidCountC");
    const QString AA_D = QStringLiteral("AminoAcidCountD");
    const QString AA_E = QStringLiteral("AminoAcidCountE");
    const QString AA_F = QStringLiteral("AminoAcidCountF");
    const QString AA_G = QStringLiteral("AminoAcidCountG");
    const QString AA_H = QStringLiteral("AminoAcidCountH");
    const QString AA_I = QStringLiteral("AminoAcidCountI");
    const QString AA_K = QStringLiteral("AminoAcidCountK");
    const QString AA_L = QStringLiteral("AminoAcidCountL");
    const QString AA_M = QStringLiteral("AminoAcidCountM");
    const QString AA_N = QStringLiteral("AminoAcidCountN");
    const QString AA_P = QStringLiteral("AminoAcidCountP");
    const QString AA_Q = QStringLiteral("AminoAcidCountQ");
    const QString AA_R = QStringLiteral("AminoAcidCountR");
    const QString AA_S = QStringLiteral("AminoAcidCountS");
    const QString AA_T = QStringLiteral("AminoAcidCountT");
    const QString AA_V = QStringLiteral("AminoAcidCountV");
    const QString AA_W = QStringLiteral("AminoAcidCountW");
    const QString AA_Y = QStringLiteral("AminoAcidCountY");
    const QString MZ_FND_STDEV_1 = QStringLiteral("MzFoundStDev1");
    const QString MZ_FND_STDEV_2 = QStringLiteral("MzFoundStDev2");
    const QString MZ_FND_STDEV_3 = QStringLiteral("MzFoundStDev3");
    const QString MZ_FND_STDEV_4 = QStringLiteral("MzFoundStDev4");
    const QString MZ_FND_STDEV_5 = QStringLiteral("MzFoundStDev5");
    const QString MZ_FND_STDEV_6 = QStringLiteral("MzFoundStDev6");
    const QString MZ_FND_STDEV_7 = QStringLiteral("MzFoundStDev7");
    const QString MZ_FND_STDEV_8 = QStringLiteral("MzFoundStDev8");
    const QString MZ_FND_STDEV_9 = QStringLiteral("MzFoundStDev9");
    const QString MZ_FND_STDEV_10 = QStringLiteral("MzFoundStDev10");
    const QString MZ_FND_STDEV_11 = QStringLiteral("MzFoundStDev11");
    const QString MZ_FND_STDEV_12 = QStringLiteral("MzFoundStDev12");
    const QString MZ_ACC_1 = QStringLiteral("MzAccuracy1");
    const QString MZ_ACC_2 = QStringLiteral("MzAccuracy2");
    const QString MZ_ACC_3 = QStringLiteral("MzAccuracy3");
    const QString MZ_ACC_4 = QStringLiteral("MzAccuracy4");
    const QString MZ_ACC_5 = QStringLiteral("MzAccuracy5");
    const QString MZ_ACC_6 = QStringLiteral("MzAccuracy6");
    const QString MZ_ACC_7 = QStringLiteral("MzAccuracy7");
    const QString MZ_ACC_8 = QStringLiteral("MzAccuracy8");
    const QString MZ_ACC_9 = QStringLiteral("MzAccuracy9");
    const QString MZ_ACC_10 = QStringLiteral("MzAccuracy10");
    const QString MZ_ACC_11 = QStringLiteral("MzAccuracy11");
    const QString MZ_ACC_12 = QStringLiteral("MzAccuracy12");
    const QString TARG_KEY = QStringLiteral("TargetKey");
    const QString PEP_STR_W_MODS = QStringLiteral("PeptideStringWithMods");
    const QString PROT_GRP = QStringLiteral("ProteinGroup");
    const QString IS_DECOY = QStringLiteral("IsDecoy");
    const QString SCAN_NUM = QStringLiteral("ScanNumber");
    const QString SCAN_TIME = QStringLiteral("ScanTime");
    const QString CLASS_SCR = QStringLiteral("ClassifierScore");
    const QString DISC_SCR = QStringLiteral("DiscriminantScore");
    const QString Q_VAL = QStringLiteral("QValue");
    const QString DECOY_RATIO = QStringLiteral("DecoyRatio");
}//namespace

struct ALGORITHMSFFLIB_EXPORTS CandidateScoresReaderRow : public CSVReaderInputBase {

    float cosineSimSum100 = -1.0;
    float allignedMaxIndexesCount = -1.0;
    float cosineSim100MS1 = -1.0;
    float cosineSimSpectrumCubed = -1.0;
    float klDivSpectrumCubeRoot = -1.0;
    float cosineSimSum45 = -1.0;
    float cosineSimSum20 = -1.0;
    float cosineSimSumTop6 = -1.0;
    float cosineSimSumBottom6 = -1.0;
    float topBottomRatio = -1.0;
    float topBottomRatioNorm = -1.0;
    float chargeNorm = -1.0;
    float mass = -1.0;
    float scanTimeDelta = -1.0;
    float scanTimeRange = -1.0;
    float scanTimePd = -1.0;
    float scanIonCount = -1.0;
    float mzNorm = -1.0;
    float klDivSpectrum = -1.0;
    float cosineSimSpectrum = -1.0;
    float cosineSim45MS1 = -1.0;
    float cosineSim20MS1 = -1.0;
    float cosineSim100MS1PreMono = -1.0;
    float cosineSim100MS1Iso1 = -1.0;
    float cosineSim100MS1Iso2 = -1.0;
    float peptideLengthNorm = -1.0;
    float scanTimePredicted = -1.0;
    float theoFragmentCount = -1.0;
    float totalIntensityLog = -1.0;
    float peakShapeRatio1 = -1.0;
    float peakShapeRatio2 = -1.0;
    float peakShapeRatio3 = -1.0;
    float shadowsCosineSimSum = -1.0;
    float iRtPredicted = -1.0;
    float cosineSimToAnchor1 = -1.0;
    float cosineSimToAnchor2 = -1.0;
    float cosineSimToAnchor3 = -1.0;
    float cosineSimToAnchor4 = -1.0;
    float cosineSimToAnchor5 = -1.0;
    float cosineSimToAnchor6 = -1.0;
    float cosineSimToAnchor7 = -1.0;
    float cosineSimToAnchor8 = -1.0;
    float cosineSimToAnchor9 = -1.0;
    float cosineSimToAnchor10 = -1.0;
    float cosineSimToAnchor11 = -1.0;
    float cosineSimToAnchor12 = -1.0;
    float cosineSimShadowsToAnchor1 = -1.0;
    float cosineSimShadowsToAnchor2 = -1.0;
    float cosineSimShadowsToAnchor3 = -1.0;
    float cosineSimShadowsToAnchor4 = -1.0;
    float cosineSimShadowsToAnchor5 = -1.0;
    float cosineSimShadowsToAnchor6 = -1.0;
    float cosineSimShadowsToAnchor7 = -1.0;
    float cosineSimShadowsToAnchor8 = -1.0;
    float cosineSimShadowsToAnchor9 = -1.0;
    float cosineSimShadowsToAnchor10 = -1.0;
    float cosineSimShadowsToAnchor11 = -1.0;
    float cosineSimShadowsToAnchor12 = -1.0;
    float shadowsIntensityRatio1 = -1.0;
    float shadowsIntensityRatio2 = -1.0;
    float shadowsIntensityRatio3 = -1.0;
    float shadowsIntensityRatio4 = -1.0;
    float shadowsIntensityRatio5 = -1.0;
    float shadowsIntensityRatio6 = -1.0;
    float shadowsIntensityRatio7 = -1.0;
    float shadowsIntensityRatio8 = -1.0;
    float shadowsIntensityRatio9 = -1.0;
    float shadowsIntensityRatio10 = -1.0;
    float shadowsIntensityRatio11 = -1.0;
    float shadowsIntensityRatio12 = -1.0;
    float mzSearched1 = -1.0;
    float mzSearched2 = -1.0;
    float mzSearched3 = -1.0;
    float mzSearched4 = -1.0;
    float mzSearched5 = -1.0;
    float mzSearched6 = -1.0;
    float mzSearched7 = -1.0;
    float mzSearched8 = -1.0;
    float mzSearched9 = -1.0;
    float mzSearched10 = -1.0;
    float mzSearched11 = -1.0;
    float mzSearched12 = -1.0;
    float theoIntensity1 = -1.0;
    float theoIntensity2 = -1.0;
    float theoIntensity3 = -1.0;
    float theoIntensity4 = -1.0;
    float theoIntensity5 = -1.0;
    float theoIntensity6 = -1.0;
    float theoIntensity7 = -1.0;
    float theoIntensity8 = -1.0;
    float theoIntensity9 = -1.0;
    float theoIntensity10 = -1.0;
    float theoIntensity11 = -1.0;
    float theoIntensity12 = -1.0;
    float mzFoundMean1 = -1.0;
    float mzFoundMean2 = -1.0;
    float mzFoundMean3 = -1.0;
    float mzFoundMean4 = -1.0;
    float mzFoundMean5 = -1.0;
    float mzFoundMean6 = -1.0;
    float mzFoundMean7 = -1.0;
    float mzFoundMean8 = -1.0;
    float mzFoundMean9 = -1.0;
    float mzFoundMean10 = -1.0;
    float mzFoundMean11 = -1.0;
    float mzFoundMean12 = -1.0;
    float intensityFoundMax1 = -1.0;
    float intensityFoundMax2 = -1.0;
    float intensityFoundMax3 = -1.0;
    float intensityFoundMax4 = -1.0;
    float intensityFoundMax5 = -1.0;
    float intensityFoundMax6 = -1.0;
    float intensityFoundMax7 = -1.0;
    float intensityFoundMax8 = -1.0;
    float intensityFoundMax9 = -1.0;
    float intensityFoundMax10 = -1.0;
    float intensityFoundMax11 = -1.0;
    float intensityFoundMax12 = -1.0;
    float mzPeakLengthsNorm1 = -1.0;
    float mzPeakLengthsNorm2 = -1.0;
    float mzPeakLengthsNorm3 = -1.0;
    float mzPeakLengthsNorm4 = -1.0;
    float mzPeakLengthsNorm5 = -1.0;
    float mzPeakLengthsNorm6 = -1.0;
    float mzPeakLengthsNorm7 = -1.0;
    float mzPeakLengthsNorm8 = -1.0;
    float mzPeakLengthsNorm9 = -1.0;
    float mzPeakLengthsNorm10 = -1.0;
    float mzPeakLengthsNorm11 = -1.0;
    float mzPeakLengthsNorm12 = -1.0;
    float columnApexIndexRatiosToAnchor1 = -1.0;
    float columnApexIndexRatiosToAnchor2 = -1.0;
    float columnApexIndexRatiosToAnchor3 = -1.0;
    float columnApexIndexRatiosToAnchor4 = -1.0;
    float columnApexIndexRatiosToAnchor5 = -1.0;
    float columnApexIndexRatiosToAnchor6 = -1.0;
    float columnApexIndexRatiosToAnchor7 = -1.0;
    float columnApexIndexRatiosToAnchor8 = -1.0;
    float columnApexIndexRatiosToAnchor9 = -1.0;
    float columnApexIndexRatiosToAnchor10 = -1.0;
    float columnApexIndexRatiosToAnchor11 = -1.0;
    float columnApexIndexRatiosToAnchor12 = -1.0;
    float aminoAcidCountA = -1.0;
    float aminoAcidCountC = -1.0;
    float aminoAcidCountD = -1.0;
    float aminoAcidCountE = -1.0;
    float aminoAcidCountF = -1.0;
    float aminoAcidCountG = -1.0;
    float aminoAcidCountH = -1.0;
    float aminoAcidCountI = -1.0;
    float aminoAcidCountK = -1.0;
    float aminoAcidCountL = -1.0;
    float aminoAcidCountM = -1.0;
    float aminoAcidCountN = -1.0;
    float aminoAcidCountP = -1.0;
    float aminoAcidCountQ = -1.0;
    float aminoAcidCountR = -1.0;
    float aminoAcidCountS = -1.0;
    float aminoAcidCountT = -1.0;
    float aminoAcidCountV = -1.0;
    float aminoAcidCountW = -1.0;
    float aminoAcidCountY = -1.0;
    float mzFoundStDev1 = -1.0;
    float mzFoundStDev2 = -1.0;
    float mzFoundStDev3 = -1.0;
    float mzFoundStDev4 = -1.0;
    float mzFoundStDev5 = -1.0;
    float mzFoundStDev6 = -1.0;
    float mzFoundStDev7 = -1.0;
    float mzFoundStDev8 = -1.0;
    float mzFoundStDev9 = -1.0;
    float mzFoundStDev10 = -1.0;
    float mzFoundStDev11 = -1.0;
    float mzFoundStDev12 = -1.0;
    float mzAccuracy1 = -1.0;
    float mzAccuracy2 = -1.0;
    float mzAccuracy3 = -1.0;
    float mzAccuracy4 = -1.0;
    float mzAccuracy5 = -1.0;
    float mzAccuracy6 = -1.0;
    float mzAccuracy7 = -1.0;
    float mzAccuracy8 = -1.0;
    float mzAccuracy9 = -1.0;
    float mzAccuracy10 = -1.0;
    float mzAccuracy11 = -1.0;
    float mzAccuracy12 = -1.0;
    QString targetKey;
    QString peptideStringWithMods;
    QString proteinGroup;
    bool isDecoy = false;
    ScanNumber scanNumber = -1;
    ScanTime scanTime = -1.0;
    double classifierScore = -1.0;
    double discriminantScore = -1.0;
    double qValue = 1.0;
    double decoyRatio = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace CandidateScoresReaderRowNamespace;

        return {
                {COS_SIM_SUM_100, QVariant(cosineSimSum100)},
                {ALL_MAX_IND_CNT, QVariant(allignedMaxIndexesCount)},
                {COS_SIM_SUM_MS1_100, QVariant(cosineSim100MS1)},
                {COS_SIM_SPEC_CUBED, QVariant(cosineSimSpectrumCubed)},
                {KL_DIV_SPEC_CUBE_RT, QVariant(klDivSpectrumCubeRoot)},
                {COS_SIM_SUM_45, QVariant(cosineSimSum45)},
                {COS_SIM_SUM_20, QVariant(cosineSimSum20)},
                {COS_SIM_SUM_TOP_6, QVariant(cosineSimSumTop6)},
                {COS_SIM_SUM_BOTTOM_6, QVariant(cosineSimSumBottom6)},
                {TOP_BOTTOM_RATIO, QVariant(topBottomRatio)},
                {TOP_BOTTOME_RATIO_NORM, QVariant(topBottomRatioNorm)},
                {CHARGE_NORM, QVariant(chargeNorm)},
                {MASS, QVariant(mass)},
                {SCAN_TIME_DELTA, QVariant(scanTimeDelta)},
                {SCAN_TIME_RANGE, QVariant(scanTimeRange)},
                {SCAN_TIME_PD, QVariant(scanTimePd)},
                {SCAN_ION_CNT, QVariant(scanIonCount)},
                {MZ_NORM, QVariant(mzNorm)},
                {KL_DIV_SPEC, QVariant(klDivSpectrum)},
                {COSINE_SIM_SPEC, QVariant(cosineSimSpectrum)},
                {COSINE_SIM_SUM_MS1_45, QVariant(cosineSim45MS1)},
                {COSINE_SIM_SUM_MS1_20, QVariant(cosineSim20MS1)},
                {COSINE_SIM_SUM_MS1_PRE_MONO, QVariant(cosineSim100MS1PreMono)},
                {COSINE_SIM_SUM_MS1_ISO_1, QVariant(cosineSim100MS1Iso1)},
                {COSINE_SIM_SUM_MS1_ISO_2, QVariant(cosineSim100MS1Iso2)},
                {PEP_LEN_NORM, QVariant(peptideLengthNorm)},
                {SCAN_TIME_PRED, QVariant(scanTimePredicted)},
                {THEO_FRAG_CNT, QVariant(theoFragmentCount)},
                {TOT_INT_LOG, QVariant(totalIntensityLog)},
                {PEAK_RATIO_1, QVariant(peakShapeRatio1)},
                {PEAK_RATIO_2, QVariant(peakShapeRatio2)},
                {PEAK_RATIO_3, QVariant(peakShapeRatio3)},
                {SHADOW_COSINE_SIM_SUM, QVariant(shadowsCosineSimSum)},
                {IRT_PRED, QVariant(iRtPredicted)},
                {COS_SIM_ANCH_1, QVariant(cosineSimToAnchor1)},
                {COS_SIM_ANCH_2, QVariant(cosineSimToAnchor2)},
                {COS_SIM_ANCH_3, QVariant(cosineSimToAnchor3)},
                {COS_SIM_ANCH_4, QVariant(cosineSimToAnchor4)},
                {COS_SIM_ANCH_5, QVariant(cosineSimToAnchor5)},
                {COS_SIM_ANCH_6, QVariant(cosineSimToAnchor6)},
                {COS_SIM_ANCH_7, QVariant(cosineSimToAnchor7)},
                {COS_SIM_ANCH_8, QVariant(cosineSimToAnchor8)},
                {COS_SIM_ANCH_9, QVariant(cosineSimToAnchor9)},
                {COS_SIM_ANCH_10, QVariant(cosineSimToAnchor10)},
                {COS_SIM_ANCH_11, QVariant(cosineSimToAnchor11)},
                {COS_SIM_ANCH_12, QVariant(cosineSimToAnchor12)},
                {COS_SIM_ANCH_SHADOW_1, QVariant(cosineSimShadowsToAnchor1)},
                {COS_SIM_ANCH_SHADOW_2, QVariant(cosineSimShadowsToAnchor2)},
                {COS_SIM_ANCH_SHADOW_3, QVariant(cosineSimShadowsToAnchor3)},
                {COS_SIM_ANCH_SHADOW_4, QVariant(cosineSimShadowsToAnchor4)},
                {COS_SIM_ANCH_SHADOW_5, QVariant(cosineSimShadowsToAnchor5)},
                {COS_SIM_ANCH_SHADOW_6, QVariant(cosineSimShadowsToAnchor6)},
                {COS_SIM_ANCH_SHADOW_7, QVariant(cosineSimShadowsToAnchor7)},
                {COS_SIM_ANCH_SHADOW_8, QVariant(cosineSimShadowsToAnchor8)},
                {COS_SIM_ANCH_SHADOW_9, QVariant(cosineSimShadowsToAnchor9)},
                {COS_SIM_ANCH_SHADOW_10, QVariant(cosineSimShadowsToAnchor10)},
                {COS_SIM_ANCH_SHADOW_11, QVariant(cosineSimShadowsToAnchor11)},
                {COS_SIM_ANCH_SHADOW_12, QVariant(cosineSimShadowsToAnchor12)},
                {SHAD_INTS_RATIO_1, QVariant(shadowsIntensityRatio1)},
                {SHAD_INTS_RATIO_2, QVariant(shadowsIntensityRatio2)},
                {SHAD_INTS_RATIO_3, QVariant(shadowsIntensityRatio3)},
                {SHAD_INTS_RATIO_4, QVariant(shadowsIntensityRatio4)},
                {SHAD_INTS_RATIO_5, QVariant(shadowsIntensityRatio5)},
                {SHAD_INTS_RATIO_6, QVariant(shadowsIntensityRatio6)},
                {SHAD_INTS_RATIO_7, QVariant(shadowsIntensityRatio7)},
                {SHAD_INTS_RATIO_8, QVariant(shadowsIntensityRatio8)},
                {SHAD_INTS_RATIO_9, QVariant(shadowsIntensityRatio9)},
                {SHAD_INTS_RATIO_10, QVariant(shadowsIntensityRatio10)},
                {SHAD_INTS_RATIO_11, QVariant(shadowsIntensityRatio11)},
                {SHAD_INTS_RATIO_12, QVariant(shadowsIntensityRatio12)},
                {MZ_SEARCHED_1, QVariant(mzSearched1)},
                {MZ_SEARCHED_2, QVariant(mzSearched2)},
                {MZ_SEARCHED_3, QVariant(mzSearched3)},
                {MZ_SEARCHED_4, QVariant(mzSearched4)},
                {MZ_SEARCHED_5, QVariant(mzSearched5)},
                {MZ_SEARCHED_6, QVariant(mzSearched6)},
                {MZ_SEARCHED_7, QVariant(mzSearched7)},
                {MZ_SEARCHED_8, QVariant(mzSearched8)},
                {MZ_SEARCHED_9, QVariant(mzSearched9)},
                {MZ_SEARCHED_10, QVariant(mzSearched10)},
                {MZ_SEARCHED_11, QVariant(mzSearched11)},
                {MZ_SEARCHED_12, QVariant(mzSearched12)},
                {THEO_INTS_1, QVariant(theoIntensity1)},
                {THEO_INTS_2, QVariant(theoIntensity2)},
                {THEO_INTS_3, QVariant(theoIntensity3)},
                {THEO_INTS_4, QVariant(theoIntensity4)},
                {THEO_INTS_5, QVariant(theoIntensity5)},
                {THEO_INTS_6, QVariant(theoIntensity6)},
                {THEO_INTS_7, QVariant(theoIntensity7)},
                {THEO_INTS_8, QVariant(theoIntensity8)},
                {THEO_INTS_9, QVariant(theoIntensity9)},
                {THEO_INTS_10, QVariant(theoIntensity10)},
                {THEO_INTS_11, QVariant(theoIntensity11)},
                {THEO_INTS_12, QVariant(theoIntensity12)},
                {MZ_FND_MEAN_1, QVariant(mzFoundMean1)},
                {MZ_FND_MEAN_2, QVariant(mzFoundMean2)},
                {MZ_FND_MEAN_3, QVariant(mzFoundMean3)},
                {MZ_FND_MEAN_4, QVariant(mzFoundMean4)},
                {MZ_FND_MEAN_5, QVariant(mzFoundMean5)},
                {MZ_FND_MEAN_6, QVariant(mzFoundMean6)},
                {MZ_FND_MEAN_7, QVariant(mzFoundMean7)},
                {MZ_FND_MEAN_8, QVariant(mzFoundMean8)},
                {MZ_FND_MEAN_9, QVariant(mzFoundMean9)},
                {MZ_FND_MEAN_10, QVariant(mzFoundMean10)},
                {MZ_FND_MEAN_11, QVariant(mzFoundMean11)},
                {MZ_FND_MEAN_12, QVariant(mzFoundMean12)},
                {INTS_FND_MAX_1, QVariant(intensityFoundMax1)},
                {INTS_FND_MAX_2, QVariant(intensityFoundMax2)},
                {INTS_FND_MAX_3, QVariant(intensityFoundMax3)},
                {INTS_FND_MAX_4, QVariant(intensityFoundMax4)},
                {INTS_FND_MAX_5, QVariant(intensityFoundMax5)},
                {INTS_FND_MAX_6, QVariant(intensityFoundMax6)},
                {INTS_FND_MAX_7, QVariant(intensityFoundMax7)},
                {INTS_FND_MAX_8, QVariant(intensityFoundMax8)},
                {INTS_FND_MAX_9, QVariant(intensityFoundMax9)},
                {INTS_FND_MAX_10, QVariant(intensityFoundMax10)},
                {INTS_FND_MAX_11, QVariant(intensityFoundMax11)},
                {INTS_FND_MAX_12, QVariant(intensityFoundMax12)},
                {MZ_PK_LEN_NORM_1, QVariant(mzPeakLengthsNorm1)},
                {MZ_PK_LEN_NORM_2, QVariant(mzPeakLengthsNorm2)},
                {MZ_PK_LEN_NORM_3, QVariant(mzPeakLengthsNorm3)},
                {MZ_PK_LEN_NORM_4, QVariant(mzPeakLengthsNorm4)},
                {MZ_PK_LEN_NORM_5, QVariant(mzPeakLengthsNorm5)},
                {MZ_PK_LEN_NORM_6, QVariant(mzPeakLengthsNorm6)},
                {MZ_PK_LEN_NORM_7, QVariant(mzPeakLengthsNorm7)},
                {MZ_PK_LEN_NORM_8, QVariant(mzPeakLengthsNorm8)},
                {MZ_PK_LEN_NORM_9, QVariant(mzPeakLengthsNorm9)},
                {MZ_PK_LEN_NORM_10, QVariant(mzPeakLengthsNorm10)},
                {MZ_PK_LEN_NORM_11, QVariant(mzPeakLengthsNorm11)},
                {MZ_PK_LEN_NORM_12, QVariant(mzPeakLengthsNorm12)},
                {COL_APX_IND_RATIO_TO_ANCH_1, QVariant(columnApexIndexRatiosToAnchor1)},
                {COL_APX_IND_RATIO_TO_ANCH_2, QVariant(columnApexIndexRatiosToAnchor2)},
                {COL_APX_IND_RATIO_TO_ANCH_3, QVariant(columnApexIndexRatiosToAnchor3)},
                {COL_APX_IND_RATIO_TO_ANCH_4, QVariant(columnApexIndexRatiosToAnchor4)},
                {COL_APX_IND_RATIO_TO_ANCH_5, QVariant(columnApexIndexRatiosToAnchor5)},
                {COL_APX_IND_RATIO_TO_ANCH_6, QVariant(columnApexIndexRatiosToAnchor6)},
                {COL_APX_IND_RATIO_TO_ANCH_7, QVariant(columnApexIndexRatiosToAnchor7)},
                {COL_APX_IND_RATIO_TO_ANCH_8, QVariant(columnApexIndexRatiosToAnchor8)},
                {COL_APX_IND_RATIO_TO_ANCH_9, QVariant(columnApexIndexRatiosToAnchor9)},
                {COL_APX_IND_RATIO_TO_ANCH_10, QVariant(columnApexIndexRatiosToAnchor10)},
                {COL_APX_IND_RATIO_TO_ANCH_11, QVariant(columnApexIndexRatiosToAnchor11)},
                {COL_APX_IND_RATIO_TO_ANCH_12, QVariant(columnApexIndexRatiosToAnchor12)},
                {AA_A, QVariant(aminoAcidCountA)},
                {AA_C, QVariant(aminoAcidCountC)},
                {AA_D, QVariant(aminoAcidCountD)},
                {AA_E, QVariant(aminoAcidCountE)},
                {AA_F, QVariant(aminoAcidCountF)},
                {AA_G, QVariant(aminoAcidCountG)},
                {AA_H, QVariant(aminoAcidCountH)},
                {AA_I, QVariant(aminoAcidCountI)},
                {AA_K, QVariant(aminoAcidCountK)},
                {AA_L, QVariant(aminoAcidCountL)},
                {AA_M, QVariant(aminoAcidCountM)},
                {AA_N, QVariant(aminoAcidCountN)},
                {AA_P, QVariant(aminoAcidCountP)},
                {AA_Q, QVariant(aminoAcidCountQ)},
                {AA_R, QVariant(aminoAcidCountR)},
                {AA_S, QVariant(aminoAcidCountS)},
                {AA_T, QVariant(aminoAcidCountT)},
                {AA_V, QVariant(aminoAcidCountV)},
                {AA_W, QVariant(aminoAcidCountW)},
                {AA_Y, QVariant(aminoAcidCountY)},
                {MZ_FND_STDEV_1, QVariant(mzFoundStDev1)},
                {MZ_FND_STDEV_2, QVariant(mzFoundStDev2)},
                {MZ_FND_STDEV_3, QVariant(mzFoundStDev3)},
                {MZ_FND_STDEV_4, QVariant(mzFoundStDev4)},
                {MZ_FND_STDEV_5, QVariant(mzFoundStDev5)},
                {MZ_FND_STDEV_6, QVariant(mzFoundStDev6)},
                {MZ_FND_STDEV_7, QVariant(mzFoundStDev7)},
                {MZ_FND_STDEV_8, QVariant(mzFoundStDev8)},
                {MZ_FND_STDEV_9, QVariant(mzFoundStDev9)},
                {MZ_FND_STDEV_10, QVariant(mzFoundStDev10)},
                {MZ_FND_STDEV_11, QVariant(mzFoundStDev11)},
                {MZ_FND_STDEV_12, QVariant(mzFoundStDev12)},
                {MZ_ACC_1, QVariant(mzAccuracy1)},
                {MZ_ACC_2, QVariant(mzAccuracy2)},
                {MZ_ACC_3, QVariant(mzAccuracy3)},
                {MZ_ACC_4, QVariant(mzAccuracy4)},
                {MZ_ACC_5, QVariant(mzAccuracy5)},
                {MZ_ACC_6, QVariant(mzAccuracy6)},
                {MZ_ACC_7, QVariant(mzAccuracy7)},
                {MZ_ACC_8, QVariant(mzAccuracy8)},
                {MZ_ACC_9, QVariant(mzAccuracy9)},
                {MZ_ACC_10, QVariant(mzAccuracy10)},
                {MZ_ACC_11, QVariant(mzAccuracy11)},
                {MZ_ACC_12, QVariant(mzAccuracy12)},
                {TARG_KEY, QVariant(targetKey)},
                {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
                {PROT_GRP, QVariant(proteinGroup)},
                {IS_DECOY, QVariant(isDecoy)},
                {SCAN_NUM, QVariant(scanNumber)},
                {SCAN_TIME, QVariant(scanTime)},
                {CLASS_SCR, QVariant(classifierScore)},
                {DISC_SCR, QVariant(discriminantScore)},
                {Q_VAL, QVariant(qValue)},
                {DECOY_RATIO, QVariant(decoyRatio)}
        };
    }

    static CandidateScoresReaderRow buildCandidateScoresReaderRow(const CandidateScores* candidateScores) {

        CandidateScoresReaderRow row;

        row.cosineSimSum100 = candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100],
        row.allignedMaxIndexesCount = candidateScores->featuresArray[CandidateScores::Features::AllignedMaxIndexesCount],
        row.cosineSim100MS1 = candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1],
        row.cosineSimSpectrumCubed = candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumCubed],
        row.klDivSpectrumCubeRoot = candidateScores->featuresArray[CandidateScores::Features::KlDivSpectrumCubeRoot],
        row.cosineSimSum45 = candidateScores->featuresArray[CandidateScores::Features::CosineSimSum45],
        row.cosineSimSum20 = candidateScores->featuresArray[CandidateScores::Features::CosineSimSum20],
        row.cosineSimSumTop6 = candidateScores->featuresArray[CandidateScores::Features::CosineSimSumTop6],
        row.cosineSimSumBottom6 = candidateScores->featuresArray[CandidateScores::Features::CosineSimSumBottom6],
        row.topBottomRatio = candidateScores->featuresArray[CandidateScores::Features::TopBottomRatio],
        row.topBottomRatioNorm = candidateScores->featuresArray[CandidateScores::Features::TopBottomRatioNorm],
        row.chargeNorm = candidateScores->featuresArray[CandidateScores::Features::ChargeNorm],
        row.mass = candidateScores->featuresArray[CandidateScores::Features::Mass],
        row.scanTimeDelta = candidateScores->featuresArray[CandidateScores::Features::ScanTimeDelta],
        row.scanTimeRange = candidateScores->featuresArray[CandidateScores::Features::ScanTimeRange],
        row.scanTimePd = candidateScores->featuresArray[CandidateScores::Features::ScanTimePd],
        row.scanIonCount = candidateScores->featuresArray[CandidateScores::Features::ScanIonCount],
        row.mzNorm = candidateScores->featuresArray[CandidateScores::Features::MzNorm],
        row.klDivSpectrum = candidateScores->featuresArray[CandidateScores::Features::KlDivSpectrum],
        row.cosineSimSpectrum = candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrum],
        row.cosineSim45MS1 = candidateScores->featuresArray[CandidateScores::Features::CosineSim45MS1],
        row.cosineSim20MS1 = candidateScores->featuresArray[CandidateScores::Features::CosineSim20MS1],
        row.cosineSim100MS1PreMono = candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1PreMono],
        row.cosineSim100MS1Iso1 = candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso1],
        row.cosineSim100MS1Iso2 = candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso2],
        row.peptideLengthNorm = candidateScores->featuresArray[CandidateScores::Features::PeptideLengthNorm],
        row.scanTimePredicted = candidateScores->featuresArray[CandidateScores::Features::ScanTimePredicted],
        row.theoFragmentCount = candidateScores->featuresArray[CandidateScores::Features::TheoFragmentCount],
        row.totalIntensityLog = candidateScores->featuresArray[CandidateScores::Features::TotalIntensityLog],
        row.peakShapeRatio1 = candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio1],
        row.peakShapeRatio2 = candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio2],
        row.peakShapeRatio3 = candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio3],
        row.shadowsCosineSimSum = candidateScores->featuresArray[CandidateScores::Features::ShadowsCosineSimSum],
        row.iRtPredicted = candidateScores->featuresArray[CandidateScores::Features::IRTPredicted],
        row.cosineSimToAnchor1 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor1],
        row.cosineSimToAnchor2 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor2],
        row.cosineSimToAnchor3 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor3],
        row.cosineSimToAnchor4 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor4],
        row.cosineSimToAnchor5 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor5],
        row.cosineSimToAnchor6 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor6],
        row.cosineSimToAnchor7 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor7],
        row.cosineSimToAnchor8 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor8],
        row.cosineSimToAnchor9 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor9],
        row.cosineSimToAnchor10 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor10],
        row.cosineSimToAnchor11 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor11],
        row.cosineSimToAnchor12 = candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor12],
        row.cosineSimShadowsToAnchor1 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor1],
        row.cosineSimShadowsToAnchor2 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor2],
        row.cosineSimShadowsToAnchor3 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor3],
        row.cosineSimShadowsToAnchor4 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor4],
        row.cosineSimShadowsToAnchor5 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor5],
        row.cosineSimShadowsToAnchor6 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor6],
        row.cosineSimShadowsToAnchor7 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor7],
        row.cosineSimShadowsToAnchor8 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor8],
        row.cosineSimShadowsToAnchor9 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor9],
        row.cosineSimShadowsToAnchor10 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor10],
        row.cosineSimShadowsToAnchor11 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor11],
        row.cosineSimShadowsToAnchor12 = candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor12],
        row.shadowsIntensityRatio1 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio1],
        row.shadowsIntensityRatio2 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio2],
        row.shadowsIntensityRatio3 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio3],
        row.shadowsIntensityRatio4 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio4],
        row.shadowsIntensityRatio5 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio5],
        row.shadowsIntensityRatio6 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio6],
        row.shadowsIntensityRatio7 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio7],
        row.shadowsIntensityRatio8 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio8],
        row.shadowsIntensityRatio9 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio9],
        row.shadowsIntensityRatio10 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio10],
        row.shadowsIntensityRatio11 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio11],
        row.shadowsIntensityRatio12 = candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio12],
        row.mzSearched1 = candidateScores->featuresArray[CandidateScores::Features::MzSearched1],
        row.mzSearched2 = candidateScores->featuresArray[CandidateScores::Features::MzSearched2],
        row.mzSearched3 = candidateScores->featuresArray[CandidateScores::Features::MzSearched3],
        row.mzSearched4 = candidateScores->featuresArray[CandidateScores::Features::MzSearched4],
        row.mzSearched5 = candidateScores->featuresArray[CandidateScores::Features::MzSearched5],
        row.mzSearched6 = candidateScores->featuresArray[CandidateScores::Features::MzSearched6],
        row.mzSearched7 = candidateScores->featuresArray[CandidateScores::Features::MzSearched7],
        row.mzSearched8 = candidateScores->featuresArray[CandidateScores::Features::MzSearched8],
        row.mzSearched9 = candidateScores->featuresArray[CandidateScores::Features::MzSearched9],
        row.mzSearched10 = candidateScores->featuresArray[CandidateScores::Features::MzSearched10],
        row.mzSearched11 = candidateScores->featuresArray[CandidateScores::Features::MzSearched11],
        row.mzSearched12 = candidateScores->featuresArray[CandidateScores::Features::MzSearched12],
        row.theoIntensity1 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity1],
        row.theoIntensity2 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity2],
        row.theoIntensity3 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity3],
        row.theoIntensity4 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity4],
        row.theoIntensity5 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity5],
        row.theoIntensity6 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity6],
        row.theoIntensity7 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity7],
        row.theoIntensity8 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity8],
        row.theoIntensity9 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity9],
        row.theoIntensity10 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity10],
        row.theoIntensity11 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity11],
        row.theoIntensity12 = candidateScores->featuresArray[CandidateScores::Features::TheoIntensity12],
        row.mzFoundMean1 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean1],
        row.mzFoundMean2 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean2],
        row.mzFoundMean3 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean3],
        row.mzFoundMean4 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean4],
        row.mzFoundMean5 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean5],
        row.mzFoundMean6 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean6],
        row.mzFoundMean7 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean7],
        row.mzFoundMean8 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean8],
        row.mzFoundMean9 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean9],
        row.mzFoundMean10 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean10],
        row.mzFoundMean11 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean11],
        row.mzFoundMean12 = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean12],
        row.intensityFoundMax1 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax1],
        row.intensityFoundMax2 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax2],
        row.intensityFoundMax3 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax3],
        row.intensityFoundMax4 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax4],
        row.intensityFoundMax5 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax5],
        row.intensityFoundMax6 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax6],
        row.intensityFoundMax7 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax7],
        row.intensityFoundMax8 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax8],
        row.intensityFoundMax9 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax9],
        row.intensityFoundMax10 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax10],
        row.intensityFoundMax11 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax11],
        row.intensityFoundMax12 = candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax12],
        row.mzPeakLengthsNorm1 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm1],
        row.mzPeakLengthsNorm2 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm2],
        row.mzPeakLengthsNorm3 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm3],
        row.mzPeakLengthsNorm4 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm4],
        row.mzPeakLengthsNorm5 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm5],
        row.mzPeakLengthsNorm6 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm6],
        row.mzPeakLengthsNorm7 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm7],
        row.mzPeakLengthsNorm8 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm8],
        row.mzPeakLengthsNorm9 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm9],
        row.mzPeakLengthsNorm10 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm10],
        row.mzPeakLengthsNorm11 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm11],
        row.mzPeakLengthsNorm12 = candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm12],
        row.columnApexIndexRatiosToAnchor1 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor1],
        row.columnApexIndexRatiosToAnchor2 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor2],
        row.columnApexIndexRatiosToAnchor3 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor3],
        row.columnApexIndexRatiosToAnchor4 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor4],
        row.columnApexIndexRatiosToAnchor5 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor5],
        row.columnApexIndexRatiosToAnchor6 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor6],
        row.columnApexIndexRatiosToAnchor7 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor7],
        row.columnApexIndexRatiosToAnchor8 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor8],
        row.columnApexIndexRatiosToAnchor9 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor9],
        row.columnApexIndexRatiosToAnchor10 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor10],
        row.columnApexIndexRatiosToAnchor11 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor11],
        row.columnApexIndexRatiosToAnchor12 = candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor12],
        row.aminoAcidCountA = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountA],
        row.aminoAcidCountC = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountC],
        row.aminoAcidCountD = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountD],
        row.aminoAcidCountE = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountE],
        row.aminoAcidCountF = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountF],
        row.aminoAcidCountG = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountG],
        row.aminoAcidCountH = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountH],
        row.aminoAcidCountI = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountI],
        row.aminoAcidCountK = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountK],
        row.aminoAcidCountL = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountL],
        row.aminoAcidCountM = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountM],
        row.aminoAcidCountN = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountN],
        row.aminoAcidCountP = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountP],
        row.aminoAcidCountQ = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountQ],
        row.aminoAcidCountR = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountR],
        row.aminoAcidCountS = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountS],
        row.aminoAcidCountT = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountT],
        row.aminoAcidCountV = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountV],
        row.aminoAcidCountW = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountW],
        row.aminoAcidCountY = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountY],
        row.mzFoundStDev1 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev1],
        row.mzFoundStDev2 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev2],
        row.mzFoundStDev3 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev3],
        row.mzFoundStDev4 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev4],
        row.mzFoundStDev5 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev5],
        row.mzFoundStDev6 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev6],
        row.mzFoundStDev7 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev7],
        row.mzFoundStDev8 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev8],
        row.mzFoundStDev9 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev9],
        row.mzFoundStDev10 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev10],
        row.mzFoundStDev11 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev11],
        row.mzFoundStDev12 = candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev12],
        row.mzAccuracy1 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy1],
        row.mzAccuracy2 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy2],
        row.mzAccuracy3 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy3],
        row.mzAccuracy4 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy4],
        row.mzAccuracy5 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy5],
        row.mzAccuracy6 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy6],
        row.mzAccuracy7 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy7],
        row.mzAccuracy8 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy8],
        row.mzAccuracy9 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy9],
        row.mzAccuracy10 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy10],
        row.mzAccuracy11 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy11],
        row.mzAccuracy12 = candidateScores->featuresArray[CandidateScores::Features::MzAccuracy12],
        row.targetKey = candidateScores->targetKey;
        row.peptideStringWithMods = candidateScores->targetDecoyCandidatePair->peptideStringWithMods();
        row.proteinGroup = candidateScores->proteinGroup;
        row.isDecoy = candidateScores->isDecoy;
        row.scanNumber = candidateScores->scanNumber;
        row.scanTime = candidateScores->scanTime;
        row.classifierScore = candidateScores->classifierScore;
        row.discriminantScore = candidateScores->discriminantScore;
        row.qValue = candidateScores->qValue;
        row.decoyRatio = candidateScores->decoyRatio;

        return row;
    }

};

#endif //PYTHIADIACPP_CANDIDATESCORES_H
