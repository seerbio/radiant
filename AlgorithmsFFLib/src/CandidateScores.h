//
// Created by anichols on 10/20/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORES_H
#define PYTHIADIACPP_CANDIDATESCORES_H


#include "AlgorithmsFFLib_Exports.h"
#include "AminoAcids.h"
#include "ParquetReader.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "PeptideStringWithMods.h"
#include "TargetDecoyCandidatePair.h"

/**
 * @brief Class representing scores for a candidate.
 *
 * This class defines a set of features and provides storage for scores related to a candidate.
 * It includes various features such as cosine similarity, peak intensities, charge, mass, scan time, etc.
 *
 */
class ALGORITHMSFFLIB_EXPORTS CandidateScores {

public:

    enum Features {
        CosineSimSum100 = 0,
        CosineSimSum100GreaterThan80,
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
        Charge,
        ScanTimeDelta,
        ScanTimeRange,
        ScanTimePd,
        ScanIonCount,
        MzNorm,
        Mass,
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
        AminoAcidCountB,
        AminoAcidCountJ,
        AminoAcidCountO,
        AminoAcidCountU,
        AminoAcidCountX,
        AminoAcidCountZ,
        MzFoundStDev1,
        MzFoundStDev2,
        MzFoundStDev3,//170
        MzFoundStDev4,
        MzFoundStDev5,
        MzFoundStDev6,
        MzFoundStDev7,
        MzFoundStDev8,
        MzFoundStDev9,
        MzFoundStDev10,
        MzFoundStDev11,
        MzFoundStDev12,
        MzAccuracy1,//180
        MzAccuracy2,
        MzAccuracy3,
        MzAccuracy4,
        MzAccuracy5,
        MzAccuracy6,
        MzAccuracy7,
        MzAccuracy8,
        MzAccuracy9,
        MzAccuracy10,
        MzAccuracy11,//190
        MzAccuracy12,
        AltTargetKeyIdCosineSimSumCharge1_OG,
        AltTargetKeyIdCosineSimSumCharge1_1,
        AltTargetKeyIdCosineSimSumCharge1_2,
        AltTargetKeyIdCosineSimSumCharge1_3,
        AltTargetKeyIdCosineSimSumCharge2_OG,
        AltTargetKeyIdCosineSimSumCharge2_1,
        AltTargetKeyIdCosineSimSumCharge2_2,
        AltTargetKeyIdCosineSimSumCharge2_3,
        AltTargetKeyIdCosineSimSumCharge3_OG,//200
        AltTargetKeyIdCosineSimSumCharge3_1,
        AltTargetKeyIdCosineSimSumCharge3_2,
        AltTargetKeyIdCosineSimSumCharge3_3,
        AltTargetKeyIdCosineSimSumCharge4_OG,
        AltTargetKeyIdCosineSimSumCharge4_1,
        AltTargetKeyIdCosineSimSumCharge4_2,
        AltTargetKeyIdCosineSimSumCharge4_3,
        AltTargetKeyIdTimeDeltaCharge1_1,
        AltTargetKeyIdTimeDeltaCharge1_2,
        AltTargetKeyIdTimeDeltaCharge1_3,//210
        AltTargetKeyIdTimeDeltaCharge2_1,
        AltTargetKeyIdTimeDeltaCharge2_2,
        AltTargetKeyIdTimeDeltaCharge2_3,
        AltTargetKeyIdTimeDeltaCharge3_1,
        AltTargetKeyIdTimeDeltaCharge3_2,
        AltTargetKeyIdTimeDeltaCharge3_3,
        AltTargetKeyIdTimeDeltaCharge4_1,
        AltTargetKeyIdTimeDeltaCharge4_2,
        AltTargetKeyIdTimeDeltaCharge4_3,
        Ms1MzMeanFound100,//220
        Ms1MzMeanFound45,
        Ms1MzMeanFound20,
        Ms1MzMeanFoundPreMono,
        Ms1MzMeanFoundIso1,
        Ms1MzMeanFoundIso2,
        Ms1MzMeanFound100PPM,
        Ms1MzMeanFound45PPM,
        Ms1MzMeanFound20PPM,
        Ms1MzMeanFoundPreMonoPPM,
        Ms1MzMeanFoundIso1PPM,//230
        Ms1MzMeanFoundIso2PPM,
        Ms1MzStDevFound100,
        Ms1MzStDevFound45,
        Ms1MzStDevFound20,
        Ms1MzStDevFoundPreMono,
        Ms1MzStDevFoundIso1,
        Ms1MzStDevFoundIso2,
        Ms1IntensityFound100,
        Ms1IntensityFound45,
        Ms1IntensityFound20,//240
        Ms1IntensityFoundPreMono,
        Ms1IntensityFoundIso1,
        Ms1IntensityFoundIso2,
        CosineSimSpectrumOverTime,
        CosineSimSpectrumOverTimeCubed,
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

    /**
    * @brief Returns a reference to the features array.
    * @return QVector<float>* A pointer to the features array.
    */
    [[nodiscard]] QVector<float>* featuresArrayRef();

    /**
    * @brief Initializes the features array with default values.
    */
    void initFeaturesArray();

    QVector<float> selectFeaturesArrayFeatures(const QVector<Features> &enumFeatures);

    static QVector<float> selectFeaturesArrayFeatures(
            const QVector<float> &featureArray,
            const QVector<Features> &enumFeatures
            );

};


namespace CandidateScoresReaderRowNamespace {

    const QString COS_SIM_SUM_100 = QStringLiteral("CosineSimSum100");
    const QString COS_SIM_SUM_100_GREATER_80 = QStringLiteral("CosineSimSum100GreaterThan80");
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
    const QString CHARGE = QStringLiteral("Charge");
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
    const QString AA_B = QStringLiteral("AminoAcidCountB");
    const QString AA_J = QStringLiteral("AminoAcidCountJ");
    const QString AA_O = QStringLiteral("AminoAcidCountO");
    const QString AA_U = QStringLiteral("AminoAcidCountU");
    const QString AA_X = QStringLiteral("AminoAcidCountX");
    const QString AA_Z = QStringLiteral("AminoAcidCountZ");
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

    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG1_OG = QStringLiteral("AltTargetKeyIdCosineSimSumCharge1_OG");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG1_1 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge1_1");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG1_2 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge1_2");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG1_3 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge1_3");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG2_OG = QStringLiteral("AltTargetKeyIdCosineSimSumCharge2_OG");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG2_1 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge2_1");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG2_2 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge2_2");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG2_3 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge2_3");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG3_OG = QStringLiteral("AltTargetKeyIdCosineSimSumCharge3_OG");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG3_1 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge3_1");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG3_2 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge3_2");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG3_3 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge3_3");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG4_OG = QStringLiteral("AltTargetKeyIdCosineSimSumCharge4_OG");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG4_1 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge4_1");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG4_2 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge4_2");
    const QString ALT_TARG_ID_COS_SIM_SUM_CHRG4_3 = QStringLiteral("AltTargetKeyIdCosineSimSumCharge4_3");

    const QString ALT_TARG_ID_TIME_DELTA_CHRG1_1 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge1_1");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG1_2 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge1_2");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG1_3 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge1_3");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG2_1 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge2_1");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG2_2 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge2_2");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG2_3 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge2_3");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG3_1 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge3_1");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG3_2 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge3_2");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG3_3 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge3_3");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG4_1 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge4_1");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG4_2 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge4_2");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG4_3 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge4_3");

    const QString MS1_MZ_MEAN_FND_100 = QStringLiteral("Ms1MzMeanFound100");
    const QString MS1_MZ_MEAN_FND_45 = QStringLiteral("Ms1MzMeanFound45");
    const QString MS1_MZ_MEAN_FND_20 = QStringLiteral("Ms1MzMeanFound20");
    const QString MS1_MZ_MEAN_FND_PRE_MONO = QStringLiteral("Ms1MzMeanFoundPreMono");
    const QString MS1_MZ_MEAN_FND_ISO1 = QStringLiteral("Ms1MzMeanFoundIso1");
    const QString MS1_MZ_MEAN_FND_ISO2 = QStringLiteral("Ms1MzMeanFoundIso2");
    const QString MS1_MZ_MEAN_FND_100_PPM = QStringLiteral("Ms1MzMeanFound100PPM");
    const QString MS1_MZ_MEAN_FND_45_PPM = QStringLiteral("Ms1MzMeanFound45PPM");
    const QString MS1_MZ_MEAN_FND_20_PPM = QStringLiteral("Ms1MzMeanFound20PPM");
    const QString MS1_MZ_MEAN_FND_PRE_MONO_PPM = QStringLiteral("Ms1MzMeanFoundPreMonoPPM");
    const QString MS1_MZ_MEAN_FND_ISO_1_PPM = QStringLiteral("Ms1MzMeanFoundIso1PPM");
    const QString MS1_MZ_MEAN_FND_ISO_2_PPM = QStringLiteral("Ms1MzMeanFoundIso2PPM");
    const QString MS1_MZ_MEAN_FND_100_STD = QStringLiteral("Ms1MzStDevFound100");
    const QString MS1_MZ_MEAN_FND_45_STD = QStringLiteral("Ms1MzStDevFound45");
    const QString MS1_MZ_MEAN_FND_20_STD = QStringLiteral("Ms1MzStDevFound20");
    const QString MS1_MZ_MEAN_FND_PRE_MONO_STD = QStringLiteral("Ms1MzStDevFoundPreMono");
    const QString MS1_MZ_MEAN_FND_ISO_1_STD = QStringLiteral("Ms1MzStDevFoundIso1");
    const QString MS1_MZ_MEAN_FND_ISO_2_STD = QStringLiteral("Ms1MzStDevFoundIso2");
    const QString MS1_INTZ_FND_100 = QStringLiteral("Ms1IntensityFound100");
    const QString MS1_INTZ_FND_45 = QStringLiteral("Ms1IntensityFound45");
    const QString MS1_INTZ_FND_20 = QStringLiteral("Ms1IntensityFound20");
    const QString MS1_INTZ_FND_PRE_MONO = QStringLiteral("Ms1IntensityFoundPreMono");
    const QString MS1_INTZ_FND_ISO_1 = QStringLiteral("Ms1IntensityFoundIso1");
    const QString MS1_INTZ_FND_ISO_2 = QStringLiteral("Ms1IntensityFoundIso2");

    const QStringList keysToCheck = {
            COS_SIM_SUM_100,
            COS_SIM_SUM_100_GREATER_80,
            ALL_MAX_IND_CNT,
            COS_SIM_SUM_MS1_100,
            COS_SIM_SPEC_CUBED,
            KL_DIV_SPEC_CUBE_RT,
            COS_SIM_SUM_45,
            COS_SIM_SUM_20,
            COS_SIM_SUM_TOP_6,
            COS_SIM_SUM_BOTTOM_6,
            TOP_BOTTOM_RATIO,
            TOP_BOTTOME_RATIO_NORM,
            CHARGE,
            MASS,
            SCAN_TIME_DELTA,
            SCAN_TIME_RANGE,
            SCAN_TIME_PD,
            SCAN_ION_CNT,
            MZ_NORM,
            KL_DIV_SPEC,
            COSINE_SIM_SPEC,
            COSINE_SIM_SUM_MS1_45,
            COSINE_SIM_SUM_MS1_20,
            COSINE_SIM_SUM_MS1_PRE_MONO,
            COSINE_SIM_SUM_MS1_ISO_1,
            COSINE_SIM_SUM_MS1_ISO_2,
            PEP_LEN_NORM,
            SCAN_TIME_PRED,
            THEO_FRAG_CNT,
            TOT_INT_LOG,
            PEAK_RATIO_1,
            PEAK_RATIO_2,
            PEAK_RATIO_3,
            SHADOW_COSINE_SIM_SUM,
            IRT_PRED,
            COS_SIM_ANCH_1,
            COS_SIM_ANCH_2,
            COS_SIM_ANCH_3,
            COS_SIM_ANCH_4,
            COS_SIM_ANCH_5,
            COS_SIM_ANCH_6,
            COS_SIM_ANCH_7,
            COS_SIM_ANCH_8,
            COS_SIM_ANCH_9,
            COS_SIM_ANCH_10,
            COS_SIM_ANCH_11,
            COS_SIM_ANCH_12,
            COS_SIM_ANCH_SHADOW_1,
            COS_SIM_ANCH_SHADOW_2,
            COS_SIM_ANCH_SHADOW_3,
            COS_SIM_ANCH_SHADOW_4,
            COS_SIM_ANCH_SHADOW_5,
            COS_SIM_ANCH_SHADOW_6,
            COS_SIM_ANCH_SHADOW_7,
            COS_SIM_ANCH_SHADOW_8,
            COS_SIM_ANCH_SHADOW_9,
            COS_SIM_ANCH_SHADOW_10,
            COS_SIM_ANCH_SHADOW_11,
            COS_SIM_ANCH_SHADOW_12,
            SHAD_INTS_RATIO_1,
            SHAD_INTS_RATIO_2,
            SHAD_INTS_RATIO_3,
            SHAD_INTS_RATIO_4,
            SHAD_INTS_RATIO_5,
            SHAD_INTS_RATIO_6,
            SHAD_INTS_RATIO_7,
            SHAD_INTS_RATIO_8,
            SHAD_INTS_RATIO_9,
            SHAD_INTS_RATIO_10,
            SHAD_INTS_RATIO_11,
            SHAD_INTS_RATIO_12,
            MZ_SEARCHED_1,
            MZ_SEARCHED_2,
            MZ_SEARCHED_3,
            MZ_SEARCHED_4,
            MZ_SEARCHED_5,
            MZ_SEARCHED_6,
            MZ_SEARCHED_7,
            MZ_SEARCHED_8,
            MZ_SEARCHED_9,
            MZ_SEARCHED_10,
            MZ_SEARCHED_11,
            MZ_SEARCHED_12,
            THEO_INTS_1,
            THEO_INTS_2,
            THEO_INTS_3,
            THEO_INTS_4,
            THEO_INTS_5,
            THEO_INTS_6,
            THEO_INTS_7,
            THEO_INTS_8,
            THEO_INTS_9,
            THEO_INTS_10,
            THEO_INTS_11,
            THEO_INTS_12,
            MZ_FND_MEAN_1,
            MZ_FND_MEAN_2,
            MZ_FND_MEAN_3,
            MZ_FND_MEAN_4,
            MZ_FND_MEAN_5,
            MZ_FND_MEAN_6,
            MZ_FND_MEAN_7,
            MZ_FND_MEAN_8,
            MZ_FND_MEAN_9,
            MZ_FND_MEAN_10,
            MZ_FND_MEAN_11,
            MZ_FND_MEAN_12,
            INTS_FND_MAX_1,
            INTS_FND_MAX_2,
            INTS_FND_MAX_3,
            INTS_FND_MAX_4,
            INTS_FND_MAX_5,
            INTS_FND_MAX_6,
            INTS_FND_MAX_7,
            INTS_FND_MAX_8,
            INTS_FND_MAX_9,
            INTS_FND_MAX_10,
            INTS_FND_MAX_11,
            INTS_FND_MAX_12,
            MZ_PK_LEN_NORM_1,
            MZ_PK_LEN_NORM_2,
            MZ_PK_LEN_NORM_3,
            MZ_PK_LEN_NORM_4,
            MZ_PK_LEN_NORM_5,
            MZ_PK_LEN_NORM_6,
            MZ_PK_LEN_NORM_7,
            MZ_PK_LEN_NORM_8,
            MZ_PK_LEN_NORM_9,
            MZ_PK_LEN_NORM_10,
            MZ_PK_LEN_NORM_11,
            MZ_PK_LEN_NORM_12,
            COL_APX_IND_RATIO_TO_ANCH_1,
            COL_APX_IND_RATIO_TO_ANCH_2,
            COL_APX_IND_RATIO_TO_ANCH_3,
            COL_APX_IND_RATIO_TO_ANCH_4,
            COL_APX_IND_RATIO_TO_ANCH_5,
            COL_APX_IND_RATIO_TO_ANCH_6,
            COL_APX_IND_RATIO_TO_ANCH_7,
            COL_APX_IND_RATIO_TO_ANCH_8,
            COL_APX_IND_RATIO_TO_ANCH_9,
            COL_APX_IND_RATIO_TO_ANCH_10,
            COL_APX_IND_RATIO_TO_ANCH_11,
            COL_APX_IND_RATIO_TO_ANCH_12,
            AA_A,
            AA_C,
            AA_D,
            AA_E,
            AA_F,
            AA_G,
            AA_H,
            AA_I,
            AA_K,
            AA_L,
            AA_M,
            AA_N,
            AA_P,
            AA_Q,
            AA_R,
            AA_S,
            AA_T,
            AA_V,
            AA_W,
            AA_Y,
            AA_B,
            AA_J,
            AA_O,
            AA_U,
            AA_X,
            AA_Z,
            MZ_FND_STDEV_1,
            MZ_FND_STDEV_2,
            MZ_FND_STDEV_3,
            MZ_FND_STDEV_4,
            MZ_FND_STDEV_5,
            MZ_FND_STDEV_6,
            MZ_FND_STDEV_7,
            MZ_FND_STDEV_8,
            MZ_FND_STDEV_9,
            MZ_FND_STDEV_10,
            MZ_FND_STDEV_11,
            MZ_FND_STDEV_12,
            MZ_ACC_1,
            MZ_ACC_2,
            MZ_ACC_3,
            MZ_ACC_4,
            MZ_ACC_5,
            MZ_ACC_6,
            MZ_ACC_7,
            MZ_ACC_8,
            MZ_ACC_9,
            MZ_ACC_10,
            MZ_ACC_11,
            MZ_ACC_12,
            TARG_KEY,
            PEP_STR_W_MODS,
            PROT_GRP,
            IS_DECOY,
            SCAN_NUM,
            SCAN_TIME,
            CLASS_SCR,
            DISC_SCR,
            Q_VAL,
            DECOY_RATIO,
            ALT_TARG_ID_COS_SIM_SUM_CHRG1_OG,
            ALT_TARG_ID_COS_SIM_SUM_CHRG1_1,
            ALT_TARG_ID_COS_SIM_SUM_CHRG1_2,
            ALT_TARG_ID_COS_SIM_SUM_CHRG1_3,
            ALT_TARG_ID_COS_SIM_SUM_CHRG2_OG,
            ALT_TARG_ID_COS_SIM_SUM_CHRG2_1,
            ALT_TARG_ID_COS_SIM_SUM_CHRG2_2,
            ALT_TARG_ID_COS_SIM_SUM_CHRG2_3,
            ALT_TARG_ID_COS_SIM_SUM_CHRG3_OG,
            ALT_TARG_ID_COS_SIM_SUM_CHRG3_1,
            ALT_TARG_ID_COS_SIM_SUM_CHRG3_2,
            ALT_TARG_ID_COS_SIM_SUM_CHRG3_3,
            ALT_TARG_ID_COS_SIM_SUM_CHRG4_OG,
            ALT_TARG_ID_COS_SIM_SUM_CHRG4_1,
            ALT_TARG_ID_COS_SIM_SUM_CHRG4_2,
            ALT_TARG_ID_COS_SIM_SUM_CHRG4_3,
            ALT_TARG_ID_TIME_DELTA_CHRG1_1,
            ALT_TARG_ID_TIME_DELTA_CHRG1_2,
            ALT_TARG_ID_TIME_DELTA_CHRG1_3,
            ALT_TARG_ID_TIME_DELTA_CHRG2_1,
            ALT_TARG_ID_TIME_DELTA_CHRG2_2,
            ALT_TARG_ID_TIME_DELTA_CHRG2_3,
            ALT_TARG_ID_TIME_DELTA_CHRG3_1,
            ALT_TARG_ID_TIME_DELTA_CHRG3_2,
            ALT_TARG_ID_TIME_DELTA_CHRG3_3,
            ALT_TARG_ID_TIME_DELTA_CHRG4_1,
            ALT_TARG_ID_TIME_DELTA_CHRG4_2,
            ALT_TARG_ID_TIME_DELTA_CHRG4_3,
            MS1_MZ_MEAN_FND_100,
            MS1_MZ_MEAN_FND_45,
            MS1_MZ_MEAN_FND_20,
            MS1_MZ_MEAN_FND_PRE_MONO,
            MS1_MZ_MEAN_FND_ISO1,
            MS1_MZ_MEAN_FND_ISO2,
            MS1_MZ_MEAN_FND_100_PPM,
            MS1_MZ_MEAN_FND_45_PPM,
            MS1_MZ_MEAN_FND_20_PPM,
            MS1_MZ_MEAN_FND_PRE_MONO_PPM,
            MS1_MZ_MEAN_FND_ISO_1_PPM,
            MS1_MZ_MEAN_FND_ISO_2_PPM,
            MS1_MZ_MEAN_FND_100_STD,
            MS1_MZ_MEAN_FND_45_STD,
            MS1_MZ_MEAN_FND_20_STD,
            MS1_MZ_MEAN_FND_PRE_MONO_STD,
            MS1_MZ_MEAN_FND_ISO_1_STD,
            MS1_MZ_MEAN_FND_ISO_2_STD,
            MS1_INTZ_FND_100,
            MS1_INTZ_FND_45,
            MS1_INTZ_FND_20,
            MS1_INTZ_FND_PRE_MONO,
            MS1_INTZ_FND_ISO_1,
            MS1_INTZ_FND_ISO_2
    };

}//namespace

struct ALGORITHMSFFLIB_EXPORTS CandidateScoresReaderRow : public ParquetReaderInputBase {

    float cosineSimSum100 = -1.0;
    float cosineSimSum100Greater80 = -1.0;
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
    float charge = -1.0;
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
    float aminoAcidCountB = -1.0;
    float aminoAcidCountJ = -1.0;
    float aminoAcidCountO = -1.0;
    float aminoAcidCountU = -1.0;
    float aminoAcidCountX = -1.0;
    float aminoAcidCountZ = -1.0;
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
    PeptideStringWithMods peptideStringWithMods;
    QString proteinGroup;
    bool isDecoy = false;
    ScanNumber scanNumber = -1;
    ScanTime scanTime = -1.0;
    double classifierScore = -1.0;
    double discriminantScore = -1.0;
    double qValue = 1.0;
    double decoyRatio = -1.0;

    float altTargetKeyIdCosineSimSumCharge1_OG = -1.0;
    float altTargetKeyIdCosineSimSumCharge1_1 = -1.0;
    float altTargetKeyIdCosineSimSumCharge1_2 = -1.0;
    float altTargetKeyIdCosineSimSumCharge1_3 = -1.0;
    float altTargetKeyIdCosineSimSumCharge2_OG = -1.0;
    float altTargetKeyIdCosineSimSumCharge2_1 = -1.0;
    float altTargetKeyIdCosineSimSumCharge2_2 = -1.0;
    float altTargetKeyIdCosineSimSumCharge2_3 = -1.0;
    float altTargetKeyIdCosineSimSumCharge3_OG = -1.0;
    float altTargetKeyIdCosineSimSumCharge3_1 = -1.0;
    float altTargetKeyIdCosineSimSumCharge3_2 = -1.0;
    float altTargetKeyIdCosineSimSumCharge3_3 = -1.0;
    float altTargetKeyIdCosineSimSumCharge4_OG = -1.0;
    float altTargetKeyIdCosineSimSumCharge4_1 = -1.0;
    float altTargetKeyIdCosineSimSumCharge4_2 = -1.0;
    float altTargetKeyIdCosineSimSumCharge4_3 = -1.0;

    float altTargetKeyIdTimeDeltaCharge1_1 = -1.0;
    float altTargetKeyIdTimeDeltaCharge1_2 = -1.0;
    float altTargetKeyIdTimeDeltaCharge1_3 = -1.0;
    float altTargetKeyIdTimeDeltaCharge2_1 = -1.0;
    float altTargetKeyIdTimeDeltaCharge2_2 = -1.0;
    float altTargetKeyIdTimeDeltaCharge2_3 = -1.0;
    float altTargetKeyIdTimeDeltaCharge3_1 = -1.0;
    float altTargetKeyIdTimeDeltaCharge3_2 = -1.0;
    float altTargetKeyIdTimeDeltaCharge3_3 = -1.0;
    float altTargetKeyIdTimeDeltaCharge4_1 = -1.0;
    float altTargetKeyIdTimeDeltaCharge4_2 = -1.0;
    float altTargetKeyIdTimeDeltaCharge4_3 = -1.0;

    float ms1MzMeanFound100 = -1.0;
    float ms1MzMeanFound45 = -1.0;
    float ms1MzMeanFound20 = -1.0;
    float ms1MzMeanFoundPreMono = -1.0;
    float ms1MzMeanFoundIso1 = -1.0;
    float ms1MzMeanFoundIso2 = -1.0;
    float ms1MzMeanFound100PPM = -1.0;
    float ms1MzMeanFound45PPM = -1.0;
    float ms1MzMeanFound20PPM = -1.0;
    float ms1MzMeanFoundPreMonoPPM = -1.0;
    float ms1MzMeanFoundIso1PPM = -1.0;
    float ms1MzMeanFoundIso2PPM = -1.0;
    float ms1MzStDevFound100 = -1.0;
    float ms1MzStDevFound45 = -1.0;
    float ms1MzStDevFound20 = -1.0;
    float ms1MzStDevFoundPreMono = -1.0;
    float ms1MzStDevFoundIso1 = -1.0;
    float ms1MzStDevFoundIso2 = -1.0;
    float ms1IntensityFound100 = -1.0;
    float ms1IntensityFound45 = -1.0;
    float ms1IntensityFound20 = -1.0;
    float ms1IntensityFoundPreMono = -1.0;
    float ms1IntensityFoundIso1 = -1.0;
    float ms1IntensityFoundIso2 = -1.0;

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace CandidateScoresReaderRowNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        cosineSimSum100 = dataMap.value(COS_SIM_SUM_100).toFloat();
        cosineSimSum100Greater80 = dataMap.value(COS_SIM_SUM_100_GREATER_80).toFloat();
        allignedMaxIndexesCount = dataMap.value(ALL_MAX_IND_CNT).toFloat();
        cosineSim100MS1 = dataMap.value(COS_SIM_SUM_MS1_100).toFloat();
        cosineSimSpectrumCubed = dataMap.value(COS_SIM_SPEC_CUBED).toFloat();
        klDivSpectrumCubeRoot = dataMap.value(KL_DIV_SPEC_CUBE_RT).toFloat();
        cosineSimSum45 = dataMap.value(COS_SIM_SUM_45).toFloat();
        cosineSimSum20 = dataMap.value(COS_SIM_SUM_20).toFloat();
        cosineSimSumTop6 = dataMap.value(COS_SIM_SUM_TOP_6).toFloat();
        cosineSimSumBottom6 = dataMap.value(COS_SIM_SUM_BOTTOM_6).toFloat();
        topBottomRatio = dataMap.value(TOP_BOTTOM_RATIO).toFloat();
        topBottomRatioNorm = dataMap.value(TOP_BOTTOME_RATIO_NORM).toFloat();
        charge = dataMap.value(CHARGE).toFloat();
        mass = dataMap.value(MASS).toFloat();
        scanTimeDelta = dataMap.value(SCAN_TIME_DELTA).toFloat();
        scanTimeRange = dataMap.value(SCAN_TIME_RANGE).toFloat();
        scanTimePd = dataMap.value(SCAN_TIME_PD).toFloat();
        scanIonCount = dataMap.value(SCAN_ION_CNT).toFloat();
        mzNorm = dataMap.value(MZ_NORM).toFloat();
        klDivSpectrum = dataMap.value(KL_DIV_SPEC).toFloat();
        cosineSimSpectrum = dataMap.value(COSINE_SIM_SPEC).toFloat();
        cosineSim45MS1 = dataMap.value(COSINE_SIM_SUM_MS1_45).toFloat();
        cosineSim20MS1 = dataMap.value(COSINE_SIM_SUM_MS1_20).toFloat();
        cosineSim100MS1PreMono = dataMap.value(COSINE_SIM_SUM_MS1_PRE_MONO).toFloat();
        cosineSim100MS1Iso1 = dataMap.value(COSINE_SIM_SUM_MS1_ISO_1).toFloat();
        cosineSim100MS1Iso2 = dataMap.value(COSINE_SIM_SUM_MS1_ISO_2).toFloat();
        peptideLengthNorm = dataMap.value(PEP_LEN_NORM).toFloat();
        scanTimePredicted = dataMap.value(SCAN_TIME_PRED).toFloat();
        theoFragmentCount = dataMap.value(THEO_FRAG_CNT).toFloat();
        totalIntensityLog = dataMap.value(TOT_INT_LOG).toFloat();
        peakShapeRatio1 = dataMap.value(PEAK_RATIO_1).toFloat();
        peakShapeRatio2 = dataMap.value(PEAK_RATIO_2).toFloat();
        peakShapeRatio3 = dataMap.value(PEAK_RATIO_3).toFloat();
        shadowsCosineSimSum = dataMap.value(SHADOW_COSINE_SIM_SUM).toFloat();
        iRtPredicted = dataMap.value(IRT_PRED).toFloat();
        cosineSimToAnchor1 = dataMap.value(COS_SIM_ANCH_1).toFloat();
        cosineSimToAnchor2 = dataMap.value(COS_SIM_ANCH_2).toFloat();
        cosineSimToAnchor3 = dataMap.value(COS_SIM_ANCH_3).toFloat();
        cosineSimToAnchor4 = dataMap.value(COS_SIM_ANCH_4).toFloat();
        cosineSimToAnchor5 = dataMap.value(COS_SIM_ANCH_5).toFloat();
        cosineSimToAnchor6 = dataMap.value(COS_SIM_ANCH_6).toFloat();
        cosineSimToAnchor7 = dataMap.value(COS_SIM_ANCH_7).toFloat();
        cosineSimToAnchor8 = dataMap.value(COS_SIM_ANCH_8).toFloat();
        cosineSimToAnchor9 = dataMap.value(COS_SIM_ANCH_9).toFloat();
        cosineSimToAnchor10 = dataMap.value(COS_SIM_ANCH_10).toFloat();
        cosineSimToAnchor11 = dataMap.value(COS_SIM_ANCH_11).toFloat();
        cosineSimToAnchor12 = dataMap.value(COS_SIM_ANCH_12).toFloat();
        cosineSimShadowsToAnchor1 = dataMap.value(COS_SIM_ANCH_SHADOW_1).toFloat();
        cosineSimShadowsToAnchor2 = dataMap.value(COS_SIM_ANCH_SHADOW_2).toFloat();
        cosineSimShadowsToAnchor3 = dataMap.value(COS_SIM_ANCH_SHADOW_3).toFloat();
        cosineSimShadowsToAnchor4 = dataMap.value(COS_SIM_ANCH_SHADOW_4).toFloat();
        cosineSimShadowsToAnchor5 = dataMap.value(COS_SIM_ANCH_SHADOW_5).toFloat();
        cosineSimShadowsToAnchor6 = dataMap.value(COS_SIM_ANCH_SHADOW_6).toFloat();
        cosineSimShadowsToAnchor7 = dataMap.value(COS_SIM_ANCH_SHADOW_7).toFloat();
        cosineSimShadowsToAnchor8 = dataMap.value(COS_SIM_ANCH_SHADOW_8).toFloat();
        cosineSimShadowsToAnchor9 = dataMap.value(COS_SIM_ANCH_SHADOW_9).toFloat();
        cosineSimShadowsToAnchor10 = dataMap.value(COS_SIM_ANCH_SHADOW_10).toFloat();
        cosineSimShadowsToAnchor11 = dataMap.value(COS_SIM_ANCH_SHADOW_11).toFloat();
        cosineSimShadowsToAnchor12 = dataMap.value(COS_SIM_ANCH_SHADOW_12).toFloat();
        shadowsIntensityRatio1 = dataMap.value(SHAD_INTS_RATIO_1).toFloat();
        shadowsIntensityRatio2 = dataMap.value(SHAD_INTS_RATIO_2).toFloat();
        shadowsIntensityRatio3 = dataMap.value(SHAD_INTS_RATIO_3).toFloat();
        shadowsIntensityRatio4 = dataMap.value(SHAD_INTS_RATIO_4).toFloat();
        shadowsIntensityRatio5 = dataMap.value(SHAD_INTS_RATIO_5).toFloat();
        shadowsIntensityRatio6 = dataMap.value(SHAD_INTS_RATIO_6).toFloat();
        shadowsIntensityRatio7 = dataMap.value(SHAD_INTS_RATIO_7).toFloat();
        shadowsIntensityRatio8 = dataMap.value(SHAD_INTS_RATIO_8).toFloat();
        shadowsIntensityRatio9 = dataMap.value(SHAD_INTS_RATIO_9).toFloat();
        shadowsIntensityRatio10 = dataMap.value(SHAD_INTS_RATIO_10).toFloat();
        shadowsIntensityRatio11 = dataMap.value(SHAD_INTS_RATIO_11).toFloat();
        shadowsIntensityRatio12 = dataMap.value(SHAD_INTS_RATIO_12).toFloat();
        mzSearched1 = dataMap.value(MZ_SEARCHED_1).toFloat();
        mzSearched2 = dataMap.value(MZ_SEARCHED_2).toFloat();
        mzSearched3 = dataMap.value(MZ_SEARCHED_3).toFloat();
        mzSearched4 = dataMap.value(MZ_SEARCHED_4).toFloat();
        mzSearched5 = dataMap.value(MZ_SEARCHED_5).toFloat();
        mzSearched6 = dataMap.value(MZ_SEARCHED_6).toFloat();
        mzSearched7 = dataMap.value(MZ_SEARCHED_7).toFloat();
        mzSearched8 = dataMap.value(MZ_SEARCHED_8).toFloat();
        mzSearched9 = dataMap.value(MZ_SEARCHED_9).toFloat();
        mzSearched10 = dataMap.value(MZ_SEARCHED_10).toFloat();
        mzSearched11 = dataMap.value(MZ_SEARCHED_11).toFloat();
        mzSearched12 = dataMap.value(MZ_SEARCHED_12).toFloat();
        theoIntensity1 = dataMap.value(THEO_INTS_1).toFloat();
        theoIntensity2 = dataMap.value(THEO_INTS_2).toFloat();
        theoIntensity3 = dataMap.value(THEO_INTS_3).toFloat();
        theoIntensity4 = dataMap.value(THEO_INTS_4).toFloat();
        theoIntensity5 = dataMap.value(THEO_INTS_5).toFloat();
        theoIntensity6 = dataMap.value(THEO_INTS_6).toFloat();
        theoIntensity7 = dataMap.value(THEO_INTS_7).toFloat();
        theoIntensity8 = dataMap.value(THEO_INTS_8).toFloat();
        theoIntensity9 = dataMap.value(THEO_INTS_9).toFloat();
        theoIntensity10 = dataMap.value(THEO_INTS_10).toFloat();
        theoIntensity11 = dataMap.value(THEO_INTS_11).toFloat();
        theoIntensity12 = dataMap.value(THEO_INTS_12).toFloat();
        mzFoundMean1 = dataMap.value(MZ_FND_MEAN_1).toFloat();
        mzFoundMean2 = dataMap.value(MZ_FND_MEAN_2).toFloat();
        mzFoundMean3 = dataMap.value(MZ_FND_MEAN_3).toFloat();
        mzFoundMean4 = dataMap.value(MZ_FND_MEAN_4).toFloat();
        mzFoundMean5 = dataMap.value(MZ_FND_MEAN_5).toFloat();
        mzFoundMean6 = dataMap.value(MZ_FND_MEAN_6).toFloat();
        mzFoundMean7 = dataMap.value(MZ_FND_MEAN_7).toFloat();
        mzFoundMean8 = dataMap.value(MZ_FND_MEAN_8).toFloat();
        mzFoundMean9 = dataMap.value(MZ_FND_MEAN_9).toFloat();
        mzFoundMean10 = dataMap.value(MZ_FND_MEAN_10).toFloat();
        mzFoundMean11 = dataMap.value(MZ_FND_MEAN_11).toFloat();
        mzFoundMean12 = dataMap.value(MZ_FND_MEAN_12).toFloat();
        intensityFoundMax1 = dataMap.value(INTS_FND_MAX_1).toFloat();
        intensityFoundMax2 = dataMap.value(INTS_FND_MAX_2).toFloat();
        intensityFoundMax3 = dataMap.value(INTS_FND_MAX_3).toFloat();
        intensityFoundMax4 = dataMap.value(INTS_FND_MAX_4).toFloat();
        intensityFoundMax5 = dataMap.value(INTS_FND_MAX_5).toFloat();
        intensityFoundMax6 = dataMap.value(INTS_FND_MAX_6).toFloat();
        intensityFoundMax7 = dataMap.value(INTS_FND_MAX_7).toFloat();
        intensityFoundMax8 = dataMap.value(INTS_FND_MAX_8).toFloat();
        intensityFoundMax9 = dataMap.value(INTS_FND_MAX_9).toFloat();
        intensityFoundMax10 = dataMap.value(INTS_FND_MAX_10).toFloat();
        intensityFoundMax11 = dataMap.value(INTS_FND_MAX_11).toFloat();
        intensityFoundMax12 = dataMap.value(INTS_FND_MAX_12).toFloat();
        mzPeakLengthsNorm1 = dataMap.value(MZ_PK_LEN_NORM_1).toFloat();
        mzPeakLengthsNorm2 = dataMap.value(MZ_PK_LEN_NORM_2).toFloat();
        mzPeakLengthsNorm3 = dataMap.value(MZ_PK_LEN_NORM_3).toFloat();
        mzPeakLengthsNorm4 = dataMap.value(MZ_PK_LEN_NORM_4).toFloat();
        mzPeakLengthsNorm5 = dataMap.value(MZ_PK_LEN_NORM_5).toFloat();
        mzPeakLengthsNorm6 = dataMap.value(MZ_PK_LEN_NORM_6).toFloat();
        mzPeakLengthsNorm7 = dataMap.value(MZ_PK_LEN_NORM_7).toFloat();
        mzPeakLengthsNorm8 = dataMap.value(MZ_PK_LEN_NORM_8).toFloat();
        mzPeakLengthsNorm9 = dataMap.value(MZ_PK_LEN_NORM_9).toFloat();
        mzPeakLengthsNorm10 = dataMap.value(MZ_PK_LEN_NORM_10).toFloat();
        mzPeakLengthsNorm11 = dataMap.value(MZ_PK_LEN_NORM_11).toFloat();
        mzPeakLengthsNorm12 = dataMap.value(MZ_PK_LEN_NORM_12).toFloat();
        columnApexIndexRatiosToAnchor1 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_1).toFloat();
        columnApexIndexRatiosToAnchor2 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_2).toFloat();
        columnApexIndexRatiosToAnchor3 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_3).toFloat();
        columnApexIndexRatiosToAnchor4 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_4).toFloat();
        columnApexIndexRatiosToAnchor5 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_5).toFloat();
        columnApexIndexRatiosToAnchor6 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_6).toFloat();
        columnApexIndexRatiosToAnchor7 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_7).toFloat();
        columnApexIndexRatiosToAnchor8 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_8).toFloat();
        columnApexIndexRatiosToAnchor9 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_9).toFloat();
        columnApexIndexRatiosToAnchor10 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_10).toFloat();
        columnApexIndexRatiosToAnchor11 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_11).toFloat();
        columnApexIndexRatiosToAnchor12 = dataMap.value(COL_APX_IND_RATIO_TO_ANCH_12).toFloat();
        aminoAcidCountA = dataMap.value(AA_A).toFloat();
        aminoAcidCountC = dataMap.value(AA_C).toFloat();
        aminoAcidCountD = dataMap.value(AA_D).toFloat();
        aminoAcidCountE = dataMap.value(AA_E).toFloat();
        aminoAcidCountF = dataMap.value(AA_F).toFloat();
        aminoAcidCountG = dataMap.value(AA_G).toFloat();
        aminoAcidCountH = dataMap.value(AA_H).toFloat();
        aminoAcidCountI = dataMap.value(AA_I).toFloat();
        aminoAcidCountK = dataMap.value(AA_K).toFloat();
        aminoAcidCountL = dataMap.value(AA_L).toFloat();
        aminoAcidCountM = dataMap.value(AA_M).toFloat();
        aminoAcidCountN = dataMap.value(AA_N).toFloat();
        aminoAcidCountP = dataMap.value(AA_P).toFloat();
        aminoAcidCountQ = dataMap.value(AA_Q).toFloat();
        aminoAcidCountR = dataMap.value(AA_R).toFloat();
        aminoAcidCountS = dataMap.value(AA_S).toFloat();
        aminoAcidCountT = dataMap.value(AA_T).toFloat();
        aminoAcidCountV = dataMap.value(AA_V).toFloat();
        aminoAcidCountW = dataMap.value(AA_W).toFloat();
        aminoAcidCountY = dataMap.value(AA_Y).toFloat();
        aminoAcidCountB = dataMap.value(AA_B).toFloat();
        aminoAcidCountJ = dataMap.value(AA_J).toFloat();
        aminoAcidCountO = dataMap.value(AA_O).toFloat();
        aminoAcidCountU = dataMap.value(AA_U).toFloat();
        aminoAcidCountX = dataMap.value(AA_X).toFloat();
        aminoAcidCountZ = dataMap.value(AA_Z).toFloat();
        mzFoundStDev1 = dataMap.value(MZ_FND_MEAN_1).toFloat();
        mzFoundStDev2 = dataMap.value(MZ_FND_MEAN_2).toFloat();
        mzFoundStDev3 = dataMap.value(MZ_FND_MEAN_3).toFloat();
        mzFoundStDev4 = dataMap.value(MZ_FND_MEAN_4).toFloat();
        mzFoundStDev5 = dataMap.value(MZ_FND_MEAN_5).toFloat();
        mzFoundStDev6 = dataMap.value(MZ_FND_MEAN_6).toFloat();
        mzFoundStDev7 = dataMap.value(MZ_FND_MEAN_8).toFloat();
        mzFoundStDev8 = dataMap.value(MZ_FND_MEAN_8).toFloat();
        mzFoundStDev9 = dataMap.value(MZ_FND_MEAN_9).toFloat();
        mzFoundStDev10 = dataMap.value(MZ_FND_MEAN_10).toFloat();
        mzFoundStDev11 = dataMap.value(MZ_FND_MEAN_11).toFloat();
        mzFoundStDev12 = dataMap.value(MZ_FND_MEAN_12).toFloat();
        mzAccuracy1 = dataMap.value(MZ_ACC_1).toFloat();
        mzAccuracy2 = dataMap.value(MZ_ACC_2).toFloat();
        mzAccuracy3 = dataMap.value(MZ_ACC_3).toFloat();
        mzAccuracy4 = dataMap.value(MZ_ACC_4).toFloat();
        mzAccuracy5 = dataMap.value(MZ_ACC_5).toFloat();
        mzAccuracy6 = dataMap.value(MZ_ACC_6).toFloat();
        mzAccuracy7 = dataMap.value(MZ_ACC_7).toFloat();
        mzAccuracy8 = dataMap.value(MZ_ACC_8).toFloat();
        mzAccuracy9 = dataMap.value(MZ_ACC_9).toFloat();
        mzAccuracy10 = dataMap.value(MZ_ACC_10).toFloat();
        mzAccuracy11 = dataMap.value(MZ_ACC_11).toFloat();
        mzAccuracy12 = dataMap.value(MZ_ACC_12).toFloat();

        targetKey = dataMap.value(TARG_KEY).toString();
        peptideStringWithMods = PeptideStringWithMods(dataMap.value(PEP_STR_W_MODS).toString());
        proteinGroup = dataMap.value(PROT_GRP).toString();
        isDecoy = dataMap.value(IS_DECOY).toBool();
        scanNumber = dataMap.value(SCAN_NUM).toInt();
        scanTime = dataMap.value(SCAN_TIME).toFloat();
        classifierScore = dataMap.value(CLASS_SCR).toDouble();
        discriminantScore = dataMap.value(DISC_SCR).toDouble();
        qValue = dataMap.value(Q_VAL).toDouble();
        decoyRatio = dataMap.value(DECOY_RATIO).toDouble();

        altTargetKeyIdCosineSimSumCharge1_OG = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG1_OG).toFloat();
        altTargetKeyIdCosineSimSumCharge1_1 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG1_1).toFloat();
        altTargetKeyIdCosineSimSumCharge1_2 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG1_2).toFloat();
        altTargetKeyIdCosineSimSumCharge1_3 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG1_3).toFloat();
        altTargetKeyIdCosineSimSumCharge2_OG = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG2_OG).toFloat();
        altTargetKeyIdCosineSimSumCharge2_1 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG2_1).toFloat();
        altTargetKeyIdCosineSimSumCharge2_2 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG2_2).toFloat();
        altTargetKeyIdCosineSimSumCharge2_3 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG2_3).toFloat();
        altTargetKeyIdCosineSimSumCharge3_OG = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG3_OG).toFloat();
        altTargetKeyIdCosineSimSumCharge3_1 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG3_1).toFloat();
        altTargetKeyIdCosineSimSumCharge3_2 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG3_2).toFloat();
        altTargetKeyIdCosineSimSumCharge3_3 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG3_3).toFloat();
        altTargetKeyIdCosineSimSumCharge4_OG = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG4_OG).toFloat();
        altTargetKeyIdCosineSimSumCharge4_1 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG4_1).toFloat();
        altTargetKeyIdCosineSimSumCharge4_2 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG4_2).toFloat();
        altTargetKeyIdCosineSimSumCharge4_3 = dataMap.value(ALT_TARG_ID_COS_SIM_SUM_CHRG4_3).toFloat();
        altTargetKeyIdTimeDeltaCharge1_1 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG1_1).toFloat();
        altTargetKeyIdTimeDeltaCharge1_2 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG1_2).toFloat();
        altTargetKeyIdTimeDeltaCharge1_3 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG1_3).toFloat();
        altTargetKeyIdTimeDeltaCharge2_1 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG2_1).toFloat();
        altTargetKeyIdTimeDeltaCharge2_2 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG2_2).toFloat();
        altTargetKeyIdTimeDeltaCharge2_3 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG2_3).toFloat();
        altTargetKeyIdTimeDeltaCharge3_1 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG3_1).toFloat();
        altTargetKeyIdTimeDeltaCharge3_2 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG3_2).toFloat();
        altTargetKeyIdTimeDeltaCharge3_3 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG3_3).toFloat();
        altTargetKeyIdTimeDeltaCharge4_1 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG4_1).toFloat();
        altTargetKeyIdTimeDeltaCharge4_2 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG4_2).toFloat();
        altTargetKeyIdTimeDeltaCharge4_3 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG4_3).toFloat();
        ms1MzMeanFound100 = dataMap.value(MS1_MZ_MEAN_FND_100).toFloat();
        ms1MzMeanFound45 = dataMap.value(MS1_MZ_MEAN_FND_45).toFloat();
        ms1MzMeanFound20 = dataMap.value(MS1_MZ_MEAN_FND_20).toFloat();
        ms1MzMeanFoundPreMono = dataMap.value(MS1_MZ_MEAN_FND_PRE_MONO).toFloat();
        ms1MzMeanFoundIso1 = dataMap.value(MS1_MZ_MEAN_FND_ISO1).toFloat();
        ms1MzMeanFoundIso2 = dataMap.value(MS1_MZ_MEAN_FND_ISO2).toFloat();
        ms1MzMeanFound100PPM = dataMap.value(MS1_MZ_MEAN_FND_100_PPM).toFloat();
        ms1MzMeanFound45PPM = dataMap.value(MS1_MZ_MEAN_FND_45_PPM).toFloat();
        ms1MzMeanFound20PPM = dataMap.value(MS1_MZ_MEAN_FND_20_PPM).toFloat();
        ms1MzMeanFoundPreMonoPPM = dataMap.value(MS1_MZ_MEAN_FND_PRE_MONO_PPM).toFloat();
        ms1MzMeanFoundIso1PPM = dataMap.value(MS1_MZ_MEAN_FND_ISO1).toFloat();
        ms1MzMeanFoundIso2PPM = dataMap.value(MS1_MZ_MEAN_FND_ISO2).toFloat();
        ms1MzStDevFound100 = dataMap.value(MS1_MZ_MEAN_FND_100_STD).toFloat();
        ms1MzStDevFound45 = dataMap.value(MS1_MZ_MEAN_FND_45_STD).toFloat();
        ms1MzStDevFound20 = dataMap.value(MS1_MZ_MEAN_FND_20_STD).toFloat();
        ms1MzStDevFoundPreMono = dataMap.value(MS1_MZ_MEAN_FND_PRE_MONO_STD).toFloat();
        ms1MzStDevFoundIso1 = dataMap.value(MS1_MZ_MEAN_FND_ISO_1_STD).toFloat();
        ms1MzStDevFoundIso2 = dataMap.value(MS1_MZ_MEAN_FND_ISO_2_STD).toFloat();
        ms1IntensityFound100 = dataMap.value(MS1_INTZ_FND_100).toFloat();
        ms1IntensityFound45 = dataMap.value(MS1_INTZ_FND_45).toFloat();
        ms1IntensityFound20 = dataMap.value(MS1_INTZ_FND_20).toFloat();
        ms1IntensityFoundPreMono = dataMap.value(MS1_INTZ_FND_PRE_MONO).toFloat();
        ms1IntensityFoundIso1 = dataMap.value(MS1_INTZ_FND_ISO_1).toFloat();
        ms1IntensityFoundIso2 = dataMap.value(MS1_INTZ_FND_ISO_2).toFloat();

        ERR_RETURN
    }

    QMap<QString, QVariant> map() override {

        using namespace CandidateScoresReaderRowNamespace;

        return {
                {COS_SIM_SUM_100, QVariant(cosineSimSum100)},
                {COS_SIM_SUM_100_GREATER_80, QVariant(cosineSimSum100Greater80)},
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
                {CHARGE, QVariant(charge)},
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
                {AA_B, QVariant(aminoAcidCountB)},
                {AA_J, QVariant(aminoAcidCountJ)},
                {AA_O, QVariant(aminoAcidCountO)},
                {AA_U, QVariant(aminoAcidCountU)},
                {AA_X, QVariant(aminoAcidCountX)},
                {AA_Z, QVariant(aminoAcidCountZ)},
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
                {DECOY_RATIO, QVariant(decoyRatio)},

                {ALT_TARG_ID_COS_SIM_SUM_CHRG1_OG, QVariant(altTargetKeyIdCosineSimSumCharge1_OG)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG1_1, QVariant(altTargetKeyIdCosineSimSumCharge1_1)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG1_2, QVariant(altTargetKeyIdCosineSimSumCharge1_2)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG1_3, QVariant(altTargetKeyIdCosineSimSumCharge1_3)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG2_OG, QVariant(altTargetKeyIdCosineSimSumCharge2_OG)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG2_1, QVariant(altTargetKeyIdCosineSimSumCharge2_1)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG2_2, QVariant(altTargetKeyIdCosineSimSumCharge2_2)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG2_3, QVariant(altTargetKeyIdCosineSimSumCharge2_3)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG3_OG, QVariant(altTargetKeyIdCosineSimSumCharge3_OG)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG3_1, QVariant(altTargetKeyIdCosineSimSumCharge3_1)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG3_2, QVariant(altTargetKeyIdCosineSimSumCharge3_2)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG3_3, QVariant(altTargetKeyIdCosineSimSumCharge3_3)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG4_OG, QVariant(altTargetKeyIdCosineSimSumCharge4_OG)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG4_1, QVariant(altTargetKeyIdCosineSimSumCharge4_1)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG4_2, QVariant(altTargetKeyIdCosineSimSumCharge4_2)},
                {ALT_TARG_ID_COS_SIM_SUM_CHRG4_3, QVariant(altTargetKeyIdCosineSimSumCharge4_3)},
                {ALT_TARG_ID_TIME_DELTA_CHRG1_1 , QVariant(altTargetKeyIdTimeDeltaCharge1_1)},
                {ALT_TARG_ID_TIME_DELTA_CHRG1_2 , QVariant(altTargetKeyIdTimeDeltaCharge1_2)},
                {ALT_TARG_ID_TIME_DELTA_CHRG1_3 , QVariant(altTargetKeyIdTimeDeltaCharge1_3)},
                {ALT_TARG_ID_TIME_DELTA_CHRG2_1 , QVariant(altTargetKeyIdTimeDeltaCharge2_1)},
                {ALT_TARG_ID_TIME_DELTA_CHRG2_2 , QVariant(altTargetKeyIdTimeDeltaCharge2_2)},
                {ALT_TARG_ID_TIME_DELTA_CHRG2_3 , QVariant(altTargetKeyIdTimeDeltaCharge2_3)},
                {ALT_TARG_ID_TIME_DELTA_CHRG3_1 , QVariant(altTargetKeyIdTimeDeltaCharge3_1)},
                {ALT_TARG_ID_TIME_DELTA_CHRG3_2 , QVariant(altTargetKeyIdTimeDeltaCharge3_2)},
                {ALT_TARG_ID_TIME_DELTA_CHRG3_3 , QVariant(altTargetKeyIdTimeDeltaCharge3_3)},
                {ALT_TARG_ID_TIME_DELTA_CHRG4_1 , QVariant(altTargetKeyIdTimeDeltaCharge4_1)},
                {ALT_TARG_ID_TIME_DELTA_CHRG4_2 , QVariant(altTargetKeyIdTimeDeltaCharge4_2)},
                {ALT_TARG_ID_TIME_DELTA_CHRG4_3 , QVariant(altTargetKeyIdTimeDeltaCharge4_3)},
                {MS1_MZ_MEAN_FND_100, QVariant(ms1MzMeanFound100)},
                {MS1_MZ_MEAN_FND_45, QVariant(ms1MzMeanFound45)},
                {MS1_MZ_MEAN_FND_20, QVariant(ms1MzMeanFound20)},
                {MS1_MZ_MEAN_FND_PRE_MONO, QVariant(ms1MzMeanFoundPreMono)},
                {MS1_MZ_MEAN_FND_ISO1, QVariant(ms1MzMeanFoundIso1)},
                {MS1_MZ_MEAN_FND_ISO2, QVariant(ms1MzMeanFoundIso2)},
                {MS1_MZ_MEAN_FND_100_PPM, QVariant(ms1MzMeanFound100PPM)},
                {MS1_MZ_MEAN_FND_45_PPM, QVariant(ms1MzMeanFound45PPM)},
                {MS1_MZ_MEAN_FND_20_PPM, QVariant(ms1MzMeanFound20PPM)},
                {MS1_MZ_MEAN_FND_PRE_MONO_PPM, QVariant(ms1MzMeanFoundPreMonoPPM)},
                {MS1_MZ_MEAN_FND_ISO_1_PPM, QVariant(ms1MzMeanFoundIso1PPM)},
                {MS1_MZ_MEAN_FND_ISO_2_PPM, QVariant(ms1MzMeanFoundIso2PPM)},
                {MS1_MZ_MEAN_FND_100_STD, QVariant(ms1MzStDevFound100)},
                {MS1_MZ_MEAN_FND_45_STD, QVariant(ms1MzStDevFound45)},
                {MS1_MZ_MEAN_FND_20_STD, QVariant(ms1MzStDevFound20)},
                {MS1_MZ_MEAN_FND_PRE_MONO_STD, QVariant(ms1MzStDevFoundPreMono)},
                {MS1_MZ_MEAN_FND_ISO_1_STD, QVariant(ms1MzStDevFoundIso1)},
                {MS1_MZ_MEAN_FND_ISO_2_STD, QVariant(ms1MzStDevFoundIso2)},
                {MS1_INTZ_FND_100, QVariant(ms1IntensityFound100)},
                {MS1_INTZ_FND_45, QVariant(ms1IntensityFound45)},
                {MS1_INTZ_FND_20, QVariant(ms1IntensityFound20)},
                {MS1_INTZ_FND_PRE_MONO, QVariant(ms1IntensityFoundPreMono)},
                {MS1_INTZ_FND_ISO_1, QVariant(ms1IntensityFoundIso1)},
                {MS1_INTZ_FND_ISO_2, QVariant(ms1IntensityFoundIso2)}
        };
    }

    static CandidateScoresReaderRow buildCandidateScoresReaderRow(const CandidateScores* candidateScores) {

        CandidateScoresReaderRow row;

        row.cosineSimSum100 = candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100],
        row.cosineSimSum100Greater80 = candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100GreaterThan80],
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
        row.charge = candidateScores->featuresArray[CandidateScores::Features::Charge],
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
        row.aminoAcidCountB = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountB],
        row.aminoAcidCountJ = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountJ],
        row.aminoAcidCountO = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountO],
        row.aminoAcidCountU = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountU],
        row.aminoAcidCountX = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountX],
        row.aminoAcidCountZ = candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountZ],
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
        row.peptideStringWithMods = candidateScores->isDecoy
                ? AminoAcids::mutatePenultimatePeptideResidues(candidateScores->targetDecoyCandidatePair->peptideStringWithMods())
                : candidateScores->targetDecoyCandidatePair->peptideStringWithMods();

        row.proteinGroup = candidateScores->proteinGroup;
        row.isDecoy = candidateScores->isDecoy;
        row.scanNumber = candidateScores->scanNumber;
        row.scanTime = candidateScores->scanTime;
        row.classifierScore = candidateScores->classifierScore;
        row.discriminantScore = candidateScores->discriminantScore;
        row.qValue = candidateScores->qValue;
        row.decoyRatio = candidateScores->decoyRatio;

        row.altTargetKeyIdCosineSimSumCharge1_OG = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG];
        row.altTargetKeyIdCosineSimSumCharge1_1 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_1];
        row.altTargetKeyIdCosineSimSumCharge1_2 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_2];
        row.altTargetKeyIdCosineSimSumCharge1_3 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_3];
        row.altTargetKeyIdCosineSimSumCharge2_OG = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG];
        row.altTargetKeyIdCosineSimSumCharge2_1 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_1];
        row.altTargetKeyIdCosineSimSumCharge2_2 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_2];
        row.altTargetKeyIdCosineSimSumCharge2_3 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_3];
        row.altTargetKeyIdCosineSimSumCharge3_OG = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG];
        row.altTargetKeyIdCosineSimSumCharge3_1 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_1];
        row.altTargetKeyIdCosineSimSumCharge3_2 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2];
        row.altTargetKeyIdCosineSimSumCharge3_3 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_3];
        row.altTargetKeyIdCosineSimSumCharge4_OG = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG];
        row.altTargetKeyIdCosineSimSumCharge4_1 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_1];
        row.altTargetKeyIdCosineSimSumCharge4_2 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_2];
        row.altTargetKeyIdCosineSimSumCharge4_3 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_3];

        row.altTargetKeyIdTimeDeltaCharge1_1 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_1];
        row.altTargetKeyIdTimeDeltaCharge1_2 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_2];
        row.altTargetKeyIdTimeDeltaCharge1_3 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_3];
        row.altTargetKeyIdTimeDeltaCharge2_1 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_1];
        row.altTargetKeyIdTimeDeltaCharge2_2 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_2];
        row.altTargetKeyIdTimeDeltaCharge2_3 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_3];
        row.altTargetKeyIdTimeDeltaCharge3_1 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_1];
        row.altTargetKeyIdTimeDeltaCharge3_2 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_2];
        row.altTargetKeyIdTimeDeltaCharge3_3 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_3];
        row.altTargetKeyIdTimeDeltaCharge4_1 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_1];
        row.altTargetKeyIdTimeDeltaCharge4_2 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_2];
        row.altTargetKeyIdTimeDeltaCharge4_3 = candidateScores->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_3];

        row.ms1MzMeanFound100 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound100];
        row.ms1MzMeanFound45 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound45];
        row.ms1MzMeanFound20 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound20];
        row.ms1MzMeanFoundPreMono = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMono];
        row.ms1MzMeanFoundIso1 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1];
        row.ms1MzMeanFoundIso2 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso2];
        row.ms1MzMeanFound100PPM = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound100PPM];
        row.ms1MzMeanFound45PPM = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound45PPM];
        row.ms1MzMeanFound20PPM = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound20PPM];
        row.ms1MzMeanFoundPreMonoPPM = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMonoPPM];
        row.ms1MzMeanFoundIso1PPM = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1PPM];
        row.ms1MzMeanFoundIso2PPM = candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso2PPM];
        row.ms1MzStDevFound100 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound100];
        row.ms1MzStDevFound45 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound45];
        row.ms1MzStDevFound20 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound20];
        row.ms1MzStDevFoundPreMono = candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundPreMono];
        row.ms1MzStDevFoundIso1 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso1];
        row.ms1MzStDevFoundIso2 = candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso2];
        row.ms1IntensityFound100 = candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound100];
        row.ms1IntensityFound45 = candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound45];
        row.ms1IntensityFound20 = candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound20];
        row.ms1IntensityFoundPreMono = candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundPreMono];
        row.ms1IntensityFoundIso1 = candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundIso1];
        row.ms1IntensityFoundIso2 = candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundIso2];

        return row;
    }

    static QVector<float> featuresArrayFromCandidateScoresReaderRow(const CandidateScoresReaderRow &candidateScoresReaderRow) {

        QVector<float> featuresArray(CandidateScores::Features::FeaturesSize);

        featuresArray[CandidateScores::Features::CosineSimSum100] = candidateScoresReaderRow.cosineSimSum100;
        featuresArray[CandidateScores::Features::CosineSimSum100GreaterThan80] = candidateScoresReaderRow.cosineSimSum100Greater80;
        featuresArray[CandidateScores::Features::AllignedMaxIndexesCount] = candidateScoresReaderRow.allignedMaxIndexesCount;
        featuresArray[CandidateScores::Features::CosineSim100MS1] = candidateScoresReaderRow.cosineSim100MS1;
        featuresArray[CandidateScores::Features::CosineSimSpectrumCubed] = candidateScoresReaderRow.cosineSimSpectrumCubed;
        featuresArray[CandidateScores::Features::KlDivSpectrumCubeRoot] = candidateScoresReaderRow.klDivSpectrumCubeRoot;
        featuresArray[CandidateScores::Features::CosineSim45MS1] = candidateScoresReaderRow.cosineSimSum45;
        featuresArray[CandidateScores::Features::CosineSim20MS1] = candidateScoresReaderRow.cosineSimSum20;
        featuresArray[CandidateScores::Features::CosineSimSumTop6] = candidateScoresReaderRow.cosineSimSumTop6;
        featuresArray[CandidateScores::Features::CosineSimSumBottom6] = candidateScoresReaderRow.cosineSimSumBottom6;
        featuresArray[CandidateScores::Features::TopBottomRatio] = candidateScoresReaderRow.topBottomRatio;
        featuresArray[CandidateScores::Features::TopBottomRatioNorm] = candidateScoresReaderRow.topBottomRatioNorm;
        featuresArray[CandidateScores::Features::Charge] = candidateScoresReaderRow.charge;
        featuresArray[CandidateScores::Features::Mass] = candidateScoresReaderRow.mass;
        featuresArray[CandidateScores::Features::ScanTimeDelta] = candidateScoresReaderRow.scanTimeDelta;
        featuresArray[CandidateScores::Features::ScanTimeRange] = candidateScoresReaderRow.scanTimeRange;
        featuresArray[CandidateScores::Features::ScanTimePd] = candidateScoresReaderRow.scanTimePd;
        featuresArray[CandidateScores::Features::ScanIonCount] = candidateScoresReaderRow.scanIonCount;
        featuresArray[CandidateScores::Features::MzNorm] = candidateScoresReaderRow.mzNorm;
        featuresArray[CandidateScores::Features::KlDivSpectrum] = candidateScoresReaderRow.klDivSpectrum;
        featuresArray[CandidateScores::Features::CosineSimSpectrum] = candidateScoresReaderRow.cosineSimSpectrum;
        featuresArray[CandidateScores::Features::CosineSim45MS1] = candidateScoresReaderRow.cosineSim45MS1;
        featuresArray[CandidateScores::Features::CosineSim20MS1] = candidateScoresReaderRow.cosineSim20MS1;
        featuresArray[CandidateScores::Features::CosineSim100MS1PreMono] = candidateScoresReaderRow.cosineSim100MS1PreMono;
        featuresArray[CandidateScores::Features::CosineSim100MS1Iso1] = candidateScoresReaderRow.cosineSim100MS1Iso1;
        featuresArray[CandidateScores::Features::CosineSim100MS1Iso2] = candidateScoresReaderRow.cosineSim100MS1Iso2;
        featuresArray[CandidateScores::Features::PeptideLengthNorm] = candidateScoresReaderRow.peptideLengthNorm;
        featuresArray[CandidateScores::Features::ScanTimePredicted] = candidateScoresReaderRow.scanTimePredicted;
        featuresArray[CandidateScores::Features::TheoFragmentCount] = candidateScoresReaderRow.theoFragmentCount;
        featuresArray[CandidateScores::Features::TotalIntensityLog] = candidateScoresReaderRow.totalIntensityLog;
        featuresArray[CandidateScores::Features::PeakShapeRatio1] = candidateScoresReaderRow.peakShapeRatio1;
        featuresArray[CandidateScores::Features::PeakShapeRatio2] = candidateScoresReaderRow.peakShapeRatio2;
        featuresArray[CandidateScores::Features::PeakShapeRatio3] = candidateScoresReaderRow.peakShapeRatio3;
        featuresArray[CandidateScores::Features::ShadowsCosineSimSum] = candidateScoresReaderRow.shadowsCosineSimSum;
        featuresArray[CandidateScores::Features::IRTPredicted] = candidateScoresReaderRow.iRtPredicted;
        featuresArray[CandidateScores::Features::CosineSimToAnchor1] = candidateScoresReaderRow.cosineSimToAnchor1;
        featuresArray[CandidateScores::Features::CosineSimToAnchor2] = candidateScoresReaderRow.cosineSimToAnchor2;
        featuresArray[CandidateScores::Features::CosineSimToAnchor3] = candidateScoresReaderRow.cosineSimToAnchor3;
        featuresArray[CandidateScores::Features::CosineSimToAnchor4] = candidateScoresReaderRow.cosineSimToAnchor4;
        featuresArray[CandidateScores::Features::CosineSimToAnchor5] = candidateScoresReaderRow.cosineSimToAnchor5;
        featuresArray[CandidateScores::Features::CosineSimToAnchor6] = candidateScoresReaderRow.cosineSimToAnchor6;
        featuresArray[CandidateScores::Features::CosineSimToAnchor7] = candidateScoresReaderRow.cosineSimToAnchor7;
        featuresArray[CandidateScores::Features::CosineSimToAnchor8] = candidateScoresReaderRow.cosineSimToAnchor8;
        featuresArray[CandidateScores::Features::CosineSimToAnchor9] = candidateScoresReaderRow.cosineSimToAnchor9;
        featuresArray[CandidateScores::Features::CosineSimToAnchor10] = candidateScoresReaderRow.cosineSimToAnchor10;
        featuresArray[CandidateScores::Features::CosineSimToAnchor11] = candidateScoresReaderRow.cosineSimToAnchor11;
        featuresArray[CandidateScores::Features::CosineSimToAnchor12] = candidateScoresReaderRow.cosineSimToAnchor12;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor1] = candidateScoresReaderRow.cosineSimShadowsToAnchor1;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor2] = candidateScoresReaderRow.cosineSimShadowsToAnchor2;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor3] = candidateScoresReaderRow.cosineSimShadowsToAnchor3;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor4] = candidateScoresReaderRow.cosineSimShadowsToAnchor4;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor5] = candidateScoresReaderRow.cosineSimShadowsToAnchor5;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor6] = candidateScoresReaderRow.cosineSimShadowsToAnchor6;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor7] = candidateScoresReaderRow.cosineSimShadowsToAnchor7;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor8] = candidateScoresReaderRow.cosineSimShadowsToAnchor8;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor9] = candidateScoresReaderRow.cosineSimShadowsToAnchor9;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor10] = candidateScoresReaderRow.cosineSimShadowsToAnchor10;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor11] = candidateScoresReaderRow.cosineSimShadowsToAnchor11;
        featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor12] = candidateScoresReaderRow.cosineSimShadowsToAnchor12;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio1] = candidateScoresReaderRow.shadowsIntensityRatio1;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio2] = candidateScoresReaderRow.shadowsIntensityRatio2;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio3] = candidateScoresReaderRow.shadowsIntensityRatio3;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio4] = candidateScoresReaderRow.shadowsIntensityRatio4;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio5] = candidateScoresReaderRow.shadowsIntensityRatio5;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio6] = candidateScoresReaderRow.shadowsIntensityRatio6;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio7] = candidateScoresReaderRow.shadowsIntensityRatio7;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio8] = candidateScoresReaderRow.shadowsIntensityRatio8;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio9] = candidateScoresReaderRow.shadowsIntensityRatio9;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio10] = candidateScoresReaderRow.shadowsIntensityRatio10;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio11] = candidateScoresReaderRow.shadowsIntensityRatio11;
        featuresArray[CandidateScores::Features::ShadowsIntensityRatio12] = candidateScoresReaderRow.shadowsIntensityRatio12;
        featuresArray[CandidateScores::Features::MzSearched1] = candidateScoresReaderRow.mzSearched1;
        featuresArray[CandidateScores::Features::MzSearched2] = candidateScoresReaderRow.mzSearched2;
        featuresArray[CandidateScores::Features::MzSearched3] = candidateScoresReaderRow.mzSearched3;
        featuresArray[CandidateScores::Features::MzSearched4] = candidateScoresReaderRow.mzSearched4;
        featuresArray[CandidateScores::Features::MzSearched5] = candidateScoresReaderRow.mzSearched5;
        featuresArray[CandidateScores::Features::MzSearched6] = candidateScoresReaderRow.mzSearched6;
        featuresArray[CandidateScores::Features::MzSearched7] = candidateScoresReaderRow.mzSearched7;
        featuresArray[CandidateScores::Features::MzSearched8] = candidateScoresReaderRow.mzSearched8;
        featuresArray[CandidateScores::Features::MzSearched9] = candidateScoresReaderRow.mzSearched9;
        featuresArray[CandidateScores::Features::MzSearched10] = candidateScoresReaderRow.mzSearched10;
        featuresArray[CandidateScores::Features::MzSearched11] = candidateScoresReaderRow.mzSearched11;
        featuresArray[CandidateScores::Features::MzSearched12] = candidateScoresReaderRow.mzSearched12;
        featuresArray[CandidateScores::Features::TheoIntensity1] = candidateScoresReaderRow.theoIntensity1;
        featuresArray[CandidateScores::Features::TheoIntensity2] = candidateScoresReaderRow.theoIntensity2;
        featuresArray[CandidateScores::Features::TheoIntensity3] = candidateScoresReaderRow.theoIntensity3;
        featuresArray[CandidateScores::Features::TheoIntensity4] = candidateScoresReaderRow.theoIntensity4;
        featuresArray[CandidateScores::Features::TheoIntensity5] = candidateScoresReaderRow.theoIntensity5;
        featuresArray[CandidateScores::Features::TheoIntensity6] = candidateScoresReaderRow.theoIntensity6;
        featuresArray[CandidateScores::Features::TheoIntensity7] = candidateScoresReaderRow.theoIntensity7;
        featuresArray[CandidateScores::Features::TheoIntensity8] = candidateScoresReaderRow.theoIntensity8;
        featuresArray[CandidateScores::Features::TheoIntensity9] = candidateScoresReaderRow.theoIntensity9;
        featuresArray[CandidateScores::Features::TheoIntensity10] = candidateScoresReaderRow.theoIntensity10;
        featuresArray[CandidateScores::Features::TheoIntensity11] = candidateScoresReaderRow.theoIntensity11;
        featuresArray[CandidateScores::Features::TheoIntensity12] = candidateScoresReaderRow.theoIntensity12;
        featuresArray[CandidateScores::Features::MzFoundMean1] = candidateScoresReaderRow.mzFoundMean1;
        featuresArray[CandidateScores::Features::MzFoundMean2] = candidateScoresReaderRow.mzFoundMean2;
        featuresArray[CandidateScores::Features::MzFoundMean3] = candidateScoresReaderRow.mzFoundMean3;
        featuresArray[CandidateScores::Features::MzFoundMean4] = candidateScoresReaderRow.mzFoundMean4;
        featuresArray[CandidateScores::Features::MzFoundMean5] = candidateScoresReaderRow.mzFoundMean5;
        featuresArray[CandidateScores::Features::MzFoundMean6] = candidateScoresReaderRow.mzFoundMean6;
        featuresArray[CandidateScores::Features::MzFoundMean7] = candidateScoresReaderRow.mzFoundMean7;
        featuresArray[CandidateScores::Features::MzFoundMean8] = candidateScoresReaderRow.mzFoundMean8;
        featuresArray[CandidateScores::Features::MzFoundMean9] = candidateScoresReaderRow.mzFoundMean9;
        featuresArray[CandidateScores::Features::MzFoundMean10] = candidateScoresReaderRow.mzFoundMean10;
        featuresArray[CandidateScores::Features::MzFoundMean11] = candidateScoresReaderRow.mzFoundMean11;
        featuresArray[CandidateScores::Features::MzFoundMean12] = candidateScoresReaderRow.mzFoundMean12;
        featuresArray[CandidateScores::Features::IntensityFoundMax1] = candidateScoresReaderRow.intensityFoundMax1;
        featuresArray[CandidateScores::Features::IntensityFoundMax2] = candidateScoresReaderRow.intensityFoundMax2;
        featuresArray[CandidateScores::Features::IntensityFoundMax3] = candidateScoresReaderRow.intensityFoundMax3;
        featuresArray[CandidateScores::Features::IntensityFoundMax4] = candidateScoresReaderRow.intensityFoundMax4;
        featuresArray[CandidateScores::Features::IntensityFoundMax5] = candidateScoresReaderRow.intensityFoundMax5;
        featuresArray[CandidateScores::Features::IntensityFoundMax6] = candidateScoresReaderRow.intensityFoundMax6;
        featuresArray[CandidateScores::Features::IntensityFoundMax7] = candidateScoresReaderRow.intensityFoundMax7;
        featuresArray[CandidateScores::Features::IntensityFoundMax8] = candidateScoresReaderRow.intensityFoundMax8;
        featuresArray[CandidateScores::Features::IntensityFoundMax9] = candidateScoresReaderRow.intensityFoundMax9;
        featuresArray[CandidateScores::Features::IntensityFoundMax10] = candidateScoresReaderRow.intensityFoundMax10;
        featuresArray[CandidateScores::Features::IntensityFoundMax11] = candidateScoresReaderRow.intensityFoundMax11;
        featuresArray[CandidateScores::Features::IntensityFoundMax12] = candidateScoresReaderRow.intensityFoundMax12;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm1] = candidateScoresReaderRow.mzPeakLengthsNorm1;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm2] = candidateScoresReaderRow.mzPeakLengthsNorm2;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm3] = candidateScoresReaderRow.mzPeakLengthsNorm3;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm4] = candidateScoresReaderRow.mzPeakLengthsNorm4;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm5] = candidateScoresReaderRow.mzPeakLengthsNorm5;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm6] = candidateScoresReaderRow.mzPeakLengthsNorm6;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm7] = candidateScoresReaderRow.mzPeakLengthsNorm7;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm8] = candidateScoresReaderRow.mzPeakLengthsNorm8;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm9] = candidateScoresReaderRow.mzPeakLengthsNorm9;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm10] = candidateScoresReaderRow.mzPeakLengthsNorm10;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm11] = candidateScoresReaderRow.mzPeakLengthsNorm11;
        featuresArray[CandidateScores::Features::MzPeakLengthsNorm12] = candidateScoresReaderRow.mzPeakLengthsNorm12;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor1] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor1;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor2] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor2;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor3] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor3;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor4] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor4;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor5] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor5;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor6] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor6;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor7] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor7;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor8] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor8;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor9] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor9;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor10] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor10;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor11] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor11;
        featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor12] = candidateScoresReaderRow.columnApexIndexRatiosToAnchor12;
        featuresArray[CandidateScores::Features::AminoAcidCountA] = candidateScoresReaderRow.aminoAcidCountA;
        featuresArray[CandidateScores::Features::AminoAcidCountC] = candidateScoresReaderRow.aminoAcidCountC;
        featuresArray[CandidateScores::Features::AminoAcidCountD] = candidateScoresReaderRow.aminoAcidCountD;
        featuresArray[CandidateScores::Features::AminoAcidCountE] = candidateScoresReaderRow.aminoAcidCountE;
        featuresArray[CandidateScores::Features::AminoAcidCountF] = candidateScoresReaderRow.aminoAcidCountF;
        featuresArray[CandidateScores::Features::AminoAcidCountG] = candidateScoresReaderRow.aminoAcidCountG;
        featuresArray[CandidateScores::Features::AminoAcidCountH] = candidateScoresReaderRow.aminoAcidCountH;
        featuresArray[CandidateScores::Features::AminoAcidCountI] = candidateScoresReaderRow.aminoAcidCountI;
        featuresArray[CandidateScores::Features::AminoAcidCountK] = candidateScoresReaderRow.aminoAcidCountK;
        featuresArray[CandidateScores::Features::AminoAcidCountL] = candidateScoresReaderRow.aminoAcidCountL;
        featuresArray[CandidateScores::Features::AminoAcidCountM] = candidateScoresReaderRow.aminoAcidCountM;
        featuresArray[CandidateScores::Features::AminoAcidCountN] = candidateScoresReaderRow.aminoAcidCountN;
        featuresArray[CandidateScores::Features::AminoAcidCountP] = candidateScoresReaderRow.aminoAcidCountP;
        featuresArray[CandidateScores::Features::AminoAcidCountQ] = candidateScoresReaderRow.aminoAcidCountQ;
        featuresArray[CandidateScores::Features::AminoAcidCountR] = candidateScoresReaderRow.aminoAcidCountR;
        featuresArray[CandidateScores::Features::AminoAcidCountS] = candidateScoresReaderRow.aminoAcidCountS;
        featuresArray[CandidateScores::Features::AminoAcidCountT] = candidateScoresReaderRow.aminoAcidCountT;
        featuresArray[CandidateScores::Features::AminoAcidCountV] = candidateScoresReaderRow.aminoAcidCountV;
        featuresArray[CandidateScores::Features::AminoAcidCountW] = candidateScoresReaderRow.aminoAcidCountW;
        featuresArray[CandidateScores::Features::AminoAcidCountY] = candidateScoresReaderRow.aminoAcidCountY;
        featuresArray[CandidateScores::Features::AminoAcidCountB] = candidateScoresReaderRow.aminoAcidCountB;
        featuresArray[CandidateScores::Features::AminoAcidCountJ] = candidateScoresReaderRow.aminoAcidCountJ;
        featuresArray[CandidateScores::Features::AminoAcidCountO] = candidateScoresReaderRow.aminoAcidCountO;
        featuresArray[CandidateScores::Features::AminoAcidCountU] = candidateScoresReaderRow.aminoAcidCountU;
        featuresArray[CandidateScores::Features::AminoAcidCountX] = candidateScoresReaderRow.aminoAcidCountX;
        featuresArray[CandidateScores::Features::AminoAcidCountZ] = candidateScoresReaderRow.aminoAcidCountZ;
        featuresArray[CandidateScores::Features::MzFoundStDev1] = candidateScoresReaderRow.mzFoundStDev1;
        featuresArray[CandidateScores::Features::MzFoundStDev2] = candidateScoresReaderRow.mzFoundStDev2;
        featuresArray[CandidateScores::Features::MzFoundStDev3] = candidateScoresReaderRow.mzFoundStDev3;
        featuresArray[CandidateScores::Features::MzFoundStDev4] = candidateScoresReaderRow.mzFoundStDev4;
        featuresArray[CandidateScores::Features::MzFoundStDev5] = candidateScoresReaderRow.mzFoundStDev5;
        featuresArray[CandidateScores::Features::MzFoundStDev6] = candidateScoresReaderRow.mzFoundStDev6;
        featuresArray[CandidateScores::Features::MzFoundStDev7] = candidateScoresReaderRow.mzFoundStDev7;
        featuresArray[CandidateScores::Features::MzFoundStDev8] = candidateScoresReaderRow.mzFoundStDev8;
        featuresArray[CandidateScores::Features::MzFoundStDev9] = candidateScoresReaderRow.mzFoundStDev9;
        featuresArray[CandidateScores::Features::MzFoundStDev10] = candidateScoresReaderRow.mzFoundStDev10;
        featuresArray[CandidateScores::Features::MzFoundStDev11] = candidateScoresReaderRow.mzFoundStDev11;
        featuresArray[CandidateScores::Features::MzFoundStDev12] = candidateScoresReaderRow.mzFoundStDev12;
        featuresArray[CandidateScores::Features::MzAccuracy1] = candidateScoresReaderRow.mzAccuracy1;
        featuresArray[CandidateScores::Features::MzAccuracy2] = candidateScoresReaderRow.mzAccuracy2;
        featuresArray[CandidateScores::Features::MzAccuracy3] = candidateScoresReaderRow.mzAccuracy3;
        featuresArray[CandidateScores::Features::MzAccuracy4] = candidateScoresReaderRow.mzAccuracy4;
        featuresArray[CandidateScores::Features::MzAccuracy5] = candidateScoresReaderRow.mzAccuracy5;
        featuresArray[CandidateScores::Features::MzAccuracy6] = candidateScoresReaderRow.mzAccuracy6;
        featuresArray[CandidateScores::Features::MzAccuracy7] = candidateScoresReaderRow.mzAccuracy7;
        featuresArray[CandidateScores::Features::MzAccuracy8] = candidateScoresReaderRow.mzAccuracy8;
        featuresArray[CandidateScores::Features::MzAccuracy9] = candidateScoresReaderRow.mzAccuracy9;
        featuresArray[CandidateScores::Features::MzAccuracy10] = candidateScoresReaderRow.mzAccuracy10;
        featuresArray[CandidateScores::Features::MzAccuracy11] = candidateScoresReaderRow.mzAccuracy11;
        featuresArray[CandidateScores::Features::MzAccuracy12] = candidateScoresReaderRow.mzAccuracy12;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge1_OG;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_1] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge1_1;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_2] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge1_2;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_3] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge1_3;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge2_OG;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_1] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge2_1;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_2] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge2_2;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_3] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge2_3;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge3_OG;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_1] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge3_1;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge3_2;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_3] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge3_3;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge4_OG;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_1] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge4_1;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge4_2;
        featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_3] = candidateScoresReaderRow.altTargetKeyIdCosineSimSumCharge4_3;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge1_1;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge1_2;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge1_3;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge2_1;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_2] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge2_2;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_3] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge2_3;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge3_1;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_2] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge3_2;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_3] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge3_3;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge4_1;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_2] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge4_2;
        featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_3] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge4_3;
        featuresArray[CandidateScores::Features::Ms1MzMeanFound100] = candidateScoresReaderRow.ms1MzMeanFound100;
        featuresArray[CandidateScores::Features::Ms1MzMeanFound45] = candidateScoresReaderRow.ms1MzMeanFound45;
        featuresArray[CandidateScores::Features::Ms1MzMeanFound20] = candidateScoresReaderRow.ms1MzMeanFound20;
        featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMono] = candidateScoresReaderRow.ms1MzMeanFoundPreMono;
        featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1] = candidateScoresReaderRow.ms1MzMeanFoundIso1;
        featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso2] = candidateScoresReaderRow.ms1MzMeanFoundIso2;
        featuresArray[CandidateScores::Features::Ms1MzMeanFound100PPM] = candidateScoresReaderRow.ms1MzMeanFound100PPM;
        featuresArray[CandidateScores::Features::Ms1MzMeanFound45PPM] = candidateScoresReaderRow.ms1MzMeanFound45PPM;
        featuresArray[CandidateScores::Features::Ms1MzMeanFound20PPM] = candidateScoresReaderRow.ms1MzMeanFound20PPM;
        featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMonoPPM] = candidateScoresReaderRow.ms1MzMeanFoundPreMonoPPM;
        featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1PPM] = candidateScoresReaderRow.ms1MzMeanFoundIso1PPM;
        featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso2PPM] = candidateScoresReaderRow.ms1MzMeanFoundIso2PPM;
        featuresArray[CandidateScores::Features::Ms1MzStDevFound100] = candidateScoresReaderRow.ms1MzStDevFound100;
        featuresArray[CandidateScores::Features::Ms1MzStDevFound45] = candidateScoresReaderRow.ms1MzStDevFound45;
        featuresArray[CandidateScores::Features::Ms1MzStDevFound20] = candidateScoresReaderRow.ms1MzStDevFound20;
        featuresArray[CandidateScores::Features::Ms1MzStDevFoundPreMono] = candidateScoresReaderRow.ms1MzStDevFoundPreMono;
        featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso1] = candidateScoresReaderRow.ms1MzStDevFoundIso1;
        featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso2] = candidateScoresReaderRow.ms1MzStDevFoundIso2;
        featuresArray[CandidateScores::Features::Ms1IntensityFound100] = candidateScoresReaderRow.ms1IntensityFound100;
        featuresArray[CandidateScores::Features::Ms1IntensityFound45] = candidateScoresReaderRow.ms1IntensityFound45;
        featuresArray[CandidateScores::Features::Ms1IntensityFound20] = candidateScoresReaderRow.ms1IntensityFound20;
        featuresArray[CandidateScores::Features::Ms1IntensityFoundPreMono] = candidateScoresReaderRow.ms1IntensityFoundPreMono;
        featuresArray[CandidateScores::Features::Ms1IntensityFoundIso1] = candidateScoresReaderRow.ms1IntensityFoundIso1;
        featuresArray[CandidateScores::Features::Ms1IntensityFoundIso2] = candidateScoresReaderRow.ms1IntensityFoundIso2;

        return featuresArray;
    }

};

#endif //PYTHIADIACPP_CANDIDATESCORES_H
