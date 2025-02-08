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

enum Features {
    CosineSimSum100 = 0,
    CosineSimSum100Top12,
    CosineSimSum100GreaterThan80,
    AllignedMaxIndexesCount,
    CosineSim100MS1, //5
    CosineSimSpectrumCubed,
    KlDivSpectrumCubeRoot,
    CosineSimSum45,
    CosineSimSumTop,
    CosineSimSumBottom, //10
    TopBottomRatio,
    TopBottomRatioNorm,
    Charge,
    ScanTimeDelta,
    ScanTimeDeltaAbs, //15
    ScanTimePd,
    ScanTimePdAbs,
    ScanIonCount,
    MzNorm,
    Mass, //20
    KlDivSpectrum,
    CosineSimSpectrum,
    CosineSim45MS1,//20
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
    IntensityFoundMaxNorm1,
    IntensityFoundMaxNorm2,
    IntensityFoundMaxNorm3,
    IntensityFoundMaxNorm4,
    IntensityFoundMaxNorm5,//110
    IntensityFoundMaxNorm6,
    IntensityFoundMaxNorm7,
    IntensityFoundMaxNorm8,
    IntensityFoundMaxNorm9,
    IntensityFoundMaxNorm10,
    IntensityFoundMaxNorm11,
    IntensityFoundMaxNorm12,
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
    AltTargetKeyIdDiscScoreChargeOG_alt,
    AltTargetKeyIdDiscScoreCharge1_1,
    AltTargetKeyIdDiscScoreCharge1_2,
    AltTargetKeyIdDiscScoreCharge2_1,
    AltTargetKeyIdDiscScoreCharge2_2,
    AltTargetKeyIdDiscScoreCharge3_1,
    AltTargetKeyIdDiscScoreCharge3_2,
    AltTargetKeyIdDiscScoreCharge4_1,
    AltTargetKeyIdDiscScoreCharge4_2,
    AltTargetKeyIdTimeDeltaCharge1_1,
    AltTargetKeyIdTimeDeltaCharge2_1,
    AltTargetKeyIdTimeDeltaCharge3_1,
    AltTargetKeyIdTimeDeltaCharge4_1,
    Ms1MzMeanFound100,//220
    Ms1MzMeanFound45,
    Ms1MzMeanFoundPreMono,
    Ms1MzMeanFoundIso1,
    Ms1MzMeanFoundIso2,
    Ms1MzMeanFound100PPM,
    Ms1MzMeanFound45PPM,
    Ms1MzMeanFoundPreMonoPPM,
    Ms1MzMeanFoundIso1PPM,//230
    Ms1MzMeanFoundIso2PPM,
    Ms1MzStDevFound100,
    Ms1MzStDevFound45,
    Ms1MzStDevFoundPreMono,
    Ms1MzStDevFoundIso1,
    Ms1MzStDevFoundIso2,
    Ms1IntensityFound100,
    Ms1IntensityFoundApex100,

    Ms1IntensityFoundApex100IM,


    Ms1IntensityFound45,
    Ms1IntensityFoundPreMono,
    Ms1IntensityFoundIso1,
    Ms1IntensityFoundIso2,
    CosineSimSpectrumOverTime,
    CosineSimSpectrumOverTimeCubed,
    CosineSimSpectrumStDev,
    CosineSimSum100MS1,
    MS1Averagine,
    CosineSimSum100Window1p5X,//250
    CosineSimSum100Window2X,
    TotalIntensityPeakHeights,
    TotalIntensityRaw,
    TargetWindowLocation,
    TargetWindowLocationAbs,
    DiscriminantScore,
    AlignmentIndexMean,
    AlignmentIndexStDev,
    AlignmentCombinedScore,
    MatrixZeroPercentage,
    MzPPMMeanAbs,
    MzPPMMean,
    MzPPMStd,
    FoundB,
    FoundY,
    FoundPercent,
    DiscScore1stRunnerUp,
    DiscScore2ndRunnerUp,
    DiscScoresCount,
    DiscScoresMean,
    DiscScoresStDev,
    YIonSeriesMax,
    YIonSeriesCount,
    YIonSeriesMean,
    YIonSeriesStd,
    YIonSeriesTheoMax,
    YIonSeriesTheoCount,
    YIonSeriesTheoMean,
    YIonSeriesTheoStd,
    YIonSeriesMaxFoundToTheoFraction,
    YIonSeriesCountRatio,
    BIonSeriesMax,
    BIonSeriesCount,
    BIonSeriesMean,
    BIonSeriesStd,
    BIonSeriesTheoMax,
    BIonSeriesTheoCount,
    BIonSeriesTheoMean,
    BIonSeriesTheoStd,
    BIonSeriesMaxFoundToTheoFraction,
    BIonSeriesCountRatio,
    CosineSimFullTheo,
    IonsFoundFractionFull,
    CosineSimFullTheoXIonsFoundFractionFull,

    // MzFoundMean1PPM,
    // MzFoundMean2PPM,
    // MzFoundMean3PPM,
    // MzFoundMean4PPM,
    // MzFoundMean5PPM,
    // MzFoundMean6PPM,
    // MzFoundMean7PPM,//100
    // MzFoundMean8PPM,
    // MzFoundMean9PPM,
    // MzFoundMean10PPM,
    // MzFoundMean11PPM,
    // MzFoundMean12PPM,
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
    // CosineSimShadowsToAnchor1,
    // CosineSimShadowsToAnchor2,
    // CosineSimShadowsToAnchor3,
    // CosineSimShadowsToAnchor4,
    // CosineSimShadowsToAnchor5,//50
    // CosineSimShadowsToAnchor6,
    // CosineSimShadowsToAnchor7,
    // CosineSimShadowsToAnchor8,
    // CosineSimShadowsToAnchor9,
    // CosineSimShadowsToAnchor10,
    // CosineSimShadowsToAnchor11,
    // CosineSimShadowsToAnchor12,
    // MzPeakLengthsNorm1,
    // MzPeakLengthsNorm2,
    // MzPeakLengthsNorm3,//120
    // MzPeakLengthsNorm4,
    // MzPeakLengthsNorm5,
    // MzPeakLengthsNorm6,
    // MzPeakLengthsNorm7,
    // MzPeakLengthsNorm8,
    // MzPeakLengthsNorm9,
    // MzPeakLengthsNorm10,
    // MzPeakLengthsNorm11,
    // MzPeakLengthsNorm12,
    // MzFoundStDev7,
    // MzFoundStDev8,
    // MzFoundStDev9,
    // MzFoundStDev10,
    // MzFoundStDev11,
    // MzFoundStDev12,


    FeaturesSize
    };


/**
 * @brief Class representing scores for a candidate.
 *
 * This class defines a set of features and provides storage for scores related to a candidate.
 * It includes various features such as cosine similarity, peak intensities, charge, mass, scan time, etc.
 *
 */
class ALGORITHMSFFLIB_EXPORTS CandidateScores {

public:

    CandidateScores() = default;
    ~CandidateScores() = default;

    TargetDecoyCandidatePair *targetDecoyCandidatePair = nullptr;

    MzTargetKey targetKey;
    QString proteinGroup;

    bool isDecoy = false;
    FrameIndex frameIndex = -1;
    FrameIndex frameIndexStart = -1;
    FrameIndex frameIndexEnd = -1;
    ScanNumber scanNumber = -1;
    ScanNumber scanNumberStart = -1;
    ScanNumber scanNumberEnd = -1;
    ScanTime scanTime = -1.0;
    ScanTime scanTimeStart = -1.0;
    ScanTime scanTimeEnd = -1.0;
    ScanTime scanTimePredicted = -1.0;

    IonMobilityIndex ionMobilityIndex = -1;
    IonMobilityIndex ionMobilityIndexStart = -1;
    IonMobilityIndex ionMobilityIndexEnd = -1;
    float imDriftTime = -1.0;

    double classifierScore = -1.0;
    double discriminantScore = -1.0;
    double qValue = 1.0;
    double decoyRatio = -1.0;

    QVector<float> featuresArray;
    QVector<float> integrations;

    PeptideSequenceWithModsChargeAndTargetKey peptideSequenceWithModsChargeAndTargetKey;
    CandidateScores *converseCandidateScore;


    /**
    * @brief Initializes the features array with default values.
    */
    void initFeaturesArray();

    QVector<float> selectFeaturesArrayFeatures(const QVector<Features> &enumFeatures);

    static QVector<float> selectFeaturesArrayFeatures(
            const QVector<float> &featureArray,
            const QVector<Features> &enumFeatures
            );

    void printFeatures(const QVector<Features> &featuresToPrint);

};


namespace CandidateScoresReaderRowNamespace {

    const QString COS_SIM_SUM_100 = QStringLiteral("CosineSimSum100");
    const QString COS_SIM_SUM_100_GREATER_80 = QStringLiteral("CosineSimSum100GreaterThan80");
    const QString ALL_MAX_IND_CNT = QStringLiteral("AllignedMaxIndexesCount");
    const QString COS_SIM_SUM_MS1_100 = QStringLiteral("CosineSim100MS1");
    const QString COS_SIM_SPEC_CUBED = QStringLiteral("CosineSimSpectrumCubed");
    const QString KL_DIV_SPEC_CUBE_RT = QStringLiteral("KlDivSpectrumCubeRoot");
    const QString COS_SIM_SUM_45 = QStringLiteral("CosineSimSum45");
    const QString COS_SIM_SUM_TOP_6 = QStringLiteral("CosineSimSumTop");
    const QString COS_SIM_SUM_BOTTOM_6 = QStringLiteral("CosineSimSumBottom");
    const QString TOP_BOTTOM_RATIO = QStringLiteral("TopBottomRatio");
    const QString TOP_BOTTOME_RATIO_NORM = QStringLiteral("TopBottomRatioNorm");
    const QString CHARGE = QStringLiteral("Charge");
    const QString MASS = QStringLiteral("Mass");
    const QString SCAN_TIME_DELTA = QStringLiteral("ScanTimeDelta");
    const QString SCAN_TIME_PD = QStringLiteral("ScanTimePd");
    const QString SCAN_ION_CNT = QStringLiteral("ScanIonCount");
    const QString MZ_NORM = QStringLiteral("MzNorm");
    const QString KL_DIV_SPEC = QStringLiteral("KlDivSpectrum");
    const QString COSINE_SIM_SPEC = QStringLiteral("CosineSimSpectrum");
    const QString COSINE_SIM_SUM_MS1_45 = QStringLiteral("CosineSim45MS1");
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
    // const QString COS_SIM_ANCH_SHADOW_1 = QStringLiteral("CosineSimShadowsToAnchor1");
    // const QString COS_SIM_ANCH_SHADOW_2 = QStringLiteral("CosineSimShadowsToAnchor2");
    // const QString COS_SIM_ANCH_SHADOW_3 = QStringLiteral("CosineSimShadowsToAnchor3");
    // const QString COS_SIM_ANCH_SHADOW_4 = QStringLiteral("CosineSimShadowsToAnchor4");
    // const QString COS_SIM_ANCH_SHADOW_5 = QStringLiteral("CosineSimShadowsToAnchor5");
    // const QString COS_SIM_ANCH_SHADOW_6 = QStringLiteral("CosineSimShadowsToAnchor6");
    // const QString COS_SIM_ANCH_SHADOW_7 = QStringLiteral("CosineSimShadowsToAnchor7");
    // const QString COS_SIM_ANCH_SHADOW_8 = QStringLiteral("CosineSimShadowsToAnchor8");
    // const QString COS_SIM_ANCH_SHADOW_9 = QStringLiteral("CosineSimShadowsToAnchor9");
    // const QString COS_SIM_ANCH_SHADOW_10 = QStringLiteral("CosineSimShadowsToAnchor10");
    // const QString COS_SIM_ANCH_SHADOW_11 = QStringLiteral("CosineSimShadowsToAnchor11");
    // const QString COS_SIM_ANCH_SHADOW_12 = QStringLiteral("CosineSimShadowsToAnchor12");
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
    // const QString MZ_PK_LEN_NORM_1 = QStringLiteral("MzPeakLengthsNorm1");
    // const QString MZ_PK_LEN_NORM_2 = QStringLiteral("MzPeakLengthsNorm2");
    // const QString MZ_PK_LEN_NORM_3 = QStringLiteral("MzPeakLengthsNorm3");
    // const QString MZ_PK_LEN_NORM_4 = QStringLiteral("MzPeakLengthsNorm4");
    // const QString MZ_PK_LEN_NORM_5 = QStringLiteral("MzPeakLengthsNorm5");
    // const QString MZ_PK_LEN_NORM_6 = QStringLiteral("MzPeakLengthsNorm6");
    // const QString MZ_PK_LEN_NORM_7 = QStringLiteral("MzPeakLengthsNorm7");
    // const QString MZ_PK_LEN_NORM_8 = QStringLiteral("MzPeakLengthsNorm8");
    // const QString MZ_PK_LEN_NORM_9 = QStringLiteral("MzPeakLengthsNorm9");
    // const QString MZ_PK_LEN_NORM_10 = QStringLiteral("MzPeakLengthsNorm10");
    // const QString MZ_PK_LEN_NORM_11 = QStringLiteral("MzPeakLengthsNorm11");
    // const QString MZ_PK_LEN_NORM_12 = QStringLiteral("MzPeakLengthsNorm12");
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
    // const QString MZ_FND_STDEV_7 = QStringLiteral("MzFoundStDev7");
    // const QString MZ_FND_STDEV_8 = QStringLiteral("MzFoundStDev8");
    // const QString MZ_FND_STDEV_9 = QStringLiteral("MzFoundStDev9");
    // const QString MZ_FND_STDEV_10 = QStringLiteral("MzFoundStDev10");
    // const QString MZ_FND_STDEV_11 = QStringLiteral("MzFoundStDev11");
    // const QString MZ_FND_STDEV_12 = QStringLiteral("MzFoundStDev12");
    const QString TARG_KEY = QStringLiteral("TargetKey");
    const QString PEP_STR_W_MODS = QStringLiteral("PeptideStringWithMods");
    const QString PROT_GRP = QStringLiteral("ProteinGroup");
    const QString IS_DECOY = QStringLiteral("IsDecoy");
    const QString SCAN_NUM = QStringLiteral("ScanNumber");
    const QString SCAN_TIME = QStringLiteral("ScanTime");
    const QString SCAN_TIME_START = QStringLiteral("ScanTimeStart");
    const QString SCAN_TIME_END = QStringLiteral("ScanTimeEnd");
    const QString CLASS_SCR = QStringLiteral("ClassifierScore");
    const QString DISC_SCR = QStringLiteral("DiscriminantScore");
    const QString Q_VAL = QStringLiteral("QValue");
    const QString DECOY_RATIO = QStringLiteral("DecoyRatio");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG_OG_ALT = QStringLiteral("AltTargetKeyIdDiscScoreChargeOG_alt");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG1_1 = QStringLiteral("AltTargetKeyIdCosineDiscScoreCharge1_1");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG1_2 = QStringLiteral("AltTargetKeyIdCosineDiscScoreCharge1_2");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG2_1 = QStringLiteral("AltTargetKeyIdCosineDiscScoreCharge2_1");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG2_2 = QStringLiteral("AltTargetKeyIdCosineDiscScoreCharge2_2");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG3_1 = QStringLiteral("AltTargetKeyIdCosineDiscScoreCharge3_1");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG3_2 = QStringLiteral("AltTargetKeyIdCosineDiscScoreCharge3_2");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG4_1 = QStringLiteral("AltTargetKeyIdCosineDiscScoreCharge4_1");
    const QString ALT_TARG_ID_DISC_SCORE_CHRG4_2 = QStringLiteral("AltTargetKeyIdCosineDiscScoreCharge4_2");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG1_1 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge1_1");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG2_1 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge2_1");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG3_1 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge3_1");
    const QString ALT_TARG_ID_TIME_DELTA_CHRG4_1 = QStringLiteral("AltTargetKeyIdTimeDeltaCharge4_1");
    const QString MS1_MZ_MEAN_FND_100 = QStringLiteral("Ms1MzMeanFound100");
    const QString MS1_MZ_MEAN_FND_45 = QStringLiteral("Ms1MzMeanFound45");
    const QString MS1_MZ_MEAN_FND_PRE_MONO = QStringLiteral("Ms1MzMeanFoundPreMono");
    const QString MS1_MZ_MEAN_FND_ISO1 = QStringLiteral("Ms1MzMeanFoundIso1");
    const QString MS1_MZ_MEAN_FND_ISO2 = QStringLiteral("Ms1MzMeanFoundIso2");
    const QString MS1_MZ_MEAN_FND_100_PPM = QStringLiteral("Ms1MzMeanFound100PPM");
    const QString MS1_MZ_MEAN_FND_45_PPM = QStringLiteral("Ms1MzMeanFound45PPM");
    const QString MS1_MZ_MEAN_FND_PRE_MONO_PPM = QStringLiteral("Ms1MzMeanFoundPreMonoPPM");
    const QString MS1_MZ_MEAN_FND_ISO_1_PPM = QStringLiteral("Ms1MzMeanFoundIso1PPM");
    const QString MS1_MZ_MEAN_FND_ISO_2_PPM = QStringLiteral("Ms1MzMeanFoundIso2PPM");
    const QString MS1_MZ_MEAN_FND_100_STD = QStringLiteral("Ms1MzStDevFound100");
    const QString MS1_MZ_MEAN_FND_45_STD = QStringLiteral("Ms1MzStDevFound45");
    const QString MS1_MZ_MEAN_FND_PRE_MONO_STD = QStringLiteral("Ms1MzStDevFoundPreMono");
    const QString MS1_MZ_MEAN_FND_ISO_1_STD = QStringLiteral("Ms1MzStDevFoundIso1");
    const QString MS1_MZ_MEAN_FND_ISO_2_STD = QStringLiteral("Ms1MzStDevFoundIso2");
    const QString MS1_INTZ_FND_100 = QStringLiteral("Ms1IntensityFound100");
    const QString MS1_INTZ_FND_45 = QStringLiteral("Ms1IntensityFound45");
    const QString MS1_INTZ_FND_PRE_MONO = QStringLiteral("Ms1IntensityFoundPreMono");
    const QString MS1_INTZ_FND_ISO_1 = QStringLiteral("Ms1IntensityFoundIso1");
    const QString MS1_INTZ_FND_ISO_2 = QStringLiteral("Ms1IntensityFoundIso2");
    const QString COS_SIM_SPEC_TIME = QStringLiteral("CosineSimSpectrumOverTime");
    const QString COS_SIM_SPEC_TIME_CUBED = QStringLiteral("CosineSimSpectrumOverTimeCubed");
    const QString COS_SIM_SPEC_STDEV = QStringLiteral("CosineSimSpectrumStDev");
    const QString COS_SIM_SUM100_MS1 = QStringLiteral("CosineSimSum100MS1");
    const QString MS1_AVERAGINE = QStringLiteral("MS1Averagine");
    const QString COS_SIM_SUM100_WIN_1p5X = QStringLiteral("CosineSimSum100Window1p5X");
    const QString COS_SIM_SUM100_WIN_2X = QStringLiteral("CosineSimSum100Window2X");
    const QString TOT_INTENSITY_PEAK_HEIGHTS = QStringLiteral("TotalIntensityPeakHeights");
    const QString TOT_INTENSITY_RAW = QStringLiteral("TotalIntensityRaw");
    const QString TARGET_WINDOW_LOCATION = QStringLiteral("TargetWindowLocation");
    // const QString TRAP_AREA_1 = QStringLiteral("TrapArea1");
    // const QString TRAP_AREA_2 = QStringLiteral("TrapArea2");
    // const QString TRAP_AREA_3 = QStringLiteral("TrapArea3");
    // const QString TRAP_AREA_4 = QStringLiteral("TrapArea4");
    // const QString TRAP_AREA_5 = QStringLiteral("TrapArea5");
    // const QString TRAP_AREA_6 = QStringLiteral("TrapArea6");
    // const QString TRAP_AREA_7 = QStringLiteral("TrapArea7");
    // const QString TRAP_AREA_8 = QStringLiteral("TrapArea8");
    // const QString TRAP_AREA_9 = QStringLiteral("TrapArea9");
    // const QString TRAP_AREA_10 = QStringLiteral("TrapArea10");
    // const QString TRAP_AREA_11 = QStringLiteral("TrapArea11");
    // const QString TRAP_AREA_12 = QStringLiteral("TrapArea12");

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


    const QStringList keysToCheck = {
            COS_SIM_SUM_100,
            COS_SIM_SUM_100_GREATER_80,
            ALL_MAX_IND_CNT,
            COS_SIM_SUM_MS1_100,
            COS_SIM_SPEC_CUBED,
            KL_DIV_SPEC_CUBE_RT,
            COS_SIM_SUM_45,
            COS_SIM_SUM_TOP_6,
            COS_SIM_SUM_BOTTOM_6,
            TOP_BOTTOM_RATIO,
            TOP_BOTTOME_RATIO_NORM,
            CHARGE,
            MASS,
            SCAN_TIME_DELTA,
            SCAN_TIME_PD,
            SCAN_ION_CNT,
            MZ_NORM,
            KL_DIV_SPEC,
            COSINE_SIM_SPEC,
            COSINE_SIM_SUM_MS1_45,
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
            // COS_SIM_ANCH_SHADOW_1,
            // COS_SIM_ANCH_SHADOW_2,
            // COS_SIM_ANCH_SHADOW_3,
            // COS_SIM_ANCH_SHADOW_4,
            // COS_SIM_ANCH_SHADOW_5,
            // COS_SIM_ANCH_SHADOW_6,
            // COS_SIM_ANCH_SHADOW_7,
            // COS_SIM_ANCH_SHADOW_8,
            // COS_SIM_ANCH_SHADOW_9,
            // COS_SIM_ANCH_SHADOW_10,
            // COS_SIM_ANCH_SHADOW_11,
            // COS_SIM_ANCH_SHADOW_12,
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
            // MZ_PK_LEN_NORM_1,
            // MZ_PK_LEN_NORM_2,
            // MZ_PK_LEN_NORM_3,
            // MZ_PK_LEN_NORM_4,
            // MZ_PK_LEN_NORM_5,
            // MZ_PK_LEN_NORM_6,
            // MZ_PK_LEN_NORM_7,
            // MZ_PK_LEN_NORM_8,
            // MZ_PK_LEN_NORM_9,
            // MZ_PK_LEN_NORM_10,
            // MZ_PK_LEN_NORM_11,
            // MZ_PK_LEN_NORM_12,
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
            // MZ_FND_STDEV_7,
            // MZ_FND_STDEV_8,
            // MZ_FND_STDEV_9,
            // MZ_FND_STDEV_10,
            // MZ_FND_STDEV_11,
            // MZ_FND_STDEV_12,
            TARG_KEY,
            PEP_STR_W_MODS,
            PROT_GRP,
            IS_DECOY,
            SCAN_NUM,
            SCAN_TIME,
            SCAN_TIME_START,
            SCAN_TIME_END,
            CLASS_SCR,
            DISC_SCR,
            Q_VAL,
            DECOY_RATIO,
            ALT_TARG_ID_DISC_SCORE_CHRG_OG_ALT,
            ALT_TARG_ID_DISC_SCORE_CHRG1_1,
            ALT_TARG_ID_DISC_SCORE_CHRG1_2,
            ALT_TARG_ID_DISC_SCORE_CHRG2_1,
            ALT_TARG_ID_DISC_SCORE_CHRG2_2,
            ALT_TARG_ID_DISC_SCORE_CHRG3_1,
            ALT_TARG_ID_DISC_SCORE_CHRG3_2,
            ALT_TARG_ID_DISC_SCORE_CHRG4_1,
            ALT_TARG_ID_DISC_SCORE_CHRG4_2,
            ALT_TARG_ID_TIME_DELTA_CHRG1_1,
            ALT_TARG_ID_TIME_DELTA_CHRG2_1,
            ALT_TARG_ID_TIME_DELTA_CHRG3_1,
            ALT_TARG_ID_TIME_DELTA_CHRG4_1,
            MS1_MZ_MEAN_FND_100,
            MS1_MZ_MEAN_FND_45,
            MS1_MZ_MEAN_FND_PRE_MONO,
            MS1_MZ_MEAN_FND_ISO1,
            MS1_MZ_MEAN_FND_ISO2,
            MS1_MZ_MEAN_FND_100_PPM,
            MS1_MZ_MEAN_FND_45_PPM,
            MS1_MZ_MEAN_FND_PRE_MONO_PPM,
            MS1_MZ_MEAN_FND_ISO_1_PPM,
            MS1_MZ_MEAN_FND_ISO_2_PPM,
            MS1_MZ_MEAN_FND_100_STD,
            MS1_MZ_MEAN_FND_45_STD,
            MS1_MZ_MEAN_FND_PRE_MONO_STD,
            MS1_MZ_MEAN_FND_ISO_1_STD,
            MS1_MZ_MEAN_FND_ISO_2_STD,
            MS1_INTZ_FND_100,
            MS1_INTZ_FND_45,
            MS1_INTZ_FND_PRE_MONO,
            MS1_INTZ_FND_ISO_1,
            MS1_INTZ_FND_ISO_2,
            COS_SIM_SPEC_TIME,
            COS_SIM_SPEC_TIME_CUBED,
            COS_SIM_SPEC_STDEV,
            COS_SIM_SUM100_MS1,
            MS1_AVERAGINE,
            COS_SIM_SUM100_WIN_1p5X,
            COS_SIM_SUM100_WIN_2X,
            TOT_INTENSITY_PEAK_HEIGHTS,
            TOT_INTENSITY_RAW,
            TARGET_WINDOW_LOCATION,
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
            MZ_SEARCHED_12
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
    float cosineSimSumTop = -1.0;
    float cosineSimSumBottom = -1.0;
    float topBottomRatio = -1.0;
    float topBottomRatioNorm = -1.0;
    float charge = -1.0;
    float mass = -1.0;
    float scanTimeDelta = -1.0;
    float scanTimePd = -1.0;
    float scanIonCount = -1.0;
    float mzNorm = -1.0;
    float klDivSpectrum = -1.0;
    float cosineSimSpectrum = -1.0;
    float cosineSim45MS1 = -1.0;
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
    // float cosineSimShadowsToAnchor1 = -1.0;
    // float cosineSimShadowsToAnchor2 = -1.0;
    // float cosineSimShadowsToAnchor3 = -1.0;
    // float cosineSimShadowsToAnchor4 = -1.0;
    // float cosineSimShadowsToAnchor5 = -1.0;
    // float cosineSimShadowsToAnchor6 = -1.0;
    // float cosineSimShadowsToAnchor7 = -1.0;
    // float cosineSimShadowsToAnchor8 = -1.0;
    // float cosineSimShadowsToAnchor9 = -1.0;
    // float cosineSimShadowsToAnchor10 = -1.0;
    // float cosineSimShadowsToAnchor11 = -1.0;
    // float cosineSimShadowsToAnchor12 = -1.0;
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
    // float mzPeakLengthsNorm1 = -1.0;
    // float mzPeakLengthsNorm2 = -1.0;
    // float mzPeakLengthsNorm3 = -1.0;
    // float mzPeakLengthsNorm4 = -1.0;
    // float mzPeakLengthsNorm5 = -1.0;
    // float mzPeakLengthsNorm6 = -1.0;
    // float mzPeakLengthsNorm7 = -1.0;
    // float mzPeakLengthsNorm8 = -1.0;
    // float mzPeakLengthsNorm9 = -1.0;
    // float mzPeakLengthsNorm10 = -1.0;
    // float mzPeakLengthsNorm11 = -1.0;
    // float mzPeakLengthsNorm12 = -1.0;
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
    // float mzFoundStDev7 = -1.0;
    // float mzFoundStDev8 = -1.0;
    // float mzFoundStDev9 = -1.0;
    // float mzFoundStDev10 = -1.0;
    // float mzFoundStDev11 = -1.0;
    // float mzFoundStDev12 = -1.0;

    QString targetKey;
    PeptideStringWithMods peptideStringWithMods;
    QString proteinGroup;
    bool isDecoy = false;
    ScanNumber scanNumber = -1;
    ScanTime scanTime = -1.0;
    ScanTime scanTimeStart = -1.0;
    ScanTime scanTimeEnd = -1.0;
    double classifierScore = -1.0;
    double discriminantScore = -1.0;
    double qValue = 1.0;
    double decoyRatio = -1.0;

    float altTargetKeyIdDiscScoreChargeOG_alt = -1.0;
    float altTargetKeyIdDiscScoreCharge1_1 = -1.0;
    float altTargetKeyIdDiscScoreCharge1_2 = -1.0;
    float altTargetKeyIdDiscScoreCharge2_1 = -1.0;
    float altTargetKeyIdDiscScoreCharge2_2 = -1.0;
    float altTargetKeyIdDiscScoreCharge3_1 = -1.0;
    float altTargetKeyIdDiscScoreCharge3_2 = -1.0;
    float altTargetKeyIdDiscScoreCharge4_1 = -1.0;
    float altTargetKeyIdDiscScoreCharge4_2 = -1.0;

    float altTargetKeyIdTimeDeltaCharge1_1 = -1.0;
    float altTargetKeyIdTimeDeltaCharge2_1 = -1.0;
    float altTargetKeyIdTimeDeltaCharge3_1 = -1.0;
    float altTargetKeyIdTimeDeltaCharge4_1 = -1.0;

    float ms1MzMeanFound100 = -1.0;
    float ms1MzMeanFound45 = -1.0;
    float ms1MzMeanFoundPreMono = -1.0;
    float ms1MzMeanFoundIso1 = -1.0;
    float ms1MzMeanFoundIso2 = -1.0;
    float ms1MzMeanFound100PPM = -1.0;
    float ms1MzMeanFound45PPM = -1.0;
    float ms1MzMeanFoundPreMonoPPM = -1.0;
    float ms1MzMeanFoundIso1PPM = -1.0;
    float ms1MzMeanFoundIso2PPM = -1.0;
    float ms1MzStDevFound100 = -1.0;
    float ms1MzStDevFound45 = -1.0;
    float ms1MzStDevFoundPreMono = -1.0;
    float ms1MzStDevFoundIso1 = -1.0;
    float ms1MzStDevFoundIso2 = -1.0;
    float ms1IntensityFound100 = -1.0;
    float ms1IntensityFound45 = -1.0;
    float ms1IntensityFoundPreMono = -1.0;
    float ms1IntensityFoundIso1 = -1.0;
    float ms1IntensityFoundIso2 = -1.0;

    float cosineSimSpectrumOverTime = -1.0;
    float cosineSimSpectrumOverTimeCubed = -1.0;
    float cosineSimSpectrumStDev = -1.0;
    float cosineSimSum100MS1 = -1.0;
    float ms1Averagine = -1.0;
    float cosineSimSum100Window1p5X = -1.0;
    float cosineSimSum100Window2X = -1.0;
    float totalIntensityPeakHeights = -1.0;
    float totalIntensityRaw = -1.0;
    float targetWindowLocation = -1.0;

    // float trapArea1 = -1.0;
    // float trapArea2 = -1.0;
    // float trapArea3 = -1.0;
    // float trapArea4 = -1.0;
    // float trapArea5 = -1.0;
    // float trapArea6 = -1.0;
    // float trapArea7 = -1.0;
    // float trapArea8 = -1.0;
    // float trapArea9 = -1.0;
    // float trapArea10 = -1.0;
    // float trapArea11 = -1.0;
    // float trapArea12 = -1.0;

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
        cosineSimSumTop = dataMap.value(COS_SIM_SUM_TOP_6).toFloat();
        cosineSimSumBottom = dataMap.value(COS_SIM_SUM_BOTTOM_6).toFloat();
        topBottomRatio = dataMap.value(TOP_BOTTOM_RATIO).toFloat();
        topBottomRatioNorm = dataMap.value(TOP_BOTTOME_RATIO_NORM).toFloat();
        charge = dataMap.value(CHARGE).toFloat();
        mass = dataMap.value(MASS).toFloat();
        scanTimeDelta = dataMap.value(SCAN_TIME_DELTA).toFloat();
        scanTimePd = dataMap.value(SCAN_TIME_PD).toFloat();
        scanIonCount = dataMap.value(SCAN_ION_CNT).toFloat();
        mzNorm = dataMap.value(MZ_NORM).toFloat();
        klDivSpectrum = dataMap.value(KL_DIV_SPEC).toFloat();
        cosineSimSpectrum = dataMap.value(COSINE_SIM_SPEC).toFloat();
        cosineSim45MS1 = dataMap.value(COSINE_SIM_SUM_MS1_45).toFloat();
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
        // cosineSimShadowsToAnchor1 = dataMap.value(COS_SIM_ANCH_SHADOW_1).toFloat();
        // cosineSimShadowsToAnchor2 = dataMap.value(COS_SIM_ANCH_SHADOW_2).toFloat();
        // cosineSimShadowsToAnchor3 = dataMap.value(COS_SIM_ANCH_SHADOW_3).toFloat();
        // cosineSimShadowsToAnchor4 = dataMap.value(COS_SIM_ANCH_SHADOW_4).toFloat();
        // cosineSimShadowsToAnchor5 = dataMap.value(COS_SIM_ANCH_SHADOW_5).toFloat();
        // cosineSimShadowsToAnchor6 = dataMap.value(COS_SIM_ANCH_SHADOW_6).toFloat();
        // cosineSimShadowsToAnchor7 = dataMap.value(COS_SIM_ANCH_SHADOW_7).toFloat();
        // cosineSimShadowsToAnchor8 = dataMap.value(COS_SIM_ANCH_SHADOW_8).toFloat();
        // cosineSimShadowsToAnchor9 = dataMap.value(COS_SIM_ANCH_SHADOW_9).toFloat();
        // cosineSimShadowsToAnchor10 = dataMap.value(COS_SIM_ANCH_SHADOW_10).toFloat();
        // cosineSimShadowsToAnchor11 = dataMap.value(COS_SIM_ANCH_SHADOW_11).toFloat();
        // cosineSimShadowsToAnchor12 = dataMap.value(COS_SIM_ANCH_SHADOW_12).toFloat();
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
        // mzPeakLengthsNorm1 = dataMap.value(MZ_PK_LEN_NORM_1).toFloat();
        // mzPeakLengthsNorm2 = dataMap.value(MZ_PK_LEN_NORM_2).toFloat();
        // mzPeakLengthsNorm3 = dataMap.value(MZ_PK_LEN_NORM_3).toFloat();
        // mzPeakLengthsNorm4 = dataMap.value(MZ_PK_LEN_NORM_4).toFloat();
        // mzPeakLengthsNorm5 = dataMap.value(MZ_PK_LEN_NORM_5).toFloat();
        // mzPeakLengthsNorm6 = dataMap.value(MZ_PK_LEN_NORM_6).toFloat();
        // mzPeakLengthsNorm7 = dataMap.value(MZ_PK_LEN_NORM_7).toFloat();
        // mzPeakLengthsNorm8 = dataMap.value(MZ_PK_LEN_NORM_8).toFloat();
        // mzPeakLengthsNorm9 = dataMap.value(MZ_PK_LEN_NORM_9).toFloat();
        // mzPeakLengthsNorm10 = dataMap.value(MZ_PK_LEN_NORM_10).toFloat();
        // mzPeakLengthsNorm11 = dataMap.value(MZ_PK_LEN_NORM_11).toFloat();
        // mzPeakLengthsNorm12 = dataMap.value(MZ_PK_LEN_NORM_12).toFloat();
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
        // mzFoundStDev7 = dataMap.value(MZ_FND_MEAN_8).toFloat();
        // mzFoundStDev8 = dataMap.value(MZ_FND_MEAN_8).toFloat();
        // mzFoundStDev9 = dataMap.value(MZ_FND_MEAN_9).toFloat();
        // mzFoundStDev10 = dataMap.value(MZ_FND_MEAN_10).toFloat();
        // mzFoundStDev11 = dataMap.value(MZ_FND_MEAN_11).toFloat();
        // mzFoundStDev12 = dataMap.value(MZ_FND_MEAN_12).toFloat();

        targetKey = dataMap.value(TARG_KEY).toString();
        peptideStringWithMods = PeptideStringWithMods(dataMap.value(PEP_STR_W_MODS).toString());
        proteinGroup = dataMap.value(PROT_GRP).toString();
        isDecoy = dataMap.value(IS_DECOY).toBool();
        scanNumber = dataMap.value(SCAN_NUM).toInt();
        scanTime = dataMap.value(SCAN_TIME).toFloat();
        scanTimeStart = dataMap.value(SCAN_TIME_START).toFloat();
        scanTimeEnd = dataMap.value(SCAN_TIME_END).toFloat();
        classifierScore = dataMap.value(CLASS_SCR).toDouble();
        discriminantScore = dataMap.value(DISC_SCR).toDouble();
        qValue = dataMap.value(Q_VAL).toDouble();
        decoyRatio = dataMap.value(DECOY_RATIO).toDouble();

        altTargetKeyIdDiscScoreChargeOG_alt = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG_OG_ALT).toFloat();
        altTargetKeyIdDiscScoreCharge1_1 = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG1_1).toFloat();
        altTargetKeyIdDiscScoreCharge1_2 = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG1_2).toFloat();
        altTargetKeyIdDiscScoreCharge2_1 = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG2_1).toFloat();
        altTargetKeyIdDiscScoreCharge2_2 = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG2_2).toFloat();
        altTargetKeyIdDiscScoreCharge3_1 = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG3_1).toFloat();
        altTargetKeyIdDiscScoreCharge3_2 = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG3_2).toFloat();
        altTargetKeyIdDiscScoreCharge4_1 = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG4_1).toFloat();
        altTargetKeyIdDiscScoreCharge4_2 = dataMap.value(ALT_TARG_ID_DISC_SCORE_CHRG4_2).toFloat();
        altTargetKeyIdTimeDeltaCharge1_1 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG1_1).toFloat();
        altTargetKeyIdTimeDeltaCharge2_1 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG2_1).toFloat();
        altTargetKeyIdTimeDeltaCharge3_1 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG3_1).toFloat();
        altTargetKeyIdTimeDeltaCharge4_1 = dataMap.value(ALT_TARG_ID_TIME_DELTA_CHRG4_1).toFloat();
        ms1MzMeanFound100 = dataMap.value(MS1_MZ_MEAN_FND_100).toFloat();
        ms1MzMeanFound45 = dataMap.value(MS1_MZ_MEAN_FND_45).toFloat();
        ms1MzMeanFoundPreMono = dataMap.value(MS1_MZ_MEAN_FND_PRE_MONO).toFloat();
        ms1MzMeanFoundIso1 = dataMap.value(MS1_MZ_MEAN_FND_ISO1).toFloat();
        ms1MzMeanFoundIso2 = dataMap.value(MS1_MZ_MEAN_FND_ISO2).toFloat();
        ms1MzMeanFound100PPM = dataMap.value(MS1_MZ_MEAN_FND_100_PPM).toFloat();
        ms1MzMeanFound45PPM = dataMap.value(MS1_MZ_MEAN_FND_45_PPM).toFloat();
        ms1MzMeanFoundPreMonoPPM = dataMap.value(MS1_MZ_MEAN_FND_PRE_MONO_PPM).toFloat();
        ms1MzMeanFoundIso1PPM = dataMap.value(MS1_MZ_MEAN_FND_ISO1).toFloat();
        ms1MzMeanFoundIso2PPM = dataMap.value(MS1_MZ_MEAN_FND_ISO2).toFloat();
        ms1MzStDevFound100 = dataMap.value(MS1_MZ_MEAN_FND_100_STD).toFloat();
        ms1MzStDevFound45 = dataMap.value(MS1_MZ_MEAN_FND_45_STD).toFloat();
        ms1MzStDevFoundPreMono = dataMap.value(MS1_MZ_MEAN_FND_PRE_MONO_STD).toFloat();
        ms1MzStDevFoundIso1 = dataMap.value(MS1_MZ_MEAN_FND_ISO_1_STD).toFloat();
        ms1MzStDevFoundIso2 = dataMap.value(MS1_MZ_MEAN_FND_ISO_2_STD).toFloat();
        ms1IntensityFound100 = dataMap.value(MS1_INTZ_FND_100).toFloat();
        ms1IntensityFound45 = dataMap.value(MS1_INTZ_FND_45).toFloat();
        ms1IntensityFoundPreMono = dataMap.value(MS1_INTZ_FND_PRE_MONO).toFloat();
        ms1IntensityFoundIso1 = dataMap.value(MS1_INTZ_FND_ISO_1).toFloat();
        ms1IntensityFoundIso2 = dataMap.value(MS1_INTZ_FND_ISO_2).toFloat();
        cosineSimSpectrumOverTime = dataMap.value(COS_SIM_SPEC_TIME).toFloat();
        cosineSimSpectrumOverTimeCubed = dataMap.value(COS_SIM_SPEC_TIME_CUBED).toFloat();
        cosineSimSpectrumStDev = dataMap.value(COS_SIM_SPEC_STDEV).toFloat();
        cosineSimSum100MS1 = dataMap.value(COS_SIM_SUM100_MS1).toFloat();
        ms1Averagine = dataMap.value(MS1_AVERAGINE).toFloat();
        cosineSimSum100Window1p5X = dataMap.value(COS_SIM_SUM100_WIN_1p5X).toFloat();
        cosineSimSum100Window2X = dataMap.value(COS_SIM_SUM100_WIN_2X).toFloat();
        totalIntensityPeakHeights = dataMap.value(TOT_INTENSITY_PEAK_HEIGHTS).toFloat();
        totalIntensityRaw = dataMap.value(TOT_INTENSITY_RAW).toFloat();
        targetWindowLocation = dataMap.value(TARGET_WINDOW_LOCATION).toFloat();

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
                {COS_SIM_SUM_TOP_6, QVariant(cosineSimSumTop)},
                {COS_SIM_SUM_BOTTOM_6, QVariant(cosineSimSumBottom)},
                {TOP_BOTTOM_RATIO, QVariant(topBottomRatio)},
                {TOP_BOTTOME_RATIO_NORM, QVariant(topBottomRatioNorm)},
                {CHARGE, QVariant(charge)},
                {MASS, QVariant(mass)},
                {SCAN_TIME_DELTA, QVariant(scanTimeDelta)},
                {SCAN_TIME_PD, QVariant(scanTimePd)},
                {SCAN_ION_CNT, QVariant(scanIonCount)},
                {MZ_NORM, QVariant(mzNorm)},
                {KL_DIV_SPEC, QVariant(klDivSpectrum)},
                {COSINE_SIM_SPEC, QVariant(cosineSimSpectrum)},
                {COSINE_SIM_SUM_MS1_45, QVariant(cosineSim45MS1)},
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
                // {COS_SIM_ANCH_SHADOW_1, QVariant(cosineSimShadowsToAnchor1)},
                // {COS_SIM_ANCH_SHADOW_2, QVariant(cosineSimShadowsToAnchor2)},
                // {COS_SIM_ANCH_SHADOW_3, QVariant(cosineSimShadowsToAnchor3)},
                // {COS_SIM_ANCH_SHADOW_4, QVariant(cosineSimShadowsToAnchor4)},
                // {COS_SIM_ANCH_SHADOW_5, QVariant(cosineSimShadowsToAnchor5)},
                // {COS_SIM_ANCH_SHADOW_6, QVariant(cosineSimShadowsToAnchor6)},
                // {COS_SIM_ANCH_SHADOW_7, QVariant(cosineSimShadowsToAnchor7)},
                // {COS_SIM_ANCH_SHADOW_8, QVariant(cosineSimShadowsToAnchor8)},
                // {COS_SIM_ANCH_SHADOW_9, QVariant(cosineSimShadowsToAnchor9)},
                // {COS_SIM_ANCH_SHADOW_10, QVariant(cosineSimShadowsToAnchor10)},
                // {COS_SIM_ANCH_SHADOW_11, QVariant(cosineSimShadowsToAnchor11)},
                // {COS_SIM_ANCH_SHADOW_12, QVariant(cosineSimShadowsToAnchor12)},
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
                // {MZ_PK_LEN_NORM_1, QVariant(mzPeakLengthsNorm1)},
                // {MZ_PK_LEN_NORM_2, QVariant(mzPeakLengthsNorm2)},
                // {MZ_PK_LEN_NORM_3, QVariant(mzPeakLengthsNorm3)},
                // {MZ_PK_LEN_NORM_4, QVariant(mzPeakLengthsNorm4)},
                // {MZ_PK_LEN_NORM_5, QVariant(mzPeakLengthsNorm5)},
                // {MZ_PK_LEN_NORM_6, QVariant(mzPeakLengthsNorm6)},
                // {MZ_PK_LEN_NORM_7, QVariant(mzPeakLengthsNorm7)},
                // {MZ_PK_LEN_NORM_8, QVariant(mzPeakLengthsNorm8)},
                // {MZ_PK_LEN_NORM_9, QVariant(mzPeakLengthsNorm9)},
                // {MZ_PK_LEN_NORM_10, QVariant(mzPeakLengthsNorm10)},
                // {MZ_PK_LEN_NORM_11, QVariant(mzPeakLengthsNorm11)},
                // {MZ_PK_LEN_NORM_12, QVariant(mzPeakLengthsNorm12)},
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
                // {MZ_FND_STDEV_7, QVariant(mzFoundStDev7)},
                // {MZ_FND_STDEV_8, QVariant(mzFoundStDev8)},
                // {MZ_FND_STDEV_9, QVariant(mzFoundStDev9)},
                // {MZ_FND_STDEV_10, QVariant(mzFoundStDev10)},
                // {MZ_FND_STDEV_11, QVariant(mzFoundStDev11)},
                // {MZ_FND_STDEV_12, QVariant(mzFoundStDev12)},
                {TARG_KEY, QVariant(targetKey)},
                {PEP_STR_W_MODS, QVariant(peptideStringWithMods)},
                {PROT_GRP, QVariant(proteinGroup)},
                {IS_DECOY, QVariant(isDecoy)},
                {SCAN_NUM, QVariant(scanNumber)},
                {SCAN_TIME, QVariant(scanTime)},
                {SCAN_TIME_START, QVariant(scanTimeStart)},
                {SCAN_TIME_END, QVariant(scanTimeEnd)},
                {CLASS_SCR, QVariant(classifierScore)},
                {DISC_SCR, QVariant(discriminantScore)},
                {Q_VAL, QVariant(qValue)},
                {DECOY_RATIO, QVariant(decoyRatio)},

                {ALT_TARG_ID_DISC_SCORE_CHRG_OG_ALT, QVariant(altTargetKeyIdDiscScoreChargeOG_alt)},
                {ALT_TARG_ID_DISC_SCORE_CHRG1_1, QVariant(altTargetKeyIdDiscScoreCharge1_1)},
                {ALT_TARG_ID_DISC_SCORE_CHRG1_2, QVariant(altTargetKeyIdDiscScoreCharge1_2)},
                {ALT_TARG_ID_DISC_SCORE_CHRG2_1, QVariant(altTargetKeyIdDiscScoreCharge2_1)},
                {ALT_TARG_ID_DISC_SCORE_CHRG2_2, QVariant(altTargetKeyIdDiscScoreCharge2_2)},
                {ALT_TARG_ID_DISC_SCORE_CHRG3_1, QVariant(altTargetKeyIdDiscScoreCharge3_1)},
                {ALT_TARG_ID_DISC_SCORE_CHRG3_2, QVariant(altTargetKeyIdDiscScoreCharge3_2)},
                {ALT_TARG_ID_DISC_SCORE_CHRG4_1, QVariant(altTargetKeyIdDiscScoreCharge4_1)},
                {ALT_TARG_ID_DISC_SCORE_CHRG4_2, QVariant(altTargetKeyIdDiscScoreCharge4_2)},
                {ALT_TARG_ID_TIME_DELTA_CHRG1_1 , QVariant(altTargetKeyIdTimeDeltaCharge1_1)},
                {ALT_TARG_ID_TIME_DELTA_CHRG2_1 , QVariant(altTargetKeyIdTimeDeltaCharge2_1)},
                {ALT_TARG_ID_TIME_DELTA_CHRG3_1 , QVariant(altTargetKeyIdTimeDeltaCharge3_1)},
                {ALT_TARG_ID_TIME_DELTA_CHRG4_1 , QVariant(altTargetKeyIdTimeDeltaCharge4_1)},
                {MS1_MZ_MEAN_FND_100, QVariant(ms1MzMeanFound100)},
                {MS1_MZ_MEAN_FND_45, QVariant(ms1MzMeanFound45)},
                {MS1_MZ_MEAN_FND_PRE_MONO, QVariant(ms1MzMeanFoundPreMono)},
                {MS1_MZ_MEAN_FND_ISO1, QVariant(ms1MzMeanFoundIso1)},
                {MS1_MZ_MEAN_FND_ISO2, QVariant(ms1MzMeanFoundIso2)},
                {MS1_MZ_MEAN_FND_100_PPM, QVariant(ms1MzMeanFound100PPM)},
                {MS1_MZ_MEAN_FND_45_PPM, QVariant(ms1MzMeanFound45PPM)},
                {MS1_MZ_MEAN_FND_PRE_MONO_PPM, QVariant(ms1MzMeanFoundPreMonoPPM)},
                {MS1_MZ_MEAN_FND_ISO_1_PPM, QVariant(ms1MzMeanFoundIso1PPM)},
                {MS1_MZ_MEAN_FND_ISO_2_PPM, QVariant(ms1MzMeanFoundIso2PPM)},
                {MS1_MZ_MEAN_FND_100_STD, QVariant(ms1MzStDevFound100)},
                {MS1_MZ_MEAN_FND_45_STD, QVariant(ms1MzStDevFound45)},
                {MS1_MZ_MEAN_FND_PRE_MONO_STD, QVariant(ms1MzStDevFoundPreMono)},
                {MS1_MZ_MEAN_FND_ISO_1_STD, QVariant(ms1MzStDevFoundIso1)},
                {MS1_MZ_MEAN_FND_ISO_2_STD, QVariant(ms1MzStDevFoundIso2)},
                {MS1_INTZ_FND_100, QVariant(ms1IntensityFound100)},
                {MS1_INTZ_FND_45, QVariant(ms1IntensityFound45)},
                {MS1_INTZ_FND_PRE_MONO, QVariant(ms1IntensityFoundPreMono)},
                {MS1_INTZ_FND_ISO_1, QVariant(ms1IntensityFoundIso1)},
                {MS1_INTZ_FND_ISO_2, QVariant(ms1IntensityFoundIso2)},

                {COS_SIM_SPEC_TIME, QVariant(cosineSimSpectrumOverTime)},
                {COS_SIM_SPEC_TIME_CUBED, QVariant(cosineSimSpectrumOverTimeCubed)},
                {COS_SIM_SPEC_STDEV, QVariant(cosineSimSpectrumStDev)},
                {COS_SIM_SUM100_MS1, QVariant(cosineSimSum100MS1)},
                {MS1_AVERAGINE, QVariant(ms1Averagine)},
                {COS_SIM_SUM100_WIN_1p5X, QVariant(cosineSimSum100Window1p5X)},
                {COS_SIM_SUM100_WIN_2X, QVariant(cosineSimSum100Window2X)},
                {TOT_INTENSITY_PEAK_HEIGHTS, QVariant(totalIntensityPeakHeights)},
                {TOT_INTENSITY_RAW, QVariant(totalIntensityRaw)},
                {TARGET_WINDOW_LOCATION, QVariant(targetWindowLocation)},

                // {TRAP_AREA_1, QVariant(trapArea1)},
                // {TRAP_AREA_2, QVariant(trapArea2)},
                // {TRAP_AREA_3, QVariant(trapArea3)},
                // {TRAP_AREA_4, QVariant(trapArea4)},
                // {TRAP_AREA_5, QVariant(trapArea5)},
                // {TRAP_AREA_6, QVariant(trapArea6)},
                // {TRAP_AREA_7, QVariant(trapArea7)},
                // {TRAP_AREA_8, QVariant(trapArea8)},
                // {TRAP_AREA_9, QVariant(trapArea9)},
                // {TRAP_AREA_10, QVariant(trapArea10)},
                // {TRAP_AREA_11, QVariant(trapArea11)},
                // {TRAP_AREA_12, QVariant(trapArea12)},

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
        };
    }

    static CandidateScoresReaderRow buildCandidateScoresReaderRow(const CandidateScores* candidateScores) {

        CandidateScoresReaderRow row;

        row.cosineSimSum100 = candidateScores->featuresArray[Features::CosineSimSum100];
        row.cosineSimSum100Greater80 = candidateScores->featuresArray[Features::CosineSimSum100GreaterThan80];
        row.allignedMaxIndexesCount = candidateScores->featuresArray[Features::AllignedMaxIndexesCount];
        row.cosineSim100MS1 = candidateScores->featuresArray[Features::CosineSim100MS1];
        row.cosineSimSpectrumCubed = candidateScores->featuresArray[Features::CosineSimSpectrumCubed];
        row.klDivSpectrumCubeRoot = candidateScores->featuresArray[Features::KlDivSpectrumCubeRoot];
        row.cosineSimSum45 = candidateScores->featuresArray[Features::CosineSimSum45];
        row.cosineSimSumTop = candidateScores->featuresArray[Features::CosineSimSumTop];
        row.cosineSimSumBottom = candidateScores->featuresArray[Features::CosineSimSumBottom];
        row.topBottomRatio = candidateScores->featuresArray[Features::TopBottomRatio];
        row.topBottomRatioNorm = candidateScores->featuresArray[Features::TopBottomRatioNorm];
        row.charge = candidateScores->featuresArray[Features::Charge];
        row.mass = candidateScores->featuresArray[Features::Mass];
        row.scanTimeDelta = candidateScores->featuresArray[Features::ScanTimeDelta];
        row.scanTimePd = candidateScores->featuresArray[Features::ScanTimePd];
        row.scanIonCount = candidateScores->featuresArray[Features::ScanIonCount];
        row.mzNorm = candidateScores->featuresArray[Features::MzNorm];
        row.klDivSpectrum = candidateScores->featuresArray[Features::KlDivSpectrum];
        row.cosineSimSpectrum = candidateScores->featuresArray[Features::CosineSimSpectrum];
        row.cosineSim45MS1 = candidateScores->featuresArray[Features::CosineSim45MS1];
        row.cosineSim100MS1PreMono = candidateScores->featuresArray[Features::CosineSim100MS1PreMono];
        row.cosineSim100MS1Iso1 = candidateScores->featuresArray[Features::CosineSim100MS1Iso1];
        row.cosineSim100MS1Iso2 = candidateScores->featuresArray[Features::CosineSim100MS1Iso2];
        row.peptideLengthNorm = candidateScores->featuresArray[Features::PeptideLengthNorm];
        row.scanTimePredicted = candidateScores->featuresArray[Features::ScanTimePredicted];
        row.theoFragmentCount = candidateScores->featuresArray[Features::TheoFragmentCount];
        row.totalIntensityLog = candidateScores->featuresArray[Features::TotalIntensityLog];
        row.peakShapeRatio1 = candidateScores->featuresArray[Features::PeakShapeRatio1];
        row.peakShapeRatio2 = candidateScores->featuresArray[Features::PeakShapeRatio2];
        row.peakShapeRatio3 = candidateScores->featuresArray[Features::PeakShapeRatio3];
        row.shadowsCosineSimSum = candidateScores->featuresArray[Features::ShadowsCosineSimSum];
        row.iRtPredicted = candidateScores->featuresArray[Features::IRTPredicted];
        row.cosineSimToAnchor1 = candidateScores->featuresArray[Features::CosineSimToAnchor1];
        row.cosineSimToAnchor2 = candidateScores->featuresArray[Features::CosineSimToAnchor2];
        row.cosineSimToAnchor3 = candidateScores->featuresArray[Features::CosineSimToAnchor3];
        row.cosineSimToAnchor4 = candidateScores->featuresArray[Features::CosineSimToAnchor4];
        row.cosineSimToAnchor5 = candidateScores->featuresArray[Features::CosineSimToAnchor5];
        row.cosineSimToAnchor6 = candidateScores->featuresArray[Features::CosineSimToAnchor6];
        row.cosineSimToAnchor7 = candidateScores->featuresArray[Features::CosineSimToAnchor7];
        row.cosineSimToAnchor8 = candidateScores->featuresArray[Features::CosineSimToAnchor8];
        row.cosineSimToAnchor9 = candidateScores->featuresArray[Features::CosineSimToAnchor9];
        row.cosineSimToAnchor10 = candidateScores->featuresArray[Features::CosineSimToAnchor10];
        row.cosineSimToAnchor11 = candidateScores->featuresArray[Features::CosineSimToAnchor11];
        row.cosineSimToAnchor12 = candidateScores->featuresArray[Features::CosineSimToAnchor12];
        // row.cosineSimShadowsToAnchor1 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor1];
        // row.cosineSimShadowsToAnchor2 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor2];
        // row.cosineSimShadowsToAnchor3 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor3];
        // row.cosineSimShadowsToAnchor4 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor4];
        // row.cosineSimShadowsToAnchor5 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor5];
        // row.cosineSimShadowsToAnchor6 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor6];
        // row.cosineSimShadowsToAnchor7 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor7];
        // row.cosineSimShadowsToAnchor8 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor8];
        // row.cosineSimShadowsToAnchor9 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor9];
        // row.cosineSimShadowsToAnchor10 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor10];
        // row.cosineSimShadowsToAnchor11 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor11];
        // row.cosineSimShadowsToAnchor12 = candidateScores->featuresArray[Features::CosineSimShadowsToAnchor12];
        row.mzFoundMean1 = candidateScores->featuresArray[Features::MzFoundMean1];
        row.mzFoundMean2 = candidateScores->featuresArray[Features::MzFoundMean2];
        row.mzFoundMean3 = candidateScores->featuresArray[Features::MzFoundMean3];
        row.mzFoundMean4 = candidateScores->featuresArray[Features::MzFoundMean4];
        row.mzFoundMean5 = candidateScores->featuresArray[Features::MzFoundMean5];
        row.mzFoundMean6 = candidateScores->featuresArray[Features::MzFoundMean6];
        row.mzFoundMean7 = candidateScores->featuresArray[Features::MzFoundMean7];
        row.mzFoundMean8 = candidateScores->featuresArray[Features::MzFoundMean8];
        row.mzFoundMean9 = candidateScores->featuresArray[Features::MzFoundMean9];
        row.mzFoundMean10 = candidateScores->featuresArray[Features::MzFoundMean10];
        row.mzFoundMean11 = candidateScores->featuresArray[Features::MzFoundMean11];
        row.mzFoundMean12 = candidateScores->featuresArray[Features::MzFoundMean12];
        row.intensityFoundMax1 = candidateScores->integrations.at(0);
        row.intensityFoundMax2 = candidateScores->integrations.at(1);
        row.intensityFoundMax3 = candidateScores->integrations.at(2);
        row.intensityFoundMax4 = candidateScores->integrations.at(3);
        row.intensityFoundMax5 = candidateScores->integrations.at(4);
        row.intensityFoundMax6 = candidateScores->integrations.at(5);
        row.intensityFoundMax7 = candidateScores->integrations.at(6);
        row.intensityFoundMax8 = candidateScores->integrations.at(7);
        row.intensityFoundMax9 = candidateScores->integrations.at(8);
        row.intensityFoundMax10 = candidateScores->integrations.at(9);
        row.intensityFoundMax11 = candidateScores->integrations.at(10);
        row.intensityFoundMax12 = candidateScores->integrations.at(11);
        // row.mzPeakLengthsNorm1 = candidateScores->featuresArray[Features::MzPeakLengthsNorm1];
        // row.mzPeakLengthsNorm2 = candidateScores->featuresArray[Features::MzPeakLengthsNorm2];
        // row.mzPeakLengthsNorm3 = candidateScores->featuresArray[Features::MzPeakLengthsNorm3];
        // row.mzPeakLengthsNorm4 = candidateScores->featuresArray[Features::MzPeakLengthsNorm4];
        // row.mzPeakLengthsNorm5 = candidateScores->featuresArray[Features::MzPeakLengthsNorm5];
        // row.mzPeakLengthsNorm6 = candidateScores->featuresArray[Features::MzPeakLengthsNorm6];
        // row.mzPeakLengthsNorm7 = candidateScores->featuresArray[Features::MzPeakLengthsNorm7];
        // row.mzPeakLengthsNorm8 = candidateScores->featuresArray[Features::MzPeakLengthsNorm8];
        // row.mzPeakLengthsNorm9 = candidateScores->featuresArray[Features::MzPeakLengthsNorm9];
        // row.mzPeakLengthsNorm10 = candidateScores->featuresArray[Features::MzPeakLengthsNorm10];
        // row.mzPeakLengthsNorm11 = candidateScores->featuresArray[Features::MzPeakLengthsNorm11];
        // row.mzPeakLengthsNorm12 = candidateScores->featuresArray[Features::MzPeakLengthsNorm12];
        row.aminoAcidCountA = candidateScores->featuresArray[Features::AminoAcidCountA];
        row.aminoAcidCountC = candidateScores->featuresArray[Features::AminoAcidCountC];
        row.aminoAcidCountD = candidateScores->featuresArray[Features::AminoAcidCountD];
        row.aminoAcidCountE = candidateScores->featuresArray[Features::AminoAcidCountE];
        row.aminoAcidCountF = candidateScores->featuresArray[Features::AminoAcidCountF];
        row.aminoAcidCountG = candidateScores->featuresArray[Features::AminoAcidCountG];
        row.aminoAcidCountH = candidateScores->featuresArray[Features::AminoAcidCountH];
        row.aminoAcidCountI = candidateScores->featuresArray[Features::AminoAcidCountI];
        row.aminoAcidCountK = candidateScores->featuresArray[Features::AminoAcidCountK];
        row.aminoAcidCountL = candidateScores->featuresArray[Features::AminoAcidCountL];
        row.aminoAcidCountM = candidateScores->featuresArray[Features::AminoAcidCountM];
        row.aminoAcidCountN = candidateScores->featuresArray[Features::AminoAcidCountN];
        row.aminoAcidCountP = candidateScores->featuresArray[Features::AminoAcidCountP];
        row.aminoAcidCountQ = candidateScores->featuresArray[Features::AminoAcidCountQ];
        row.aminoAcidCountR = candidateScores->featuresArray[Features::AminoAcidCountR];
        row.aminoAcidCountS = candidateScores->featuresArray[Features::AminoAcidCountS];
        row.aminoAcidCountT = candidateScores->featuresArray[Features::AminoAcidCountT];
        row.aminoAcidCountV = candidateScores->featuresArray[Features::AminoAcidCountV];
        row.aminoAcidCountW = candidateScores->featuresArray[Features::AminoAcidCountW];
        row.aminoAcidCountY = candidateScores->featuresArray[Features::AminoAcidCountY];
        row.aminoAcidCountB = candidateScores->featuresArray[Features::AminoAcidCountB];
        row.aminoAcidCountJ = candidateScores->featuresArray[Features::AminoAcidCountJ];
        row.aminoAcidCountO = candidateScores->featuresArray[Features::AminoAcidCountO];
        row.aminoAcidCountU = candidateScores->featuresArray[Features::AminoAcidCountU];
        row.aminoAcidCountX = candidateScores->featuresArray[Features::AminoAcidCountX];
        row.aminoAcidCountZ = candidateScores->featuresArray[Features::AminoAcidCountZ];
        row.mzFoundStDev1 = candidateScores->featuresArray[Features::MzFoundStDev1];
        row.mzFoundStDev2 = candidateScores->featuresArray[Features::MzFoundStDev2];
        row.mzFoundStDev3 = candidateScores->featuresArray[Features::MzFoundStDev3];
        row.mzFoundStDev4 = candidateScores->featuresArray[Features::MzFoundStDev4];
        row.mzFoundStDev5 = candidateScores->featuresArray[Features::MzFoundStDev5];
        row.mzFoundStDev6 = candidateScores->featuresArray[Features::MzFoundStDev6];
        // row.mzFoundStDev7 = candidateScores->featuresArray[Features::MzFoundStDev7];
        // row.mzFoundStDev8 = candidateScores->featuresArray[Features::MzFoundStDev8];
        // row.mzFoundStDev9 = candidateScores->featuresArray[Features::MzFoundStDev9];
        // row.mzFoundStDev10 = candidateScores->featuresArray[Features::MzFoundStDev10];
        // row.mzFoundStDev11 = candidateScores->featuresArray[Features::MzFoundStDev11];
        // row.mzFoundStDev12 = candidateScores->featuresArray[Features::MzFoundStDev12];
        row.targetKey = candidateScores->targetKey;
        row.peptideStringWithMods = candidateScores->isDecoy
                ? AminoAcids::mutatePenultimatePeptideResidues(candidateScores->targetDecoyCandidatePair->peptideStringWithMods())
                : candidateScores->targetDecoyCandidatePair->peptideStringWithMods();

        row.proteinGroup = candidateScores->proteinGroup;
        row.isDecoy = candidateScores->isDecoy;
        row.scanNumber = candidateScores->scanNumber;
        row.scanTime = candidateScores->scanTime;
        row.scanTimeStart = candidateScores->scanTimeStart;
        row.scanTimeEnd = candidateScores->scanTimeEnd;
        row.classifierScore = candidateScores->classifierScore;
        row.discriminantScore = candidateScores->discriminantScore;
        row.qValue = candidateScores->qValue;
        row.decoyRatio = candidateScores->decoyRatio;

        row.altTargetKeyIdDiscScoreChargeOG_alt = candidateScores->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt];
        row.altTargetKeyIdDiscScoreCharge1_1 = candidateScores->featuresArray[AltTargetKeyIdDiscScoreCharge1_1];
        row.altTargetKeyIdDiscScoreCharge1_2 = candidateScores->featuresArray[AltTargetKeyIdDiscScoreCharge1_2];
        row.altTargetKeyIdDiscScoreCharge2_1 = candidateScores->featuresArray[AltTargetKeyIdDiscScoreCharge2_1];
        row.altTargetKeyIdDiscScoreCharge2_2 = candidateScores->featuresArray[AltTargetKeyIdDiscScoreCharge2_2];
        row.altTargetKeyIdDiscScoreCharge3_1 = candidateScores->featuresArray[AltTargetKeyIdDiscScoreCharge3_1];
        row.altTargetKeyIdDiscScoreCharge3_2 = candidateScores->featuresArray[AltTargetKeyIdDiscScoreCharge3_2];
        row.altTargetKeyIdDiscScoreCharge4_1 = candidateScores->featuresArray[AltTargetKeyIdDiscScoreCharge4_1];
        row.altTargetKeyIdDiscScoreCharge4_2 = candidateScores->featuresArray[AltTargetKeyIdDiscScoreCharge4_2];

        row.altTargetKeyIdTimeDeltaCharge1_1 = candidateScores->featuresArray[AltTargetKeyIdTimeDeltaCharge1_1];
        row.altTargetKeyIdTimeDeltaCharge2_1 = candidateScores->featuresArray[AltTargetKeyIdTimeDeltaCharge2_1];
        row.altTargetKeyIdTimeDeltaCharge3_1 = candidateScores->featuresArray[AltTargetKeyIdTimeDeltaCharge3_1];
        row.altTargetKeyIdTimeDeltaCharge4_1 = candidateScores->featuresArray[AltTargetKeyIdTimeDeltaCharge4_1];

        row.ms1MzMeanFound100 = candidateScores->featuresArray[Ms1MzMeanFound100];
        row.ms1MzMeanFound45 = candidateScores->featuresArray[Ms1MzMeanFound45];
        row.ms1MzMeanFoundPreMono = candidateScores->featuresArray[Ms1MzMeanFoundPreMono];
        row.ms1MzMeanFoundIso1 = candidateScores->featuresArray[Ms1MzMeanFoundIso1];
        row.ms1MzMeanFoundIso2 = candidateScores->featuresArray[Ms1MzMeanFoundIso2];
        row.ms1MzMeanFound100PPM = candidateScores->featuresArray[Ms1MzMeanFound100PPM];
        row.ms1MzMeanFound45PPM = candidateScores->featuresArray[Ms1MzMeanFound45PPM];
        row.ms1MzMeanFoundPreMonoPPM = candidateScores->featuresArray[Ms1MzMeanFoundPreMonoPPM];
        row.ms1MzMeanFoundIso1PPM = candidateScores->featuresArray[Ms1MzMeanFoundIso1PPM];
        row.ms1MzMeanFoundIso2PPM = candidateScores->featuresArray[Ms1MzMeanFoundIso2PPM];
        row.ms1MzStDevFound100 = candidateScores->featuresArray[Ms1MzStDevFound100];
        row.ms1MzStDevFound45 = candidateScores->featuresArray[Ms1MzStDevFound45];
        row.ms1MzStDevFoundPreMono = candidateScores->featuresArray[Ms1MzStDevFoundPreMono];
        row.ms1MzStDevFoundIso1 = candidateScores->featuresArray[Ms1MzStDevFoundIso1];
        row.ms1MzStDevFoundIso2 = candidateScores->featuresArray[Ms1MzStDevFoundIso2];
        row.ms1IntensityFound100 = candidateScores->featuresArray[Ms1IntensityFound100];
        row.ms1IntensityFound45 = candidateScores->featuresArray[Ms1IntensityFound45];
        row.ms1IntensityFoundPreMono = candidateScores->featuresArray[Ms1IntensityFoundPreMono];
        row.ms1IntensityFoundIso1 = candidateScores->featuresArray[Ms1IntensityFoundIso1];
        row.ms1IntensityFoundIso2 = candidateScores->featuresArray[Ms1IntensityFoundIso2];

        row.cosineSimSpectrumOverTime = candidateScores->featuresArray[CosineSimSpectrumOverTime];
        row.cosineSimSpectrumOverTimeCubed = candidateScores->featuresArray[CosineSimSpectrumOverTimeCubed];
        row.cosineSimSpectrumStDev = candidateScores->featuresArray[CosineSimSpectrumStDev];
        row.cosineSimSum100MS1 = candidateScores->featuresArray[CosineSimSum100MS1];
        row.ms1Averagine = candidateScores->featuresArray[MS1Averagine];
        row.cosineSimSum100Window1p5X = candidateScores->featuresArray[CosineSimSum100Window1p5X];
        row.cosineSimSum100Window2X = candidateScores->featuresArray[CosineSimSum100Window2X];
        row.totalIntensityPeakHeights = candidateScores->featuresArray[TotalIntensityPeakHeights];
        row.totalIntensityRaw = std::accumulate(candidateScores->integrations.begin(), candidateScores->integrations.end(), 0.0f);
        row.targetWindowLocation = candidateScores->featuresArray[TargetWindowLocation];

        const QVector<MS2Ion> &ms2Ions = candidateScores->isDecoy
                                       ? candidateScores->targetDecoyCandidatePair->ms2IonsDecoy()
                                       : candidateScores->targetDecoyCandidatePair->ms2IonsTarget();

        for (int i = 0; i < ms2Ions.size(); ++i) {

            switch (i) {
                case 0:
                    row.mzSearched1 = ms2Ions.at(0).mz;
                    break;
                case 1:
                    row.mzSearched2 = ms2Ions.at(1).mz;
                    break;
                case 2:
                    row.mzSearched3 = ms2Ions.at(2).mz;
                    break;
                case 3:
                    row.mzSearched4 = ms2Ions.at(3).mz;
                    break;
                case 4:
                    row.mzSearched5 = ms2Ions.at(4).mz;
                    break;
                case 5:
                    row.mzSearched6 = ms2Ions.at(5).mz;
                    break;
                case 6:
                    row.mzSearched7 = ms2Ions.at(6).mz;
                    break;
                case 7:
                    row.mzSearched8 = ms2Ions.at(7).mz;
                    break;
                case 8:
                    row.mzSearched9 = ms2Ions.at(8).mz;
                    break;
                case 9:
                    row.mzSearched10 = ms2Ions.at(9).mz;
                    break;
                case 10:
                    row.mzSearched11 = ms2Ions.at(10).mz;
                    break;
                case 11:
                    row.mzSearched12 = ms2Ions.at(11).mz;
                    break;
                default:
                    break;
            }

        }

        return row;
    }

    static QVector<float> featuresArrayFromCandidateScoresReaderRow(const CandidateScoresReaderRow &candidateScoresReaderRow) {

        QVector<float> featuresArray(FeaturesSize);

        featuresArray[Features::CosineSimSum100] = candidateScoresReaderRow.cosineSimSum100;
        featuresArray[Features::CosineSimSum100GreaterThan80] = candidateScoresReaderRow.cosineSimSum100Greater80;
        featuresArray[Features::AllignedMaxIndexesCount] = candidateScoresReaderRow.allignedMaxIndexesCount;
        featuresArray[Features::CosineSim100MS1] = candidateScoresReaderRow.cosineSim100MS1;
        featuresArray[Features::CosineSimSpectrumCubed] = candidateScoresReaderRow.cosineSimSpectrumCubed;
        featuresArray[Features::KlDivSpectrumCubeRoot] = candidateScoresReaderRow.klDivSpectrumCubeRoot;
        featuresArray[Features::CosineSim45MS1] = candidateScoresReaderRow.cosineSimSum45;
        featuresArray[Features::CosineSimSumTop] = candidateScoresReaderRow.cosineSimSumTop;
        featuresArray[Features::CosineSimSumBottom] = candidateScoresReaderRow.cosineSimSumBottom;
        featuresArray[Features::TopBottomRatio] = candidateScoresReaderRow.topBottomRatio;
        featuresArray[Features::TopBottomRatioNorm] = candidateScoresReaderRow.topBottomRatioNorm;
        featuresArray[Features::Charge] = candidateScoresReaderRow.charge;
        featuresArray[Features::Mass] = candidateScoresReaderRow.mass;
        featuresArray[Features::ScanTimeDelta] = candidateScoresReaderRow.scanTimeDelta;
        featuresArray[Features::ScanTimePd] = candidateScoresReaderRow.scanTimePd;
        featuresArray[Features::ScanIonCount] = candidateScoresReaderRow.scanIonCount;
        featuresArray[Features::MzNorm] = candidateScoresReaderRow.mzNorm;
        featuresArray[Features::KlDivSpectrum] = candidateScoresReaderRow.klDivSpectrum;
        featuresArray[Features::CosineSimSpectrum] = candidateScoresReaderRow.cosineSimSpectrum;
        featuresArray[Features::CosineSim45MS1] = candidateScoresReaderRow.cosineSim45MS1;
        featuresArray[Features::CosineSim100MS1PreMono] = candidateScoresReaderRow.cosineSim100MS1PreMono;
        featuresArray[Features::CosineSim100MS1Iso1] = candidateScoresReaderRow.cosineSim100MS1Iso1;
        featuresArray[Features::CosineSim100MS1Iso2] = candidateScoresReaderRow.cosineSim100MS1Iso2;
        featuresArray[Features::PeptideLengthNorm] = candidateScoresReaderRow.peptideLengthNorm;
        featuresArray[Features::ScanTimePredicted] = candidateScoresReaderRow.scanTimePredicted;
        featuresArray[Features::TheoFragmentCount] = candidateScoresReaderRow.theoFragmentCount;
        featuresArray[Features::TotalIntensityLog] = candidateScoresReaderRow.totalIntensityLog;
        featuresArray[Features::PeakShapeRatio1] = candidateScoresReaderRow.peakShapeRatio1;
        featuresArray[Features::PeakShapeRatio2] = candidateScoresReaderRow.peakShapeRatio2;
        featuresArray[Features::PeakShapeRatio3] = candidateScoresReaderRow.peakShapeRatio3;
        featuresArray[Features::ShadowsCosineSimSum] = candidateScoresReaderRow.shadowsCosineSimSum;
        featuresArray[Features::IRTPredicted] = candidateScoresReaderRow.iRtPredicted;
        featuresArray[Features::CosineSimToAnchor1] = candidateScoresReaderRow.cosineSimToAnchor1;
        featuresArray[Features::CosineSimToAnchor2] = candidateScoresReaderRow.cosineSimToAnchor2;
        featuresArray[Features::CosineSimToAnchor3] = candidateScoresReaderRow.cosineSimToAnchor3;
        featuresArray[Features::CosineSimToAnchor4] = candidateScoresReaderRow.cosineSimToAnchor4;
        featuresArray[Features::CosineSimToAnchor5] = candidateScoresReaderRow.cosineSimToAnchor5;
        featuresArray[Features::CosineSimToAnchor6] = candidateScoresReaderRow.cosineSimToAnchor6;
        featuresArray[Features::CosineSimToAnchor7] = candidateScoresReaderRow.cosineSimToAnchor7;
        featuresArray[Features::CosineSimToAnchor8] = candidateScoresReaderRow.cosineSimToAnchor8;
        featuresArray[Features::CosineSimToAnchor9] = candidateScoresReaderRow.cosineSimToAnchor9;
        featuresArray[Features::CosineSimToAnchor10] = candidateScoresReaderRow.cosineSimToAnchor10;
        featuresArray[Features::CosineSimToAnchor11] = candidateScoresReaderRow.cosineSimToAnchor11;
        featuresArray[Features::CosineSimToAnchor12] = candidateScoresReaderRow.cosineSimToAnchor12;
        // featuresArray[Features::CosineSimShadowsToAnchor1] = candidateScoresReaderRow.cosineSimShadowsToAnchor1;
        // featuresArray[Features::CosineSimShadowsToAnchor2] = candidateScoresReaderRow.cosineSimShadowsToAnchor2;
        // featuresArray[Features::CosineSimShadowsToAnchor3] = candidateScoresReaderRow.cosineSimShadowsToAnchor3;
        // featuresArray[Features::CosineSimShadowsToAnchor4] = candidateScoresReaderRow.cosineSimShadowsToAnchor4;
        // featuresArray[Features::CosineSimShadowsToAnchor5] = candidateScoresReaderRow.cosineSimShadowsToAnchor5;
        // featuresArray[Features::CosineSimShadowsToAnchor6] = candidateScoresReaderRow.cosineSimShadowsToAnchor6;
        // featuresArray[Features::CosineSimShadowsToAnchor7] = candidateScoresReaderRow.cosineSimShadowsToAnchor7;
        // featuresArray[Features::CosineSimShadowsToAnchor8] = candidateScoresReaderRow.cosineSimShadowsToAnchor8;
        // featuresArray[Features::CosineSimShadowsToAnchor9] = candidateScoresReaderRow.cosineSimShadowsToAnchor9;
        // featuresArray[Features::CosineSimShadowsToAnchor10] = candidateScoresReaderRow.cosineSimShadowsToAnchor10;
        // featuresArray[Features::CosineSimShadowsToAnchor11] = candidateScoresReaderRow.cosineSimShadowsToAnchor11;
        // featuresArray[Features::CosineSimShadowsToAnchor12] = candidateScoresReaderRow.cosineSimShadowsToAnchor12;
        featuresArray[Features::MzFoundMean1] = candidateScoresReaderRow.mzFoundMean1;
        featuresArray[Features::MzFoundMean2] = candidateScoresReaderRow.mzFoundMean2;
        featuresArray[Features::MzFoundMean3] = candidateScoresReaderRow.mzFoundMean3;
        featuresArray[Features::MzFoundMean4] = candidateScoresReaderRow.mzFoundMean4;
        featuresArray[Features::MzFoundMean5] = candidateScoresReaderRow.mzFoundMean5;
        featuresArray[Features::MzFoundMean6] = candidateScoresReaderRow.mzFoundMean6;
        featuresArray[Features::MzFoundMean7] = candidateScoresReaderRow.mzFoundMean7;
        featuresArray[Features::MzFoundMean8] = candidateScoresReaderRow.mzFoundMean8;
        featuresArray[Features::MzFoundMean9] = candidateScoresReaderRow.mzFoundMean9;
        featuresArray[Features::MzFoundMean10] = candidateScoresReaderRow.mzFoundMean10;
        featuresArray[Features::MzFoundMean11] = candidateScoresReaderRow.mzFoundMean11;
        featuresArray[Features::MzFoundMean12] = candidateScoresReaderRow.mzFoundMean12;
        featuresArray[Features::IntensityFoundMax1] = candidateScoresReaderRow.intensityFoundMax1;
        featuresArray[Features::IntensityFoundMax2] = candidateScoresReaderRow.intensityFoundMax2;
        featuresArray[Features::IntensityFoundMax3] = candidateScoresReaderRow.intensityFoundMax3;
        featuresArray[Features::IntensityFoundMax4] = candidateScoresReaderRow.intensityFoundMax4;
        featuresArray[Features::IntensityFoundMax5] = candidateScoresReaderRow.intensityFoundMax5;
        featuresArray[Features::IntensityFoundMax6] = candidateScoresReaderRow.intensityFoundMax6;
        featuresArray[Features::IntensityFoundMax7] = candidateScoresReaderRow.intensityFoundMax7;
        featuresArray[Features::IntensityFoundMax8] = candidateScoresReaderRow.intensityFoundMax8;
        featuresArray[Features::IntensityFoundMax9] = candidateScoresReaderRow.intensityFoundMax9;
        featuresArray[Features::IntensityFoundMax10] = candidateScoresReaderRow.intensityFoundMax10;
        featuresArray[Features::IntensityFoundMax11] = candidateScoresReaderRow.intensityFoundMax11;
        featuresArray[Features::IntensityFoundMax12] = candidateScoresReaderRow.intensityFoundMax12;
        // featuresArray[Features::MzPeakLengthsNorm1] = candidateScoresReaderRow.mzPeakLengthsNorm1;
        // featuresArray[Features::MzPeakLengthsNorm2] = candidateScoresReaderRow.mzPeakLengthsNorm2;
        // featuresArray[Features::MzPeakLengthsNorm3] = candidateScoresReaderRow.mzPeakLengthsNorm3;
        // featuresArray[Features::MzPeakLengthsNorm4] = candidateScoresReaderRow.mzPeakLengthsNorm4;
        // featuresArray[Features::MzPeakLengthsNorm5] = candidateScoresReaderRow.mzPeakLengthsNorm5;
        // featuresArray[Features::MzPeakLengthsNorm6] = candidateScoresReaderRow.mzPeakLengthsNorm6;
        // featuresArray[Features::MzPeakLengthsNorm7] = candidateScoresReaderRow.mzPeakLengthsNorm7;
        // featuresArray[Features::MzPeakLengthsNorm8] = candidateScoresReaderRow.mzPeakLengthsNorm8;
        // featuresArray[Features::MzPeakLengthsNorm9] = candidateScoresReaderRow.mzPeakLengthsNorm9;
        // featuresArray[Features::MzPeakLengthsNorm10] = candidateScoresReaderRow.mzPeakLengthsNorm10;
        // featuresArray[Features::MzPeakLengthsNorm11] = candidateScoresReaderRow.mzPeakLengthsNorm11;
        // featuresArray[Features::MzPeakLengthsNorm12] = candidateScoresReaderRow.mzPeakLengthsNorm12;
        featuresArray[Features::AminoAcidCountA] = candidateScoresReaderRow.aminoAcidCountA;
        featuresArray[Features::AminoAcidCountC] = candidateScoresReaderRow.aminoAcidCountC;
        featuresArray[Features::AminoAcidCountD] = candidateScoresReaderRow.aminoAcidCountD;
        featuresArray[Features::AminoAcidCountE] = candidateScoresReaderRow.aminoAcidCountE;
        featuresArray[Features::AminoAcidCountF] = candidateScoresReaderRow.aminoAcidCountF;
        featuresArray[Features::AminoAcidCountG] = candidateScoresReaderRow.aminoAcidCountG;
        featuresArray[Features::AminoAcidCountH] = candidateScoresReaderRow.aminoAcidCountH;
        featuresArray[Features::AminoAcidCountI] = candidateScoresReaderRow.aminoAcidCountI;
        featuresArray[Features::AminoAcidCountK] = candidateScoresReaderRow.aminoAcidCountK;
        featuresArray[Features::AminoAcidCountL] = candidateScoresReaderRow.aminoAcidCountL;
        featuresArray[Features::AminoAcidCountM] = candidateScoresReaderRow.aminoAcidCountM;
        featuresArray[Features::AminoAcidCountN] = candidateScoresReaderRow.aminoAcidCountN;
        featuresArray[Features::AminoAcidCountP] = candidateScoresReaderRow.aminoAcidCountP;
        featuresArray[Features::AminoAcidCountQ] = candidateScoresReaderRow.aminoAcidCountQ;
        featuresArray[Features::AminoAcidCountR] = candidateScoresReaderRow.aminoAcidCountR;
        featuresArray[Features::AminoAcidCountS] = candidateScoresReaderRow.aminoAcidCountS;
        featuresArray[Features::AminoAcidCountT] = candidateScoresReaderRow.aminoAcidCountT;
        featuresArray[Features::AminoAcidCountV] = candidateScoresReaderRow.aminoAcidCountV;
        featuresArray[Features::AminoAcidCountW] = candidateScoresReaderRow.aminoAcidCountW;
        featuresArray[Features::AminoAcidCountY] = candidateScoresReaderRow.aminoAcidCountY;
        featuresArray[Features::AminoAcidCountB] = candidateScoresReaderRow.aminoAcidCountB;
        featuresArray[Features::AminoAcidCountJ] = candidateScoresReaderRow.aminoAcidCountJ;
        featuresArray[Features::AminoAcidCountO] = candidateScoresReaderRow.aminoAcidCountO;
        featuresArray[Features::AminoAcidCountU] = candidateScoresReaderRow.aminoAcidCountU;
        featuresArray[Features::AminoAcidCountX] = candidateScoresReaderRow.aminoAcidCountX;
        featuresArray[Features::AminoAcidCountZ] = candidateScoresReaderRow.aminoAcidCountZ;
        featuresArray[Features::MzFoundStDev1] = candidateScoresReaderRow.mzFoundStDev1;
        featuresArray[Features::MzFoundStDev2] = candidateScoresReaderRow.mzFoundStDev2;
        featuresArray[Features::MzFoundStDev3] = candidateScoresReaderRow.mzFoundStDev3;
        featuresArray[Features::MzFoundStDev4] = candidateScoresReaderRow.mzFoundStDev4;
        featuresArray[Features::MzFoundStDev5] = candidateScoresReaderRow.mzFoundStDev5;
        featuresArray[Features::MzFoundStDev6] = candidateScoresReaderRow.mzFoundStDev6;
        // featuresArray[Features::MzFoundStDev7] = candidateScoresReaderRow.mzFoundStDev7;
        // featuresArray[Features::MzFoundStDev8] = candidateScoresReaderRow.mzFoundStDev8;
        // featuresArray[Features::MzFoundStDev9] = candidateScoresReaderRow.mzFoundStDev9;
        // featuresArray[Features::MzFoundStDev10] = candidateScoresReaderRow.mzFoundStDev10;
        // featuresArray[Features::MzFoundStDev11] = candidateScoresReaderRow.mzFoundStDev11;
        // featuresArray[Features::MzFoundStDev12] = candidateScoresReaderRow.mzFoundStDev12;
        featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt] = candidateScoresReaderRow.altTargetKeyIdDiscScoreChargeOG_alt;
        featuresArray[Features::AltTargetKeyIdDiscScoreCharge1_1] = candidateScoresReaderRow.altTargetKeyIdDiscScoreCharge1_1;
        featuresArray[Features::AltTargetKeyIdDiscScoreCharge1_2] = candidateScoresReaderRow.altTargetKeyIdDiscScoreCharge1_2;
        featuresArray[Features::AltTargetKeyIdDiscScoreCharge2_1] = candidateScoresReaderRow.altTargetKeyIdDiscScoreCharge2_1;
        featuresArray[Features::AltTargetKeyIdDiscScoreCharge2_2] = candidateScoresReaderRow.altTargetKeyIdDiscScoreCharge2_2;
        featuresArray[Features::AltTargetKeyIdDiscScoreCharge3_1] = candidateScoresReaderRow.altTargetKeyIdDiscScoreCharge3_1;
        featuresArray[Features::AltTargetKeyIdDiscScoreCharge3_2] = candidateScoresReaderRow.altTargetKeyIdDiscScoreCharge3_2;
        featuresArray[Features::AltTargetKeyIdDiscScoreCharge3_1] = candidateScoresReaderRow.altTargetKeyIdDiscScoreCharge4_1;
        featuresArray[Features::AltTargetKeyIdDiscScoreCharge3_2] = candidateScoresReaderRow.altTargetKeyIdDiscScoreCharge4_2;
        featuresArray[Features::AltTargetKeyIdTimeDeltaCharge1_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge1_1;
        featuresArray[Features::AltTargetKeyIdTimeDeltaCharge2_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge2_1;
        featuresArray[Features::AltTargetKeyIdTimeDeltaCharge3_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge3_1;
        featuresArray[Features::AltTargetKeyIdTimeDeltaCharge4_1] = candidateScoresReaderRow.altTargetKeyIdTimeDeltaCharge4_1;
        featuresArray[Features::Ms1MzMeanFound100] = candidateScoresReaderRow.ms1MzMeanFound100;
        featuresArray[Features::Ms1MzMeanFound45] = candidateScoresReaderRow.ms1MzMeanFound45;
        featuresArray[Features::Ms1MzMeanFoundPreMono] = candidateScoresReaderRow.ms1MzMeanFoundPreMono;
        featuresArray[Features::Ms1MzMeanFoundIso1] = candidateScoresReaderRow.ms1MzMeanFoundIso1;
        featuresArray[Features::Ms1MzMeanFoundIso2] = candidateScoresReaderRow.ms1MzMeanFoundIso2;
        featuresArray[Features::Ms1MzMeanFound100PPM] = candidateScoresReaderRow.ms1MzMeanFound100PPM;
        featuresArray[Features::Ms1MzMeanFound45PPM] = candidateScoresReaderRow.ms1MzMeanFound45PPM;
        featuresArray[Features::Ms1MzMeanFoundPreMonoPPM] = candidateScoresReaderRow.ms1MzMeanFoundPreMonoPPM;
        featuresArray[Features::Ms1MzMeanFoundIso1PPM] = candidateScoresReaderRow.ms1MzMeanFoundIso1PPM;
        featuresArray[Features::Ms1MzMeanFoundIso2PPM] = candidateScoresReaderRow.ms1MzMeanFoundIso2PPM;
        featuresArray[Features::Ms1MzStDevFound100] = candidateScoresReaderRow.ms1MzStDevFound100;
        featuresArray[Features::Ms1MzStDevFound45] = candidateScoresReaderRow.ms1MzStDevFound45;
        featuresArray[Features::Ms1MzStDevFoundPreMono] = candidateScoresReaderRow.ms1MzStDevFoundPreMono;
        featuresArray[Features::Ms1MzStDevFoundIso1] = candidateScoresReaderRow.ms1MzStDevFoundIso1;
        featuresArray[Features::Ms1MzStDevFoundIso2] = candidateScoresReaderRow.ms1MzStDevFoundIso2;
        featuresArray[Features::Ms1IntensityFound100] = candidateScoresReaderRow.ms1IntensityFound100;
        featuresArray[Features::Ms1IntensityFound45] = candidateScoresReaderRow.ms1IntensityFound45;
        featuresArray[Features::Ms1IntensityFoundPreMono] = candidateScoresReaderRow.ms1IntensityFoundPreMono;
        featuresArray[Features::Ms1IntensityFoundIso1] = candidateScoresReaderRow.ms1IntensityFoundIso1;
        featuresArray[Features::Ms1IntensityFoundIso2] = candidateScoresReaderRow.ms1IntensityFoundIso2;

        featuresArray[CosineSimSpectrumOverTime] = candidateScoresReaderRow.cosineSimSpectrumOverTime;
        featuresArray[CosineSimSpectrumOverTimeCubed] = candidateScoresReaderRow.cosineSimSpectrumOverTimeCubed;
        featuresArray[CosineSimSpectrumStDev] = candidateScoresReaderRow.cosineSimSpectrumStDev;
        featuresArray[CosineSimSum100MS1] = candidateScoresReaderRow.cosineSimSum100MS1;
        featuresArray[MS1Averagine] = candidateScoresReaderRow.ms1Averagine;
        featuresArray[CosineSimSum100Window1p5X] = candidateScoresReaderRow.cosineSimSum100Window1p5X;
        featuresArray[CosineSimSum100Window2X] = candidateScoresReaderRow.cosineSimSum100Window2X;
        featuresArray[TotalIntensityPeakHeights] = candidateScoresReaderRow.totalIntensityPeakHeights;
        featuresArray[TotalIntensityRaw] = candidateScoresReaderRow.totalIntensityRaw;
        featuresArray[TargetWindowLocation] = candidateScoresReaderRow.targetWindowLocation;

        return featuresArray;
    }

};

using CandidateScoresTarget = CandidateScores;
using CandidateScoresDecoy = CandidateScores;


#endif //PYTHIADIACPP_CANDIDATESCORES_H
