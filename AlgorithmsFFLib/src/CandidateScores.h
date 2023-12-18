//
// Created by anichols on 10/20/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORES_H
#define PYTHIADIACPP_CANDIDATESCORES_H


#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "TargetDecoyCandidatePair.h"
#include "ParquetReader.h"


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
        TopBottomRatioNorm,
        ChargeNorm,
        Mass,
        ScanTimeDelta,
        ScanTimeRange,
        ScanTimePd,
        ScanIonCount,
        MzNorm,
//        KlDivSum,
        KlDivSpectrum,
        CosineSimSpectrum,
//        CosineSim45MS1,
//        CosineSim20MS1,
        CosineSim100MS1PreMono,
        CosineSim100MS1Iso1,
        CosineSim100MS1Iso2,
        PeptideLengthNorm,
//        ScanTimePredicted,
        TheoFragmentCount,
        TotalIntensityLog,
//        DiscriminateScore,
//        ScanNumberCandidateCount,
        PeakShapeRatio1,
        PeakShapeRatio2,
        PeakShapeRatio3,
        UnfragPrecursorCosineSim,
        ShadowsCosineSimSum,
        IRTPredicted,
        CosineSimToAnchor1,
        CosineSimToAnchor2,
        CosineSimToAnchor3,
        CosineSimToAnchor4,
        CosineSimToAnchor5,
        CosineSimToAnchor6,
        CosineSimToAnchor7,
        CosineSimToAnchor8,
        CosineSimToAnchor9,
        CosineSimToAnchor10,
        CosineSimToAnchor11,
        CosineSimToAnchor12,
        CosineSimShadowsToAnchor1,
        CosineSimShadowsToAnchor2,
        CosineSimShadowsToAnchor3,
        CosineSimShadowsToAnchor4,
        CosineSimShadowsToAnchor5,
        CosineSimShadowsToAnchor6,
        CosineSimShadowsToAnchor7,
        CosineSimShadowsToAnchor8,
        CosineSimShadowsToAnchor9,
        CosineSimShadowsToAnchor10,
        CosineSimShadowsToAnchor11,
        CosineSimShadowsToAnchor12,
        ShadowsIntensityRatio1,
        ShadowsIntensityRatio2,
        ShadowsIntensityRatio3,
        ShadowsIntensityRatio4,
        ShadowsIntensityRatio5,
        ShadowsIntensityRatio6,
        ShadowsIntensityRatio7,
        ShadowsIntensityRatio8,
        ShadowsIntensityRatio9,
        ShadowsIntensityRatio10,
        ShadowsIntensityRatio11,
        ShadowsIntensityRatio12,
        MzSearched1,
        MzSearched2,
        MzSearched3,
        MzSearched4,
        MzSearched5,
        MzSearched6,
        MzSearched7,
        MzSearched8,
        MzSearched9,
        MzSearched10,
        MzSearched11,
        MzSearched12,
        TheoIntensity1,
        TheoIntensity2,
        TheoIntensity3,
        TheoIntensity4,
        TheoIntensity5,
        TheoIntensity6,
        TheoIntensity7,
        TheoIntensity8,
        TheoIntensity9,
        TheoIntensity10,
        TheoIntensity11,
        TheoIntensity12,
        MzFoundMean1,
        MzFoundMean2,
        MzFoundMean3,
        MzFoundMean4,
        MzFoundMean5,
        MzFoundMean6,
        MzFoundMean7,
        MzFoundMean8,
        MzFoundMean9,
        MzFoundMean10,
        MzFoundMean11,
        MzFoundMean12,
        IntensityFoundMax1,
        IntensityFoundMax2,
        IntensityFoundMax3,
        IntensityFoundMax4,
        IntensityFoundMax5,
        IntensityFoundMax6,
        IntensityFoundMax7,
        IntensityFoundMax8,
        IntensityFoundMax9,
        IntensityFoundMax10,
        IntensityFoundMax11,
        IntensityFoundMax12,
        MzPeakLengthsNorm1,
        MzPeakLengthsNorm2,
        MzPeakLengthsNorm3,
        MzPeakLengthsNorm4,
        MzPeakLengthsNorm5,
        MzPeakLengthsNorm6,
        MzPeakLengthsNorm7,
        MzPeakLengthsNorm8,
        MzPeakLengthsNorm9,
        MzPeakLengthsNorm10,
        MzPeakLengthsNorm11,
        MzPeakLengthsNorm12,
//        ColumnApexIndexRatiosToAnchor1,
//        ColumnApexIndexRatiosToAnchor2,
//        ColumnApexIndexRatiosToAnchor3,
//        ColumnApexIndexRatiosToAnchor4,
//        ColumnApexIndexRatiosToAnchor5,
//        ColumnApexIndexRatiosToAnchor6,
//        ColumnApexIndexRatiosToAnchor7,
//        ColumnApexIndexRatiosToAnchor8,
//        ColumnApexIndexRatiosToAnchor9,
//        ColumnApexIndexRatiosToAnchor10,
//        ColumnApexIndexRatiosToAnchor11,
//        ColumnApexIndexRatiosToAnchor12,
        AminoAcidCountA,
        AminoAcidCountC,
        AminoAcidCountD,
        AminoAcidCountE,
        AminoAcidCountF,
        AminoAcidCountG,
        AminoAcidCountH,
        AminoAcidCountI,
        AminoAcidCountK,
        AminoAcidCountL,
        AminoAcidCountM,
        AminoAcidCountN,
        AminoAcidCountP,
        AminoAcidCountQ,
        AminoAcidCountR,
        AminoAcidCountS,
        AminoAcidCountT,
        AminoAcidCountV,
        AminoAcidCountW,
        AminoAcidCountY,
        MzFoundStDev1,
        MzFoundStDev2,
        MzFoundStDev3,
        MzFoundStDev4,
        MzFoundStDev5,
        MzFoundStDev6,
        MzFoundStDev7,
        MzFoundStDev8,
        MzFoundStDev9,
        MzFoundStDev10,
        MzFoundStDev11,
        MzFoundStDev12,
        FeaturesSize
    };

    CandidateScores() = default;
    ~CandidateScores() = default;

    TargetDecoyCandidatePair *targetDecoyCandidatePair;
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

//    [[nodiscard]] Err featuresArrayEssentials(QVector<float> *vecOutput) const;

    [[nodiscard]] QVector<float>* featuresArrayRef();

//    [[nodiscard]] Err featuresArrayNeuralNet(QVector<float> **vecOutput) const;

    void initFeaturesArray();

};


#endif //PYTHIADIACPP_CANDIDATESCORES_H
