//
// Created by anichols on 1/9/24.
//

#include "DiscriminantScoretron.h"

#include "ClassifierWeightsManager.h"
#include "ParallelUtils.h"

#include <QtConcurrent/QtConcurrent>

#include "EigenUtils.h"


Err DiscriminantScoretron::trainLDAClassifier(
        const QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> &targetDecoyCandidateScoresPair,
        int threadCount,
        int verbosity,
        QVector<float> *weights
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targetDecoyCandidateScoresPair); ree;

    QElapsedTimer et;
    et.start();

    weights->clear();

    const QPair<QVector<FeaturesArrayTargets*>, QVector<FeaturesArrayDecoys*>> unzipped
                                                        = ParallelUtils::unZip(targetDecoyCandidateScoresPair);

    QVector<QVector<float>> A;
    QVector<float> b;
    e = ClassifierWeightsManager::buildDataClassifier2(
        unzipped.first,
        unzipped.second,
        &A,
        &b
        ); ree;

    e = ClassifierWeightsManager::fitWeights(A, b, weights); ree;

    if (verbosity > 0) {
        qDebug() << "fit weights" << et.restart() << "mSec";
        qDebug() << "Weights:" << *weights;
        qDebug() << "b:" << b;
    }

    ERR_RETURN
}

namespace {

    QPair<Err, QVector<float>> applyWeightsLogic(
        const QVector<float> &weights,
        const QVector<FeaturesArray*> &featuresArray
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(weights); rree;
        e = ErrorUtils::isNotEmpty(featuresArray); rree;

        QVector<float> discriminantScores;
        e = ClassifierWeightsManager::applyWeights(
            featuresArray,
            weights,
            &discriminantScores
        ); rree;

        e = ErrorUtils::isEqual(
            discriminantScores.size(),
            featuresArray.size()
            ); rree;

        return {e, discriminantScores};
    }

}//namespace
Err DiscriminantScoretron::applyWeights(
    const QVector<float>& weights,
    int threadCount,
    const QVector<FeaturesArray*> &featuresArray,
    QVector<float> *discriminantScores
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(weights); ree;
    e = ErrorUtils::isNotEmpty(featuresArray); ree;

    QVector<QVector<FeaturesArray*>> featuresArrayPntrsTranched;
    e = ParallelUtils::trancheVectorForParallelizationInOrder(
        featuresArray,
        threadCount,
        0,
        &featuresArrayPntrsTranched
    ); ree;

    const auto applyWeightsBinder = std::bind(
        applyWeightsLogic,
        weights,
        std::placeholders::_1
        );

    QFuture<QPair<Err, QVector<float>>> result = QtConcurrent::mapped(
        featuresArrayPntrsTranched,
        applyWeightsBinder
        );
    result.waitForFinished();

    for (const QPair<Err, QVector<float>> &res : result) {
        e = res.first; ree;
        discriminantScores->append(res.second);
    }

    ERR_RETURN
}

QVector<float> DiscriminantScoretron::scoreVectorLogic(
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            CandidateScores* candidateScores
            ) {

        ERR_INIT

        const QVector<CandidateScores::Features> baseFeatures = {

            CandidateScores::Features::CosineSimSum100GreaterThan80,
            CandidateScores::Features::AllignedMaxIndexesCount,
            CandidateScores::Features::CosineSimSpectrumOverTimeCubed,
            CandidateScores::Features::CosineSimSpectrumStDev,
            CandidateScores::Features::CosineSim100MS1,
            CandidateScores::Features::CosineSim100MS1Iso1,
            CandidateScores::Features::CosineSim100MS1Iso2,
            CandidateScores::Features::CosineSim100MS1PreMono,
            CandidateScores::Features::CosineSimSpectrumCubed,
            CandidateScores::Features::KlDivSpectrumCubeRoot,
            CandidateScores::Features::CosineSimSum45,
            // CandidateScores::Features::CosineSimSum20,//TODO delete
            CandidateScores::Features::CosineSimSumTop,
            CandidateScores::Features::Charge,
            CandidateScores::Features::CosineSimSumBottom,
            CandidateScores::Features::TopBottomRatio,
            CandidateScores::Features::TopBottomRatioNorm
        };

        if (useNeuralNetworkScores) {

#define PRODUCTION
#ifdef PRODUCTION
            QVector<float> &vec = candidateScores->featuresArray;
            return vec;
#else
            const QVector<CandidateScores::Features> nnFeatures = {
                CandidateScores::Features::CosineSimSum100,
                CandidateScores::Features::CosineSimSum100Top12,
                CandidateScores::Features::CosineSimSum100GreaterThan80,
                CandidateScores::Features::AllignedMaxIndexesCount,
                CandidateScores::Features::CosineSim100MS1,
                CandidateScores::Features::CosineSimSpectrumCubed,
                CandidateScores::Features::KlDivSpectrumCubeRoot,
                CandidateScores::Features::CosineSimSum45,
                // CandidateScores::Features::CosineSimSum20,//TODO delete
                CandidateScores::Features::CosineSimSumTop,
                CandidateScores::Features::CosineSimSumBottom,
                CandidateScores::Features::TopBottomRatio,
                CandidateScores::Features::TopBottomRatioNorm,
                CandidateScores::Features::Charge,
                CandidateScores::Features::ScanTimeDelta,
                CandidateScores::Features::ScanTimePd,
                CandidateScores::Features::ScanIonCount,
                CandidateScores::Features::MzNorm,
                CandidateScores::Features::Mass,
                CandidateScores::Features::KlDivSpectrum,
                CandidateScores::Features::CosineSimSpectrum,
                CandidateScores::Features::CosineSim45MS1,
                // CandidateScores::Features::CosineSim20MS1,//TODO delete
                CandidateScores::Features::CosineSim100MS1PreMono,
                CandidateScores::Features::CosineSim100MS1Iso1,
                CandidateScores::Features::CosineSim100MS1Iso2,
                CandidateScores::Features::PeptideLengthNorm,
                CandidateScores::Features::ScanTimePredicted,
                CandidateScores::Features::TheoFragmentCount,
                CandidateScores::Features::TotalIntensityLog,
                CandidateScores::Features::PeakShapeRatio1,
                CandidateScores::Features::PeakShapeRatio2,
                CandidateScores::Features::PeakShapeRatio3,
                CandidateScores::Features::ShadowsCosineSimSum,
                CandidateScores::Features::IRTPredicted,
                CandidateScores::Features::CosineSimToAnchor1,
                CandidateScores::Features::CosineSimToAnchor2,
                CandidateScores::Features::CosineSimToAnchor3,
                CandidateScores::Features::CosineSimToAnchor4,
                CandidateScores::Features::CosineSimToAnchor5,
                CandidateScores::Features::CosineSimToAnchor6,
                CandidateScores::Features::CosineSimToAnchor7,
                CandidateScores::Features::CosineSimToAnchor8,
                CandidateScores::Features::CosineSimToAnchor9,
                CandidateScores::Features::CosineSimToAnchor10,
                CandidateScores::Features::CosineSimToAnchor11,
                CandidateScores::Features::CosineSimToAnchor12,
                CandidateScores::Features::CosineSimShadowsToAnchor1,
                CandidateScores::Features::CosineSimShadowsToAnchor2,
                CandidateScores::Features::CosineSimShadowsToAnchor3,
                CandidateScores::Features::CosineSimShadowsToAnchor4,
                CandidateScores::Features::CosineSimShadowsToAnchor5,
                CandidateScores::Features::CosineSimShadowsToAnchor6,
                CandidateScores::Features::CosineSimShadowsToAnchor7,
                CandidateScores::Features::CosineSimShadowsToAnchor8,
                CandidateScores::Features::CosineSimShadowsToAnchor9,
                CandidateScores::Features::CosineSimShadowsToAnchor10,
                CandidateScores::Features::CosineSimShadowsToAnchor11,
                CandidateScores::Features::CosineSimShadowsToAnchor12,
                CandidateScores::Features::ShadowsIntensityRatio1,
                CandidateScores::Features::ShadowsIntensityRatio2,
                CandidateScores::Features::ShadowsIntensityRatio3,
                CandidateScores::Features::ShadowsIntensityRatio4,
                CandidateScores::Features::ShadowsIntensityRatio5,
                CandidateScores::Features::ShadowsIntensityRatio6,
                CandidateScores::Features::ShadowsIntensityRatio7,
                CandidateScores::Features::ShadowsIntensityRatio8,
                CandidateScores::Features::ShadowsIntensityRatio9,
                CandidateScores::Features::ShadowsIntensityRatio10,
                CandidateScores::Features::ShadowsIntensityRatio11,
                CandidateScores::Features::ShadowsIntensityRatio12,
                CandidateScores::Features::MzSearched1,
                CandidateScores::Features::MzSearched2,
                CandidateScores::Features::MzSearched3,
                CandidateScores::Features::MzSearched4,
                CandidateScores::Features::MzSearched5,
                CandidateScores::Features::MzSearched6,
                CandidateScores::Features::MzSearched7,
                CandidateScores::Features::MzSearched8,
                CandidateScores::Features::MzSearched9,
                CandidateScores::Features::MzSearched10,
                CandidateScores::Features::MzSearched11,
                CandidateScores::Features::MzSearched12,
                CandidateScores::Features::TheoIntensity1,
                CandidateScores::Features::TheoIntensity2,
                CandidateScores::Features::TheoIntensity3,
                CandidateScores::Features::TheoIntensity4,
                CandidateScores::Features::TheoIntensity5,
                CandidateScores::Features::TheoIntensity6,
                CandidateScores::Features::TheoIntensity7,
                CandidateScores::Features::TheoIntensity8,
                CandidateScores::Features::TheoIntensity9,
                CandidateScores::Features::TheoIntensity10,
                CandidateScores::Features::TheoIntensity11,
                CandidateScores::Features::TheoIntensity12,
                CandidateScores::Features::MzFoundMean1,
                CandidateScores::Features::MzFoundMean2,
                CandidateScores::Features::MzFoundMean3,
                CandidateScores::Features::MzFoundMean4,
                CandidateScores::Features::MzFoundMean5,
                CandidateScores::Features::MzFoundMean6,
                CandidateScores::Features::MzFoundMean7,
                CandidateScores::Features::MzFoundMean8,
                CandidateScores::Features::MzFoundMean9,
                CandidateScores::Features::MzFoundMean10,
                CandidateScores::Features::MzFoundMean11,
                CandidateScores::Features::MzFoundMean12,
                CandidateScores::Features::IntensityFoundMax1,
                CandidateScores::Features::IntensityFoundMax2,
                CandidateScores::Features::IntensityFoundMax3,
                CandidateScores::Features::IntensityFoundMax4,
                CandidateScores::Features::IntensityFoundMax5,
                CandidateScores::Features::IntensityFoundMax6,
                CandidateScores::Features::IntensityFoundMax7,
                CandidateScores::Features::IntensityFoundMax8,
                CandidateScores::Features::IntensityFoundMax9,
                CandidateScores::Features::IntensityFoundMax10,
                CandidateScores::Features::IntensityFoundMax11,
                CandidateScores::Features::IntensityFoundMax12,
                CandidateScores::Features::MzPeakLengthsNorm1,
                CandidateScores::Features::MzPeakLengthsNorm2,
                CandidateScores::Features::MzPeakLengthsNorm3,
                CandidateScores::Features::MzPeakLengthsNorm4,
                CandidateScores::Features::MzPeakLengthsNorm5,
                CandidateScores::Features::MzPeakLengthsNorm6,
                CandidateScores::Features::MzPeakLengthsNorm7,
                CandidateScores::Features::MzPeakLengthsNorm8,
                CandidateScores::Features::MzPeakLengthsNorm9,
                CandidateScores::Features::MzPeakLengthsNorm10,
                CandidateScores::Features::MzPeakLengthsNorm11,
                CandidateScores::Features::MzPeakLengthsNorm12,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor1,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor2,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor3,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor4,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor5,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor6,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor7,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor8,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor9,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor10,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor11,
                CandidateScores::Features::ColumnApexIndexRatiosToAnchor12,
                CandidateScores::Features::AminoAcidCountA,
                CandidateScores::Features::AminoAcidCountC,
                CandidateScores::Features::AminoAcidCountD,
                CandidateScores::Features::AminoAcidCountE,
                CandidateScores::Features::AminoAcidCountF,
                CandidateScores::Features::AminoAcidCountG,
                CandidateScores::Features::AminoAcidCountH,
                CandidateScores::Features::AminoAcidCountI,
                CandidateScores::Features::AminoAcidCountK,
                CandidateScores::Features::AminoAcidCountL,
                CandidateScores::Features::AminoAcidCountM,
                CandidateScores::Features::AminoAcidCountN,
                CandidateScores::Features::AminoAcidCountP,
                CandidateScores::Features::AminoAcidCountQ,
                CandidateScores::Features::AminoAcidCountR,
                CandidateScores::Features::AminoAcidCountS,
                CandidateScores::Features::AminoAcidCountT,
                CandidateScores::Features::AminoAcidCountV,
                CandidateScores::Features::AminoAcidCountW,
                CandidateScores::Features::AminoAcidCountY,
                CandidateScores::Features::AminoAcidCountB,
                CandidateScores::Features::AminoAcidCountJ,
                CandidateScores::Features::AminoAcidCountO,
                CandidateScores::Features::AminoAcidCountU,
                CandidateScores::Features::AminoAcidCountX,
                CandidateScores::Features::AminoAcidCountZ,
                CandidateScores::Features::MzFoundStDev1,
                CandidateScores::Features::MzFoundStDev2,
                CandidateScores::Features::MzFoundStDev3,
                CandidateScores::Features::MzFoundStDev4,
                CandidateScores::Features::MzFoundStDev5,
                CandidateScores::Features::MzFoundStDev6,
                CandidateScores::Features::MzFoundStDev7,
                CandidateScores::Features::MzFoundStDev8,
                CandidateScores::Features::MzFoundStDev9,
                CandidateScores::Features::MzFoundStDev10,
                CandidateScores::Features::MzFoundStDev11,
                CandidateScores::Features::MzFoundStDev12,
                CandidateScores::Features::MzAccuracy1,
                CandidateScores::Features::MzAccuracy2,
                CandidateScores::Features::MzAccuracy3,
                CandidateScores::Features::MzAccuracy4,
                CandidateScores::Features::MzAccuracy5,
                CandidateScores::Features::MzAccuracy6,
                CandidateScores::Features::MzAccuracy7,
                CandidateScores::Features::MzAccuracy8,
                CandidateScores::Features::MzAccuracy9,
                CandidateScores::Features::MzAccuracy10,
                CandidateScores::Features::MzAccuracy11,
                CandidateScores::Features::MzAccuracy12,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_1,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_2,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_3,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_1,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_2,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_3,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_1,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_3,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_1,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_2,
                CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_3,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_1,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_2,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_3,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_1,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_2,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_3,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_1,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_2,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_3,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_1,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_2,
                CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_3,
                CandidateScores::Features::Ms1MzMeanFound100,
                CandidateScores::Features::Ms1MzMeanFound45,
                // CandidateScores::Features::Ms1MzMeanFound20,//TODO delete
                CandidateScores::Features::Ms1MzMeanFoundPreMono,
                CandidateScores::Features::Ms1MzMeanFoundIso1,
                CandidateScores::Features::Ms1MzMeanFoundIso2,
                CandidateScores::Features::Ms1MzMeanFound100PPM,
                CandidateScores::Features::Ms1MzMeanFound45PPM,
                // CandidateScores::Features::Ms1MzMeanFound20PPM,//TODO delete
                CandidateScores::Features::Ms1MzMeanFoundPreMonoPPM,
                CandidateScores::Features::Ms1MzMeanFoundIso1PPM,
                CandidateScores::Features::Ms1MzMeanFoundIso2PPM,
                CandidateScores::Features::Ms1MzStDevFound100,
                CandidateScores::Features::Ms1MzStDevFound45,
                // CandidateScores::Features::Ms1MzStDevFound20,//TODO delete
                CandidateScores::Features::Ms1MzStDevFoundPreMono,
                CandidateScores::Features::Ms1MzStDevFoundIso1,
                CandidateScores::Features::Ms1MzStDevFoundIso2,
                CandidateScores::Features::Ms1IntensityFound100,
                CandidateScores::Features::Ms1IntensityFound45,
                // CandidateScores::Features::Ms1IntensityFound20,//TODO delete
                CandidateScores::Features::Ms1IntensityFoundPreMono,
                CandidateScores::Features::Ms1IntensityFoundIso1,
                CandidateScores::Features::Ms1IntensityFoundIso2,
                CandidateScores::Features::CosineSimSpectrumOverTime,
                CandidateScores::Features::CosineSimSpectrumOverTimeCubed,
                CandidateScores::Features::CosineSimSpectrumStDev,
                CandidateScores::Features::CosineSimSum100MS1,
                CandidateScores::Features::MS1Averagine,
                CandidateScores::Features::CosineSimSum100Window1p5X,
                CandidateScores::Features::CosineSimSum100Window2X
            };
            const QVector<float> nnVec = candidateScores->selectFeaturesArrayFeatures(nnFeatures);
            return nnVec;
#endif

        }
        else if (useExtendedScores) {
            QVector<float> vec = candidateScores->selectFeaturesArrayFeatures(baseFeatures);
            vec.append(candidateScores->selectFeaturesArrayFeatures({
                            CandidateScores::Features::CosineSimSum100,
                            CandidateScores::Features::ScanTimeDeltaAbs,
                            CandidateScores::Features::ScanTimePdAbs,
                            CandidateScores::Features::PeakShapeRatio1,
                            CandidateScores::Features::PeakShapeRatio2,
                            CandidateScores::Features::PeakShapeRatio3,
                            CandidateScores::Features::ShadowsCosineSimSum,
                            CandidateScores::Features::CosineSimToAnchor1,
                            CandidateScores::Features::CosineSimToAnchor2,
                            CandidateScores::Features::CosineSimToAnchor3,
                            CandidateScores::Features::CosineSimToAnchor4,
                            CandidateScores::Features::CosineSimToAnchor5,
                            CandidateScores::Features::CosineSimToAnchor6,
                            CandidateScores::Features::CosineSimToAnchor7,
                            CandidateScores::Features::CosineSimToAnchor8,
                            // CandidateScores::Features::CosineSimToAnchor9,
                            // CandidateScores::Features::CosineSimToAnchor10,
                            // CandidateScores::Features::CosineSimToAnchor11,
                            // CandidateScores::Features::CosineSimToAnchor12,
                            CandidateScores::Features::CosineSimShadowsToAnchor1,
                            CandidateScores::Features::CosineSimShadowsToAnchor2,
                            CandidateScores::Features::CosineSimShadowsToAnchor3,
                            CandidateScores::Features::CosineSimShadowsToAnchor4,
                            CandidateScores::Features::CosineSimShadowsToAnchor5,
                            CandidateScores::Features::CosineSimShadowsToAnchor6,
                            CandidateScores::Features::CosineSimShadowsToAnchor7,
                            CandidateScores::Features::CosineSimShadowsToAnchor8,
                            CandidateScores::Features::CosineSimShadowsToAnchor9,
                            CandidateScores::Features::CosineSimShadowsToAnchor10,
                            CandidateScores::Features::CosineSimShadowsToAnchor11,
                            CandidateScores::Features::CosineSimShadowsToAnchor12,
                            CandidateScores::Features::MzPeakLengthsNorm1,
                            CandidateScores::Features::MzPeakLengthsNorm2,
                            CandidateScores::Features::MzPeakLengthsNorm3,
                            CandidateScores::Features::MzPeakLengthsNorm4,
                            CandidateScores::Features::MzPeakLengthsNorm5,
                            CandidateScores::Features::MzPeakLengthsNorm6,
                            CandidateScores::Features::MzPeakLengthsNorm7,
                            CandidateScores::Features::MzPeakLengthsNorm8,
                            CandidateScores::Features::MzPeakLengthsNorm9,
                            CandidateScores::Features::MzPeakLengthsNorm10,
                            CandidateScores::Features::MzPeakLengthsNorm11,
                            CandidateScores::Features::MzPeakLengthsNorm12,
                            CandidateScores::Features::CosineSimSpectrumOverTime,
                            CandidateScores::Features::TotalIntensityLog,
                            CandidateScores::CosineSimSum100Window1p5X,
                            CandidateScores::CosineSimSum100Window2X,
                
                            CandidateScores::TargetWindowLocationAbs
                            }));
            return vec;
        }

        const QVector<float> vec = candidateScores->selectFeaturesArrayFeatures(baseFeatures);

        return vec;
    }

Err DiscriminantScoretron::convertScoreCandidatesToFeaturesArrays(
    const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>>& candidateScoresTargetVsDecoy,
    bool useExtendedScores,
    bool useNeuralNetworkScores,
    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *featuresArrayTargetVsDecoy
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScoresTargetVsDecoy); ree;

    featuresArrayTargetVsDecoy->clear();
    featuresArrayTargetVsDecoy->reserve(candidateScoresTargetVsDecoy.size());

    for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : candidateScoresTargetVsDecoy) {
        const FeaturesArray featuresArrayTarget = scoreVectorLogic(
            useExtendedScores,
            useNeuralNetworkScores,
            pr.first
            );

        const FeaturesArray featuresArrayDecoy = scoreVectorLogic(
            useExtendedScores,
            useNeuralNetworkScores,
            pr.second
            );

        featuresArrayTargetVsDecoy->push_back({featuresArrayTarget, featuresArrayDecoy});
    }

    ERR_RETURN
}
