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

    e = ClassifierWeightsManager::fitWeights(
        A,
        b,
        weights
        ); ree;

    if (verbosity > 0) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "fit weights" << et.restart() << "mSec";
        qDebug() << "Weights:" << *weights;
        qDebug() << "b:" << b;
        for (const QVector<float> &v : A) {
            qDebug() << v;
        }
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

    if (threadCount == 1) {
        const QPair<Err, QVector<float>> result = applyWeightsLogic(weights, featuresArray);
        e = result.first; ree;
        discriminantScores->append(result.second);
        ERR_RETURN
    }

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

        const QVector<Features> baseFeatures = {

            Features::CosineSimSum100GreaterThan80,
            Features::CosineSimSpectrumOverTimeCubed,
            Features::CosineSimSpectrumStDev,
            Features::CosineSim100MS1,
            Features::CosineSim100MS1Iso1,
            Features::CosineSim100MS1Iso2,
            Features::CosineSim100MS1PreMono,
            Features::CosineSimSpectrumCubed,
            Features::CosineSimSum45,
            Features::CosineSimSumTop,
            Features::CosineSimSumBottom,
            Features::TopBottomRatio,
            Features::TopBottomRatioNorm,
            Features::PeakShapeRatio1,
            Features::PeakShapeRatio2,
            Features::PeakShapeRatio3,
            Features::MatrixZeroPercentage,
            Features::MzPPMMeanAbs,
            Features::FoundB,
            Features::FoundY,
            Features::FoundPercent,
            DiscScoresCount,
            DiscScoresMean,
            DiscScoresStDev,
            DiscScore1stRunnerUp,
            DiscScore2ndRunnerUp
        };

        if (useNeuralNetworkScores) {

            const QVector<Features> nnFeatures = {
                Features::CosineSimSum100,
                Features::CosineSimSum100Top12,
                Features::CosineSimSum100GreaterThan80,
                Features::AllignedMaxIndexesCount,
                Features::CosineSim100MS1,
                Features::CosineSimSpectrumCubed,
                Features::KlDivSpectrumCubeRoot,
                Features::CosineSimSum45,
                Features::CosineSimSumTop,
                Features::CosineSimSumBottom,
                Features::TopBottomRatio,
                Features::TopBottomRatioNorm,
                Features::Charge,
                Features::ScanTimeDelta,
                Features::ScanTimeDeltaAbs,
                Features::ScanTimePd,
                Features::ScanTimePdAbs,
                Features::ScanIonCount,
                Features::MzNorm,
                Features::Mass,
                Features::KlDivSpectrum,
                Features::CosineSimSpectrum,
                Features::CosineSim45MS1,
                Features::CosineSim100MS1PreMono,
                Features::CosineSim100MS1Iso1,
                Features::CosineSim100MS1Iso2,
                Features::PeptideLengthNorm,
                Features::ScanTimePredicted,
                Features::TheoFragmentCount,
                Features::TotalIntensityLog,
                Features::PeakShapeRatio1,
                Features::PeakShapeRatio2,
                Features::PeakShapeRatio3,
                Features::ShadowsCosineSimSum,
                Features::IRTPredicted,
                Features::IntensityFoundMaxNorm1,
                Features::IntensityFoundMaxNorm2,
                Features::IntensityFoundMaxNorm3,
                Features::IntensityFoundMaxNorm4,
                Features::IntensityFoundMaxNorm5,
                Features::IntensityFoundMaxNorm6,
                Features::IntensityFoundMaxNorm7,
                Features::IntensityFoundMaxNorm8,
                Features::IntensityFoundMaxNorm9,
                Features::IntensityFoundMaxNorm10,
                Features::IntensityFoundMaxNorm11,
                Features::IntensityFoundMaxNorm12,
                Features::AminoAcidCountA,
                Features::AminoAcidCountC,
                Features::AminoAcidCountD,
                Features::AminoAcidCountE,
                Features::AminoAcidCountF,
                Features::AminoAcidCountG,
                Features::AminoAcidCountH,
                Features::AminoAcidCountI,
                Features::AminoAcidCountK,
                Features::AminoAcidCountL,
                Features::AminoAcidCountM,
                Features::AminoAcidCountN,
                Features::AminoAcidCountP,
                Features::AminoAcidCountQ,
                Features::AminoAcidCountR,
                Features::AminoAcidCountS,
                Features::AminoAcidCountT,
                Features::AminoAcidCountV,
                Features::AminoAcidCountW,
                Features::AminoAcidCountY,
                Features::AminoAcidCountB,
                Features::AminoAcidCountJ,
                Features::AminoAcidCountO,
                Features::AminoAcidCountU,
                Features::AminoAcidCountX,
                Features::AminoAcidCountZ,
                Features::AltTargetKeyIdDiscScoreChargeOG_alt,
                Features::AltTargetKeyIdDiscScoreCharge1_1,
                Features::AltTargetKeyIdDiscScoreCharge1_2,
                Features::AltTargetKeyIdDiscScoreCharge2_1,
                Features::AltTargetKeyIdDiscScoreCharge2_2,
                Features::AltTargetKeyIdDiscScoreCharge3_1,
                Features::AltTargetKeyIdDiscScoreCharge3_2,
                Features::AltTargetKeyIdDiscScoreCharge4_1,
                Features::AltTargetKeyIdDiscScoreCharge4_2,
                Features::AltTargetKeyIdTimeDeltaCharge1_1,
                Features::AltTargetKeyIdTimeDeltaCharge2_1,
                Features::AltTargetKeyIdTimeDeltaCharge3_1,
                Features::AltTargetKeyIdTimeDeltaCharge4_1,
                Features::Ms1MzMeanFound100,
                Features::Ms1MzMeanFound45,
                Features::Ms1MzMeanFoundPreMono,
                Features::Ms1MzMeanFoundIso1,
                Features::Ms1MzMeanFoundIso2,
                Features::Ms1MzMeanFound100PPM,
                Features::Ms1MzMeanFound45PPM,
                Features::Ms1MzMeanFoundPreMonoPPM,
                Features::Ms1MzMeanFoundIso1PPM,
                Features::Ms1MzMeanFoundIso2PPM,
                Features::Ms1MzStDevFound100,
                Features::Ms1MzStDevFound45,
                Features::Ms1MzStDevFoundPreMono,
                Features::Ms1MzStDevFoundIso1,
                Features::Ms1MzStDevFoundIso2,
                Features::Ms1IntensityFound100,
                Features::Ms1IntensityFound45,
                Features::Ms1IntensityFoundPreMono,
                Features::Ms1IntensityFoundIso1,
                Features::Ms1IntensityFoundIso2,
                Features::CosineSimSpectrumOverTime,
                Features::CosineSimSpectrumOverTimeCubed,
                Features::CosineSimSpectrumStDev,
                Features::CosineSimSum100MS1,
                Features::MS1Averagine,
                Features::CosineSimSum100Window1p5X,
                Features::CosineSimSum100Window2X,
                Features::TotalIntensityRaw,
                Features::TargetWindowLocation,
                Features::DiscriminantScore,
                Features::MatrixZeroPercentage,

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
                BIonSeriesMax,
                BIonSeriesCount,
                BIonSeriesMean,
                BIonSeriesStd,
                BIonSeriesTheoMax,
                BIonSeriesTheoCount,
                BIonSeriesTheoMean,
                BIonSeriesTheoStd,
                // BIonSeriesMaxFoundToTheoFraction,
                // YIonSeriesCountRatio,
                // BIonSeriesCountRatio,

                // CosineSimFullTheo,
                IonsFoundFractionFull,
                CosineSimFullTheoXIonsFoundFractionFull,
                MzPPMMean,
                FoundB,
                FoundY,
                FoundPercent,
                MzPPMStd

            };
            const QVector<float> nnVec = candidateScores->selectFeaturesArrayFeatures(nnFeatures);
            return nnVec;

        }
        else if (useExtendedScores) {
            const QVector<Features> vec = {
                            Features::CosineSimSum100GreaterThan80,
                            Features::CosineSimSpectrumOverTimeCubed,
                            Features::CosineSimSpectrumStDev,
                            Features::CosineSim100MS1,
                            Features::CosineSim100MS1Iso1, //5
                            Features::CosineSim100MS1Iso2,
                            Features::CosineSim100MS1PreMono,
                            Features::CosineSimSpectrumCubed,
                            Features::CosineSimSum45,
                            Features::CosineSimSumTop, //10
                            Features::CosineSimSumBottom,
                            Features::TopBottomRatio,
                            Features::TopBottomRatioNorm,
                            Features::PeakShapeRatio1,
                            Features::PeakShapeRatio2,//15
                            Features::PeakShapeRatio3,
                            Features::MatrixZeroPercentage,
                            Features::MzPPMMeanAbs,
                            Features::FoundY,
                            Features::FoundPercent, //20
                            Features::CosineSimSum100Top12,
                            Features::ScanTimeDeltaAbs,
                            Features::ScanTimePdAbs,
                            Features::ShadowsCosineSimSum,
                            Features::TotalIntensityLog, //25
                            Features::CosineSimSum100Window1p5X,
                            Features::CosineSimSum100Window2X,
                            Features::TargetWindowLocationAbs
                            };

            const QVector<float> vecScores = candidateScores->selectFeaturesArrayFeatures(vec);
            return vecScores;
        }

        const QVector<float> vec = candidateScores->selectFeaturesArrayFeatures(baseFeatures);
        return vec;
    }

QVector<float> DiscriminantScoretron::defaultWeights(
    bool useExtendedScores,
    bool useNeuralNetworkScores
    ) {

    CandidateScores cs;
    cs.initFeaturesArray();

    cs.featuresArray[Features::CosineSimSum100GreaterThan80] = 1.0f;
    cs.featuresArray[Features::CosineSimSum100Top12] = 1.0f;
    cs.featuresArray[Features::CosineSimSpectrumOverTimeCubed] = 1.0f;
    // cs.featuresArray[Features::KlDivSpectrumCubeRoot] = -1.0f;
    // cs.featuresArray[Features::CosineSim100MS1] = 1.0f;
    // cs.featuresArray[Features::CosineSimSpectrumStDev] = -1.0f;
    cs.featuresArray[Features::ScanTimeDeltaAbs] = -0.5f;
    // cs.featuresArray[Features::ShadowsCosineSimSum] = -1.0f;

    return scoreVectorLogic(useExtendedScores, useNeuralNetworkScores, &cs);
}

namespace {

    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> extractFeaturesLogic(
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &pairs,
        bool useExtendedScores,
        bool useNeuralNetworkScores
        ) {

        QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> output;
        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : pairs) {
            const FeaturesArray featuresArrayTarget = DiscriminantScoretron::scoreVectorLogic(
                useExtendedScores,
                useNeuralNetworkScores,
                pr.first
                );

            const FeaturesArray featuresArrayDecoy = DiscriminantScoretron::scoreVectorLogic(
                useExtendedScores,
                useNeuralNetworkScores,
                pr.second
                );

            output.push_back({featuresArrayTarget, featuresArrayDecoy});
        }

        return output;
    }


}//namespace
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

    const auto trainingLogicBinder = std::bind(
        extractFeaturesLogic,
        std::placeholders::_1,
        useExtendedScores,
        useNeuralNetworkScores
    );

    QVector<QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>>> candidateScoresTargetVsDecoyTranched;
    e = ParallelUtils::trancheVectorForParallelizationInOrder(
        candidateScoresTargetVsDecoy,
        ParallelUtils::numberOfAvailableSystemProcessors(),
        0,
        &candidateScoresTargetVsDecoyTranched
        ); ree;

    QFuture<QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>>> futures = QtConcurrent::mapped(
        candidateScoresTargetVsDecoyTranched,
        trainingLogicBinder
        );
    futures.waitForFinished();

    for (const QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> &res : futures) {
        featuresArrayTargetVsDecoy->append(res);
    }

    ERR_RETURN
}
