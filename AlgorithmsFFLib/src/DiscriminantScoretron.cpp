//
// Created by anichols on 1/9/24.
//

#include "DiscriminantScoretron.h"

#include "ClassifierWeightsManager.h"
#include "ParallelUtils.h"

#include <QtConcurrent/QtConcurrent>

#include "EigenUtils.h"


QVector<Features> DiscriminantScoretron::featuresCalibration() {
    const QVector<Features> baseFeatures = {
            CosineSimSpectrumOverTimeCubed, // 1 1
            YIonSeriesMax, // 0.963 0.95
            BIonSeriesCountRatio, // 0.929 1
            DiscScore1stRunnerUp, // 0.906 0.6
            // TargetWindowLocationAbs, // 0.903 0.55
            CosineSimSpectrum, // 0.899 0.9
            // YIonSeriesTheoCount, // 0.891 0.6
            // ScanIonCount, // 0.873 0.8
            CosineSimSum100GreaterThan80, // 0.846 0.75
            CosineSim45MS1, // 0.843 0.65
            // MzPPMMeanAbs, // 0.839 0.6
            CosineSimSum45, // 0.824 0.7
            // ScanTimePredicted, // 0.824 0.65
            CosineSim100MS1Iso2, // 0.82 0.85
            // CosineSimSpectrumStDev, // 0.801 0.55
            // IRTPredicted, // 0.798 0.25
            FoundB, // 0.798 0.3
            // KlDivSpectrumCubeRoot, // 0.79 0.15
            TopBottomRatio, // 0.779 0.6
            // YIonSeriesCount, // 0.772 0.3
            // CosineSimFullTheo, // 0.772 0.6
            // ScanTimePdAbs, // 0.749 0.5
            CosineSimSpectrumOverTime, // 0.749 0.9
            BIonSeriesTheoStd, // 0.749 0.95
            MS1Averagine, // 0.479 0.75
            // DiscScore2ndRunnerUp, // 0.745 0.35
            // KlDivSpectrum, // 0.738 0.55
            // YIonSeriesTheoStd, // 0.689 0.95
            // CosineSimSum100Top12, // 0.685 0.25
            // AllignedMaxIndexesCount, // 0.682 0.45
            // IonsFoundFractionFull, // 0.682 0.5
            // BIonSeriesTheoCount, // 0.674 0.6
            // FoundY, // 0.67 0.7
            // MatrixZeroPercentage, // 0.659 0.35
            // YIonSeriesTheoMax, // 0.659 0.45
            // YIonSeriesCountRatio, // 0.655 0.25
            // CosineSim100MS1, // 0.637 0.6
            // BIonSeriesTheoMean, // 0.637 0.85
            // PeptideLengthNorm, // 0.633 0.65
            // ScanTimeDeltaAbs, // 0.625 0.85
            // DiscScoresStDev, // 0.625 0.75
            // MonoPreMonoRatio, // 0.603 0.45
            // CosineSimSumTop, // 0.599 0.65
            // CosineSim100MS1PreMono, // 0.599 0.45
            // YIonSeriesTheoMean, // 0.592 0.15
            // BIonSeriesMax, // 0.547 0.5
            // YIonSeriesMaxFoundToTheoFraction, // 0.487 0.5
            // CosineSimSum100, // 0.468 0.35
            // TopBottomRatioNorm, // 0.438 0.15
            // TotalIntensityLog, // 0.404 0.65
            // YIonSeriesMean, // 0.401 0.35
            // BIonSeriesTheoMax, // 0.39 0.35
            // MzNorm, // 0.382 0.35
            // CosineSim100MS1Iso1, // 0.371 0.7
            // CosineSimFullTheoXIonsFoundFractionFull, // 0.311 0.5
            // DiscScoresMean, // 0.285 0.5
            // BIonSeriesMean, // 0.262 0.55
            // FoundPercent, // 0.258 0.5
            // BIonSeriesStd, // 0.247 0.25
            // PeakShapeRatio2, // 0.232 0.45
            // CosineSimSumBottom, // 0.206 0.4
            // CosineSimSum100MS1, // 0.206 0.5
            // DiscScoresCount, // 0.206 0.25
            // YIonSeriesStd, // 0.184 0.3
            // CosineSimSum100Window2X, // 0.176 0.5
            // PeakShapeRatio3, // 0.165 0.45
            // TheoFragmentCount, // 0.146 0.1
            // MzPPMStd, // 0.112 0.3
            // CosineSimSpectrumCubed, // 0.079 0.15
            // BIonSeriesMaxFoundToTheoFraction, // 0.079 0.15
            // CosineSimSum100Window1p5X, // 0.067 0.35
            // BIonSeriesCount, // 0.052 0.1
            // ShadowsCosineSimSum, // 0.049 0.05
            // PeakShapeRatio1, // 0.019 0

            // CosineSimSpectrumOverTimeCubed,// 1
            // YIonSeriesMax,// 0.963
            // BIonSeriesCountRatio,// 0.929
            // DiscScore1stRunnerUp,// 0.906
            // TargetWindowLocationAbs,// 0.903
            // CosineSimSpectrum,// 0.899
            // YIonSeriesTheoCount,// 0.891
            // // ScanIonCount,// 0.873
            // CosineSimSum100GreaterThan80,// 0.846
            // CosineSim45MS1,// 0.843
            // // MzPPMMeanAbs,// 0.839
            // CosineSimSum45,// 0.824
            // // ScanTimePredicted,// 0.824
            // CosineSim100MS1Iso2,// 0.82
            // CosineSimSpectrumStDev,// 0.801
            // // IRTPredicted,// 0.798
            // FoundB,// 0.798
            // TopBottomRatio,// 0.779
            // YIonSeriesCount,// 0.772
            // CosineSimFullTheo,// 0.772
            // // ScanTimePdAbs,// 0.749
            // // CosineSimSpectrumOverTime,// 0.749
            // // BIonSeriesTheoStd,// 0.749
            // // DiscScore2ndRunnerUp,// 0.745
            // // KlDivSpectrum,// 0.738
            // // YIonSeriesTheoStd,// 0.689
            // CosineSimSum100Top12,// 0.685
            // // AllignedMaxIndexesCount,// 0.682
            // // IonsFoundFractionFull,// 0.682
            // // BIonSeriesTheoCount,// 0.674
            // // FoundY,// 0.67
            // // MatrixZeroPercentage,// 0.659
            // // YIonSeriesTheoMax,// 0.659
            // // YIonSeriesCountRatio,// 0.655
            // // CosineSim100MS1,// 0.637
            // // BIonSeriesTheoMean,// 0.637
            // // PeptideLengthNorm,// 0.633
            // // ScanTimeDeltaAbs,// 0.625
            // // DiscScoresStDev,// 0.625
            // // MonoPreMonoRatio,// 0.603
            // // CosineSimSumTop,// 0.599
            // // CosineSim100MS1PreMono,// 0.599
            // // YIonSeriesTheoMean,// 0.592
            // // BIonSeriesMax,// 0.547
            // // YIonSeriesMaxFoundToTheoFraction,// 0.487
            // // MS1Averagine,// 0.479
            // // CosineSimSum100,// 0.468
            // // TopBottomRatioNorm,// 0.438
            // // TotalIntensityLog,// 0.404
            // // YIonSeriesMean,// 0.401
            // // BIonSeriesTheoMax,// 0.39
            // // MzNorm,// 0.382
            // // CosineSim100MS1Iso1,// 0.371
            // // CosineSimFullTheoXIonsFoundFractionFull,// 0.311
            // // DiscScoresMean,// 0.285
            // // BIonSeriesMean,// 0.262
            // // FoundPercent,// 0.258
            // // BIonSeriesStd,// 0.247
            // // PeakShapeRatio2,// 0.232
            // // CosineSimSumBottom,// 0.206
            // // CosineSimSum100MS1,// 0.206
            // // DiscScoresCount,// 0.206
            // // YIonSeriesStd,// 0.184
            // // CosineSimSum100Window2X,// 0.176
            // // PeakShapeRatio3,// 0.165
            // // TheoFragmentCount,// 0.146
            // // MzPPMStd,// 0.112
            // // CosineSimSpectrumCubed,// 0.079
            // // BIonSeriesMaxFoundToTheoFraction,// 0.079
            // // CosineSimSum100Window1p5X,// 0.067
            // // BIonSeriesCount,// 0.052
            // // ShadowsCosineSimSum,// 0.049
            // // PeakShapeRatio1,// 0.019



        // CosineSimSum100GreaterThan80,
        // CosineSimSpectrumOverTimeCubed,
        // CosineSimSpectrumStDev,
        // CosineSim100MS1,
        // Ms1MzStDevFound100,
        // CosineSim100MS1Iso1,
        // CosineSim100MS1Iso2,
        // CosineSimSpectrumCubed,
        // CosineSimSum45,
        // CosineSimSumTop,
        // CosineSimSumBottom,
        // TopBottomRatio,
        // TopBottomRatioNorm,
        // PeakShapeRatio1,
        // PeakShapeRatio2,
        // PeakShapeRatio3,
        // MatrixZeroPercentage,
        // MzPPMMeanAbs,
        // FoundB,
        // FoundY,
        // FoundPercent,
        // DiscScoresCount,
        // DiscScoresMean,
        // DiscScoresStDev,
        // DiscScore1stRunnerUp,
        // DiscScore2ndRunnerUp
    };

    return baseFeatures;
}

QVector<Features> DiscriminantScoretron::featuresOptimization() {

    const QVector<Features> vec = {
        CosineSimSum100GreaterThan80,
        CosineSimSpectrumOverTimeCubed,
        CosineSimSpectrumStDev,
        CosineSimSum100MS1,
        CosineSim100MS1,
        CosineSim100MS1Iso1,
        CosineSim100MS1Iso2,
        MonoPreMonoRatio,
        CosineSimSpectrumCubed,
        CosineSimSum45,
        CosineSimSumTop,
        CosineSimSumBottom,
        TopBottomRatio,
        TopBottomRatioNorm,
        PeakShapeRatio1,
        PeakShapeRatio2,
        PeakShapeRatio3,
        MatrixZeroPercentage,
        MzPPMMeanAbs,
        FoundY,
        FoundPercent,
        CosineSimSum100Top12,
        ScanTimeDeltaAbs,
        ScanTimePdAbs,
        ShadowsCosineSimSum,
        TotalIntensityLog,
        CosineSimSum100Window1p5X,
        CosineSimSum100Window2X,
        TargetWindowLocationAbs
        };

    return vec;
}

QVector<Features> DiscriminantScoretron::featuresNeuralNetwork() {
    const QVector<Features> nnFeatures = {
                CosineSimSum100,
                CosineSimSum100Top12,
                CosineSimSum100GreaterThan80,
                AllignedMaxIndexesCount,
                CosineSim100MS1,
                CosineSimSpectrumCubed,
                KlDivSpectrumCubeRoot,
                CosineSimSum45,
                CosineSimSumTop,
                CosineSimSumBottom,
                TopBottomRatio,
                TopBottomRatioNorm,
                Charge,
                ScanTimeDelta,
                ScanTimeDeltaAbs,
                ScanTimePd,
                ScanTimePdAbs,
                ScanIonCount,
                MzNorm,
                Mass,
                KlDivSpectrum,
                CosineSimSpectrum,
                CosineSim45MS1,
                CosineSim100MS1PreMono,
                CosineSim100MS1Iso1,
                CosineSim100MS1Iso2,
                PeptideLengthNorm,
                ScanTimePredicted,
                TheoFragmentCount,
                TotalIntensityLog,
                PeakShapeRatio1,
                PeakShapeRatio2,
                PeakShapeRatio3,
                ShadowsCosineSimSum,
                IRTPredicted,
                IntensityFoundMaxNorm1,
                IntensityFoundMaxNorm2,
                IntensityFoundMaxNorm3,
                IntensityFoundMaxNorm4,
                IntensityFoundMaxNorm5,
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
                AminoAcidCountB,
                AminoAcidCountJ,
                AminoAcidCountO,
                AminoAcidCountU,
                AminoAcidCountX,
                AminoAcidCountZ,
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
                Ms1MzMeanFound100,
                Ms1MzMeanFoundPreMono,
                Ms1MzMeanFound100PPM,
                Ms1MzMeanFoundPreMonoPPM,
                Ms1MzStDevFound100,
                Ms1MzStDevFoundPreMono,
                Ms1IntensityFound100,
                Ms1IntensityFound45,
                Ms1IntensityFoundPreMono,
                Ms1IntensityFoundIso1,
                Ms1IntensityFoundIso2,
                CosineSimSpectrumOverTime,
                CosineSimSpectrumOverTimeCubed,
                CosineSimSpectrumStDev,
                CosineSimSum100MS1,
                MS1Averagine,
                CosineSimSum100Window1p5X,
                CosineSimSum100Window2X,
                TotalIntensityRaw,
                TargetWindowLocation,
                DiscriminantScore,
                MatrixZeroPercentage,
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
                MonoPreMonoRatio,
                IonsFoundFractionFull,
                CosineSimFullTheoXIonsFoundFractionFull,
                MzPPMMean,
                FoundB,
                FoundY,
                FoundPercent,
                MzPPMStd
            };

    return nnFeatures;
}

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

QVector<float> DiscriminantScoretron::defaultWeights(const QVector<Features> &features) {

    CandidateScores cs;
    cs.initFeaturesArray();

    cs.featuresArray[CosineSimSum100GreaterThan80] = 1.0f;
    cs.featuresArray[CosineSimSum100Top12] = 1.0f;
    cs.featuresArray[CosineSimSpectrumOverTimeCubed] = 1.0f;
    // cs.featuresArray[KlDivSpectrumCubeRoot] = -1.0f;
    // cs.featuresArray[CosineSim100MS1] = 1.0f;
    // cs.featuresArray[CosineSimSpectrumStDev] = -1.0f;
    cs.featuresArray[ScanTimeDeltaAbs] = -0.5f;
    // cs.featuresArray[ShadowsCosineSimSum] = -1.0f;

    return CandidateScores::selectFeaturesArrayFeatures(cs.featuresArray, features);
}

namespace {

    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> extractFeaturesLogic(
        const QVector<Features> &features,
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &pairs
        ) {

        QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> output;
        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : pairs) {

            const FeaturesArray featuresArrayTarget = CandidateScores::selectFeaturesArrayFeatures(
                pr.first->featuresArray,
                features
                );

            const FeaturesArray featuresArrayDecoy = CandidateScores::selectFeaturesArrayFeatures(
               pr.second->featuresArray,
               features
               );

            output.push_back({featuresArrayTarget, featuresArrayDecoy});
        }

        return output;
    }


}//namespace
Err DiscriminantScoretron::convertScoreCandidatesToFeaturesArrays(
    const QVector<Features> &features,
    const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>>& candidateScoresTargetVsDecoy,
    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *featuresArrayTargetVsDecoy
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScoresTargetVsDecoy); ree;

    featuresArrayTargetVsDecoy->clear();
    featuresArrayTargetVsDecoy->reserve(candidateScoresTargetVsDecoy.size());

    const auto trainingLogicBinder = std::bind(
        extractFeaturesLogic,
        features,
        std::placeholders::_1
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
