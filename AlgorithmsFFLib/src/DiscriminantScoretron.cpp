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
        CosineSimSum100GreaterThan80,
        CosineSimSpectrumOverTimeCubed,
        CosineSimSpectrumStDev,
        CosineSim100MS1,
        Ms1MzStDevFound100,
        CosineSim100MS1Iso1,
        CosineSim100MS1Iso2,
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
        FoundB,
        FoundY,
        FoundPercent,
        DiscScoresCount,
        DiscScoresMean,
        DiscScoresStDev,
        DiscScore1stRunnerUpDiff,
        DiscScore2ndRunnerUpDiff,
    	MzFoundOverCount650,
		MzFoundUnderCount650,
    	MzPeakLengthsMean,
    	MzPeakLengthsStd
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
    	MzPPMStd,
        FoundY,
        FoundPercent,
        CosineSimSum100Top12,
        ScanTimeDeltaAbs,
        ScanTimePdAbs,
        ShadowsCosineSimSum,
        TotalIntensityLog,
        CosineSimSum100Window1p5X,
        CosineSimSum100Window2X,
        TargetWindowLocationAbs,
    	MzPeakLengthsStd
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
                // Ms1MzMeanFound100,
                // Ms1MzMeanFoundPreMono,
                // Ms1MzMeanFound100PPM,
                // Ms1MzMeanFoundPreMonoPPM,
                // Ms1MzStDevFound100,
                // Ms1MzStDevFoundPreMono,
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
                DiscScore1stRunnerUpDiff,
                DiscScore2ndRunnerUpDiff,
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

                MzPeakLengthsMean,
                MzPeakLengthsStd,
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

                MzPPMStd,

    			MzFoundOverCount650,
				MzFoundUnderCount650,
            };

    return nnFeatures;
}

Err DiscriminantScoretron::trainLDAClassifier(
        const QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> &_targetDecoyCandidateScoresPair,
        int verbosity,
        QVector<float> *weights
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(_targetDecoyCandidateScoresPair); ree;

    QElapsedTimer et;
    et.start();

    weights->clear();

    QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> targetDecoyCandidateScoresPair = _targetDecoyCandidateScoresPair;
    // targetDecoyCandidateScoresPair.resize(
    //     std::max(targetDecoyCandidateScoresPair.size() / 2, std::min(50, _targetDecoyCandidateScoresPair.size()))
    //     ); ree;

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
