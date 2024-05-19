//
// Created by anichols on 1/9/24.
//

#include "DiscriminantScoretron.h"

#include "CandidateScores.h"
#include "ClassifierWeightsManager.h"
#include "ParallelUtils.h"

#include <QtConcurrent/QtConcurrent>

namespace {

    struct BuildClassiferParallelInput {
        QVector<ScoresTargets*> scoresTargets;
        QVector<ScoresDecoys*> scoresDecoys;
    };

    struct BuildClassifierParallelOutput {
        QVector<QVector<float>> A;
        QVector<float> b;
    };

    Err buildParallelInput(
            const QVector<ScoresTargets*> &scoresTargets,
            const QVector<ScoresDecoys*> &scoresDecoys,
            QVector<BuildClassiferParallelInput> *inputs
    ) {
        ERR_INIT

        e = ErrorUtils::isEqual(scoresTargets.size(), scoresDecoys.size()); ree;

        QVector<QVector<ScoresTargets*>> scoresTargetsTranched;
        QVector<QVector<ScoresDecoys*>> scoresDecoysTranched;

        const int tranches = std::min(ParallelUtils::numberOfAvailableSystemProcessors(), scoresTargets.size());

        e = ParallelUtils::trancheVectorForParallelization(
                scoresTargets,
                tranches,
                &scoresTargetsTranched
        ); ree;

        e = ParallelUtils::trancheVectorForParallelization(
                scoresDecoys,
                tranches,
                &scoresDecoysTranched
        ); ree;

        e = ErrorUtils::isEqual(scoresTargetsTranched.size(), scoresDecoysTranched.size()); ree;
        for (int i = 0; i < scoresDecoysTranched.size(); i++) {
            BuildClassiferParallelInput input;
            input.scoresTargets = scoresTargetsTranched.at(i);
            input.scoresDecoys = scoresDecoysTranched.at(i);
            inputs->push_back(input);
        }

        ERR_RETURN
    }

    QPair<Err, BuildClassifierParallelOutput> buildClassifierDataParallel2(
            const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &pairs,
            bool useExtendedScores,
            bool useNeuralNetworkScores
            ) {

        ERR_INIT

        BuildClassifierParallelOutput output;

        QVector<QVector<float>> scoresTargets;
        scoresTargets.reserve(pairs.size());
        QVector<QVector<float>> scoresDecoys;
        scoresDecoys.reserve(pairs.size());

        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : pairs) {
            const QVector<float> candScoresTarget = DiscriminantScoretron::scoreVectorLogic(
                    useExtendedScores,
                    useNeuralNetworkScores,
                    pr.first
                    );

            const QVector<float> candScoresDecoy = DiscriminantScoretron::scoreVectorLogic(
                    useExtendedScores,
                    useNeuralNetworkScores,
                    pr.second
            );

            scoresTargets.push_back(candScoresTarget);
            scoresDecoys.push_back(candScoresDecoy);
        }

        QVector<QVector<float>*> scoresTargetsPntrs;
        scoresTargetsPntrs.reserve(pairs.size());
        std::transform(
                scoresTargets.begin(),
                scoresTargets.end(),
                std::back_inserter(scoresTargetsPntrs),
                [](QVector<float> &v){return &v;}
                );

        QVector<QVector<float>*> scoresDecoysPntrs;
        scoresDecoysPntrs.reserve(pairs.size());
        std::transform(
                scoresDecoys.begin(),
                scoresDecoys.end(),
                std::back_inserter(scoresDecoysPntrs),
                [](QVector<float> &v){return &v;}
        );

        e = ClassifierWeightsManager::buildDataClassifier1(
                scoresTargetsPntrs,
                scoresDecoysPntrs,
                &output.A,
                &output.b
        ); rree;

        return {e, output};
    }

    Err processParallelResults(
            const QPair<Err, BuildClassifierParallelOutput> &result,
            int inputsSize,
            QVector<QVector<float>> *A,
            QVector<float> *b
    ) {

        ERR_INIT
        e = result.first; ree;

        const QVector<QVector<float>> &ALocal = result.second.A;
        const QVector<float> bLocal = result.second.b;

        if (A->isEmpty() || b->isEmpty()) {
            A->resize(ALocal.size());
            for(int row = 0; row < A->size(); row++) {
                QVector<float> v(ALocal.front().size(), 0.0);
                (*A)[row] = v;
            }
            b->resize(bLocal.size());
        }

        for (int row = 0; row < A->size(); row++) {
            for (int col = 0; col < A->front().size(); col++) {
                (*A)[row][col] += ALocal[row][col] / static_cast<float>(inputsSize);
            }
        }

        for (int row = 0; row < bLocal.size(); row++) {
            (*b)[row] += bLocal[row] / static_cast<float>(inputsSize);
        }

        ERR_RETURN
    }

}//namespace


Err DiscriminantScoretron::trainLDAClassifier(
        const QList<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &targetDecoyCandidateScoresPair,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        QVector<float> *weights
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targetDecoyCandidateScoresPair); ree;

    QElapsedTimer et;
    et.start();

    weights->clear();

    QVector<QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>>> targetDecoyCandidateScoresPairTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            targetDecoyCandidateScoresPair.toVector(),
            ParallelUtils::numberOfAvailableSystemProcessors(),
            &targetDecoyCandidateScoresPairTranched
            ); ree;

    const auto classifierBuilder = std::bind(
            buildClassifierDataParallel2,
            std::placeholders::_1,
            useExtendedScores,
            useNeuralNetworkScores
    );

    QFuture<QPair<Err, BuildClassifierParallelOutput>> futures = QtConcurrent::mapped(
            targetDecoyCandidateScoresPairTranched,
            classifierBuilder
    );
    futures.waitForFinished();

    QVector<QVector<float>> A;
    QVector<float> b;
    for (const QPair<Err, BuildClassifierParallelOutput> &result : futures) {
        e = processParallelResults(
                result,
                targetDecoyCandidateScoresPairTranched.size(),
                &A,
                &b
        ); ree;
    }

    e = ClassifierWeightsManager::fitWeights(A, b, weights); ree;

    qDebug() << "fit weights" << et.restart() << "mSec";
    qDebug() << "Weights:" << *weights;
    qDebug() << "b:" << b;

    ERR_RETURN
}

namespace {

    Err applyWeightsLogic(
        const QVector<float>& weights,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        const QVector<CandidateScores*> &candidateScoresPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(weights); ree;
        e = ErrorUtils::isNotEmpty(candidateScoresPntrs); ree;

        QVector<QVector<float>> candidateScoresMatrixVec;
        for (CandidateScores *cs : candidateScoresPntrs) {

            const QVector<float> vec = DiscriminantScoretron::scoreVectorLogic(
                useExtendedScores,
                useNeuralNetworkScores,
                cs
            );
            candidateScoresMatrixVec.push_back(vec);
        }

        QVector<QVector<float>*> candidateScoresMatrixVecPntrs;
        for (QVector<float> &f : candidateScoresMatrixVec) {
            candidateScoresMatrixVecPntrs.push_back(&f);
        }

        QVector<float> discriminantScores;
        e = ClassifierWeightsManager::applyWeights(
            candidateScoresMatrixVecPntrs,
            weights,
            &discriminantScores
        ); ree;

        e = ErrorUtils::isEqual(
            discriminantScores.size(),
            candidateScoresPntrs.size()
            ); ree;

        for (int i = 0; i < discriminantScores.size(); i++) {
            CandidateScores *cs = candidateScoresPntrs.at(i);
            cs->discriminantScore = discriminantScores.at(i);
        }

        ERR_RETURN
    }

}//namespace
Err DiscriminantScoretron::applyWeights(
    const QVector<float>& weights,
    bool useExtendedScores,
    bool useNeuralNetworkScores,
    QVector<CandidateScores*> *candidateScoresPntrs
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(weights); ree;
    e = ErrorUtils::isFalse(candidateScoresPntrs->isEmpty()); ree;

    QVector<QVector<CandidateScores*>> candidateScoresPntrsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
        *candidateScoresPntrs,
        ParallelUtils::numberOfAvailableSystemProcessors(),
        &candidateScoresPntrsTranched
    ); ree;

    const auto applyWeightsBinder = std::bind(
        applyWeightsLogic,
        weights,
        useExtendedScores,
        useNeuralNetworkScores,
        std::placeholders::_1
        );

    QFuture<Err> result = QtConcurrent::mapped(
        candidateScoresPntrsTranched,
        applyWeightsBinder
        );
    result.waitForFinished();

    for (Err res : result) {
        e = res; ree;
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

            // CandidateScores::Features::CosineSimSum100,
            CandidateScores::Features::CosineSimSum100GreaterThan80,
            CandidateScores::Features::AllignedMaxIndexesCount,
            CandidateScores::Features::CosineSimSpectrumOverTimeCubed,
            CandidateScores::Features::CosineSimSpectrumStDev,
            CandidateScores::Features::CosineSim100MS1,
            CandidateScores::Features::CosineSim100MS1Iso1,
            CandidateScores::Features::CosineSim100MS1Iso2,
            CandidateScores::Features::CosineSim100MS1PreMono,

            // CandidateScores::Features::CosineSimSum100,
            // CandidateScores::Features::CosineSimSum100GreaterThan80,
            // CandidateScores::Features::AllignedMaxIndexesCount,
            // CandidateScores::Features::CosineSim100MS1,
            // CandidateScores::Features::CosineSimSpectrumCubed,
            // CandidateScores::Features::KlDivSpectrumCubeRoot,
            CandidateScores::Features::CosineSimSum45,
            CandidateScores::Features::CosineSimSum20,
            CandidateScores::Features::CosineSimSumTop,
            CandidateScores::Features::CosineSimSumBottom,
            CandidateScores::Features::TopBottomRatio,
            CandidateScores::Features::TopBottomRatioNorm,
            CandidateScores::Features::Charge,
            CandidateScores::Features::ScanTimeDelta,
            CandidateScores::Features::ScanTimePd
        };

        if (useNeuralNetworkScores) {
            QVector<float> &vec = candidateScores->featuresArray;
            return vec;
        }
        else if (useExtendedScores) {
            QVector<float> vec = candidateScores->selectFeaturesArrayFeatures(baseFeatures);
            vec.append(candidateScores->selectFeaturesArrayFeatures({
//                            CandidateScores::Features::CosineSimSum100Frequencies,
                            // CandidateScores::Features::CosineSimSum100MS1,
                            // CandidateScores::Features::CosineSim100MS1PreMono,
                            // CandidateScores::Features::CosineSim100MS1Iso1,
                            // CandidateScores::Features::CosineSim100MS1Iso2,
                            // CandidateScores::Features::PeakShapeRatio1,
                            // CandidateScores::Features::PeakShapeRatio2,
                            // CandidateScores::Features::PeakShapeRatio3,
                            // CandidateScores::Features::ShadowsCosineSimSum,
                            // CandidateScores::Features::CosineSimToAnchor1,
                            // CandidateScores::Features::CosineSimToAnchor2,
                            // CandidateScores::Features::CosineSimToAnchor3,
                            // CandidateScores::Features::CosineSimToAnchor4,
                            // CandidateScores::Features::MzAccuracy1,
                            // CandidateScores::Features::MzAccuracy2,
                            // CandidateScores::Features::MzAccuracy3,
                            // CandidateScores::Features::MzAccuracy4,
                            // CandidateScores::Features::CosineSimShadowsToAnchor1,
                            // CandidateScores::Features::CosineSimShadowsToAnchor2,
                            // CandidateScores::Features::CosineSimShadowsToAnchor3,
                            // CandidateScores::Features::CosineSimShadowsToAnchor4,
                            // CandidateScores::Features::CosineSimSpectrumOverTime,
                            // CandidateScores::Features::TheoFragmentCount,
                            // CandidateScores::Features::TotalIntensityLog
                    }));
            return vec;
        }

        const QVector<float> vec = candidateScores->selectFeaturesArrayFeatures(baseFeatures);

        return vec;
    }
