//
// Created by anichols on 1/9/24.
//

#include "DiscriminantScoretron.h"

#include "CandidateScores.h"
#include "ClassifierWeightsManager.h"
#include "ParallelUtils.h"

namespace {

    struct TargetDecoyCandidateScores {
        CandidateScores *candidateScoresTarget;
        CandidateScores *candidateScoresDecoy;
        ScoresTargets *scoresTarget;
        ScoresDecoys  *scoresDecoy;
    };

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

    QPair<Err, BuildClassifierParallelOutput> buildClassifierDataParallel(const BuildClassiferParallelInput &input) {

        ERR_INIT

        BuildClassifierParallelOutput output;

        e = ClassifierWeightsManager::buildDataClassifier1(
                input.scoresTargets,
                input.scoresDecoys,
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

    QPair<Err, QVector<float>> scoreVectorParallel(
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            CandidateScores* candidateScores
    ) {

        ERR_INIT

        if (useNeuralNetworkScores) {
            QVector<float> &vec = candidateScores->featuresArray;
            return {e, vec};
        }
        else if (useExtendedScores) {
            QVector<float> vec = candidateScores->featuresArray.mid(0, CandidateScores::Features::CosineSimToAnchor1);
            return {e, vec};
        }

        const QVector<float> vec = candidateScores->featuresArray.mid(0, CandidateScores::Features::ScanIonCount);
        return {e, vec};

    }

}//namespace
Err DiscriminantScoretron::setDiscriminantScoreForCandidates(
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        QVector<CandidateScores*> *candidateScoresPntrs
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScoresPntrs->isEmpty()); ree;

    QElapsedTimer et;
    et.start();

    const auto scoreBuildBinder = std::bind(
            scoreVectorParallel,
            useExtendedScores,
            useNeuralNetworkScores,
            std::placeholders::_1
    );

    QFuture<QPair<Err, QVector<float>>> futuresScoreBuilder = QtConcurrent::mapped(
            *candidateScoresPntrs,
            scoreBuildBinder
    );
    futuresScoreBuilder.waitForFinished();

    QVector<QVector<float>> scoreVectors;
    for (const QPair<Err, QVector<float>> &res : futuresScoreBuilder) {
        e = res.first; ree;
        scoreVectors.push_back(res.second);
    }

    e = ErrorUtils::isEqual(candidateScoresPntrs->size(), scoreVectors.size()); ree;

    QMap<QString, TargetDecoyCandidateScores> keyVsTargetDecoyCandidateScores;

    for(int i = 0; i < scoreVectors.size(); i++) {

        CandidateScores *cs = (*candidateScoresPntrs)[i];

        const QString key = cs->targetDecoyCandidatePair->peptideStringWithMods()
                            + QString::number(cs->targetDecoyCandidatePair->charge()) + cs->targetKey;
        QVector<float> *sv = &scoreVectors[i];

        if (cs->isDecoy) {
            keyVsTargetDecoyCandidateScores[key].scoresDecoy = sv;
            keyVsTargetDecoyCandidateScores[key].candidateScoresDecoy = cs;
            continue;
        }

        keyVsTargetDecoyCandidateScores[key].candidateScoresTarget = cs;
        keyVsTargetDecoyCandidateScores[key].scoresTarget = sv;
    }
    qDebug() << "build key vs scores" << et.restart() << "mSec";

    QVector<ScoresTargets*> scoresTargets;
    QVector<ScoresDecoys*> scoresDecoys;
    QVector<CandidateScores*> candidateScoresTargetsPtrs;
    QVector<CandidateScores*> candidateScoresDecoysPtrs;

    for (auto it = keyVsTargetDecoyCandidateScores.begin(); it != keyVsTargetDecoyCandidateScores.end(); it++) {

        const QString &key = it.key();
        const TargetDecoyCandidateScores &tdcs = it.value();

        e = ErrorUtils::isEqual(tdcs.scoresTarget->size(), tdcs.scoresDecoy->size());
        if (e != eNoError) {
            qDebug() << "target decoys not paired for key" << key;
            rrr(eValueError);
        }

        scoresTargets.push_back(tdcs.scoresTarget);
        scoresDecoys.push_back(tdcs.scoresDecoy);
        candidateScoresTargetsPtrs.push_back(tdcs.candidateScoresTarget);
        candidateScoresDecoysPtrs.push_back(tdcs.candidateScoresDecoy);
    }

    QVector<BuildClassiferParallelInput> inputs;
    e = buildParallelInput(
            scoresTargets,
            scoresDecoys,
            &inputs
    ); ree;

    QFuture<QPair<Err, BuildClassifierParallelOutput>> futures = QtConcurrent::mapped(
            inputs,
            buildClassifierDataParallel
    );
    futures.waitForFinished();

    QVector<QVector<float>> A;
    QVector<float> b;
    for (const QPair<Err, BuildClassifierParallelOutput> &result : futures) {
        e = processParallelResults(
                result,
                inputs.size(),
                &A,
                &b
        ); ree;
    }

    QVector<float> weights;
    e = ClassifierWeightsManager::fitWeights(A, b, &weights); ree;
    qDebug() << "fit weights" << et.restart() << "mSec";

    qDebug() << "Weights:" << weights;
    qDebug() << "b:" << b;

//#define ENUMERATE_B
#ifdef ENUMERATE_B
    for (int i = 0; i < b.size(); i++) {
        qDebug() << i << b.at(i);
    }
    einfo;
#endif

    QVector<float> discScoreTargets;
    e = ClassifierWeightsManager::applyWeights(scoresTargets, weights, &discScoreTargets); ree;

    QVector<float> discScoreDecoys;
    e = ClassifierWeightsManager::applyWeights(scoresDecoys, weights, &discScoreDecoys); ree;

    qDebug() << "apply weights scores" << et.restart() << "mSec";

    e = ErrorUtils::isEqual(scoresTargets.size(), scoresDecoys.size()); ree;
    e = ErrorUtils::isEqual(scoresTargets.size(), candidateScoresTargetsPtrs.size()); ree;
    e = ErrorUtils::isEqual(scoresDecoys.size(), candidateScoresDecoysPtrs.size()); ree;
    e = ErrorUtils::isEqual(discScoreTargets.size(), scoresTargets.size());
    e = ErrorUtils::isEqual(discScoreDecoys.size(), scoresDecoys.size());

    for (int i = 0; i < scoresDecoys.size(); i++) {

        CandidateScores* csTarget = candidateScoresTargetsPtrs.at(i);
        csTarget->discriminantScore = discScoreTargets.at(i);

        CandidateScores* csDecoy = candidateScoresDecoysPtrs.at(i);
        csDecoy->discriminantScore = discScoreDecoys.at(i);
    }

    ERR_RETURN

}
