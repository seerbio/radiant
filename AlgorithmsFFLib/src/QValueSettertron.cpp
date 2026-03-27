//
// Created by anichols on 1/9/24.
//

#include "QValueSettertron.h"

#include <boost/math/constants/constants.hpp>

#include "CandidateScores.h"
#include "GlobalSettings.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"

namespace {

    QString buildCandidateKey(const CandidateScores &candidateScores) {
        // Note: we still use target/decoy keying here even if the candidate pair is originally a decoy!
        // This is necessary to key the pair elements separately -- target/decoy status decisions ignore this string.
        const QString decoyToString = candidateScores.isDecoy ? "_1" : "_0";
        return candidateScores.targetDecoyCandidatePair->peptideStringWithMods()
               + QString::number(candidateScores.targetDecoyCandidatePair->charge()) + candidateScores.targetKey + decoyToString;
    }

    Err buildSetQValueForCandidateScoresInputs(
            const QValueSettertron::QValueScoreType &qValueScoreType,
            QVector<CandidateScores*> *candidateScores,
            QHash<PeptideSequenceChargeKey, double> *identifierVsTargets,
            QHash<PeptideSequenceChargeKey, double> *identifierVsDecoys

    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
        identifierVsTargets->clear();
        identifierVsDecoys->clear();

        for (const CandidateScores *cs : *candidateScores) {


            const PeptideSequenceChargeKey peptideSequenceChargeKey = buildCandidateKey(*cs);

            double classifierScore = -1.0;
            if (qValueScoreType == QValueSettertron::QValueScoreType::NNClassifierScore) {
                classifierScore = 1.0 / cs->classifierScore;
            }
            else {
                classifierScore = cs->discriminantScore;
            }

            if (cs->isDecoy || cs->targetDecoyCandidatePair->isDecoy()) {
                identifierVsDecoys->insert(peptideSequenceChargeKey, classifierScore);
                continue;
            }

            identifierVsTargets->insert(peptideSequenceChargeKey, classifierScore);

        }

        ERR_RETURN
    }

    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairsParallelLogic(
            const QHash<PeptideSequenceChargeKey, QPair<double, double>> &identifierVsQValueVsDecoyRatio,
            CandidateScores* cs
            ) {

        ERR_INIT

        if (cs->isDecoy || cs->targetDecoyCandidatePair->isDecoy()) {
            ERR_RETURN
        }

        const PeptideSequenceChargeKey peptideSequenceChargeKey = buildCandidateKey(*cs);

        e = ErrorUtils::isTrue(identifierVsQValueVsDecoyRatio.contains(peptideSequenceChargeKey)); ree;

        const QPair<double, double> qValVsDecoyRatio = identifierVsQValueVsDecoyRatio.value(peptideSequenceChargeKey);

        cs->qValue = qValVsDecoyRatio.first;
        cs->decoyRatio = qValVsDecoyRatio.second;

        ERR_RETURN
    }

    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            const QHash<PeptideSequenceChargeKey, QPair<double, double>> &identifierVsQValueVsDecoyRatio,
            QVector<CandidateScores*> *candidateScores
            ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

        const auto binderLogic = std::bind(
            setQValueAndDecoyRatioToTargetDecoyCandidatePairsParallelLogic,
            identifierVsQValueVsDecoyRatio,
            std::placeholders::_1
            );

        QFuture<Err> futures = QtConcurrent::mapped(*candidateScores, binderLogic);
        futures.waitForFinished();

        for (Err result : futures) {
            e = result; ree;
        }

        ERR_RETURN
    }

}//namespace
Err QValueSettertron::setQValueForCandidates(
        const QValueSettertron::QValueScoreType &qValueScoreType,
        QVector<CandidateScores*> *candidateScores
        ){

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    QElapsedTimer et;
    et.start();

    QHash<PeptideSequenceChargeKey, double> identifierVsTargets;
    QHash<PeptideSequenceChargeKey, double> identifierVsDecoys;
    e = buildSetQValueForCandidateScoresInputs(
            qValueScoreType,
            candidateScores,
            &identifierVsTargets,
            &identifierVsDecoys
            ); ree;

    const QVector<QPair<PeptideSequenceChargeKey, double>> identifierVsTargetsScores
                                            = ParallelUtils::convertHashToVectorPairs(identifierVsTargets);

    QHash<PeptideSequenceChargeKey, QPair<double, double>> identifierVsQValueVsDecoyRatio;
    e = MathUtils::calculateQValue(
            identifierVsTargetsScores,
            identifierVsDecoys,
            &identifierVsQValueVsDecoyRatio
    ); ree;

    e = setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            identifierVsQValueVsDecoyRatio,
            candidateScores
    ); ree;

    ERR_RETURN
}

namespace {

    Err setQValueLogic(
        const QValueSettertron::QValueScoreType &qValueScoreType,
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &targetDecoyCandidateScorePairsPntrs,
        const QVector<double> &targetScores,
        const QVector<double> &decoyScores
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyCandidateScorePairsPntrs); ree;
        e = ErrorUtils::isNotEmpty(targetScores); ree;
        e = ErrorUtils::isNotEmpty(decoyScores); ree;

        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : targetDecoyCandidateScorePairsPntrs) {
            if (pr.first->targetDecoyCandidatePair->isDecoy()) {
                continue;
            }

            float score = QValueSettertron::QValueScoreType::DiscriminantScore == qValueScoreType
                    ? pr.first->discriminantScore
                    : 1.0 / pr.first->classifierScore;

            auto decoyIndexLowest = std::lower_bound(decoyScores.begin(), decoyScores.end(), score);
            if (decoyIndexLowest > decoyScores.begin()) {
                if (*(decoyIndexLowest - 1) > decoyScores.front()) {
                    score = *(--decoyIndexLowest);
                }
            }

            auto targetIndexLowest = std::lower_bound(targetScores.begin(), targetScores.end(), score);

            const long targetCount = std::distance(targetIndexLowest, targetScores.end());
            const long decoyCount = std::distance(decoyIndexLowest, decoyScores.end());

            const double qvalue = std::min(
                    1.0,
                    (static_cast<double>(std::max(static_cast<long>(0), decoyCount))) / static_cast<double>(std::max(static_cast<long>(1), targetCount))
                    );

            const double decoyRatio
                = std::min(1.0, (static_cast<double>(std::max(static_cast<long>(0), decoyCount)) / static_cast<double>(std::max(1, targetScores.size()))));

            pr.first->qValue = qvalue;
            pr.first->decoyRatio = decoyRatio;
        }

        ERR_RETURN
    }

}//namespace
Err QValueSettertron::setQValueForCandidates(
    const QValueScoreType &qValueScoreType,
    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairsPntrs
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(targetDecoyCandidateScorePairsPntrs->isEmpty()); ree;

    const int n = targetDecoyCandidateScorePairsPntrs->size();

    // Count decoy entries from the library -- both sides of these pairs should be treated as decoys!
    const int num_original_decoys = std::count_if(
        targetDecoyCandidateScorePairsPntrs->begin(),
        targetDecoyCandidateScorePairsPntrs->end(),
        [](const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr) {
            return pr.first->targetDecoyCandidatePair->isDecoy();
        });

    QVector<double> targetScores, decoyScores;
    targetScores.reserve(n - num_original_decoys);
    decoyScores.reserve(n + num_original_decoys);

    for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : *targetDecoyCandidateScorePairsPntrs) {
        const float scoreTarget = QValueScoreType::DiscriminantScore == qValueScoreType
                                ? pr.first->discriminantScore
                                : 1.0 / pr.first->classifierScore;;
        const float scoreDecoy = QValueScoreType::DiscriminantScore == qValueScoreType
                               ? pr.second->discriminantScore
                               : 1.0 / pr.second->classifierScore;;
        const bool isOriginalDecoy = pr.first->targetDecoyCandidatePair->isDecoy();
        (isOriginalDecoy ? decoyScores : targetScores).push_back(static_cast<double>(scoreTarget));
        decoyScores.push_back(static_cast<double>(scoreDecoy));
    }

    std::sort(targetScores.begin(), targetScores.end());
    std::sort(decoyScores.begin(), decoyScores.end());

    const auto setQValBinder = std::bind(
        setQValueLogic,
        qValueScoreType,
        std::placeholders::_1,
        targetScores,
        decoyScores
        );

    QVector<QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>>> targetDecoyCandidateScorePairsPntrsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
        *targetDecoyCandidateScorePairsPntrs,
        std::min(ParallelUtils::numberOfAvailableSystemProcessors(), targetDecoyCandidateScorePairsPntrs->size()),
        &targetDecoyCandidateScorePairsPntrsTranched
        );

    QFuture<Err> futures = QtConcurrent::mapped(
        targetDecoyCandidateScorePairsPntrsTranched,
        setQValBinder
        );
    futures.waitForFinished();

    for (Err result : futures) {
        e = result; ree;
    }

    ERR_RETURN
}




