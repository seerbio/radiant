//
// Created by anichols on 1/9/24.
//

#include "QValueSettertron.h"

#include <boost/math/constants/constants.hpp>

#include "CandidateScores.h"
#include "GlobalSettings.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"

#include <algorithm>
#include <cmath>
#include <functional>

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

    if (identifierVsTargetsScores.isEmpty() || identifierVsDecoys.isEmpty()) {
        for (CandidateScores *cs : *candidateScores) {
            if (cs == nullptr) {
                continue;
            }

            cs->qValue = 1.0;
            cs->decoyRatio = 1.0;
        }

        ERR_RETURN
    }

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

    struct TargetQValueEntry {
        CandidateScoresTarget *candidateScores = nullptr;
        double score = 0.0;
        double rawFdr = 1.0;
        double decoyRatio = 1.0;
    };

    double qValueScore(
        const QValueSettertron::QValueScoreType &qValueScoreType,
        const CandidateScores *candidateScores
        ) {

        return QValueSettertron::QValueScoreType::DiscriminantScore == qValueScoreType
            ? candidateScores->discriminantScore
            : 1.0 / candidateScores->classifierScore;
    }

    Err setLegacyQValueLogic(
        const QValueSettertron::QValueScoreType &qValueScoreType,
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &targetDecoyCandidateScorePairsPntrs,
        const QVector<double> &targetScores,
        const QVector<double> &decoyScores
        ) {

        ERR_INIT

        if (targetScores.isEmpty()) {
            ERR_RETURN
        }
        e = ErrorUtils::isNotEmpty(decoyScores); ree;

        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : targetDecoyCandidateScorePairsPntrs) {

            // Candidate pairs that originated from library decoys contribute to the decoy distribution only.
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

    Err setLegacyQValuesForCandidatePairs(
        const QValueSettertron::QValueScoreType &qValueScoreType,
        QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairsPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(targetDecoyCandidateScorePairsPntrs->isEmpty()); ree;

        const int n = targetDecoyCandidateScorePairsPntrs->size();

        // Count decoy entries from the library -- both sides of these pairs should be treated as decoys.
        const int num_original_decoys = std::count_if(
            targetDecoyCandidateScorePairsPntrs->begin(),
            targetDecoyCandidateScorePairsPntrs->end(),
            [](const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr) {
                return pr.first->targetDecoyCandidatePair->isDecoy();
            });

        QVector<double> targetScores;
        QVector<double> decoyScores;
        targetScores.reserve(n - num_original_decoys);
        decoyScores.reserve(n + num_original_decoys);

        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : *targetDecoyCandidateScorePairsPntrs) {
            const float scoreTarget = QValueSettertron::QValueScoreType::DiscriminantScore == qValueScoreType
                                    ? pr.first->discriminantScore
                                    : 1.0 / pr.first->classifierScore;;
            const float scoreDecoy = QValueSettertron::QValueScoreType::DiscriminantScore == qValueScoreType
                                   ? pr.second->discriminantScore
                                   : 1.0 / pr.second->classifierScore;;
            const bool isOriginalDecoy = pr.first->targetDecoyCandidatePair->isDecoy();
            (isOriginalDecoy ? decoyScores : targetScores).push_back(static_cast<double>(scoreTarget));
            decoyScores.push_back(static_cast<double>(scoreDecoy));
        }

        if (targetScores.isEmpty()) {
            ERR_RETURN
        }
        e = ErrorUtils::isNotEmpty(decoyScores); ree;

        std::sort(targetScores.begin(), targetScores.end());
        std::sort(decoyScores.begin(), decoyScores.end());

        QVector<QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>>> targetDecoyCandidateScorePairsPntrsTranches;
        e = ParallelUtils::trancheVectorForParallelization(
            *targetDecoyCandidateScorePairsPntrs,
            std::min(ParallelUtils::numberOfAvailableSystemProcessors(), targetDecoyCandidateScorePairsPntrs->size()),
            &targetDecoyCandidateScorePairsPntrsTranches
            );

        const auto binderLogic = std::bind(
            setLegacyQValueLogic,
            qValueScoreType,
            std::placeholders::_1,
            targetScores,
            decoyScores
            );

        QFuture<Err> futures = QtConcurrent::mapped(targetDecoyCandidateScorePairsPntrsTranches, binderLogic);
        futures.waitForFinished();

        for (Err result : futures) {
            e = result; ree;
        }

        ERR_RETURN
    }

    Err setMonotoneQValuesForCandidatePairs(
        const QValueSettertron::QValueScoreType &qValueScoreType,
        QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairsPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(targetDecoyCandidateScorePairsPntrs->isEmpty()); ree;

        const int n = targetDecoyCandidateScorePairsPntrs->size();

        QVector<TargetQValueEntry> targetScores;
        QVector<double> decoyScores;
        targetScores.reserve(n);
        decoyScores.reserve(n);

        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : *targetDecoyCandidateScorePairsPntrs) {
            pr.first->qValue = 1.0;
            pr.first->decoyRatio = 1.0;
            pr.second->qValue = 1.0;
            pr.second->decoyRatio = 1.0;

            const double scoreTarget = qValueScore(qValueScoreType, pr.first);
            const double scoreDecoy = qValueScore(qValueScoreType, pr.second);
            const bool isOriginalDecoy = pr.first->targetDecoyCandidatePair->isDecoy();
            if (isOriginalDecoy) {
                decoyScores.push_back(std::max(scoreTarget, scoreDecoy));
                continue;
            }

            if (scoreTarget >= scoreDecoy) {
                targetScores.push_back({pr.first, scoreTarget, 1.0, 1.0});
            }
            else {
                decoyScores.push_back(scoreDecoy);
            }
        }

        if (targetScores.isEmpty()) {
            ERR_RETURN
        }

        std::sort(
            targetScores.begin(),
            targetScores.end(),
            [](const TargetQValueEntry &l, const TargetQValueEntry &r) {
                return l.score > r.score;
            });
        std::sort(decoyScores.begin(), decoyScores.end(), std::greater<>());

        int decoyCount = 0;
        for (int targetIndex = 0; targetIndex < targetScores.size(); targetIndex++) {
            const double score = targetScores.at(targetIndex).score;
            while (decoyCount < decoyScores.size() && decoyScores.at(decoyCount) >= score) {
                decoyCount++;
            }

            const int targetCount = targetIndex + 1;
            targetScores[targetIndex].rawFdr = std::min(
                1.0,
                static_cast<double>(decoyCount) / static_cast<double>(std::max(1, targetCount))
                );
            targetScores[targetIndex].decoyRatio = std::min(
                1.0,
                static_cast<double>(decoyCount) / static_cast<double>(std::max(1, targetScores.size()))
                );
        }

        double runningQValue = 1.0;
        for (int targetIndex = targetScores.size() - 1; targetIndex >= 0; targetIndex--) {
            runningQValue = std::min(runningQValue, targetScores.at(targetIndex).rawFdr);
            targetScores[targetIndex].candidateScores->qValue = runningQValue;
            targetScores[targetIndex].candidateScores->decoyRatio = targetScores.at(targetIndex).decoyRatio;
        }

        ERR_RETURN
    }

    QString targetKeyChargeStratum(const CandidateScoresTarget *candidateScores, const double localRtBinSeconds) {
        if (candidateScores == nullptr || candidateScores->targetDecoyCandidatePair == nullptr) {
            return QString();
        }

        QString stratum = candidateScores->targetKey
                          + QStringLiteral("_z")
                          + QString::number(candidateScores->targetDecoyCandidatePair->charge());
        if (localRtBinSeconds > 0.0 && candidateScores->scanTime >= 0.0) {
            const int rtBin = static_cast<int>(std::floor(candidateScores->scanTime / localRtBinSeconds));
            stratum += QStringLiteral("_rt")
                       + QString::number(rtBin);
        }

        return stratum;
    }

    Err setMonotoneQValuesForCandidatePairsByTargetKeyAndCharge(
        const QValueSettertron::QValueScoreType &qValueScoreType,
        QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairsPntrs,
        const double localRtBinSeconds
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(targetDecoyCandidateScorePairsPntrs->isEmpty()); ree;

        QHash<QString, QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>>> stratumVsPairs;
        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr : *targetDecoyCandidateScorePairsPntrs) {
            if (pr.first == nullptr) {
                continue;
            }
            stratumVsPairs[targetKeyChargeStratum(pr.first, localRtBinSeconds)].push_back(pr);
        }

        e = ErrorUtils::isFalse(stratumVsPairs.isEmpty()); ree;

        for (auto it = stratumVsPairs.begin(); it != stratumVsPairs.end(); ++it) {
            QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> stratumPairs = it.value();
            if (stratumPairs.isEmpty()) {
                continue;
            }

            e = setMonotoneQValuesForCandidatePairs(
                qValueScoreType,
                &stratumPairs
                ); ree;
        }

        ERR_RETURN
    }

}//namespace
Err QValueSettertron::setQValueForCandidates(
    const QValueScoreType &qValueScoreType,
    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairsPntrs,
    bool useMonotoneQValues,
    bool stratifyByTargetKey,
    double localRtBinSeconds
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(targetDecoyCandidateScorePairsPntrs->isEmpty()); ree;

    if (useMonotoneQValues && stratifyByTargetKey) {
        e = setMonotoneQValuesForCandidatePairsByTargetKeyAndCharge(
            qValueScoreType,
            targetDecoyCandidateScorePairsPntrs,
            localRtBinSeconds
            ); ree;
    }
    else if (useMonotoneQValues) {
        e = setMonotoneQValuesForCandidatePairs(qValueScoreType, targetDecoyCandidateScorePairsPntrs); ree;
    }
    else {
        e = setLegacyQValuesForCandidatePairs(qValueScoreType, targetDecoyCandidateScorePairsPntrs); ree;
    }

    ERR_RETURN
}
