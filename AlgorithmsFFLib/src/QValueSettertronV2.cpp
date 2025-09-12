//
// Created by anichols on 1/9/24.
//

#include "QValueSettertronV2.h"

#include <boost/math/constants/constants.hpp>

#include "CandidateScoresV2.h"
#include "GlobalSettings.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"

namespace {

    QString buildCandidateKey(const CandidateScoresV2 &candidateScores) {
        const QString decoyToString = candidateScores.isDecoy ? "_1" : "_0";
        return candidateScores.targetDecoyCandidatePair->peptideStringWithMods()
               + QString::number(candidateScores.targetDecoyCandidatePair->charge()) + candidateScores.targetKey + decoyToString;
    }

    Err buildSetQValueForCandidateScoresInputs(
            const QValueSettertronV2::QValueScoreTypeV2 &qValueScoreType,
            QVector<CandidateScoresV2*> *candidateScores,
            QHash<PeptideSequenceChargeKey, double> *identifierVsTargets,
            QHash<PeptideSequenceChargeKey, double> *identifierVsDecoys

    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
        identifierVsTargets->clear();
        identifierVsDecoys->clear();

        for (const CandidateScoresV2 *cs : *candidateScores) {


            const PeptideSequenceChargeKey peptideSequenceChargeKey = buildCandidateKey(*cs);

            double classifierScore = -1.0;
            if (qValueScoreType == QValueSettertronV2::QValueScoreTypeV2::NNClassifierScore) {
                classifierScore = 1.0 / cs->classifierScore;
            }
            else if (qValueScoreType == QValueSettertronV2::QValueScoreTypeV2::DiscriminantScore) {
                classifierScore = cs->discriminantScore;
            }
            else if (qValueScoreType == QValueSettertronV2::QValueScoreTypeV2::CosineSimSumMeanCorrelation) {
            	classifierScore = cs->featuresArray[FTR::CosineSimSumMeanCorrelation];
            }
        	else {
        		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "QValueScoreTypeV2 score type not implemented";
        		rrr(eValueError);
        	}

            if (cs->isDecoy) {
                identifierVsDecoys->insert(peptideSequenceChargeKey, classifierScore);
                continue;
            }

            identifierVsTargets->insert(peptideSequenceChargeKey, classifierScore);

        }

        ERR_RETURN
    }

    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairsParallelLogic(
            const QHash<PeptideSequenceChargeKey, QPair<double, double>> &identifierVsQValueVsDecoyRatio,
            CandidateScoresV2* cs
            ) {

        ERR_INIT

        if (cs->isDecoy) {
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
            QVector<CandidateScoresV2*> *candidateScores
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
Err QValueSettertronV2::setQValueForCandidates(
        const QValueSettertronV2::QValueScoreTypeV2 &qValueScoreType,
        QVector<CandidateScoresV2*> *candidateScores
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
        const QValueSettertronV2::QValueScoreTypeV2 &qValueScoreType,
        const QVector<QPair<CandidateScoresV2Target*, CandidateScoresV2Decoy*>> &targetDecoyCandidateScorePairsPntrs,
        const QVector<double> &targetScores,
        const QVector<double> &decoyScores
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyCandidateScorePairsPntrs); ree;
        e = ErrorUtils::isNotEmpty(targetScores); ree;
        e = ErrorUtils::isNotEmpty(decoyScores); ree;

        for (const QPair<CandidateScoresV2Target*, CandidateScoresV2Decoy*> &pr : targetDecoyCandidateScorePairsPntrs) {
            float score = -1.0;
            if (QValueSettertronV2::QValueScoreTypeV2::DiscriminantScore == qValueScoreType) {
            	score = pr.first->discriminantScore;
            }
        	else if (QValueSettertronV2::QValueScoreTypeV2::NNClassifierScore == qValueScoreType) {
        		score = 1.0f / pr.first->classifierScore;
        	}
        	else if (QValueSettertronV2::QValueScoreTypeV2::CosineSimSumMeanCorrelation == qValueScoreType) {
	            score = pr.first->featuresArray[FTR::CosineSimSumMeanCorrelation];
            }
        	else {
        		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "QValueScoreTypeV2 score type not implemented";
        		rrr(eValueError);
        	}

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
Err QValueSettertronV2::setQValueForCandidates(
    const QValueScoreTypeV2 &qValueScoreType,
    QVector<QPair<CandidateScoresV2Target*, CandidateScoresV2Decoy*>> *targetDecoyCandidateScorePairsPntrs
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(targetDecoyCandidateScorePairsPntrs->isEmpty()); ree;

    QVector<double> targetScores;
    targetScores.reserve(targetDecoyCandidateScorePairsPntrs->size());
    QVector<double> decoyScores;
    decoyScores.reserve(targetDecoyCandidateScorePairsPntrs->size());
    for (const QPair<CandidateScoresV2Target*, CandidateScoresV2Decoy*> &pr : *targetDecoyCandidateScorePairsPntrs) {
        const float scoreTarget = QValueScoreTypeV2::DiscriminantScore == qValueScoreType
                                ? pr.first->discriminantScore
                                : 1.0 / pr.first->classifierScore;;
        const float scoreDecoy = QValueScoreTypeV2::DiscriminantScore == qValueScoreType
                               ? pr.second->discriminantScore
                               : 1.0 / pr.second->classifierScore;;
        targetScores.push_back(static_cast<double>(scoreTarget));
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

    QVector<QVector<QPair<CandidateScoresV2Target*, CandidateScoresV2Decoy*>>> targetDecoyCandidateScorePairsPntrsTranched;
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




