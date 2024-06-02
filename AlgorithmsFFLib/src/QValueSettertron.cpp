//
// Created by anichols on 1/9/24.
//

#include "QValueSettertron.h"

#include "CandidateScores.h"
#include "GlobalSettings.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"

namespace {

    QString buildCandidateKey(const CandidateScores &candidateScores) {
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
            CandidateScores* cs
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
        QVector<CandidateScores *> *candidateScores
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
