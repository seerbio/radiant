//
// Created by anichols on 1/9/24.
//

#include "QValueSettertron.h"

#include "CandidateScores.h"
#include "GlobalSettings.h"
#include "ErrorUtils.h"

namespace {

    QString buildCandidateKey(const CandidateScores &candidateScores) {
        const QString decoyToString = candidateScores.isDecoy ? "_1" : "_0";
        return candidateScores.targetDecoyCandidatePair->peptideStringWithMods()
               + QString::number(candidateScores.targetDecoyCandidatePair->charge()) + candidateScores.targetKey + decoyToString;
    }

    Err buildsetQValueForCandidateScoresInputs(
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

    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            const QHash<PeptideSequenceChargeKey, double> &identifierVsQValue,
            const QHash<PeptideSequenceChargeKey, double> &identifierVsDecoyRatio,
            QVector<CandidateScores*> *candidateScores
    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

        for (CandidateScores *cs : *candidateScores) {

            if (cs->isDecoy) {
                continue;
            }

            const PeptideSequenceChargeKey peptideSequenceChargeKey = buildCandidateKey(*cs);

            e = ErrorUtils::isTrue(identifierVsQValue.contains(peptideSequenceChargeKey)); ree;
            e = ErrorUtils::isTrue(identifierVsDecoyRatio.contains(peptideSequenceChargeKey)); ree;

            cs->qValue = identifierVsQValue.value(peptideSequenceChargeKey);
            cs->decoyRatio = identifierVsDecoyRatio.value(peptideSequenceChargeKey);
        }

        ERR_RETURN
    }

}//namespace
Err QValueSettertron::setQValueForCandidates(
        const QValueSettertron::QValueScoreType &qValueScoreType,
        QVector<CandidateScores *> *candidateScores)
{

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    QElapsedTimer et;
    et.start();

    QHash<PeptideSequenceChargeKey, double> identifierVsTargets;
    QHash<PeptideSequenceChargeKey, double> identifierVsDecoys;
    e = buildsetQValueForCandidateScoresInputs(
            qValueScoreType,
            candidateScores,
            &identifierVsTargets,
            &identifierVsDecoys
    ); ree;

    QHash<PeptideSequenceChargeKey, double> identifierVsQValue;
    QHash<PeptideSequenceChargeKey, double> identifierVsDecoyRatio;
    e = MathUtils::calculateQValue(
            identifierVsTargets,
            identifierVsDecoys,
            &identifierVsQValue,
            &identifierVsDecoyRatio
    ); ree;

    e = setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            identifierVsQValue,
            identifierVsDecoyRatio,
            candidateScores
    ); ree;

    ERR_RETURN


}
