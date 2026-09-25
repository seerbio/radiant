#include "IdLevelQValueAnnotator.h"

#include "CandidateScores.h"
#include "ErrorUtils.h"
#include "MathUtils.h"
#include "TargetDecoyCandidatePair.h"

#include <algorithm>
#include <cmath>

namespace {

    enum class IdLevel {
        Precursor,
        Peptide,
        Protein
    };

    bool isDecoy(const CandidateScores *candidateScores) {
        return candidateScores->isDecoy || candidateScores->targetDecoyCandidatePair->isDecoy();
    }

    double scoreForFdr(const CandidateScores *candidateScores, const bool useClassifierScore) {
        if (useClassifierScore && candidateScores->classifierScore > 0.0) {
            return 1.0 / candidateScores->classifierScore;
        }
        return candidateScores->discriminantScore;
    }

    QString cleanedProteinGroup(const QString &proteinGroup) {
        QString key = proteinGroup.trimmed();
        key.replace(",", " ");
        return key;
    }

    QString keyForLevel(const CandidateScores *candidateScores, const IdLevel idLevel) {
        const QString peptide = candidateScores->targetDecoyCandidatePair->peptideStringWithMods();
        switch (idLevel) {
            case IdLevel::Precursor:
                return peptide + "|" + QString::number(candidateScores->targetDecoyCandidatePair->charge());
            case IdLevel::Peptide:
                return peptide;
            case IdLevel::Protein:
                return cleanedProteinGroup(candidateScores->proteinGroup);
        }
        return {};
    }

    void resetLevelFields(QVector<CandidateScores*> *candidateScores) {
        for (CandidateScores *candidateScore : *candidateScores) {
            candidateScore->precursorQValue = 1.0;
            candidateScore->peptideQValue = 1.0;
            candidateScore->proteinQValue = 1.0;
            candidateScore->isBestPrecursorCandidate = 0;
            candidateScore->isBestPeptideCandidate = 0;
            candidateScore->isBestProteinCandidate = 0;
        }
    }

    void setQValue(CandidateScores *candidateScores, const IdLevel idLevel, const double qValue) {
        switch (idLevel) {
            case IdLevel::Precursor:
                candidateScores->precursorQValue = qValue;
                break;
            case IdLevel::Peptide:
                candidateScores->peptideQValue = qValue;
                break;
            case IdLevel::Protein:
                candidateScores->proteinQValue = qValue;
                break;
        }
    }

    void setBestFlag(CandidateScores *candidateScores, const IdLevel idLevel) {
        switch (idLevel) {
            case IdLevel::Precursor:
                candidateScores->isBestPrecursorCandidate = 1;
                break;
            case IdLevel::Peptide:
                candidateScores->isBestPeptideCandidate = 1;
                break;
            case IdLevel::Protein:
                candidateScores->isBestProteinCandidate = 1;
                break;
        }
    }

    Err annotateLevel(
        const IdLevel idLevel,
        QVector<CandidateScores*> *candidateScores,
        const bool useClassifierScore
        ) {
        ERR_INIT

        QHash<QString, double> targetKeyVsBestScore;
        QHash<QString, double> decoyKeyVsBestScore;
        QHash<QString, CandidateScores*> scopedKeyVsBestCandidate;

        for (CandidateScores *candidateScore : *candidateScores) {
            const QString key = keyForLevel(candidateScore, idLevel);
            if (key.isEmpty()) {
                continue;
            }

            const bool decoy = isDecoy(candidateScore);
            const QString scopedKey = (decoy ? QStringLiteral("decoy|") : QStringLiteral("target|")) + key;
            const double score = scoreForFdr(candidateScore, useClassifierScore);

            CandidateScores *bestCandidate = scopedKeyVsBestCandidate.value(scopedKey, nullptr);
            if (bestCandidate == nullptr || score > scoreForFdr(bestCandidate, useClassifierScore)) {
                scopedKeyVsBestCandidate.insert(scopedKey, candidateScore);
            }

            QHash<QString, double> &keyVsBestScore = decoy ? decoyKeyVsBestScore : targetKeyVsBestScore;
            if (!keyVsBestScore.contains(key) || score > keyVsBestScore.value(key)) {
                keyVsBestScore.insert(key, score);
            }
        }

        for (CandidateScores *bestCandidate : scopedKeyVsBestCandidate) {
            setBestFlag(bestCandidate, idLevel);
        }

        if (targetKeyVsBestScore.isEmpty() || decoyKeyVsBestScore.isEmpty()) {
            ERR_RETURN
        }

        QVector<QPair<QString, double>> targetKeyVsBestScores;
        targetKeyVsBestScores.reserve(targetKeyVsBestScore.size());
        for (auto it = targetKeyVsBestScore.begin(); it != targetKeyVsBestScore.end(); ++it) {
            targetKeyVsBestScores.push_back({it.key(), it.value()});
        }

        QHash<QString, QPair<double, double>> keyVsQValueVsDecoyRatio;
        e = MathUtils::calculateQValue(
            targetKeyVsBestScores,
            decoyKeyVsBestScore,
            &keyVsQValueVsDecoyRatio
            ); ree;

        for (CandidateScores *candidateScore : *candidateScores) {
            if (isDecoy(candidateScore)) {
                continue;
            }

            const QString key = keyForLevel(candidateScore, idLevel);
            if (key.isEmpty() || !keyVsQValueVsDecoyRatio.contains(key)) {
                continue;
            }

            setQValue(candidateScore, idLevel, keyVsQValueVsDecoyRatio.value(key).first);
        }

        ERR_RETURN
    }

    QString localRtStratum(const CandidateScores *candidateScores, const double localRtBinSeconds) {
        if (localRtBinSeconds <= 0.0 || candidateScores == nullptr || candidateScores->scanTime < 0.0) {
            return QStringLiteral("global");
        }

        return QString::number(static_cast<int>(std::floor(candidateScores->scanTime / localRtBinSeconds)));
    }

    Err annotateLevelLocalRt(
        const IdLevel idLevel,
        QVector<CandidateScores*> *candidateScores,
        const bool useClassifierScore,
        const double localRtBinSeconds
        ) {

        ERR_INIT

        if (localRtBinSeconds <= 0.0) {
            e = annotateLevel(idLevel, candidateScores, useClassifierScore); ree;
            ERR_RETURN
        }

        QHash<QString, QVector<CandidateScores*>> stratumVsCandidates;
        for (CandidateScores *candidateScore : *candidateScores) {
            stratumVsCandidates[localRtStratum(candidateScore, localRtBinSeconds)].push_back(candidateScore);
        }

        for (auto it = stratumVsCandidates.begin(); it != stratumVsCandidates.end(); ++it) {
            QVector<CandidateScores*> stratumCandidates = it.value();
            e = annotateLevel(idLevel, &stratumCandidates, useClassifierScore); ree;
        }

        ERR_RETURN
    }

} // namespace

Err IdLevelQValueAnnotator::annotate(
    QVector<CandidateScores*> *candidateScores,
    const bool useClassifierScore,
    const double localRtBinSeconds
    ) {
    ERR_INIT

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    resetLevelFields(candidateScores);

    e = annotateLevelLocalRt(IdLevel::Precursor, candidateScores, useClassifierScore, localRtBinSeconds); ree;
    e = annotateLevelLocalRt(IdLevel::Peptide, candidateScores, useClassifierScore, localRtBinSeconds); ree;
    e = annotateLevelLocalRt(IdLevel::Protein, candidateScores, useClassifierScore, localRtBinSeconds); ree;

    ERR_RETURN
}
