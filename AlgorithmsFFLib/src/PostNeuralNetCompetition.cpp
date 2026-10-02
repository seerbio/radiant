#include "PostNeuralNetCompetition.h"

#include "AminoAcids.h"
#include "CandidateScores.h"
#include "FragmentCompetition.h"
#include "TargetDecoyCandidatePair.h"

#include <QHash>
#include <QSet>

#include <algorithm>
#include <cmath>
#include <vector>

namespace {

bool decoy(const CandidateScores *candidate) {
    return candidate->isDecoy || candidate->targetDecoyCandidatePair->isDecoy();
}

QString precursorKey(const CandidateScores *candidate) {
    const auto origin = candidate->targetDecoyCandidatePair->peptideStringWithMods();
    const QString reported = candidate->isDecoy
        ? AminoAcids::mutatePenultimatePeptideResidues(origin) : origin;
    // Keep origin and reported sequence separate so distinct decoy origins that
    // share a synthetic display sequence remain distinct.
    return (decoy(candidate) ? QStringLiteral("D|") : QStringLiteral("T|"))
        + origin + "|" + reported + "|"
        + QString::number(candidate->targetDecoyCandidatePair->charge());
}

} // namespace

Error::Err PostNeuralNetCompetition::apply(
    int minimumSharedFragments, QVector<CandidateScores*> *candidates) {
    if (candidates == nullptr) return Error::eValueError;
    if (minimumSharedFragments == 0) return Error::eNoError;
    if (minimumSharedFragments < 2 || minimumSharedFragments > 12) return Error::eValueError;

    QSet<const CandidateScores*> seen;
    for (const CandidateScores *candidate : *candidates) {
        if (candidate == nullptr || candidate->targetDecoyCandidatePair == nullptr
            || seen.contains(candidate)
            || !std::isfinite(candidate->classifierScore)
            || candidate->classifierScore < 0.0 || candidate->classifierScore > 1.0) {
            return Error::eValueError;
        }
        seen.insert(candidate);
    }
    if (candidates->isEmpty()) return Error::eNoError;

    // Use the held-out NN score only for equal-evidence physical ties; the
    // competition routine handles its lower-is-better direction directly.
    const auto error = FragmentCompetition::removeCompetingCandidates(
        true,
        candidates,
        minimumSharedFragments,
        FragmentCompetition::TieBreakMetric::ClassifierScore);
    if (error != Error::eNoError) return error;

    QHash<QString, CandidateScores*> best;
    best.reserve(candidates->size());
    QVector<QString> keys;
    keys.reserve(candidates->size());
    for (CandidateScores *candidate : *candidates) {
        const QString key = precursorKey(candidate);
        keys.push_back(key);
        CandidateScores *previous = best.value(key, nullptr);
        if (previous == nullptr || candidate->classifierScore < previous->classifierScore) {
            best.insert(key, candidate);
        }
    }

    QVector<CandidateScores*> retained;
    retained.reserve(best.size());
    for (int index = 0; index < candidates->size(); ++index) {
        if (best.value(keys.at(index)) == candidates->at(index)) {
            retained.push_back(candidates->at(index));
        }
    }

    // Exact probability ordering is required after physical competition can
    // remove rows from an approximately sorted NN/LDA sequence.
    std::vector<CandidateScores*> ranked(retained.begin(), retained.end());
    std::stable_sort(ranked.begin(), ranked.end(), [](const auto *left, const auto *right) {
        if (left->classifierScore != right->classifierScore) {
            return left->classifierScore < right->classifierScore;
        }
        return decoy(left) > decoy(right);
    });

    struct Block { std::size_t begin, end; double q; };
    std::vector<Block> blocks;
    int targets = 0;
    int decoys = 0;
    for (std::size_t begin = 0; begin < ranked.size();) {
        std::size_t end = begin;
        while (end < ranked.size()
               && ranked[end]->classifierScore == ranked[begin]->classifierScore) {
            decoy(ranked[end]) ? ++decoys : ++targets;
            ++end;
        }
        blocks.push_back({begin, end,
            std::min(1.0, (decoys + 1.0) / std::max(targets, 1))});
        begin = end;
    }

    double q = 1.0;
    for (auto block = blocks.rbegin(); block != blocks.rend(); ++block) {
        q = std::min(q, block->q);
        for (std::size_t index = block->begin; index < block->end; ++index) {
            CandidateScores *candidate = ranked[index];
            candidate->qValue = q;
        }
    }

    *candidates = QVector<CandidateScores*>(ranked.begin(), ranked.end());
    return Error::eNoError;
}
