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
    // Distinct origins can share the same synthetic display sequence. Retain
    // both parts so library and synthetic decoys of one origin also remain
    // distinct when their reported sequences differ.
    return (decoy(candidate) ? QStringLiteral("D|") : QStringLiteral("T|"))
        + origin + "|" + reported + "|" + QString::number(candidate->targetDecoyCandidatePair->charge());
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

    // The existing physical competition uses the LDA score only to break
    // equal-evidence ties. Use the fixed held-out NN score for that tie here,
    // then restore every original LDA score, including removed candidates.
    std::vector<std::pair<CandidateScores*, double>> originalScores;
    originalScores.reserve(candidates->size());
    for (CandidateScores *candidate : *candidates) {
        originalScores.emplace_back(candidate, candidate->discriminantScore);
        candidate->discriminantScore = -candidate->classifierScore;
    }
    const auto error = FragmentCompetition::removeCompetingCandidates(
        true, candidates, minimumSharedFragments);
    for (const auto &[candidate, score] : originalScores) candidate->discriminantScore = score;
    if (error != Error::eNoError) return error;

    QHash<QString, CandidateScores*> best;
    for (CandidateScores *candidate : *candidates) {
        const QString key = precursorKey(candidate);
        CandidateScores *previous = best.value(key, nullptr);
        if (previous == nullptr || candidate->classifierScore < previous->classifierScore) {
            best.insert(key, candidate);
        }
    }
    QVector<CandidateScores*> retained;
    retained.reserve(best.size());
    for (CandidateScores *candidate : *candidates) {
        if (best.value(precursorKey(candidate)) == candidate) retained.push_back(candidate);
    }
    // Use exact probability order for confidence and output. The earlier
    // approximate NN/LDA comparator is not transitive: removing an intermediate
    // candidate can make a previously adjacent-sorted sequence appear unsorted.
    std::vector<CandidateScores*> ranked(retained.begin(), retained.end());
    std::stable_sort(ranked.begin(), ranked.end(), [](const auto *left, const auto *right) {
        if (left->classifierScore != right->classifierScore)
            return left->classifierScore < right->classifierScore;
        return decoy(left) > decoy(right);
    });
    struct Block { std::size_t begin, end; double q; };
    std::vector<Block> blocks;
    int targets = 0, decoys = 0;
    for (std::size_t begin = 0; begin < ranked.size();) {
        std::size_t end = begin;
        while (end < ranked.size()
               && ranked[end]->classifierScore == ranked[begin]->classifierScore) {
            decoy(ranked[end]) ? ++decoys : ++targets;
            ++end;
        }
        blocks.push_back({begin, end, std::min(1.0, (decoys + 1.0) / std::max(targets, 1))});
        begin = end;
    }
    double q = 1.0;
    for (auto block = blocks.rbegin(); block != blocks.rend(); ++block) {
        q = std::min(q, block->q);
        for (std::size_t index = block->begin; index < block->end; ++index) {
            CandidateScores *candidate = ranked[index];
            candidate->qValue = q;
            candidate->precursorQValue = q;
            candidate->isBestPrecursorCandidate = 1;
        }
    }
    *candidates = QVector<CandidateScores*>(ranked.begin(), ranked.end());
    return Error::eNoError;
}
