#include "CandidatePoolSelection.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <tuple>

namespace {
QString unflanked(QString sequence) {
    int begin = 0, end = sequence.size();
    while (begin < end && sequence.at(begin) == '_') ++begin;
    while (end > begin && sequence.at(end - 1) == '_') --end;
    return sequence.mid(begin, end - begin);
}
}

Error::Err CandidatePoolSelection::rankedUnique(
    const QVector<Identity> &rows, int maximumCandidates, QVector<int> *indices) {
    if (indices == nullptr || maximumCandidates < 1) return Error::eValueError;
    for (const auto &row : rows) {
        if (row.reportedPeptide.isEmpty() || row.originPeptide.isEmpty()
            || row.targetKey.isEmpty() || row.charge <= 0 || row.view < 0
            || row.viewRank < 0 || !std::isfinite(row.apex)) return Error::eValueError;
    }
    QVector<int> ordered(rows.size());
    std::iota(ordered.begin(), ordered.end(), 0);
    std::stable_sort(ordered.begin(), ordered.end(), [&rows](int a, int b) {
        return std::tie(rows[a].viewRank, rows[a].view)
             < std::tie(rows[b].viewRank, rows[b].view);
    });
    using Key = std::tuple<QString, QString, int, bool, QString, double>;
    std::set<Key> seen;
    QVector<int> result;
    result.reserve(std::min(maximumCandidates, rows.size()));
    for (int index : ordered) {
        const auto &row = rows[index];
        if (seen.emplace(row.reportedPeptide, row.originPeptide, row.charge,
                         row.isDecoy, row.targetKey, row.apex).second) {
            result.push_back(index);
            if (result.size() == maximumCandidates) break;
        }
    }
    *indices = std::move(result);
    return Error::eNoError;
}

Error::Err CandidatePoolSelection::finalize(
    const QVector<FragmentCompetition::Evidence> &evidence,
    const QVector<double> &probabilities, int minimumSharedFragments,
    Confidence *confidence) {
    QVector<QString> origins;
    origins.reserve(evidence.size());
    for (const auto &row : evidence) origins.push_back(row.reportedPeptide);
    return finalize(evidence, origins, probabilities, minimumSharedFragments, confidence);
}

Error::Err CandidatePoolSelection::finalize(
    const QVector<FragmentCompetition::Evidence> &evidence,
    const QVector<QString> &originPeptides,
    const QVector<double> &probabilities, int minimumSharedFragments,
    Confidence *confidence) {
    if (confidence == nullptr || evidence.size() != probabilities.size()
        || evidence.size() != originPeptides.size()
        || minimumSharedFragments < 2 || minimumSharedFragments > 12) return Error::eValueError;
    for (int row = 0; row < evidence.size(); ++row) {
        const auto &e = evidence[row];
        if (!std::isfinite(probabilities[row]) || probabilities[row] < 0.0
            || probabilities[row] > 1.0 || unflanked(e.reportedPeptide).isEmpty()
            || unflanked(originPeptides[row]).isEmpty()
            || !std::isfinite(e.charge) || e.charge <= 0.0 || e.charge != std::floor(e.charge))
            return Error::eValueError;
    }
    QVector<int> ranked(evidence.size());
    std::iota(ranked.begin(), ranked.end(), 0);
    std::stable_sort(ranked.begin(), ranked.end(), [&probabilities](int a, int b) {
        return probabilities[a] < probabilities[b];
    });
    std::set<std::tuple<QString, QString, double, bool>> seen;
    QVector<int> source;
    QVector<FragmentCompetition::Evidence> unique;
    source.reserve(evidence.size());
    unique.reserve(evidence.size());
    for (int index : ranked) {
        const auto &e = evidence[index];
        const auto origin = e.isDecoy ? unflanked(originPeptides[index]) : QString();
        if (seen.emplace(unflanked(e.reportedPeptide), origin, e.charge, e.isDecoy).second) {
            source.push_back(index);
            unique.push_back(e);
            unique.last().priority = -probabilities[index];
        }
    }
    QVector<int> kept;
    const auto error = FragmentCompetition::retainEvidence(unique, minimumSharedFragments, &kept);
    if (error != Error::eNoError) return error;
    Confidence result;
    result.inputIndices.reserve(kept.size());
    for (int index : kept) result.inputIndices.push_back(source[index]);
    std::stable_sort(result.inputIndices.begin(), result.inputIndices.end(),
        [&probabilities, &evidence](int a, int b) {
            if (probabilities[a] != probabilities[b]) return probabilities[a] < probabilities[b];
            return evidence[a].isDecoy > evidence[b].isDecoy;
        });
    result.qValues.resize(kept.size());
    struct Block { int begin, end; double q; };
    QVector<Block> blocks;
    int targets = 0, decoys = 0;
    for (int begin = 0; begin < result.inputIndices.size();) {
        int end = begin;
        while (end < result.inputIndices.size()
               && probabilities[result.inputIndices[end]] == probabilities[result.inputIndices[begin]]) {
            evidence[result.inputIndices[end]].isDecoy ? ++decoys : ++targets;
            ++end;
        }
        blocks.push_back({begin, end, std::min(1.0, (decoys + 1.0) / std::max(targets, 1))});
        begin = end;
    }
    double q = 1.0;
    for (auto block = blocks.crbegin(); block != blocks.crend(); ++block) {
        q = std::min(q, block->q);
        for (int row = block->begin; row < block->end; ++row) result.qValues[row] = q;
    }
    *confidence = std::move(result);
    return Error::eNoError;
}
