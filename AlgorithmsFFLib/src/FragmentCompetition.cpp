#include "FragmentCompetition.h"

#include "CandidateScores.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

namespace {

constexpr int ionCount = 12;
constexpr double precursorTolerance = 5.0 * 1.0e-6;
constexpr double fragmentTolerance = 20.0 * 1.0e-6;
constexpr double minimumTraceCosine = 0.5;
constexpr double maximumApexWidthFraction = 0.5;

struct Candidate {
    double priority;
    QString peptide;
    int inputIndex;
    double mass;
    double charge;
    double apex;
    double width;
    bool isDecoy;
    std::array<double, ionCount> searched{};
    std::array<double, ionCount> weighted{};
    std::array<bool, ionCount> supported{};
    int supportedCount = 0;
};

bool isNoPeakPlaceholder(const CandidateScores &scores) {
    if (scores.scanNumber != -1 || !scores.integrations.isEmpty()
        || scores.featuresArray.size() != FeaturesSize
        || scores.featuresArray.at(CosineSimSum100) != 0.0f) {
        return false;
    }
    for (int ion = 0; ion < ionCount; ++ion) {
        if (scores.featuresArray.at(MzFoundMean1 + ion) != 0.0f) {
            return false;
        }
    }
    return true;
}

bool fragmentsMatch(double first, double second) {
    return std::abs(first - second) <= fragmentTolerance * std::max(first, second);
}

bool readEvidence(const FragmentCompetition::Evidence &e, int inputIndex, Candidate *candidate) {
    candidate->inputIndex = inputIndex;
    candidate->mass = e.mass;
    candidate->charge = e.charge;
    candidate->apex = e.apex;
    candidate->width = e.width;
    candidate->priority = e.priority;
    candidate->peptide = e.reportedPeptide;
    candidate->isDecoy = e.isDecoy;
    if (!std::isfinite(e.mass) || e.mass <= 0.0
        || !std::isfinite(e.charge) || e.charge <= 0.0
        || e.charge != std::floor(e.charge)
        || !std::isfinite(e.apex)
        || !std::isfinite(e.width) || e.width <= 0.0
        || !std::isfinite(e.priority)) return false;
    candidate->searched = e.searched;
    for (int ion = 0; ion < ionCount; ++ion) {
        candidate->supported[ion] = std::isfinite(e.searched[ion])
            && e.searched[ion] > 0.0
            && std::isfinite(e.intensity[ion]) && e.intensity[ion] > 0.0
            && std::isfinite(e.cosine[ion]) && e.cosine[ion] >= minimumTraceCosine;
        if (candidate->supported[ion]) {
            const double clippedCosine = std::min(e.cosine[ion], 1.0);
            candidate->weighted[ion] = e.intensity[ion] * (clippedCosine * clippedCosine);
        }
    }

    std::array<int, ionCount> ordered;
    std::iota(ordered.begin(), ordered.end(), 0);
    std::stable_sort(ordered.begin(), ordered.end(), [candidate](int left, int right) {
        return candidate->weighted[left] > candidate->weighted[right];
    });
    std::array<int, ionCount> retained{};
    for (int ion : ordered) {
        if (!candidate->supported[ion]) {
            continue;
        }
        bool duplicate = false;
        for (int previous = 0; previous < candidate->supportedCount; ++previous) {
            if (fragmentsMatch(candidate->searched[ion], candidate->searched[retained[previous]])) {
                duplicate = true;
                break;
            }
        }
        if (duplicate) {
            candidate->supported[ion] = false;
            candidate->weighted[ion] = 0.0;
        } else {
            retained[candidate->supportedCount++] = ion;
        }
    }
    return true;
}

bool shareEnoughFragments(
    const Candidate &first, const Candidate &second, int minimumSharedFragments) {
    if (first.supportedCount < minimumSharedFragments
        || second.supportedCount < minimumSharedFragments) {
        return false;
    }
    std::array<bool, ionCount> matchedFirst{};
    std::array<bool, ionCount> matchedSecond{};
    for (int left = 0; left < ionCount; ++left) {
        if (!first.supported[left]) continue;
        for (int right = 0; right < ionCount; ++right) {
            if (second.supported[right]
                && fragmentsMatch(first.searched[left], second.searched[right])) {
                matchedFirst[left] = true;
                matchedSecond[right] = true;
            }
        }
    }
    return std::count(matchedFirst.begin(), matchedFirst.end(), true) >= minimumSharedFragments
        && std::count(matchedSecond.begin(), matchedSecond.end(), true) >= minimumSharedFragments;
}

int findRoot(int index, std::vector<int> *parents) {
    while ((*parents)[index] != index) {
        (*parents)[index] = (*parents)[(*parents)[index]];
        index = (*parents)[index];
    }
    return index;
}

bool winsTie(const Candidate &first, const Candidate &second) {
    if (first.priority != second.priority) {
        return first.priority > second.priority;
    }
    if (first.isDecoy != second.isDecoy) {
        return first.isDecoy;
    }
    const int peptideOrder = QString::compare(first.peptide, second.peptide);
    if (peptideOrder != 0) return peptideOrder < 0;
    if (first.apex != second.apex) return first.apex < second.apex;
    if (first.mass != second.mass) return first.mass < second.mass;
    return first.inputIndex < second.inputIndex;
}

void resolve(
    std::vector<Candidate> &observed, int minimumSharedFragments,
    std::vector<bool> *retained) {
    auto &keep = *retained;
    const int count = static_cast<int>(observed.size());
    std::vector<int> ordered(count), parents(count);
    std::iota(ordered.begin(), ordered.end(), 0);
    std::iota(parents.begin(), parents.end(), 0);
    std::stable_sort(ordered.begin(), ordered.end(), [&observed](int left, int right) {
        if (observed[left].charge != observed[right].charge) {
            return observed[left].charge < observed[right].charge;
        }
        return observed[left].mass < observed[right].mass;
    });
    for (int position = 0; position < count; ++position) {
        const int firstIndex = ordered[position];
        const Candidate &first = observed[firstIndex];
        const double maximumMass = first.mass / (1.0 - precursorTolerance);
        for (int next = position + 1; next < count; ++next) {
            const int secondIndex = ordered[next];
            const Candidate &second = observed[secondIndex];
            if (first.charge != second.charge || second.mass > maximumMass) break;
            if (std::abs(first.apex - second.apex)
                > maximumApexWidthFraction * std::min(first.width, second.width)) continue;
            if (!shareEnoughFragments(first, second, minimumSharedFragments)) continue;
            const int firstRoot = findRoot(firstIndex, &parents);
            const int secondRoot = findRoot(secondIndex, &parents);
            if (firstRoot != secondRoot) parents[secondRoot] = firstRoot;
        }
    }

    // Linked lists avoid allocating a vector for every singleton component.
    std::vector<int> heads(count, -1), nextMember(count, -1), groupSizes(count, 0);
    for (int index = count - 1; index >= 0; --index) {
        const int root = findRoot(index, &parents);
        nextMember[index] = heads[root];
        heads[root] = index;
        ++groupSizes[root];
    }
    for (int root = 0; root < count; ++root) {
        if (groupSizes[root] < 2) continue;
        int winner = -1;
        double bestEvidence = 0.0;
        for (int index = heads[root]; index >= 0; index = nextMember[index]) {
            const Candidate &candidate = observed[index];
            keep[candidate.inputIndex] = false;
            double evidence = 0.0;
            for (int ion = 0; ion < ionCount; ++ion) {
                if (!candidate.supported[ion]) continue;
                bool shared = false;
                for (int other = heads[root]; other >= 0 && !shared; other = nextMember[other]) {
                    if (other == index) continue;
                    // Any searched fragment in another member prevents uniqueness,
                    // even if that other member has no support for that fragment.
                    for (double mz : observed[other].searched) {
                        if (std::isfinite(mz) && mz > 0.0
                            && fragmentsMatch(candidate.searched[ion], mz)) {
                            shared = true;
                            break;
                        }
                    }
                }
                if (!shared) evidence += candidate.weighted[ion];
            }
            if (winner < 0 || evidence > bestEvidence
                || (evidence == bestEvidence && winsTie(candidate, observed[winner]))) {
                winner = index;
                bestEvidence = evidence;
            }
        }
        if (winner >= 0 && bestEvidence > 0.0) keep[observed[winner].inputIndex] = true;
    }

}

bool readCandidate(const CandidateScores &scores, int inputIndex, Candidate *candidate) {
    if (scores.targetDecoyCandidatePair == nullptr || scores.featuresArray.size() != FeaturesSize
        || scores.integrations.size() < ionCount) return false;
    FragmentCompetition::Evidence e;
    e.mass = scores.featuresArray.at(Mass);
    e.charge = scores.featuresArray.at(Charge);
    e.apex = scores.scanTime;
    e.width = static_cast<double>(scores.scanTimeEnd) - scores.scanTimeStart;
    e.priority = scores.discriminantScore;
    e.isDecoy = scores.isDecoy || scores.targetDecoyCandidatePair->isDecoy();
    const auto origin = scores.targetDecoyCandidatePair->peptideStringWithMods();
    e.reportedPeptide = scores.isDecoy
        ? AminoAcids::mutatePenultimatePeptideResidues(origin) : origin;
    const auto ions = scores.isDecoy
        ? scores.targetDecoyCandidatePair->ms2IonsDecoy()
        : scores.targetDecoyCandidatePair->ms2IonsTarget();
    e.searched.fill(-1.0);
    for (int ion = 0; ion < std::min(ionCount, ions.size()); ++ion) {
        e.searched[ion] = ions.at(ion).mz;
        e.intensity[ion] = scores.integrations.at(ion);
        e.cosine[ion] = scores.featuresArray.at(CosineSimToAnchor1 + ion);
    }
    return readEvidence(e, inputIndex, candidate);
}

} // namespace

Error::Err FragmentCompetition::retainEvidence(
    const QVector<Evidence> &evidence, int minimumSharedFragments, QVector<int> *retained) {
    if (retained == nullptr || minimumSharedFragments < 2 || minimumSharedFragments > ionCount)
        return Error::eValueError;
    std::vector<Candidate> observed;
    observed.reserve(evidence.size());
    for (int index = 0; index < evidence.size(); ++index) {
        Candidate candidate;
        if (!readEvidence(evidence.at(index), index, &candidate)) return Error::eValueError;
        observed.push_back(std::move(candidate));
    }
    std::vector<bool> keep(evidence.size(), true);
    resolve(observed, minimumSharedFragments, &keep);
    QVector<int> output;
    output.reserve(evidence.size());
    for (int index = 0; index < evidence.size(); ++index) {
        if (keep[index]) output.push_back(index);
    }
    *retained = std::move(output);
    return Error::eNoError;
}

Error::Err FragmentCompetition::removeCompetingCandidates(
    bool enabled, QVector<CandidateScores*> *candidates, int minimumSharedFragments) {
    if (candidates == nullptr) return Error::eValueError;
    if (!enabled || candidates->isEmpty()) return Error::eNoError;
    if (minimumSharedFragments < 2 || minimumSharedFragments > ionCount)
        return Error::eValueError;
    std::vector<Candidate> observed;
    std::vector<bool> keep(candidates->size(), true);
    for (int input = 0; input < candidates->size(); ++input) {
        const CandidateScores *scores = candidates->at(input);
        if (scores == nullptr) {
            keep[input] = false;
            continue;
        }
        if (isNoPeakPlaceholder(*scores)) continue;
        Candidate candidate;
        if (!readCandidate(*scores, input, &candidate)) {
            qWarning() << "Invalid candidate for fragment competition at row" << input;
            return Error::eValueError;
        }
        observed.push_back(std::move(candidate));
    }
    resolve(observed, minimumSharedFragments, &keep);
    QVector<CandidateScores*> retained;
    retained.reserve(candidates->size());
    for (int index = 0; index < candidates->size(); ++index) {
        if (keep[index]) retained.push_back(candidates->at(index));
    }
    *candidates = std::move(retained);
    return Error::eNoError;
}
