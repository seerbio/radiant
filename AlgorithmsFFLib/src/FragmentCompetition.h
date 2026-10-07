#ifndef RADIANT_FRAGMENTCOMPETITION_H
#define RADIANT_FRAGMENTCOMPETITION_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"

#include <QString>
#include <QVector>

#include <array>

class CandidateScores;

class ALGORITHMSFFLIB_EXPORTS FragmentCompetition {
public:
    enum class TieBreakMetric {
        DiscriminantScore,
        ClassifierScore
    };

    // Already searched fragment coordinates; no library filtering or decoy
    // remutation. Priority is used only for equal unshared-evidence ties.
    struct Evidence {
        double mass = 0.0, charge = 0.0, apex = 0.0, width = 0.0;
        double priority = 0.0;
        bool isDecoy = false;
        QString reportedPeptide;
        std::array<double, 12> searched{};
        std::array<double, 12> intensity{};
        std::array<double, 12> cosine{};
    };

    // Returns surviving input indices in input order. Output changes only
    // on success. Uses exactly the same physical competition as native rows.
    static Error::Err retainEvidence(
        const QVector<Evidence> &evidence, int minimumSharedFragments,
        QVector<int> *retained);

    /**
     * Resolve coeluting isobaric assignments using unshared fragment evidence.
     *
     * Equal-charge candidates link at 5 ppm precursor mass, 20 ppm fragment
     * mass, at least minimumSharedFragments distinct supported fragments
     * (four by default), trace cosine >= 0.5,
     * and apex separation <= half the narrower peak width. Each connected
     * group keeps the greatest unshared intensity * cosine^2 evidence.
     * Groups with no unshared evidence are rejected; singletons are retained.
     * Ties prefer the selected metric, then decoys, peptide, apex, mass, and
     * input order. Classifier scores are interpreted as probabilities, so
     * lower values have higher priority.
     *
     * Preserves surviving pointers in input order and never edits their scores.
     * Unscored no-peak placeholders pass through to the normal NN input filter.
     * Null rows are discarded when enabled. Disabled is an exact no-op.
     */
    static Error::Err removeCompetingCandidates(
        bool enabled,
        QVector<CandidateScores*> *candidates,
        int minimumSharedFragments = 4,
        TieBreakMetric tieBreakMetric = TieBreakMetric::DiscriminantScore);
};

#endif
