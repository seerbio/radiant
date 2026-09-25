#ifndef RADIANT_FRAGMENTCOMPETITION_H
#define RADIANT_FRAGMENTCOMPETITION_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"

#include <QVector>

class CandidateScores;

class ALGORITHMSFFLIB_EXPORTS FragmentCompetition {
public:
    /**
     * Resolve coeluting isobaric assignments using unshared fragment evidence.
     *
     * Equal-charge candidates link at 5 ppm precursor mass, 20 ppm fragment
     * mass, four distinct shared supported fragments, trace cosine >= 0.5,
     * and apex separation <= half the narrower peak width. Each connected
     * group keeps the greatest unshared intensity * cosine^2 evidence.
     * Groups with no unshared evidence are rejected; singletons are retained.
     * Ties prefer higher LDA, then decoys, peptide, apex, mass, and input order.
     *
     * Preserves surviving pointers in input order and never edits their scores.
     * Unscored no-peak placeholders pass through to the normal NN input filter.
     * Null rows are discarded when enabled. Disabled is an exact no-op.
     */
    static Error::Err removeCompetingCandidates(
        bool enabled,
        QVector<CandidateScores*> *candidates);
};

#endif
