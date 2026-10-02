#ifndef RADIANT_CANDIDATE_POOL_SELECTION_H
#define RADIANT_CANDIDATE_POOL_SELECTION_H

#include "FragmentCompetition.h"

class ALGORITHMSFFLIB_EXPORTS CandidatePoolSelection {
public:
    struct Identity {
        QString reportedPeptide, originPeptide, targetKey;
        int charge = 0, view = 0, viewRank = 0;
        bool isDecoy = false;
        double apex = 0.0;
    };

    struct Confidence {
        QVector<int> inputIndices;
        QVector<double> qValues;
    };

    // Round-robin by within-view rank; lower view number breaks rank ties.
    // Exact peak identity includes reported/origin peptide, charge, decoy
    // class, target key and apex. The cap applies after deduplication.
    static Error::Err rankedUnique(
        const QVector<Identity> &rows, int maximumCandidates, QVector<int> *indices);

    // Select the best probability per precursor and decoy class. Targets use
    // the reported modified sequence; decoys additionally use the origin
    // modified sequence, since different origins can share a synthetic display
    // sequence. Keeping the reported sequence also distinguishes library and
    // synthetic decoys of one origin. Charge is part of both identities.
    // Then compete using unchanged physical
    // evidence and recompute exact-score tied +1 q-values over the survivors.
    // Outputs change only on success.
    static Error::Err finalize(
        const QVector<FragmentCompetition::Evidence> &evidence,
        const QVector<QString> &originPeptides,
        const QVector<double> &decoyProbabilities, int minimumSharedFragments,
        Confidence *confidence);

    // Compatibility entry point for evidence whose reported sequences already
    // uniquely identify each target/decoy precursor. Select the best peak,
    // then compete using unchanged physical evidence. Recompute exact-score
    // tied +1 q-values over the survivors. Outputs change only on success.
    static Error::Err finalize(
        const QVector<FragmentCompetition::Evidence> &evidence,
        const QVector<double> &decoyProbabilities, int minimumSharedFragments,
        Confidence *confidence);
};

#endif
