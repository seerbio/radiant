#ifndef POSTNEURALNETCOMPETITION_H
#define POSTNEURALNETCOMPETITION_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"

#include <QVector>

class CandidateScores;

class ALGORITHMSFFLIB_EXPORTS PostNeuralNetCompetition {
public:
    // Zero disables the entire operation. Otherwise retain one best candidate
    // per origin/reported precursor/decoy state, assign conservative q-values,
    // and emit exact probability order (decoys first for equal probabilities).
    static Error::Err apply(int minimumSharedFragments, QVector<CandidateScores*> *candidates);
};

#endif
