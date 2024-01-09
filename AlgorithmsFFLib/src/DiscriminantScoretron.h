//
// Created by anichols on 1/9/24.
//

#ifndef PYTHIADIACPP_DISCRIMINANTSCORETRON_H
#define PYTHIADIACPP_DISCRIMINANTSCORETRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;

class CandidateScores;

class ALGORITHMSFFLIB_EXPORTS DiscriminantScoretron {

public:

    static Err setDiscriminantScoreForCandidates(
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            QVector<CandidateScores*> *candidateScoresPntrs
    );


};


#endif //PYTHIADIACPP_DISCRIMINANTSCORETRON_H
