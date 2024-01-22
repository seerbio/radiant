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

    /**
    * @brief Sets discriminant scores for candidate pairs using parallel computation.
    *
    * This function calculates discriminant scores for target and decoy candidate pairs based on various inputs,
    * including extended scores and neural network scores. It utilizes parallel computation for efficiency.
    *
    * @param useExtendedScores A boolean indicating whether extended scores are used in the calculation.
    * @param useNeuralNetworkScores A boolean indicating whether neural network scores are used in the calculation.
    * @param candidateScoresPntrs A QVector of CandidateScores pointers representing candidate pairs.
    * @return An Err enum indicating the success or failure of the operation.
    */
    static Err setDiscriminantScoreForCandidates(
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            QVector<CandidateScores*> *candidateScoresPntrs
    );


};


#endif //PYTHIADIACPP_DISCRIMINANTSCORETRON_H
