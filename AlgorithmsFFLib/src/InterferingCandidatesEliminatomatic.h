//
// Created by anichols on 1/9/24.
//

#ifndef PYTHIADIACPP_INTERFERINGCANDIDATESELIMINATOMATIC_H
#define PYTHIADIACPP_INTERFERINGCANDIDATESELIMINATOMATIC_H

#include "Error.h"

#include "AlgorithmsFFLib_Exports.h"

using namespace Error;

class CandidateScores;

class ALGORITHMSFFLIB_EXPORTS InterferingCandidatesEliminatomatic {

public:

    /**
    * @brief Removes interfering candidates based on specified criteria.
    *
    * This function removes interfering candidates from the provided list based on criteria such as shared ions,
    * mz range, and discriminant scores.
    *
    * @param ionsSharedToReject An integer indicating the minimum number of shared ions to reject a candidate.
    * @param mzMinMS2 A double representing the minimum mz value in MS2.
    * @param mzMaxMS2 A double representing the maximum mz value in MS2.
    * @param candidates Input and output parameter, a QVector of CandidateScores* representing the candidate scores.
    * @return An Err enum indicating the success or failure of the operation.
    */
    static Err removeInterferingCandidates(
            int ionsSharedToReject,
            double mzMinMS2,
            double mzMaxMS2,
            QVector<CandidateScores*> *candidates
    );


};


#endif //PYTHIADIACPP_INTERFERINGCANDIDATESELIMINATOMATIC_H
