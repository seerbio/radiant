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

    static Err removeInterferingCandidates(
            int ionsSharedToReject,
            double mzMinMS2,
            double mzMaxMS2,
            QVector<CandidateScores*> *candidates
    );


};


#endif //PYTHIADIACPP_INTERFERINGCANDIDATESELIMINATOMATIC_H
