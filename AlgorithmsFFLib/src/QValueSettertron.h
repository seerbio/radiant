//
// Created by anichols on 1/9/24.
//

#ifndef PYTHIADIACPP_QVALUESETTERTRON_H
#define PYTHIADIACPP_QVALUESETTERTRON_H

#include "Error.h"

#include "AlgorithmsFFLib_Exports.h"

using namespace Error;

class CandidateScores;


class ALGORITHMSFFLIB_EXPORTS QValueSettertron {

public:

    enum class QValueScoreType {
        DiscriminantScore,
        NNClassifierScore
    };

    static Err setQValueForCandidates(
            const QValueScoreType &qValueScoreType,
            QVector<CandidateScores*> *candidateScores
    );

};


#endif //PYTHIADIACPP_QVALUESETTERTRON_H
