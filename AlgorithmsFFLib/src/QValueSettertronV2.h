//
// Created by anichols on 1/9/24.
//

#ifndef PYTHIADIACPP_QVALUESETTERTRON_H
#define PYTHIADIACPP_QVALUESETTERTRON_H

#include "Error.h"

#include "AlgorithmsFFLib_Exports.h"

using namespace Error;

class CandidateScoresV2;
using CandidateScoresV2Target = CandidateScoresV2;
using CandidateScoresV2Decoy = CandidateScoresV2;

class ALGORITHMSFFLIB_EXPORTS QValueSettertronV2 {

public:

    enum class QValueScoreTypeV2 {
        DiscriminantScore,
        NNClassifierScore,
    	CosineSimSumMeanCorrelation
    };

    /**
    * @brief Sets q-values for candidate scores based on specified q-value score type.
    *
    * This function calculates and sets q-values for candidate scores based on the specified q-value score type.
    * It takes a QVector of CandidateScores and updates their q-values accordingly.
    *
    * @param qValueScoreType The q-value score type to use for calculation.
    * @param candidateScores A pointer to the QVector of CandidateScores to set q-values for.
    * @return An error code indicating the success or failure of the q-value setting process.
    */
    static Err setQValueForCandidates(
            const QValueScoreTypeV2 &qValueScoreType,
            QVector<CandidateScoresV2*> *candidateScores
    );

    static Err setQValueForCandidates(
        const QValueScoreTypeV2 &qValueScoreType,
        QVector<QPair<CandidateScoresV2Target*, CandidateScoresV2Decoy*>> *targetDecoyCandidateScorePairsPntrs
    );

};


#endif //PYTHIADIACPP_QVALUESETTERTRON_H
