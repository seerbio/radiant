//
// Created by anichols on 1/9/24.
//

#ifndef PYTHIADIACPP_QVALUESETTERTRON_H
#define PYTHIADIACPP_QVALUESETTERTRON_H

#include "Error.h"

#include "AlgorithmsFFLib_Exports.h"

using namespace Error;

class CandidateScores;
class CandidateScoresDDA;
using CandidateScoresTarget = CandidateScores;
using CandidateScoresDecoy = CandidateScores;

class ALGORITHMSFFLIB_EXPORTS QValueSettertron {

public:

    enum class QValueScoreType {
        DiscriminantScore,
        NNClassifierScore
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
            const QValueScoreType &qValueScoreType,
            QVector<CandidateScores*> *candidateScores
    );

    static Err setQValueForCandidates(
        const QValueScoreType &qValueScoreType,
        QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairsPntrs
    );

	static Err setQValueForCandidates(
		const QValueScoreType &qValueScoreType,
		QVector<QPair<CandidateScoresDDA*, CandidateScoresDDA*>> *targetDecoyCandidateScorePairsPntrs
	);

};


#endif //PYTHIADIACPP_QVALUESETTERTRON_H
