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

using CandidateScoresTarget = CandidateScores;
using CandidateScoresDecoy = CandidateScores;

class ALGORITHMSFFLIB_EXPORTS DiscriminantScoretron {

public:

    static Err trainLDAClassifier(
            const QList<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &targetDecoyCandidateScoresPair,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            QVector<float> *weights
            );

    static Err applyWeights(
        const QVector<float> &weights,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        QVector<CandidateScores*> *candidateScoresPntrs
    );

    static QVector<float> scoreVectorLogic(
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            CandidateScores* candidateScores
            );


};


#endif //PYTHIADIACPP_DISCRIMINANTSCORETRON_H
