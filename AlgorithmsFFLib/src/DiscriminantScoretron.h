//
// Created by anichols on 1/9/24.
//

#ifndef PYTHIADIACPP_DISCRIMINANTSCORETRON_H
#define PYTHIADIACPP_DISCRIMINANTSCORETRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "CandidateScores.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;

class ALGORITHMSFFLIB_EXPORTS DiscriminantScoretron {

public:

    static Err trainLDAClassifier(
            const QVector<QPair<FeaturesArrayDecoys*, FeaturesArrayDecoys*>> &targetDecoyCandidateScoresPair,
            int threadCount,
            int verbosity,
            QVector<float> *weights
            );

    static Err applyWeights(
        const QVector<float> &weights,
        int threadCount,
        const QVector<FeaturesArray*> &candidateScoresPntrs,
        QVector<float> *discriminantScores
    );

    static QVector<float> scoreVectorLogic(
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            CandidateScores* candidateScores
            );

    static Err convertScoreCandidatesToFeaturesArrays(
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &candidateScoresTargetVsDecoy,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *featuresArrayTargetVsDecoy
        );

};


#endif //PYTHIADIACPP_DISCRIMINANTSCORETRON_H
