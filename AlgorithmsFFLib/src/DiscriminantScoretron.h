//
// Created by anichols on 1/9/24.
//

#ifndef PYTHIADIACPP_DISCRIMINANTSCORETRON_H
#define PYTHIADIACPP_DISCRIMINANTSCORETRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "CandidateScores.h"
#include "CandidateScoresDDA.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;

class ALGORITHMSFFLIB_EXPORTS DiscriminantScoretron {

public:

    static QVector<Features> featuresCalibration();

    static QVector<FeaturesDDA> featuresCalibrationDDA();

    static QVector<Features> featuresOptimization();

    static QVector<Features> featuresNeuralNetwork();

    static QVector<FeaturesDDA> featuresNeuralNetworkDDA();

    static Err trainLDAClassifier(
            const QVector<QPair<FeaturesArrayDecoys*, FeaturesArrayDecoys*>> &targetDecoyCandidateScoresPair,
            int verbosity,
            QVector<float> *weights
            );

    static Err applyWeights(
        const QVector<float> &weights,
        int threadCount,
        const QVector<FeaturesArray*> &featureArrayPntrs,
        QVector<float> *discriminantScores
    );

    static QVector<float> defaultWeights(const QVector<Features> &features);
    static QVector<float> defaultWeights(const QVector<FeaturesDDA> &features);

    static Err convertScoreCandidatesToFeaturesArrays(
        const QVector<Features> &features,
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &candidateScoresTargetVsDecoy,
        QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *featuresArrayTargetVsDecoy
        );

};


#endif //PYTHIADIACPP_DISCRIMINANTSCORETRON_H
