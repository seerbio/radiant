//
// Created by anichols on 11/13/23.
//

#ifndef PYTHIADIACPP_KARNNNEURALNET_H
#define PYTHIADIACPP_KARNNNEURALNET_H

#include "Error.h"
#include "MachineLrnLib_Exports.h"

using namespace Error;

class NN;

class MACHINELRNLIB_EXPORTS KarnnNeuralNet {

public:

    KarnnNeuralNet();
    ~KarnnNeuralNet() = default;

    Err init(
            int baggingCount,
            int hiddenLayerCount,
            float learningRate
            );

    Err run(
            const QVector<QVector<double>> &trainingData,
            const QVector<bool> &labels,
            int epochs,
            std::vector<float> *nnScores
            );

private:


private:

    float m_regularization;
    int m_threads;

    int m_bagginCount;
    int m_hiddenLayerCount;
    float m_learningRate;


};


#endif //PYTHIADIACPP_KARNNNEURALNET_H
