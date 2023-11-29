//
// Created by anichols on 11/13/23.
//

#ifndef PYTHIADIACPP_KARNNNEURALNET_H
#define PYTHIADIACPP_KARNNNEURALNET_H

#include "Error.h"
#include "MachineLrnLib_Exports.h"

using namespace Error;

class NN;

struct KarnnNNTarget {
    QString seq;
    float nnScore = 0.0;
    bool isDecoy = false;
    QVector<double> scoreVec;
    int index = -1;
};

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
            int epochs,
            int seed,
            QVector<KarnnNNTarget> *karnnNNTargets
            );

private:

    Err initNeuralNets(
            int rows,
            int cols,
            int epochs,
            int seed,
            std::vector<float*> &dataVecPointers,
            std::vector<float*> &labelVec,
            QVector<NN> *nets
            );

    Err calculatedNNScores(
            int cols,
            QVector<NN> *net,
            QVector<KarnnNNTarget> *karnnNNTargets
            );

private:

    float m_regularization;

    int m_bagginCount;
    int m_hiddenLayerCount;
    float m_learningRate;


};


#endif //PYTHIADIACPP_KARNNNEURALNET_H
