//
// Created by anichols on 9/2/23.
//

#ifndef PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H
#define PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H

#include "KarnnLib_Exports.h"

#include "Error.h"

using namespace Error;

class KARNNLIB_EXPORTS ClassifierWeightsManager {

public:

    ClassifierWeightsManager() = default;
    ~ClassifierWeightsManager() = default;

    static Err buildDataClassifier1(
            const QVector<QVector<double>> &targets,
            const QVector<QVector<double>> &decoys,
            QVector<QVector<double>> *A,
            QVector<double> *b
            );

    static Err fitWeights(
            const QVector<QVector<double>> &matA,
            const QVector<double> &vecB,
            QVector<double> *x
            );

private:

    QVector<double> m_weights;
    QVector<double> m_weightsBest;
    QVector<double> m_guideWeights;
    QVector<double> m_guideWeightsBest;
    QVector<double> m_selectionWeights;
    QVector<double> m_selectionWeightsBest;

};


#endif //PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H
