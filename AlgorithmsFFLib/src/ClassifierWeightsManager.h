//
// Created by anichols on 9/2/23.
//

#ifndef PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H
#define PYTHIADIACPP_CLASSIFIERWEIGHTSMANAGER_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"

using ScoresTargets = QVector<float>;
using ScoresDecoys = QVector<float>;
using ScoresTargetsAndDecoys = QVector<float>;

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS ClassifierWeightsManager {

public:

    ClassifierWeightsManager() = default;
    ~ClassifierWeightsManager() = default;

    /**
    * @brief Builds a data classifier using targets and decoys for classification.
    *
    * This function takes two sets of vectors (targets and decoys) and constructs a data classifier
    * based on a specific algorithm. The resulting classifier is represented by matrices A and b.
    *
    * @param targets A QVector of QVector<float>* representing the target vectors.
    * @param decoys A QVector of QVector<float>* representing the decoy vectors.
    * @param A Output parameter, a QVector of QVector<float>*, representing the matrix A of the classifier.
    * @param b Output parameter, a QVector<float>* representing the vector b of the classifier.
    * @return An Err enum indicating the success or failure of the operation.
    */
    static Err buildDataClassifier1(
        const QVector<QVector<float>*> &targets,
        const QVector<QVector<float>*> &decoys,
        QVector<QVector<float>> *A,
        QVector<float> *b
        );

    /**
    * @brief Builds a data classifier using targets and decoys for classification.
    *
    * This function implements a different algorithm compared to buildDataClassifier1.
    * It takes two sets of vectors (targets and decoys) and constructs a data classifier
    * based on a specific algorithm. The resulting classifier is represented by matrices A and b.
    *
    * @param targets A QVector of QVector<float>* representing the target vectors.
    * @param decoys A QVector of QVector<float>* representing the decoy vectors.
    * @param A Output parameter, a QVector of QVector<float>*, representing the matrix A of the classifier.
    * @param b Output parameter, a QVector<float>* representing the vector b of the classifier.
    * @return An Err enum indicating the success or failure of the operation.
    */
    static Err buildDataClassifier2(
        const QVector<QVector<float>*> &targets,
        const QVector<QVector<float>*> &decoys,
        QVector<QVector<float>> *A,
        QVector<float> *b
    );

    /**
    * @brief Fits weights for the data classifier using regularization w/ fullPivHouseholderQr.
    *
    * This function takes matrices A and b representing the data classifier and fits weights using
    * regularized fullPivHouseholderQr. The resulting weights are stored in the 'weights' vector.
    *
    * @param matA A QVector of QVector<float> representing the matrix A of the data classifier.
    * @param vecB A QVector<float> representing the vector b of the data classifier.
    * @param weights Output parameter, a QVector<float>* representing the fitted weights.
    * @return An Err enum indicating the success or failure of the operation.
    */
    static Err fitWeights(
        const QVector<QVector<float>> &matA,
        const QVector<float> &vecB,
        QVector<float> *weights
        );

    /**
    * @brief Applies weights to input data for classification.
    *
    * This function takes input data represented by matrices A and weights,
    * and applies the weights to produce classification results.
    *
    * @param matA Input data represented as a QVector of QVector<float>*.
    * @param weights Weights to be applied for classification.
    * @param results Output parameter, a QVector<float>* storing the classification results.
    * @return An Err enum indicating the success or failure of the operation.
    */
    static Err applyWeights(
        const QVector<QVector<float>*> &matA,
        const QVector<float> &weights,
        QVector<float> *results
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
