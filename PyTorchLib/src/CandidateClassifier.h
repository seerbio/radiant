//
// Created by anichols on 9/17/23.
//

#ifndef PYTHIADIACPP_CANDIDATECLASSIFIER_H
#define PYTHIADIACPP_CANDIDATECLASSIFIER_H

#include "PyTorchLib_Exports.h"

#include <QScopedPointer>
#include <QVector>


class PYTORCHLIB_EXPORTS CandidateClassifier {

public:

    CandidateClassifier();
    ~CandidateClassifier();

    /*!
    * @brief   Principal method for setting up and training the candidate classifier using input data.
    * @param   xData: The input data for training the classifier, each item in the outer vector corresponds to a single data point in the high-dimension data space.
    * @param   yData: The actual outcomes (or labels) corresponding to the inputs in "xData". Usually 1 indicates the desired outcome (e.g. presence of a certain feature), while 0 indicates the opposite.
    * @param   epochsMax: The maximum number of iterations for training the classifier.
    * @param   batchSize: The size of the subset of xData to use for each cycle of optimization.
    * @param   learningRate: The step-size for the optimizer during training.
    * @param   seed: The seed for generating random numbers, which is needed for initializing the weights in the classifier, for instance.
    * @return  bool: Returns true if the training process is successful and false otherwise. Unsuccessful training could result from invalid input (e.g. "xData" and "yData" not of the same length) or due to not achieving the desired accuracy.
    *
    * Note: This method uses the PyTorch library for building and training the neural network classifier.
    * Note2: Weights are initialized using kaiming_uniform_, (Tanh), for optimal gradient descent performance.
    */
    [[nodiscard]] bool trainCandidateClassifier(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData,
            int epochsMax,
            int batchSize,
            double learningRate,
            int seed,
            double nodeFraction,
            int verbosity
            ) const;

    /*!
    * @brief   Principal method for using the trained classifier to predict outcomes for given inputs.
    * @param   xData: The input data for making predictions, where each item in the outer vector corresponds to a single data point in the high-dimension data space.
    * @param   predictions: A pointer to the vector where the predicted outcomes will be stored.
    * @return  bool: Returns true if the prediction process is successful and false otherwise. An unsuccessful prediction could result from the classifier not being trained yet.
    *
    * Note: The caller must ensure that the neural network is trained (i.e., the "trainCandidateClassifier" method has been called) before making predictions, otherwise this method will return "false". An unsuccessful prediction could also result from "xData" not having the same structure (i.e., inner vector size) as the data it was trained with.
    */
    bool predict(
            const QVector<QVector<float>> &xData,
            QVector<float> *predictions
            ) const;

private:

    Q_DISABLE_COPY(CandidateClassifier) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_CANDIDATECLASSIFIER_H
