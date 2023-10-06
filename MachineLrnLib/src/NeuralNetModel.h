//
// Created by Drucifer on 3/6/2022.
//

#ifndef PYTHIACPP_NEURALNETMODEL_H
#define PYTHIACPP_NEURALNETMODEL_H

#include "Error.h"
#include "MachineLrnLib_Exports.h"

#include "Eigen/Dense"


using namespace Error;


class MACHINELRNLIB_EXPORTS NeuralNetModel {

public:

    NeuralNetModel();
    ~NeuralNetModel();

    Err init(const QString &modelFilePath);

    QVector<QVector<float>> batchPredict(const QVector<Eigen::MatrixXd> &mats);

    QVector<float> predict(const Eigen::VectorX<double> &vec);

private:

    Q_DISABLE_COPY(NeuralNetModel) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_NEURALNETMODEL_H
