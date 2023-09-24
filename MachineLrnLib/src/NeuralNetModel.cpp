//
// Created by Drucifer on 3/6/2022.
//

#include "NeuralNetModel.h"

#include "ErrorUtils.h"
#include "ParallelUtils.h"

#include "fdeep/fdeep.hpp"

#include "QFileInfo"
#include <QtConcurrent/QtConcurrent>


class Q_DECL_HIDDEN NeuralNetModel::Private
{
public:

    Private();
    ~Private();

    Err init(const QString &modelFilePath);

    QVector<QVector<float>> batchPredict(const QVector<Eigen::MatrixXd> &mats);

    QVector<float> predict(const Eigen::VectorX<double> &vec);


private:

    QString m_modelFilePath;

    fdeep::model *m_model;

};

NeuralNetModel::Private::Private()
: m_modelFilePath("")
, m_model(nullptr)
{}

NeuralNetModel::Private::~Private() {
    delete m_model;
}


Err NeuralNetModel::Private::init(const QString &modelFilePath) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(modelFilePath); ree;

    const QFileInfo fi(modelFilePath);
    e = ErrorUtils::isTrue(fi.exists() && fi.isFile()); ree;

    m_modelFilePath = modelFilePath;

    m_model = new fdeep::model(fdeep::load_model(m_modelFilePath.toStdString(), true, fdeep::dev_null_logger));

    ERR_RETURN
}

namespace {

    fdeep::tensor eigenMatrixToFDeepTensor(const Eigen::MatrixXd &mat) {

        const int channels = 0;
        const int rows = static_cast<int>(mat.rows());
        const int cols = static_cast<int>(mat.cols());

        fdeep::tensor_shape tensorShape(rows, cols);
        fdeep::tensor matrixTensor(tensorShape, 0.0f);

        for (int i = 0; i < rows; ++i) {

            for (int j = 0; j < cols; ++j) {
                matrixTensor.set(fdeep::tensor_pos(i, j), static_cast<float>(mat.coeff(i, j)));
            }
        }

        return matrixTensor;
    }

    fdeep::tensor eigenVecToFDeepTensor(const Eigen::VectorX<double> vec) {

        const int channels = 0;
        const int rows = static_cast<int>(vec.rows());

        fdeep::tensor_shape tensorShape(rows);
        fdeep::tensor matrixTensor(tensorShape, 0.0f);

        for (int i = 0; i < rows; ++i) {
            matrixTensor.set(fdeep::tensor_pos(i), static_cast<float>(vec.coeff(i)));
        }

        return matrixTensor;
    }

    struct ParallelInput {
        std::string modelFilePath;
        QVector<Eigen::MatrixXd> matsToPredict;
    };

    QVector<QVector<float>> parallelLogic(const ParallelInput &parallelInput) {

        std::vector<fdeep::tensors> tensorsToPredict;
        for (const Eigen::MatrixXd &mat : parallelInput.matsToPredict) {
            const fdeep::tensor t = eigenMatrixToFDeepTensor(mat);
            tensorsToPredict.push_back({t});
        }


//#define VERBOSE_OUTPUT
#ifdef VERBOSE_OUTPUT
        fdeep::model model = fdeep::load_model(parallelInput.modelFilePath, true);
#else
        fdeep::model model = fdeep::load_model(parallelInput.modelFilePath, true, fdeep::dev_null_logger);
#endif

        const std::vector<fdeep::tensors> predictions = model.predict_multi(tensorsToPredict, false);

        QVector<QVector<float>> returnPredictions;
        for (const auto &pred : predictions){
            const QVector<float> vec = QVector<float>::fromStdVector(pred.front().to_vector());
            returnPredictions.push_back(vec);
        }

        return returnPredictions;
    }


}// NAMESPACE
QVector<QVector<float>> NeuralNetModel::Private::batchPredict(const QVector<Eigen::MatrixXd> &mats) {

    const int trancheBuffer = 0;

    QVector<QVector<Eigen::MatrixXd>> matsTranched;
    ParallelUtils::trancheVectorForParallelizationInOrder(
            mats,
            ParallelUtils::numberOfAvailableSystemProcessors(),
            trancheBuffer,
            &matsTranched
            );

    QVector<QFuture<QVector<QVector<float>>>> futures;
    for(const QVector<Eigen::MatrixXd> &matTranch : matsTranched){

        ParallelInput pi;
        pi.modelFilePath = m_modelFilePath.toStdString();
        pi.matsToPredict = matTranch;

        QFuture<QVector<QVector<float>>> thread = QtConcurrent::run(parallelLogic, pi);
        futures.push_back(thread);
    }

    for(QFuture<QVector<QVector<float>>> &future : futures){
        future.waitForFinished();
    }

    QVector<QVector<float>> results;
    for(const QFuture<QVector<QVector<float>>> &future : futures){
        results.append(future.result());
    }

    return results;
}

QVector<float> NeuralNetModel::Private::predict(const Eigen::VectorX<double> &vec) {

    const fdeep::tensor t = eigenVecToFDeepTensor(vec);

    const fdeep::tensors pred = m_model->predict({t});

    const QVector<float> vecRet = QVector<float>::fromStdVector(pred.front().to_vector());

    return vecRet;
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

NeuralNetModel::NeuralNetModel()  : d_ptr(new Private()) {}

NeuralNetModel::~NeuralNetModel() {}

Err NeuralNetModel::init(const QString &modelFilePath) {
    ERR_INIT
    e = d_ptr->init(modelFilePath); ree;
    ERR_RETURN
}

// TODO change this to handle errors.
QVector<QVector<float>> NeuralNetModel::batchPredict(const QVector<Eigen::MatrixXd> &mats) {
    return d_ptr->batchPredict(mats);
}

QVector<float> NeuralNetModel::predict(const Eigen::VectorX<double> &mat) {
    return d_ptr->predict(mat);
}
