//
// Created by anichols on 9/2/23.
//

#include "ClassifierWeightsManager.h"

#include "ErrorUtils.h"
#include "EigenUtils.h"
#include "ParallelUtils.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"


Err ClassifierWeightsManager::buildDataClassifier1(
        const QVector<QVector<float>*> &targets,
        const QVector<QVector<float>*> &decoys,
        QVector<QVector<float>> *A,
        QVector<float> *b
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targets); ree;
    e = ErrorUtils::isEqual(targets.size(), decoys.size()); ree;
    e = ErrorUtils::isEqual(targets.front()->size(), decoys.front()->size()); ree;

    const Eigen::MatrixX<float> matTargets = EigenUtils::convertQVectorsToEigenMatrix(targets);
    const Eigen::MatrixX<float> matDecoys = EigenUtils::convertQVectorsToEigenMatrix(decoys);

    const Eigen::MatrixX<float> subtractionMat = matTargets - matDecoys;
    const Eigen::VectorX<float> subtractionMatSum = subtractionMat.colwise().sum();
    const Eigen::VectorX<float> meanMat = subtractionMatSum / matTargets.rows();

    *b = EigenUtils::convertEigenVectorToQVector(meanMat);

    const int rows = targets.size();
    const int cols = targets.front()->size();

    Eigen::MatrixX<float> matA(cols, cols);
    matA.setZero();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = j; k < cols; k++) {
                matA.coeffRef(j, k) += (targets[i]->at(j) - decoys[i]->at(j) - b->at(j)) * (targets[i]->at(k) - decoys[i]->at(k) - b->at(k));
                matA.coeffRef(k, j) = matA.coeffRef(j, k);
            }
        }
    }

    matA /= static_cast<float>(rows - 1);
    for (int i = 0; i < cols; i++) {
        matA.coeffRef(i, i) += std::numeric_limits<float>::min();
    }

    *A = EigenUtils::convertEigenMatrixToQVectors(matA);

    ERR_RETURN
}

Err ClassifierWeightsManager::buildDataClassifier2(
        const QVector<QVector<float>*> &targets,
        const QVector<QVector<float>*> &decoys,
        QVector<QVector<float>> *A,
        QVector<float> *b
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targets); ree;
    e = ErrorUtils::isEqual(targets.size(), decoys.size()); ree;
    e = ErrorUtils::isEqual(targets.front()->size(), decoys.front()->size()); ree;

    const Eigen::MatrixX<float> matTargets = EigenUtils::convertQVectorsToEigenMatrix(targets);
    const Eigen::MatrixX<float> matDecoys = EigenUtils::convertQVectorsToEigenMatrix(decoys);

    const Eigen::VectorX<float> matTargetsSum = matTargets.colwise().sum();
    const Eigen::VectorX<float> matDecoysSum = matDecoys.colwise().sum();

    const Eigen::VectorX<float> matTargetsSumMean = matTargetsSum / matTargets.rows();
    const Eigen::VectorX<float> matDecoysSumMean = matDecoysSum / matTargets.rows();

    const Eigen::MatrixX<float> subtractionMat = matTargets - matDecoys;
    const Eigen::VectorX<float> subtractionMatSum = subtractionMat.colwise().sum();
    const Eigen::VectorX<float> meanMat = subtractionMatSum / matTargets.rows();

    *b = EigenUtils::convertEigenVectorToQVector(meanMat);

    const int rows = targets.size();
    const int cols = targets.front()->size();

    Eigen::MatrixX<float> matA(cols, cols);
    matA.setZero();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = j; k < cols; k++) {
                matA.coeffRef(j, k) += 0.5 *
                                       ((decoys[i]->at(j) - matDecoysSumMean[j]) * (decoys[i]->at(k) - matDecoysSumMean[k]) +
                                        (targets[i]->at(j) - matTargetsSumMean[j]) * (targets[i]->at(k) - matTargetsSumMean[k]));
                matA.coeffRef(k, j) = matA.coeffRef(j, k);
            }
        }
    }

    matA /= rows - 1;
    for (int i = 0; i < cols; i++) {
        matA.coeffRef(i, i) += std::numeric_limits<float>::min();
    }

    *A = EigenUtils::convertEigenMatrixToQVectors(matA);

    ERR_RETURN
}

namespace {

    bool allMatAInputsAreSameSize(const QVector<QVector<float>> &matA) {

        const int matSize = matA.front().size();
        const auto allOfLogic = [matSize](const QVector<float> &v){return v.size() == matSize;};

        return std::all_of(matA.begin(), matA.end(), allOfLogic);
    }


}//namespace
Err ClassifierWeightsManager::fitWeights(
        const QVector<QVector<float>> &matA,
        const QVector<float> &vecB,
        QVector<float> *weights
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(matA); ree;
    e = ErrorUtils::isNotEmpty(vecB); ree;
    e = ErrorUtils::isTrue(allMatAInputsAreSameSize(matA)); ree;
    e = ErrorUtils::isEqual(matA.front().size(), vecB.size()); ree;

    const Eigen::MatrixX<float> A = EigenUtils::convertQVectorsToEigenMatrix(matA);
    const Eigen::VectorX<float> b = EigenUtils::convertQVectorToEigenVector(vecB);

    double lambda = 0.001; //TODO auto set this
    Eigen::MatrixX<float> Areg = A + lambda * Eigen::MatrixX<float>::Identity(A.rows(), A.cols());

    Eigen::VectorX<float> x = Areg.fullPivHouseholderQr().solve(b);
    *weights = EigenUtils::convertEigenVectorToQVector(x);

    ERR_RETURN
}

Err ClassifierWeightsManager::applyWeights(
        const QVector<QVector<float>*> &matA,
        const QVector<float> &weights,
        QVector<float> *results
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(matA); ree;
    e = ErrorUtils::isFalse(matA.front()->isEmpty()); ree;
    e = ErrorUtils::isEqual(matA.front()->size(), weights.size()); ree;

    const Eigen::MatrixX<float> A = EigenUtils::convertQVectorsToEigenMatrix(matA);
    const Eigen::VectorX<float> x = EigenUtils::convertQVectorToEigenVector(weights);

    const Eigen::VectorX<float> b = A * x;

    *results = EigenUtils::convertEigenVectorToQVector(b);

    ERR_RETURN
}
