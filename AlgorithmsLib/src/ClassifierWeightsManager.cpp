//
// Created by anichols on 9/2/23.
//

#include "ClassifierWeightsManager.h"

#include "ErrorUtils.h"
#include "EigenUtils.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <iostream>


const QVector<double> &ClassifierWeightsManager::getWeights() const {
    return m_weights;
}

void ClassifierWeightsManager::setWeights(const QVector<double> &weights) {
    m_weights = weights;
}

const QVector<double> &ClassifierWeightsManager::getWeightsBest() const {
    return m_weightsBest;
}

void ClassifierWeightsManager::setWeightsBest(const QVector<double> &weightsBest) {
    m_weightsBest = weightsBest;
}

const QVector<double> &ClassifierWeightsManager::getGuideWeights() const {
    return m_guideWeights;
}

void ClassifierWeightsManager::setGuideWeights(const QVector<double> &guideWeights) {
    m_guideWeights = guideWeights;
}

const QVector<double> &ClassifierWeightsManager::getGuideWeightsBest() const {
    return m_guideWeightsBest;
}

void ClassifierWeightsManager::setGuideWeightsBest(const QVector<double> &guideWeightsBest) {
    m_guideWeightsBest = guideWeightsBest;
}

const QVector<double> &ClassifierWeightsManager::getSelectionWeights() const {
    return m_selectionWeights;
}

void ClassifierWeightsManager::setSelectionWeights(const QVector<double> &selectionWeights) {
    m_selectionWeights = selectionWeights;
}

const QVector<double> &ClassifierWeightsManager::getSelectionWeightsBest() const {
    return m_selectionWeightsBest;
}

void ClassifierWeightsManager::setSelectionWeightsBest(const QVector<double> &selectionWeightsBest) {
    m_selectionWeightsBest = selectionWeightsBest;
}

namespace {

    bool allMatAInputsAreSameSize(const QVector<QVector<double>> &matA) {

        const int matSize = matA.front().size();
        const auto allOfLogic = [matSize](const QVector<double> &v){return v.size() == matSize;};

        return std::all_of(matA.begin(), matA.end(), allOfLogic);
    }


}//namespace
Err ClassifierWeightsManager::fitWeights(
        const QVector<QVector<double>> &matA,
        const QVector<double> &vecB,
        QVector<double> *weights
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(matA); ree;
    e = ErrorUtils::isNotEmpty(vecB); ree;
    e = ErrorUtils::isTrue(allMatAInputsAreSameSize(matA)); ree;
    e = ErrorUtils::isEqual(matA.front().size(), vecB.size()); ree;

    const Eigen::MatrixX<double> A = EigenUtils::convertQVectorsToEigenMatrix(matA);
    const Eigen::VectorX<double> b = EigenUtils::convertQVectorToEigenVector(vecB);

    Eigen::VectorXd x = A.fullPivHouseholderQr().solve(b);
    *weights = EigenUtils::convertEigenVectorToQVector(x);

    ERR_RETURN
}

Err ClassifierWeightsManager::buildDataClassifier1(
        const QVector<QVector<double>> &targets,
        const QVector<QVector<double>> &decoys,
        QVector<QVector<double>> *A,
        QVector<double> *b
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targets); ree;
    e = ErrorUtils::isEqual(targets.size(), decoys.size()); ree;
    e = ErrorUtils::isEqual(targets.front().size(), decoys.front().size()); ree;

    const Eigen::MatrixX<double> matTargets = EigenUtils::convertQVectorsToEigenMatrix(targets);
    const Eigen::MatrixX<double> matDecoys = EigenUtils::convertQVectorsToEigenMatrix(decoys);

    const Eigen::MatrixX<double> subtractionMat = matTargets - matDecoys;
    const Eigen::VectorX<double> subtractionMatSum = subtractionMat.colwise().sum();
    const Eigen::VectorX<double> meanMat = subtractionMatSum / matTargets.rows();

    *b = EigenUtils::convertEigenVectorToQVector(meanMat);

    const int rows = targets.size();
    const int cols = targets.front().size();

    Eigen::MatrixX<double> matA(cols, cols);
    matA.setZero();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = j; k < cols; k++) {
                matA.coeffRef(j, k) += (targets[i][j] - decoys[i][j] - b->at(j)) * (targets[i][k] - decoys[i][k] - b->at(k));
                matA.coeffRef(k, j) = matA.coeffRef(j, k);
            }
        }
    }

    matA /= rows - 1;
    for (int i = 0; i < cols; i++) {
        matA.coeffRef(i, i) += std::numeric_limits<double>::min();
    }

    *A = EigenUtils::convertEigenMatrixToQVectors(matA);

    ERR_RETURN
}

Err ClassifierWeightsManager::buildDataClassifier2(
        const QVector<QVector<double>> &targets,
        const QVector<QVector<double>> &decoys,
        QVector<QVector<double>> *A,
        QVector<double> *b
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targets); ree;
    e = ErrorUtils::isEqual(targets.size(), decoys.size()); ree;
    e = ErrorUtils::isEqual(targets.front().size(), decoys.front().size()); ree;

    const Eigen::MatrixX<double> matTargets = EigenUtils::convertQVectorsToEigenMatrix(targets);
    const Eigen::MatrixX<double> matDecoys = EigenUtils::convertQVectorsToEigenMatrix(decoys);

    const Eigen::VectorX<double> matTargetsSum = matTargets.colwise().sum();
    const Eigen::VectorX<double> matDecoysSum = matDecoys.colwise().sum();

    const Eigen::VectorX<double> matTargetsSumMean = matTargetsSum / matTargets.rows();
    const Eigen::VectorX<double> matDecoysSumMean = matDecoysSum / matTargets.rows();

    const Eigen::MatrixX<double> subtractionMat = matTargets - matDecoys;
    const Eigen::VectorX<double> subtractionMatSum = subtractionMat.colwise().sum();
    const Eigen::VectorX<double> meanMat = subtractionMatSum / matTargets.rows();

    *b = EigenUtils::convertEigenVectorToQVector(meanMat);

    const int rows = targets.size();
    const int cols = targets.front().size();

    Eigen::MatrixX<double> matA(cols, cols);
    matA.setZero();

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = j; k < cols; k++) {
                matA.coeffRef(j, k) += 0.5 *
                                       ((decoys[i][j] - matDecoysSumMean[j]) * (decoys[i][k] - matDecoysSumMean[k]) +
                                        (targets[i][j] - matTargetsSumMean[j]) * (targets[i][k] - matTargetsSumMean[k]));
                matA.coeffRef(k, j) = matA.coeffRef(j, k);
            }
        }
    }

    matA /= rows - 1;
    for (int i = 0; i < cols; i++) {
        matA.coeffRef(i, i) += std::numeric_limits<double>::min();
    }

    *A = EigenUtils::convertEigenMatrixToQVectors(matA);

    ERR_RETURN
}

Err ClassifierWeightsManager::applyWeights(
        const QVector<QVector<double>> &matA,
        const QVector<double> &weights,
        QVector<double> *results
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(matA); ree;
    e = ErrorUtils::isNotEmpty(matA.front()); ree;
    e = ErrorUtils::isEqual(matA.front().size(), weights.size()); ree;

    const Eigen::MatrixX<double> A = EigenUtils::convertQVectorsToEigenMatrix(matA);
    const Eigen::VectorX<double> x = EigenUtils::convertQVectorToEigenVector(weights);

    const Eigen::VectorX<double> b = A * x;

    *results = EigenUtils::convertEigenVectorToQVector(b);

    ERR_RETURN
}
