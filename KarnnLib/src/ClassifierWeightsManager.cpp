//
// Created by anichols on 9/2/23.
//

#include "ClassifierWeightsManager.h"

#include "ErrorUtils.h"

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <iostream>

namespace {

    bool allMatAInputsAreSameSize(const QVector<QVector<double>> &matA) {

        const int matSize = matA.front().size();
        const auto allOfLogic = [matSize](const QVector<double> &v){return v.size() == matSize;};

        return std::all_of(matA.begin(), matA.end(), allOfLogic);
    }

    Eigen::MatrixX<double> convertMatFromQVectors(const QVector<QVector<double>> &matA) {

        const int rows = matA.size();
        const int columns = matA.front().size();
        Eigen::MatrixX<double> mat(rows, columns);
        mat.setZero();

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                mat.coeffRef(i, j) = matA[i][j];
            }
        }

        return mat;
    }

    QVector<QVector<double>> convertQVectorsFromMat(const Eigen::MatrixX<double> &mat) {

        QVector<QVector<double>> vecs;
        for (int row = 0; row < mat.rows(); row++) {

            Eigen::VectorX<double> v = mat.row(row);
            const std::vector<double> vecReturn(v.data(), v.data() + v.size());
            vecs.push_back(QVector<double>::fromStdVector(vecReturn));
        }

        return vecs;
    }

}//namespace
Err ClassifierWeightsManager::fitWeights(
        const QVector<QVector<double>> &matA,
        const QVector<double> &vecB,
        QVector<double> *x
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(matA); ree;
    e = ErrorUtils::isNotEmpty(vecB); ree;
    e = ErrorUtils::isTrue(allMatAInputsAreSameSize(matA)); ree;
    e = ErrorUtils::isEqual(matA.front().size(), vecB.size()); ree;

    const Eigen::MatrixX<double> A = convertMatFromQVectors(matA);

    QVector<double> _vecB = vecB;
    const Eigen::VectorX<double> b = Eigen::Map<Eigen::VectorX<double>, Eigen::Unaligned>(_vecB.data(), _vecB.size());

    Eigen::VectorXd X = A.fullPivHouseholderQr().solve(b);
    const std::vector<double> vecReturn(X.data(), X.data() + X.size());
    *x = QVector<double>::fromStdVector(vecReturn);

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

    const Eigen::MatrixX<double> matTargets = convertMatFromQVectors(targets);
    const Eigen::MatrixX<double> matDecoys = convertMatFromQVectors(decoys);

    const Eigen::MatrixX<double> subtractionMat = matTargets - matDecoys;
    const Eigen::VectorX<double> subtractionMatSum = subtractionMat.colwise().sum();
    const Eigen::VectorX<double> meanMat = subtractionMatSum / matTargets.rows();

    const std::vector<double> vecReturn(meanMat.data(), meanMat.data() + meanMat.size());
    *b = QVector<double>::fromStdVector(vecReturn);

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

    *A = convertQVectorsFromMat(matA);

    ERR_RETURN
}
