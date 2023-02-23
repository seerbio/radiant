//
// Created by Drucifer on 11/7/2021.
//

#ifndef EIGENUTILS_H
#define EIGENUTILS_H

#include "EigenLib_Exports.h"
#include "MathUtils.h"

#include <QVector>

#include <Eigen/Dense>

#include <fstream>
#include <random>
#include <iostream>

class EIGENLIB_EXPORTS EigenUtils {

public:


    /*!
    * @brief Returns the index of a vector nearest to a given value
    */
    static QVector<int> returnIndexNearestToCutoff(const Eigen::RowVectorXd &vec,
                                                   double value,
                                                   int topX = 1) {

        topX = std::min(topX, static_cast<int>(vec.size()));
        const int sortLimit = std::min( (topX + 1), static_cast<int>(vec.size() - 1) );

        Eigen::RowVectorXd diff = (vec.array() - value).cwiseAbs();
        std::vector<double> diffVec(diff.data(), diff.data() + diff.size());

        QVector<QPair<int, double>> ogIndexPoints;
        for (int i = 0; i < diffVec.size(); i++) {
            ogIndexPoints.push_back({i, diffVec.at(i)});
        }

        const auto sortLogic = [](const QPair<int, double> &l, const QPair<int, double> &r){return l.second < r.second;};
        std::partial_sort(ogIndexPoints.begin(), ogIndexPoints.begin() + sortLimit, ogIndexPoints.end(), sortLogic);

        QVector<int> returnVec;
        for (int i = 0; i < sortLimit; i++) {
            returnVec.push_back(ogIndexPoints.at(i).first);
        }

        return returnVec;
    }


    /*!
    * @brief Returns the index of a vector nearest to a given value
    */
    template <typename T>
    static QVector<int> returnIndexNearestToCutoff(const QVector<T> &_container,
                                                   T value,
                                                   int topX = 1) {

        topX = std::min(topX, static_cast<int>(_container.size()));
        //TODO replace w/ EigenUtils::convertQvectorToEigenVector
        QVector<T> container = _container;
        Eigen::Map<Eigen::VectorX<T>> vec(container.data(), container.size());

        return returnIndexNearestToCutoff(vec.template cast<double>(), static_cast<double>(value), topX);
    }


    /*!
    * @brief Calculates root mean squared error given two Eigen vectors
    */
    static double calculateRMSE(const Eigen::RowVectorXd &vec1, const Eigen::RowVectorXd &vec2) {
        Eigen::RowVectorXd vec = (vec1 - vec2);
        return std::sqrt(vec.cwiseProduct(vec).sum() / vec.cols());
    }

    /*!
    * @brief Calculates the mean of an Eigen vector
    */
    template <typename EigenVector>
    static double calculateMeanOfVector(const EigenVector &vec) {
        return vec.sum() / vec.size();
    }

    /*!
    * @brief Calculates the Standard Deviation of an Eigen Vector
    */
    template <typename EigenVector>
    static double calculateStDevOfVector(const EigenVector &vec) {
        double mean = calculateMeanOfVector(vec);
        Eigen::RowVectorXd diff = vec.array() - mean;
        return std::sqrt(diff.cwiseProduct(diff).sum() / vec.size());
    }

    /*!
    * @brief Translates the points of an eigen vector given a rollDistance.  Mimics
    * numpy's roll function.  Note: that this wraps translated values
    */
    static Eigen::RowVectorXd roll(const Eigen::RowVectorXd &vec, int rollDistance) {

        Eigen::RowVectorXd returnVec(vec.cols());
        if (rollDistance < 0) {
            Eigen::RowVectorXd front = vec.middleCols(0, std::abs(rollDistance));
            Eigen::RowVectorXd back
                    = vec.middleCols(std::abs(rollDistance), vec.cols() - std::abs(rollDistance));
            returnVec << back, front;
            return returnVec;
        }

        else if (rollDistance > 0) {
            Eigen::RowVectorXd front = vec.middleCols(0, vec.cols() - std::abs(rollDistance));
            Eigen::RowVectorXd back
                    = vec.middleCols(vec.cols() - std::abs(rollDistance), std::abs(rollDistance));
            returnVec << back, front;
            return returnVec;
        }

        return vec;
    }


    /*!
    * @brief Fits a polynomial to an eigen vector given an order.
    */
    static void fitPolynomialQRDecomposition(const Eigen::MatrixXd points, int order,
                                             QVector<double> *coeff) {

        Eigen::MatrixXd A(points.rows(), order + 1);
        Eigen::VectorXd yv_mapped = points.col(1);

        // create decomposition matrix
        for (int i = 0; i < points.rows(); i++) {
            double pointValue = points(i, 0);
            for (size_t j = 0; j < order + 1; j++) {
                A(i, j) = pow(pointValue, j);
            }
        }

        // solve for linear least squares fit
        const Eigen::VectorXd result = A.householderQr().solve(yv_mapped);

        coeff->resize(order + 1);
        for (int i = 0; i <= order; ++i) {
            (*coeff)[i] = result[i];
        }

    }

    /*!
    * @brief Normalizes vector by dividing all elements by max element
    */
    template <typename EigenMatrix, typename T>
    static void normalizeVector(EigenMatrix *vec) {
        const T normalizerDemonator
            = MathUtils::tZero<T>(vec->maxCoeff()) ? std::numeric_limits<T>::min() : vec->maxCoeff();
        *vec /= normalizerDemonator;
    }

    /*!
    * @brief Adds noise to vector using eigen's randomizing methods.
    */
    template <typename EigenMatrix>
    static void addNoiseToVector(double scaleFactor, EigenMatrix *vec) {
        EigenMatrix noiseVector = EigenMatrix::Random(vec->size());
        noiseVector *= scaleFactor;

        *vec += noiseVector;
    }

    /*!
    * @brief  turns an array of integers into a one hot vector matrix where each column
     * corresponds to the respective position of each element in the input vec.
     *
     * Args:
     *  vec: a vector of integer values greater than or equal to 0.
     *  elementLength: The maximimum length of each one hot column vector in the matrix
    */
    static Eigen::MatrixXd vectorToOneHotMatrix(const QVector<int> &vec, int elementLength) {

        Eigen::MatrixXd mat(elementLength, vec.size());
        mat.setZero();

        for (int i =0; i < vec.size(); i++) {
            const int vecIndex = vec.at(i);
            mat.coeffRef(vecIndex, i) = 1.0;
        }

        return mat;
    }


    template <typename T>
    static void thresholdMatrix(double thresholdValue,  Eigen::MatrixX<T> *mat) {
        *mat = (mat->array() < thresholdValue).select(0.0, *mat);
    }


    template <typename T>
    static void thresholdVector(double thresholdValue,  Eigen::VectorX<T> *vec) {
        *vec = (vec->array() < thresholdValue).select(0.0, *vec);
    }


    template <typename EigenMatrix>
    static double cosineSimilarity(const EigenMatrix &v1,  const EigenMatrix &v2) {

        const double nearZero = 1e-5;
        const double dotProduct = v1.dot(v2);
        const double v1Magnitude = std::max(std::sqrt(v1.array().square().sum()), nearZero);
        const double v2Magnitude = std::max(std::sqrt(v2.array().square().sum()), nearZero);

        return MathUtils::pRound(dotProduct / (v1Magnitude * v2Magnitude), 4);
    }


    template <typename EigenMatrix>
    static double klDivergence(EigenMatrix v1, EigenMatrix v2) {

        const double v1Sum = v1.sum();
        const double v2Sum = v2.sum();

        if (MathUtils::tZero(v1Sum) || MathUtils::tZero(v2Sum)) {;
            return 1e4;
        }

        v1 /= v1Sum;
        v2 /= v2Sum;

        v1 = (v1.array() < 1e-5).select(1e-5, v1);
        v2 = (v2.array() < 1e-5).select(1e-5, v2);

        return MathUtils::pRound((v1.array() * Eigen::log2(v1.array() / v2.array())).sum(), 4);
    }


    template <typename EigenMatrix>
    static int nonZeros(EigenMatrix v1) {

        int nonZeros = 0;
        for (int i = 0; i < v1.size(); i++) {

            if (!MathUtils::tZero(v1.coeff(i))) {
                nonZeros++;
            }

        }

        return nonZeros;
    }


    template <typename T>
    static QVector<QPair<int, T>> returnTopXIndexAndValues(const Eigen::VectorX<T> &vec, int topX) {

        topX = std::min(topX, static_cast<int>(vec.size()));

        std::vector<double> stdVec(vec.data(), vec.data() + vec.size());

        QVector<QPair<int, double>> ogIndexPoints;
        for (int i = 0; i < stdVec.size(); i++) {
            ogIndexPoints.push_back({i, stdVec.at(i)});
        }

        const auto sortLogic = [](const QPair<int, double> &l, const QPair<int, double> &r){return l.second < r.second;};
        std::sort(ogIndexPoints.rbegin(), ogIndexPoints.rend(), sortLogic);

        ogIndexPoints.resize(topX);

        return ogIndexPoints;
    }


    template <typename T>
    static Eigen::VectorX<T> convertQVectorToEigenVector(const QVector<T> &_vec) {

        std::vector<T> vec = _vec.toStdVector();
        Eigen::VectorX<T> ev
                = Eigen::Map<Eigen::VectorX<T>, Eigen::Unaligned>(vec.data(), vec.size());

        return ev;
    }


    template <typename T>
    static QVector<T> convertEigenVectorToQVector(const Eigen::VectorX<T> &vec) {
        std::vector<T> vecReturn(vec.data(), vec.data() + vec.size());
        return QVector<T>::fromStdVector(vecReturn);
    }

    /*!
    * \brief Finds all apexes in a eigen vector.
     *
     *  Returns a QMap where the key is the index of the apex and the value is the value of
     *  the apex.
    */
    template <typename T>
    static QMap<int, T> apexes(
            const Eigen::VectorX<T> &vec,
            int precision = 1e4
    ){

        QMap<int, T> apexIndicies;

        const int vecSize = static_cast<int>(vec.size());
        for (int index = 0; index < vecSize; index++) {

            const auto centerPointValue = static_cast<int>(std::round(vec.coeff(index) * precision));

            if(index < 1){

                const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

                if (centerPointValue >= rightPointValue) {
                    apexIndicies.insert(index,  vec.coeff(index));
                }

                continue;
            }

            if (index >= vecSize - 1) {

                const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));

                if (centerPointValue > leftPointValue) {
                    apexIndicies.insert(index,  vec.coeff(index));
                }

                continue;
            }

            const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));
            const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

            if(centerPointValue > leftPointValue && centerPointValue >= rightPointValue){
                apexIndicies.insert(index,  vec.coeff(index));
            }
        }

        return apexIndicies;
    }


    /*!
    * \brief Finds all troughs in a eigen vector.
     *
     *  Returns a QMap where the key is the index of the apex and the value is the value of
     *  the apex.
    */
    template <typename T>
    static QMap<int, T> troughs(
            const Eigen::VectorX<T> &vec,
            int precision = 1e4
    ){

        QMap<int, T> troughtIndicies;

        const int vecSize = static_cast<int>(vec.size());
        for (int index = 0; index < vecSize; index++) {

            const auto centerPointValue = static_cast<int>(std::round(vec.coeff(index) * precision));

            if(index < 1){

                const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

                if (centerPointValue <= rightPointValue) {
                    troughtIndicies.insert(index,  vec.coeff(index));
                }

                continue;
            }

            if (index >= vecSize - 1) {

                const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));

                if (centerPointValue < leftPointValue) {
                    troughtIndicies.insert(index,  vec.coeff(index));
                }

                continue;
            }

            const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));
            const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

            if(centerPointValue < leftPointValue && centerPointValue <= rightPointValue){
                troughtIndicies.insert(index,  vec.coeff(index));
            }
        }

        return troughtIndicies;
    }


};


#endif //EIGENUTILS_H
