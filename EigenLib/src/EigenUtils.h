//
// Created by Drucifer on 11/7/2021.
//

#ifndef EIGENUTILS_H
#define EIGENUTILS_H

#include "CSVReader.h"
#include "EigenLib_Exports.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "MathUtils.h"

#include <QFile>
#include <QVector>
#include <QMap>

#include <Eigen/Dense>

#include <fstream>
#include <random>
#include <iostream>

class EIGENLIB_EXPORTS EigenUtils {

public:

    /*!
     * @brief  Replaces all NaN values in the vector with a provided value
     * @tparam T: The datatype of the elements in the Eigen::Vector. Could be any numeric type such as int, float, double etc.
     * @param replaceVal: The value which will replace all the NaNs in the vector
     * @param vec: Pointer to an Eigen::Vector. This function modifies the vector in-place.
    */
    template<typename T>
    static void replaceNaN(T replaceVal, Eigen::VectorX<T> *vec){
        *vec = vec->array().isNaN().select(replaceVal, vec->array());
    }

    /*!
     * @brief  Replaces all NaN values in the matrix with a provided value
     * @tparam T: The datatype of the elements in the Eigen::Matrix. Could be any numeric type such as int, float, double etc.
     * @param replaceVal: The value which will replace all the NaNs in the matrix
     * @param mat: Pointer to an Eigen::Matrix. This function modifies the matrix in-place.
    */
    template<typename T>
    static void replaceNaN(T replaceVal, Eigen::MatrixX<T> *mat){
        *mat = mat->array().isNaN().select(replaceVal, mat->array());
    }

    /*!
    * @brief  Replaces all Infinity (inf) values in the matrix with a provided value
    * @tparam T: The datatype of the elements in the Eigen::Matrix. Could be any numeric type such as int, float, double etc.
    * @param replaceVal: The value which will replace all the infinity values in the matrix
    * @param mat: Pointer to an Eigen::Matrix. This function modifies the matrix in-place.
    */
    template<typename T>
    static void replaceInf(T replaceVal, Eigen::MatrixX<T> *mat){
        *mat = mat->unaryExpr([replaceVal](double v) { return std::isinf(v) ? replaceVal : v; });
    }

    /*!
    * @brief  Returns the indexes of 'topX' elements in 'vec' which are closest to 'value'
    * @tparam T: The datatype of the elements in the Eigen::RowVectorX. Could be any numeric type such as int, float, double etc.
    * @param vec: An Eigen::RowVectorX from which closest elements' indexes are returned
    * @param value: The value with which closeness (absolute difference) is evaluated
    * @param topX: Number of top elements' indexes to return which are closest to 'value'. Default value is 1
    * @return QVector<int>: QVector containing indexes of 'topX' closest elements to 'value'
    */
    template<typename T>
    static QVector<int> returnIndexNearestToCutoff(
            const Eigen::RowVectorX<T> &vec,
            T value,
            int topX = 1
                    ) {

        topX = std::min(topX, static_cast<int>(vec.size()));
        const int sortLimit = std::min(topX, static_cast<int>(vec.size() - 1) );

        Eigen::RowVectorX<T> diff = (vec.array() - value).cwiseAbs();
        std::vector<T> diffVec(diff.data(), diff.data() + diff.size());

        QVector<QPair<int, T>> ogIndexPoints;
        for (int i = 0; i < diffVec.size(); i++) {
            ogIndexPoints.push_back({i, diffVec.at(i)});
        }

        const auto sortLogic = [](const QPair<int, T> &l, const QPair<int, T> &r){return l.second < r.second;};
        std::partial_sort(ogIndexPoints.begin(), ogIndexPoints.begin() + sortLimit, ogIndexPoints.end(), sortLogic);

        QVector<int> returnVec;
        for (int i = 0; i < sortLimit; i++) {
            returnVec.push_back(ogIndexPoints.at(i).first);
        }

        return returnVec;
    }

    /*!
    * @brief  Calculates the Root Mean Squared Error (RMSE) between two vectors
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param vec1: First Eigen::VectorX
    * @param vec2: Second Eigen::VectorX
    * @param rmseOutput: Pointer to a double where the result (RMSE value) will be stored
    * @return Err: Error status after attempt of operation. Error could occur if vectors' sizes do not match.
    * This error status should follow your internally defined error schemas/codes
    */
    template <typename T>
    static Err calculateRMSE(
            const Eigen::VectorX<T> &vec1,
            const Eigen::VectorX<T> &vec2,
            double *rmseOutput
            ) {
        ERR_INIT
        e = ErrorUtils::isEqual(vec1.size(), vec2.size()); ree;

        const Eigen::VectorX<T> diff = vec1 - vec2;
        T squaredNorm = diff.squaredNorm();
        *rmseOutput = std::sqrt(squaredNorm / static_cast<double>(vec1.size()));

        ERR_RETURN
    }

    /*!
    * @brief  Calculates the standard deviation of the elements in an Eigen::VectorX
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param vec: Eigen::VectorX whose standard deviation needs to be calculated
    * @return double: The calculated standard deviation of the elements in 'vec'
    */
    template <typename T>
    static double calculateStDevOfVector(const Eigen::VectorX<T> &vec) {
        Eigen::VectorX<T> diff = vec.array() - vec.mean();
        T variance = diff.squaredNorm() / static_cast<double>(vec.size());
        return std::sqrt(variance);
    }

    /*!
    * @brief  Rolls (rotates) the elements in an Eigen::RowVectorXd by a specified distance. Translated values at edges will be wrapped.
    * @param vec: Eigen::RowVectorXd which needs to be rolled
    * @param rollDistance: The distance by which vec needs to be rolled. If rollDistance is positive, it rolls the elements to the right.
    *   If it's negative, elements are rolled to the left. If rollDistance is 0, the original vec is returned.
    * @return Eigen::RowVectorXd: The rolled vector
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
    * @brief  Fits a polynomial to data points via QR Decomposition 
    * @param points: Eigen::MatrixXd with each row representing a point (x value in column 0 and y value in column 1)
    * @param order: The order of the polynomial to fit to the data points
    * @param coeffs: Pointer to a QVector<double> where the coefficients of the fitted polynomial will be stored. coeffs[i] represents the coefficient of the power i term.
    * The coefficients are in increasing order, which means that the first element is the zero order coefficient, the second element is the first order coefficient, and so on.
    */
    static void fitPolynomialQRDecomposition(
            const Eigen::MatrixXd &points,
            int order,
            QVector<double> *coeffs
    ) {

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

        coeffs->resize(order + 1);
        for (int i = 0; i <= order; ++i) {
            (*coeffs)[i] = result[i];
        }

    }

    /*!
    * @brief  Normalizes an Eigen::VectorX by dividing all elements by the maximum coefficient in the vector.
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param vec: Pointer to an Eigen::VectorX. This function modifies the vector in-place.
    */
    template <typename T>
    static void normalizeVector(Eigen::VectorX<T> *vec) {
        const double normalizerDemonator = std::max(std::numeric_limits<double>::min(), vec->maxCoeff());
        *vec /= normalizerDemonator;
    }

    /*!
    * @brief  Normalizes an Eigen::MatrixX by dividing all elements by the maximum coefficient in the Matrix.
    * @tparam T: The datatype of the elements in the Eigen::MatrixX. Could be any numeric type such as int, float, double etc.
    * @param vec: Pointer to an Eigen::MatrixX. This function modifies the Matrix in-place.
    */
    template <typename T>
    static void normalizeMatrix(Eigen::MatrixX<T> *mat) {
        const double normalizerDemonator = std::max(std::numeric_limits<double>::min(), mat->maxCoeff());
        *mat /= normalizerDemonator;
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

    /**
     * @brief   Thresholds the given matrix by setting all elements that are less than the given threshold value to zero.
     *
     * This function applies a threshold operation on an Eigen matrix. If an element of the matrix is less than the
     * specified threshold value, the function replaces that element with zero. Elements equal to or greater
     * than the threshold are left unchanged.
     *
     * @tparam  T   The datatype of the elements in the Eigen::Matrix. Could be any numeric type such as int, float, double etc.
     *
     * @param   thresholdValue  The threshold value. Any elements in the matrix that are less than this value will be set to zero.
     * @param   mat             Pointer to an Eigen::Matrix object. The matrix to be thresholded.
     *                          This function modifies the matrix in-place.
     */
    template <typename T>
    static void thresholdMatrix(T thresholdValue, Eigen::MatrixX<T> *mat) {
        *mat = (mat->array() < thresholdValue).select(0.0, *mat);
    }

    /*!
    * @brief  Applies a threshold function on an Eigen::MatrixX. If an element is less than 'thresholdValue', that element is replaced by 'fillVal'.
    * @tparam T: The datatype of the elements in the Eigen::MatrixX. Could be any numeric type such as int, float, double etc.
    * @param thresholdValue: The value below which the matrix values are replaced by 'fillVal'
    * @param fillVal: The value with which the matrix values (which are less than 'thresholdValue') are replaced
    * @param mat: Pointer to an Eigen::MatrixX. This function modifies the matrix in-place.
    */
    template <typename T>
    static void thresholdMatrix(T thresholdValue, T fillVal, Eigen::MatrixX<T> *mat) {
        *mat = (mat->array() < thresholdValue).select(fillVal, *mat);
    }

    /*!
    * @brief  Applies a threshold function on an Eigen::VectorX. If an element is less than 'thresholdValue', that element is replaced by 0.0.
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param thresholdValue: The value below which the vector values are replaced by 0.0
    * @param vec: Pointer to an Eigen::VectorX. This function modifies the vector in-place.
    */
    template <typename T>
    static void thresholdVector(T thresholdValue, Eigen::VectorX<T> *vec) {
        *vec = (vec->array() < thresholdValue).select(0.0, *vec);
    }

    /*!
    * @brief  Applies a threshold function on an Eigen::VectorX. If an element is less than 'thresholdValue', that element is replaced by 'fillVal'.
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param thresholdValue: The value below which the vector values are replaced by 'fillVal'
    * @param fillVal: The value with which the vector values (which are less than 'thresholdValue') are replaced
    * @param vec: Pointer to an Eigen::VectorX. This function modifies the vector in-place.
    */
    template <typename T>
    static void thresholdVector(T thresholdValue, T fillVal, Eigen::VectorX<T> *vec) {
        *vec = (vec->array() < thresholdValue).select(fillVal, *vec);
    }

    template <typename T>
    static T cosineSimilarity(const Eigen::VectorX<T> &v1, const Eigen::VectorX<T> &v2) {

        const T nearZero = 1e-5;
        const T dotProduct = v1.dot(v2);
        const T v1Magnitude = std::max(std::sqrt(v1.array().square().sum()), nearZero);
        const T v2Magnitude = std::max(std::sqrt(v2.array().square().sum()), nearZero);

        return MathUtils::pRound(dotProduct / (v1Magnitude * v2Magnitude), 4);
    }

    template <typename T>
    static T klDivergence(Eigen::VectorX<T> v1, Eigen::VectorX<T> v2) {

        const T nearZero = 1e-5;
        const T threshold = 0.0;
        thresholdVector(threshold, nearZero, &v1);
        thresholdVector(threshold, nearZero, &v2);

        v1 = v1.array() + nearZero;
        v2 = v2.array() + nearZero;

        v1 = v1.array() / v1.sum();
        v2 = v2.array() / v2.sum();

        Eigen::VectorX<T> log_ratio = v1.array().log() - v2.array().log();
        return v1.dot(log_ratio);
    }

    template<typename T>
    static Eigen::VectorX<T> rowWiseKLDivergenceOfMatrices(
            const Eigen::MatrixX<T> &mat1,
            const Eigen::MatrixX<T> &mat2
    ) {

        const double nearZero = 1e-5;

        Eigen::MatrixX<double> mat1Sum = mat1.array().rowwise().sum();
        thresholdMatrix(0.0, nearZero, &mat1Sum);

        Eigen::MatrixX<double> mat2Sum = mat2.array().rowwise().sum();
        thresholdMatrix(0.0, nearZero, &mat2Sum);

        mat1 /= mat1Sum;
        mat2 /= mat2Sum;

        mat1 = (mat1.array() < nearZero).select(nearZero, mat1);
        mat2 = (mat2.array() < nearZero).select(nearZero, mat2);

        const Eigen::MatrixX<double> mat1mat2Quotient = mat1.array() / mat2.array();
        const Eigen::MatrixX<double> mat1mat2QuotientLog2 = mat1mat2Quotient.array().log2();
        const double mat1mat2QuotientLog2Sum = mat1mat2QuotientLog2.array().sum();

        return mat1.array() * mat1mat2QuotientLog2Sum;
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
    static QPair<int, T> returnTopIndexAndValue(const Eigen::VectorX<T> &vec) {

        QPair<int, T> maxIndex = {-1, std::numeric_limits<T>::min()};

        for (int i = 0; i < vec.size(); i++) {

            const T currentVal = vec.coeff(i);

            if (currentVal > maxIndex.second) {
                maxIndex = {i, currentVal};
            }
        }

        return maxIndex;
    }

    template <typename T>
    static QVector<QPair<int, T>> returnTopXIndexAndValues(const Eigen::VectorX<T> &vec, int topX) {

        topX = std::min( (topX), static_cast<int>(vec.size() - 1) );

        std::vector<T> stdVec(vec.data(), vec.data() + vec.size());

        QVector<QPair<int, T>> ogIndexPoints;
        for (int i = 0; i < stdVec.size(); i++) {
            ogIndexPoints.push_back({i, stdVec.at(i)});
        }

        const auto sortLogic = [](const QPair<int, T> &l, const QPair<int, T> &r){return l.second > r.second;};
        std::partial_sort(ogIndexPoints.begin(), ogIndexPoints.begin() + topX, ogIndexPoints.end(), sortLogic);

        ogIndexPoints.resize(topX);

        return ogIndexPoints;
    }

    template <typename T>
    static Eigen::VectorX<T> convertQVectorToEigenVector(const QVector<T> &_vec) {

        std::vector<T> vec(_vec.begin(), _vec.end());
        Eigen::VectorX<T> ev
                = Eigen::Map<Eigen::VectorX<T>, Eigen::Unaligned>(vec.data(), vec.size());

        return ev;
    }

    template <typename T>
    static QVector<T> convertEigenVectorToQVector(const Eigen::VectorX<T> &vec) {
        std::vector<T> vecReturn(vec.data(), vec.data() + vec.size());
        return QVector<T>(vecReturn.begin(), vecReturn.end());
    }

    template <typename T>
    static Eigen::VectorX<T> convertQMapToEigenVector(
            const QMap<int, T> &m,
            int pointCount
            ) {

        Eigen::VectorX<T> vec(pointCount);
        vec.setZero();

        for (auto it = m.begin(); it != m.end(); it++) {

            const int vecIndex = it.key();
            const double val = it.value();

            if (vecIndex < 0 || vecIndex >= pointCount) {
                continue;
            }

            vec.coeffRef(vecIndex) = val;
        }

        return vec;
    }

    /*!
    * \brief Finds all apexes in a eigen vector.
     *
     *  Returns a QMap where the key is the index of the apex and the value is the value of
     *  the apex.
    */
    template <typename T>
    static QMap<int, T> apexes(
            const Eigen::VectorX<T> &_vec,
            int precision = 1e4
    ){

        QMap<int, T> apexIndicies;

        Eigen::VectorX<T> vec = _vec;
        if (MathUtils::tZero(vec.sum())) {
            return {};
        }

        vec /= vec.maxCoeff();

        const int vecSize = static_cast<int>(vec.size());
        for (int index = 0; index < vecSize; index++) {

            const auto centerPointValue = static_cast<int>(std::round(vec.coeff(index) * precision));

            if(index < 1){

                const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

                if (centerPointValue >= rightPointValue) {
                    apexIndicies.insert(index, _vec.coeff(index));
                }

                continue;
            }

            if (index >= vecSize - 1) {

                const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));

                if (centerPointValue > leftPointValue) {
                    apexIndicies.insert(index, _vec.coeff(index));
                }

                continue;
            }

            const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));
            const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

            if(centerPointValue > leftPointValue && centerPointValue >= rightPointValue){
                apexIndicies.insert(index, _vec.coeff(index));
            }
        }

        return apexIndicies;
    }

    /*!
    * \brief Finds all apexes in a eigen vector and returns the apexes only w/ no regard to order.
     *
     *  Returns a QVector of all found apexes
    */
    template <typename T>
    static QVector<int> apexesIndexesOnly(
            const Eigen::VectorX<T> &_vec,
            int precision = 1e4
    ){

        QVector<int> apexIndicies;

        Eigen::VectorX<T> vec = _vec;
        if (MathUtils::tZero(vec.sum())) {
            return {};
        }

        vec /= vec.maxCoeff();

        const int vecSize = static_cast<int>(vec.size());
        for (int index = 0; index < vecSize; index++) {

            const auto centerPointValue = static_cast<int>(std::round(vec.coeff(index) * precision));

            if(index < 1){

                const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

                if (centerPointValue >= rightPointValue) {
                    apexIndicies.push_back(index);
                }

                continue;
            }

            if (index >= vecSize - 1) {

                const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));

                if (centerPointValue > leftPointValue) {
                    apexIndicies.push_back(index);
                }

                continue;
            }

            const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));
            const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

            if(centerPointValue > leftPointValue && centerPointValue >= rightPointValue){
                apexIndicies.push_back(index);
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
            const Eigen::VectorX<T> &_vec,
            int precision = 1e4
                    ){

        QMap<int, T> troughtIndicies;

        Eigen::VectorX<T> vec = _vec;
        if (MathUtils::tZero(vec.sum())) {
            return {};
        }

        vec /= vec.maxCoeff();

        const int vecSize = static_cast<int>(vec.size());
        for (int index = 0; index < vecSize; index++) {

            const auto centerPointValue = static_cast<int>(std::round(vec.coeff(index) * precision));

            if(index < 1){

                const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

                if (centerPointValue <= rightPointValue) {
                    troughtIndicies.insert(index, _vec.coeff(index));
                }

                continue;
            }

            if (index >= vecSize - 1) {

                const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));

                if (centerPointValue < leftPointValue) {
                    troughtIndicies.insert(index, _vec.coeff(index));
                }

                continue;
            }

            const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));
            const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

            if(centerPointValue < leftPointValue && centerPointValue <= rightPointValue){
                troughtIndicies.insert(index, _vec.coeff(index));
            }
        }

        return troughtIndicies;
    }

    template<typename T>
    static Eigen::VectorX<T> rowWiseCosineSimilarOfMatrices(
            const Eigen::MatrixX<T> &mat1,
            const Eigen::MatrixX<T> &mat2
            ) {

        const Eigen::MatrixX<double> matElementWiseProd = mat1.cwiseProduct(mat2);
        const Eigen::MatrixX<double> matElementWiseProdSum = matElementWiseProd.rowwise().sum();

        const Eigen::MatrixX<double> mat1Norm = mat1.array().rowwise().norm();
        const Eigen::MatrixX<double> mat2Norm = mat2.array().rowwise().norm();

        Eigen::MatrixX<double> mat1NormMat2NormProd = mat1Norm.array() * mat2Norm.array();

        const double nearZero = 0.000001;
        EigenUtils::thresholdMatrix(nearZero, nearZero, &mat1NormMat2NormProd);

        const Eigen::MatrixX<double> rowWiseCosineSimilarity = matElementWiseProdSum.array() / mat1NormMat2NormProd.array();

        return rowWiseCosineSimilarity;
    }

    template <typename T>
    static Eigen::VectorX<T> rowWiseKLDivergence(
            const Eigen::MatrixX<T> &mat1,
            const Eigen::MatrixX<T> &mat2
            ) {

        Eigen::VectorX<T> mat1Sum = mat1.rowwise().sum();
        Eigen::VectorX<T> mat2Sum = mat2.rowwise().sum();

        const double nearZero = 0.000001;
        thresholdVector(nearZero, nearZero, &mat1Sum);
        thresholdVector(nearZero, nearZero, &mat2Sum);

        Eigen::MatrixX<T> m1 = mat1.array().colwise() / mat1Sum.array();
        Eigen::MatrixX<T> m2 = mat2.array().colwise() / mat2Sum.array();

        thresholdMatrix(nearZero, nearZero, &m1);
        thresholdMatrix(nearZero, nearZero, &m1);

        Eigen::MatrixX<double> klDivMat
                = (m1.array() * Eigen::log2(m1.array() / m2.array())).rowwise().sum();

        double infReplaceVal = 10000.0;
        replaceNaN(nearZero, &klDivMat);
        replaceInf(infReplaceVal, &klDivMat);

        return klDivMat;
    }

    static Eigen::VectorX<double> convertQPointFVecToEigen(
            const QVector<QPointF> &points,
            int precision,
            double valMax
            );

    enum class ThresholderDirection {
        Above = 0
        , Below
    };

    template<typename T>
    static void removeRowsAboveOrBelowThreshold(
            T threshold,
            int columnIdxToFilter,
            ThresholderDirection thresholderDirection,
            Eigen::MatrixX<T> *mat
            ) {

        Eigen::VectorXd thresholderVec;

        if (thresholderDirection == ThresholderDirection::Below) {
            const Eigen::VectorXd belowThreshold
                    = (mat->col(columnIdxToFilter).array() < threshold).rowwise().any().template cast<double>();
            thresholderVec = belowThreshold;
        }

        else {
            const Eigen::VectorXd aboveThreshold
                    = (mat->col(columnIdxToFilter).array() > threshold).rowwise().any().template cast<double>();
            thresholderVec = aboveThreshold;
        }

        const int numRows = mat->rows();
        const int numCols = mat->cols();

        const int newNumRows = static_cast<int>((thresholderVec.array() == 0).count());
        Eigen::MatrixX<T> newMatrix(newNumRows, numCols);
        int newRowIdx = 0;
        for (int i = 0; i < numRows; ++i) {
            if (thresholderVec(i) == 0) {
                newMatrix.row(newRowIdx) = mat->row(i);
                ++newRowIdx;
            }
        }
        *mat = newMatrix;
    }

    template<typename T>
    static void minMaxScaleVector(Eigen::VectorX<T> *vec) {

        const T maxVal = vec->maxCoeff();
        const T minVal = vec->minCoeff();
        const T minMaxDiff = std::max(maxVal - minVal, std::numeric_limits<T>::min());

        *vec = (vec->array() - minVal) / minMaxDiff;
    }

    template<typename T>
    static void minMaxScaleMatrix(Eigen::MatrixX<T> *mat) {
        for (int col = 0; col < mat->cols(); col++) {

            Eigen::VectorX<T> vec = mat->col(col);
            minMaxScaleVector(&vec);
            mat->col(col) = vec;
        }

        EigenUtils::replaceNaN(static_cast<T>(0.0), mat);
    }

    template<typename T>
    static Eigen::MatrixX<T> convertQVectorsToEigenMatrix(const QVector<QVector<T>> &matA) {

        const int rows = matA.size();
        const int columns = matA.front().size();
        Eigen::MatrixX<T> mat(rows, columns);
        mat.setZero();

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                mat.coeffRef(i, j) = matA[i][j];
            }
        }

        return mat;
    }

    template<typename T>
    static Eigen::MatrixX<T> convertQVectorsToEigenMatrix(const QVector<QVector<T>*> &matA) {

        const int rows = matA.size();
        const int columns = matA.front()->size();
        Eigen::MatrixX<T> mat(rows, columns);
        mat.setZero();

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < columns; j++) {
                mat.coeffRef(i, j) = matA[i]->at(j);
            }
        }

        return mat;
    }

    template<typename T>
    static QVector<QVector<T>> convertEigenMatrixToQVectors(const Eigen::MatrixX<T> &mat) {

        QVector<QVector<T>> vecs;
        for (int row = 0; row < mat.rows(); row++) {

            Eigen::VectorX<T> v = mat.row(row);
            vecs.push_back(EigenUtils::convertEigenVectorToQVector<T>(v));
        }

        return vecs;
    }

    template<typename T>
    static Eigen::VectorX<T> trimVector(const Eigen::VectorX<T> &vec) {

        const int center = EigenUtils::returnTopIndexAndValue(vec).first;

        int currentIndex = center - 1;
        int leftStart = currentIndex;
        while(currentIndex >= 0) {
            const double currentValue = vec.coeff(currentIndex);
            if (MathUtils::tZero(currentValue)) {
                break;
            }
            leftStart = currentIndex;
            currentIndex--;
        }

        currentIndex = center + 1;
        int rightStart = currentIndex;
        while(currentIndex < vec.size()) {
            const double currentValue = vec.coeff(currentIndex);
            if (MathUtils::tZero(currentValue)) {
                break;
            }
            rightStart = currentIndex;
            currentIndex++;
        }

        return vec.segment(leftStart, rightStart - leftStart + 1);
    }

};

#endif //EIGENUTILS_H
