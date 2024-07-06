//
// Created by Drucifer on 11/7/2021.
//

#ifndef EIGENUTILS_H
#define EIGENUTILS_H

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
    * @param coeffs: Pointer to a QVector<double> where the coefficients of the fitted polynomial will be stored. coeffs[i]
    *      represents the coefficient of the power i term. The coefficients are in increasing order, which means that the
    *      first element is the zero order coefficient, the second element is the first order coefficient, and so on.
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
    * @brief  Adds random noise to an Eigen::VectorX or Eigen::MatrixXd object. The object could be a vector or a matrix.
    * @tparam EigenMatrix: Type of Eigen object (could be Eigen::VectorX<T> or Eigen::MatrixXd)
    * @param scaleFactor: Factor by which the random noise is scaled before being added to the object
    * @param vec: Pointer to an Eigen object (vector or matrix). This function modifies the object in-place by adding noise.
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

    /*!
    * @brief  Calculates the cosine similarity between two Eigen::VectorX and stores the result in the passed pointer
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param v1: First Eigen::VectorX
    * @param v2: Second Eigen::VectorX
    * @param cosineSimilarity: Pointer to a variable of type T where the result (cosine similarity) will be stored. The cosine similarity is rounded to 4 decimal places.
    * @return Err: Error status after attempt of operation. Error could occur if vectors' sizes do not match.
    * This error status should follow your internally defined error schemas/codes
    */
    template <typename T>
    static Err cosineSimilarity(
            const Eigen::VectorX<T> &v1,
            const Eigen::VectorX<T> &v2,
            T *cosineSimilarity
            ) {

        ERR_INIT
        e = ErrorUtils::isEqual(v1.size(), v2.size()); ree;

        const T nearZero = 1e-5;
        const T dotProduct = v1.dot(v2);
        const T v1Magnitude = std::max(std::sqrt(v1.array().square().sum()), nearZero);
        const T v2Magnitude = std::max(std::sqrt(v2.array().square().sum()), nearZero);

        *cosineSimilarity = MathUtils::pRound(dotProduct / (v1Magnitude * v2Magnitude), 4);

        ERR_RETURN
    }

    /*!
    * @brief  Calculates the Kullback-Leibler (KL) Divergence between two Eigen::VectorX and stores the result in the passed pointer
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param v1: First Eigen::VectorX
    * @param v2: Second Eigen::VectorX
    * @param klDiv: Pointer to a variable of type T where the result (KL Divergence) will be stored.
    * @return Err: Error status after attempt of operation. Error could occur if vectors' sizes do not match.
    * This error status should follow your internally defined error schemas/codes
    * Before calculating the KL Divergence, any values less than zero in both vectors are replaced by a small arbitrary value (nearZero), and then both vectors are normalized.
    */
    template <typename T>
    static Err klDivergence(
            Eigen::VectorX<T> v1,
            Eigen::VectorX<T> v2,
            T *klDiv
            ) {

        ERR_INIT
        e = ErrorUtils::isEqual(v1.size(), v2.size()); ree;

        const T nearZero = 1e-5;
        const T threshold = 0.0;
        thresholdVector(threshold, nearZero, &v1);
        thresholdVector(threshold, nearZero, &v2);

        v1 = v1.array() + nearZero;
        v2 = v2.array() + nearZero;

        v1 = v1.array() / v1.sum();
        v2 = v2.array() / v2.sum();

        Eigen::VectorX<T> log_ratio = v1.array().log() - v2.array().log();
        *klDiv = v1.dot(log_ratio);

        ERR_RETURN
    }

    /*!
    * @brief  Calculates the row-wise cosine similarity between two Eigen::MatrixX and stores the result in the passed pointer
    * @tparam T: The datatype of the elements in the Eigen::MatrixX. Could be any numeric type such as int, float, double etc.
    * @param mat1: The first Eigen::MatrixX
    * @param mat2: The second Eigen::MatrixX
    * @param cosineSimilarities: Pointer to a Eigen::VectorX where the result (cosine similarities for each row) will be stored.
    * @return Err: Error status after attempt of operation. Error could occur if matrices' dimensions do not match.
    * This error status should follow your internally defined error schemas/codes
    * During the calculation, original matrices are not modified. If the norm of any row in either matrix is below a small threshold, it is replaced by the threshold to prevent division by zero.
    */
    template<typename T>
    static  Err rowWiseCosineSimilarOfMatrices(
            const Eigen::MatrixX<T> &mat1,
            const Eigen::MatrixX<T> &mat2,
            Eigen::VectorX<T> *cosineSimilarities
    ) {

        ERR_INIT

        e = ErrorUtils::isEqual(mat1.rows(), mat2.rows()); ree;
        e = ErrorUtils::isEqual(mat1.cols(), mat2.cols()); ree;

        const Eigen::MatrixX<T> matElementWiseProd = mat1.cwiseProduct(mat2);
        const Eigen::MatrixX<T> matElementWiseProdSum = matElementWiseProd.rowwise().sum();

        const Eigen::MatrixX<T> mat1Norm = mat1.array().rowwise().norm();
        const Eigen::MatrixX<T> mat2Norm = mat2.array().rowwise().norm();

        Eigen::MatrixX<T> mat1NormMat2NormProd = mat1Norm.array() * mat2Norm.array();

        const double nearZero = 0.000001;
        EigenUtils::thresholdMatrix<T>(nearZero, nearZero, &mat1NormMat2NormProd);

        const Eigen::MatrixX<T> rowWiseCosineSimilarity = matElementWiseProdSum.array() / mat1NormMat2NormProd.array();

        *cosineSimilarities = rowWiseCosineSimilarity;

        ERR_RETURN
    }

    /*!
    * @brief  Calculates the row-wise Kullback-Leibler (KL) Divergence between two Eigen::MatrixX and stores the result in the passed pointer
    * @tparam T: The datatype of the elements in the Eigen::MatrixX. Could be any numeric type such as int, float, double etc.
    * @param mat1: The first Eigen::MatrixX
    * @param mat2: The second Eigen::MatrixX
    * @param klDivs: Pointer to an Eigen::VectorX where the result (KL Divergences for each row) will be stored.
    * @return Err: Error status after attempt of operation. Error could occur if matrices' dimensions do not match.
    * This error status should follow your internally defined error schemas/codes
    * During the calculation, each row in the matrices are normalized, small values are replaced by a small arbitrary value (nearZero), and then KL Divergence for each row is calculated.
    * The calculated KL Divergences for each row are stored in the 'klDivs' vector.
    */
    template<typename T>
    static Err rowWiseKLDivergenceOfMatrices(
            const Eigen::MatrixX<T> &_mat1,
            const Eigen::MatrixX<T> &_mat2,
            Eigen::VectorX<T> *klDivs
    ) {

        ERR_INIT

        Eigen::MatrixX<T> mat1 = _mat1;
        Eigen::MatrixX<T> mat2 = _mat2;

        e = ErrorUtils::isEqual(mat1.rows(), mat2.rows()); ree;
        e = ErrorUtils::isEqual(mat1.cols(), mat2.cols()); ree;

        const T nearZero = 1e-5;

        Eigen::VectorX<T> mat1Sum = mat1.array().rowwise().sum();
        thresholdVector(static_cast<T>(0.0), nearZero, &mat1Sum);

        Eigen::VectorX<T> mat2Sum = mat2.array().rowwise().sum();
        thresholdVector(static_cast<T>(0.0), nearZero, &mat2Sum);

        mat1 = mat1.array() / mat1Sum.array();
        mat2 = mat2.array() / mat2Sum.array();

        mat1 = (mat1.array() < nearZero).select(nearZero, mat1);
        mat2 = (mat2.array() < nearZero).select(nearZero, mat2);

        const Eigen::MatrixX<T> mat1mat2Quotient = mat1.array() / mat2.array();
        const Eigen::MatrixX<T> mat1mat2QuotientLog2 = mat1mat2Quotient.array().log2();
        const double mat1mat2QuotientLog2Sum = mat1mat2QuotientLog2.array().sum();

        *klDivs = mat1.array() * mat1mat2QuotientLog2Sum;
        // std::cout << Eigen::RowVectorX<float>(_mat1) << std::endl;
        // std::cout << Eigen::RowVectorX<float>(_mat2) << std::endl;
        // std::cout << Eigen::RowVectorX<float>(*klDivs) << std::endl;
        // std::cout << "***********" << std::endl;

        ERR_RETURN
    }

    /*!
    * @brief  Calculates the count of non-zero elements in an Eigen object (could be Eigen::MatrixX or Eigen::VectorX)
    * @tparam EigenMatrix: Type of the Eigen object (could be Eigen::MatrixX<T>, Eigen::VectorX<T> etc.)
    * @param v1: Eigen object
    * @return int: The count of non-zero elements in 'v1'
    */
    template <typename EigenMatrix>
    static int nonZeros(EigenMatrix v1) {
        return (v1.array() != 0).count();;
    }

    /*!
    * @brief  Returns the index and maximum value in an Eigen::VectorX
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param vec: Eigen::VectorX from which the maximum value and its index are to be determined
    * @return QPair<int, T>: A QPair in which the first element is the index of the maximum value in the vector, and the second element is the maximum value itself
    */
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

    /*!
    * @brief  Returns the indices and values of the top 'topX' elements in an Eigen::VectorX
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param vec: Eigen::VectorX from which the top 'topX' values and their indices are to be determined
    * @param topX: The number of top elements to return. If 'topX' is larger than the size of vec, then all elements are returned.
    * @return QVector<QPair<int, T>>: A QVector of QPairs, each QPair's first element is an index of a top element and the second element is the corresponding value from the input vector.
    */
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

    /*!
    * @brief  Converts a QVector to an Eigen::VectorX
    * @tparam T: The datatype of the elements in the QVector and the created Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param _vec: QVector that is to be converted to Eigen::VectorX
    * @return Eigen::VectorX<T>: Resulting Eigen::VectorX which is created from the input QVector
    */
    template <typename T>
    static Eigen::VectorX<T> convertQVectorToEigenVector(const QVector<T> &_vec) {

        std::vector<T> vec(_vec.begin(), _vec.end());
        Eigen::VectorX<T> ev
                = Eigen::Map<Eigen::VectorX<T>, Eigen::Unaligned>(vec.data(), vec.size());

        return ev;
    }

    /*!
    * @brief  Converts an Eigen::VectorX to a QVector
    * @tparam T: The datatype of the elements in the Eigen::VectorX and the created QVector. Could be any numeric type such as int, float, double etc.
    * @param vec: Eigen::VectorX that is to be converted to QVector
    * @return QVector<T>: Resulting QVector which is created from the input Eigen::VectorX
    */
    template <typename T>
    static QVector<T> convertEigenVectorToQVector(const Eigen::VectorX<T> &vec) {
        std::vector<T> vecReturn(vec.data(), vec.data() + vec.size());
        return QVector<T>(vecReturn.begin(), vecReturn.end());
    }

    /*!
    * @brief  Converts a QMap to an Eigen::VectorX
    * @tparam T: The datatype of the elements in the QMap and the created Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param m: QMap that is to be converted to Eigen::VectorX. The keys in the QMap represent the indices in the vector, and the values represent the corresponding values in the vector.
    * @param pointCount: The size of the created Eigen::VectorX. Any key in the QMap that is outside of the range [0, pointCount-1] is ignored.
    * @return Eigen::VectorX<T>: Resulting Eigen::VectorX which is created from the input QMap. The vector is initialized to size 'pointCount' with all elements set to zero, and then values are filled from the QMap.
    */
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
    * @brief  Identifies and returns the apexes in an Eigen::VectorX as a QMap
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param _vec: Eigen::VectorX in which the apexes are to be identified
    * @param precision: The desired level of precision when comparing values from 'vec'. The default value is 1e4.
    * @return QMap<int, T>: QMap where the keys represent the indices of the apexes in the vector, and the corresponding values are the values of the apexes in the vector.
    * An apex is defined as a point which is greater than or equal to the points around it (on either side) in the vector.
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
    * @brief  Identifies and returns the indices of the apexes in an Eigen::VectorX as a QVector
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param _vec: Eigen::VectorX in which the apexes are to be identified
    * @param precision: The desired level of precision when comparing values from 'vec'. The default value is 1e4.
    * @return QVector<int>: QVector containing the indices of the apexes in the vector.
    * An apex is defined as a point which is greater than or equal to the points around it (on either side) in the vector.
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
    * @brief  Identifies and returns the troughs in an Eigen::VectorX as a QMap
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param _vec: Eigen::VectorX in which the troughs (local minima) are to be identified
    * @param precision: The desired level of precision when comparing values from 'vec'. The default value is 1e4.
    * @return QMap<int, T>: QMap where the keys represent the indices of the troughs in the vector, and the corresponding values are the values of the troughs in the vector.
    * A trough is defined as a point which is less than or equal to the points around it (on either side) in the vector.
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

    enum class ThresholderDirection {
        Above = 0
        , Below
    };

    /*!
    * @brief  Removes rows from an Eigen::MatrixX that are above or below a specified threshold in a specified column
    * @tparam T: The datatype of the elements in the Eigen::MatrixX. Could be any numeric type such as int, float, double etc.
    * @param threshold: The value to be used as a threshold for filtering rows
    * @param columnIdxToFilter: The index of column based on which rows should be filtered
    * @param thresholderDirection: A ThresholderDirection enum value indicating whether to filter rows that are above or below the threshold
    * @param mat: Pointer to an Eigen::MatrixX from which rows are to be filtered. The matrix will be updated in place to exclude filtered out rows.
    */
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

    /*!
    * @brief  Scales an Eigen::VectorX to the range [0, 1] using Min-Max scaling
    * @tparam T: The datatype of the elements in the Eigen::VectorX. Could be any numeric type such as int, float, double etc.
    * @param vec: Pointer to an Eigen::VectorX which is to be scaled. The vector will be updated in place.
    * The scaling is done using the formula x_sc = (x - x_min) / (x_max - x_min), where x_max and x_min are the maximum and minimum values of the elements in the vector.
    * If the vector has the same value for all elements (i.e., max = min), then all elements of the scaled vector will be zero due to the equation becoming 0/0. If so, they are set as zero to avoid division by zero.
    */
    template<typename T>
    static void minMaxScaleVector(Eigen::VectorX<T> *vec) {

        const T maxVal = vec->maxCoeff();
        const T minVal = vec->minCoeff();
        const T minMaxDiff = std::max(maxVal - minVal, std::numeric_limits<T>::min());

        *vec = (vec->array() - minVal) / minMaxDiff;
    }

    /*!
    * @brief  Scales each column of an Eigen::MatrixX to the range [0, 1] using Min-Max scaling
    * @tparam T: The datatype of the elements in the Eigen::MatrixX. Could be any numeric type such as int, float, double etc.
    * @param mat: Pointer to an Eigen::MatrixX which is to be scaled. The matrix will be updated in place.
    * Each column in the matrix is treated and scaled individually,
    * After the scaling operation, replaceNaN() is called to replace any NaN values in the matrix with 0.
    */
    template<typename T>
    static void minMaxScaleMatrix(Eigen::MatrixX<T> *mat) {
        for (int col = 0; col < mat->cols(); col++) {

            Eigen::VectorX<T> vec = mat->col(col);
            minMaxScaleVector(&vec);
            mat->col(col) = vec;
        }

        EigenUtils::replaceNaN(static_cast<T>(0.0), mat);
    }

    /*!
    * @brief  Converts a QVector of QVectors to an Eigen::MatrixX
    * @tparam T: The datatype of the elements in the QVectors and the created Eigen::MatrixX. Could be any numeric type such as int, float, double etc.
    * @param matA: QVector of QVectors that is to be converted to Eigen::MatrixX
    * @return Eigen::MatrixX<T>: Resulting Eigen::MatrixX which is created from the input QVector of QVectors
    */
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

    /*!
    * @brief  Converts a QVector of pointers to QVectors to an Eigen::MatrixX
    * @tparam T: The datatype of the elements in the QVector and the created Eigen::MatrixX. Could be any numeric type such as int, float, double etc.
    * @param matA: QVector of pointers to QVectors that is to be converted to Eigen::MatrixX
    * @return Eigen::MatrixX<T>: Resulting Eigen::MatrixX which is created from the input QVector of QVectors
    * Note: This function assumes the QVector contains non-null pointers. Best practice would include validating the pointers before calling this function.
    */
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

    /*!
    * @brief  Converts an Eigen::MatrixX to a QVector of QVectors
    * @tparam T: The datatype of the elements in the Eigen::MatrixX and the created QVector. Could be any numeric type such as int, float, double etc.
    * @param mat: The Eigen::MatrixX that is to be converted to QVectors
    * @return QVector<QVector<T>>: Resulting QVector of QVectors which is created from the input Eigen::MatrixX
    */
    template<typename T>
    static QVector<QVector<T>> convertEigenMatrixToQVectors(const Eigen::MatrixX<T> &mat) {

        QVector<QVector<T>> vecs;
        for (int row = 0; row < mat.rows(); row++) {

            Eigen::VectorX<T> v = mat.row(row);
            vecs.push_back(EigenUtils::convertEigenVectorToQVector<T>(v));
        }

        return vecs;
    }

    /*!
    * @brief  Trims an Eigen::VectorX to remove zero-value elements from start and end
    * @tparam T: The datatype of the elements in the Eigen::VectorX. The datatype should be suitable for comparison operations.
    * @param vec: Eigen::VectorX that is to be trimmed
    * @return Eigen::VectorX<T>: Resulting Eigen::VectorX which is a segment of the input vector, with initial and trailing zero-value elements removed.
    * The function works by first finding the 'center' position which is the top (maximum) value in the vector. It then moves towards both ends from the center,
    * stopping when it encounters a zero-value element or reaches the end of the vector. The segment between the two end positions (excluding zero-value elements) is returned.
    */
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

    // TODO add test
    template<typename T>
    static Err simpleIntegrator(
        const Eigen::VectorX<T> &vec,
        T stopThresholdFraction,
        QVector<QPair<PeakIntegrationIndexes, T>> *peakIntegrationIndexesVsIntensity
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(vec.size() > 0); ree;

        Eigen::VectorX<T> eVec = vec;
        EigenUtils::thresholdVector(static_cast<float>(1.01), &eVec);

        const QMap<int, T> vecApexs = EigenUtils::apexes(eVec);
        if (vecApexs.isEmpty()) {
            ERR_RETURN
        }

        Eigen::VectorX<T> apexes =EigenUtils::convertQMapToEigenVector(vecApexs, vecApexs.lastKey() + 1);
        QVector<QPair<int, T>> apexPairs = EigenUtils::returnTopXIndexAndValues(apexes, vecApexs.size());

        for (const QPair<int, T> &pr : apexPairs) {

            const int apexIndex = pr.first;
            const T apexValue = pr.second;

            if (MathUtils::tZero(apexValue)) {
                continue;
            }

            const T stopThreshold = apexValue * stopThresholdFraction;

            T rightStopVal = apexValue;
            int rightStopIndex = apexIndex;

            int rightCurrentIndex = apexIndex;
            while (rightCurrentIndex < eVec.size()) {

                const T currentValue = eVec(rightCurrentIndex);
                if (currentValue < stopThreshold) {
                    rightStopIndex = rightCurrentIndex;
                    break;
                }

                if (currentValue < rightStopVal || MathUtils::tSame(currentValue, rightStopVal)) {
                    rightStopVal = currentValue;
                    rightStopIndex = rightCurrentIndex;
                    ++rightCurrentIndex;
                    continue;
                }

                break;
            }

            T leftStopVal = apexValue;
            int leftStopIndex = apexIndex;

            int leftCurrentIndex = apexIndex;
            while (leftCurrentIndex < eVec.size()) {

                const T currentValue = eVec(leftCurrentIndex);
                if (currentValue < stopThreshold) {
                    leftStopIndex = leftCurrentIndex;
                    break;
                }

                if (currentValue < leftStopVal || MathUtils::tSame(currentValue, leftStopVal)) {
                    leftStopVal = currentValue;
                    leftStopIndex = leftCurrentIndex;
                    --leftCurrentIndex;
                    continue;
                }

                break;
            }

            peakIntegrationIndexesVsIntensity->push_back({
                 {std::max(leftStopIndex, 0), std::min(rightStopIndex, static_cast<int>(vec.size() - 1))},
                 apexValue
            }
            );

            for (int i = leftStopIndex; i <= rightStopIndex; i++) {
                eVec.coeffRef(i) = 0.0;
            }
        }

        std::sort(
            peakIntegrationIndexesVsIntensity->rbegin(),
            peakIntegrationIndexesVsIntensity->rend(),
            [](const QPair<PeakIntegrationIndexes, T> &l, const QPair<PeakIntegrationIndexes, T> &r) {
                return l.second < r.second;
            }
            );

        ERR_RETURN
    }

    template<typename T>
    static Eigen::VectorXf removeZeroElements(const Eigen::VectorX<T> &vec){

        const int count = (vec.array() != 0).count();
        Eigen::VectorX<T> result(count);

        int index = 0;
        for (int i = 0; i < vec.size(); i++) {

            const T v = vec.coeff(i);
            if (MathUtils::tZero(v)) {
                continue;
            }

            result[index++] = v;
        }

        return result;
    }

};

#endif //EIGENUTILS_H
