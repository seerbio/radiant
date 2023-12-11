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

    template<typename T>
    static void replaceNaN(T replaceVal, Eigen::VectorX<T> *vec){
        *vec = vec->array().isNaN().select(replaceVal, vec->array());
    }

    template<typename T>
    static void replaceNaN(T replaceVal, Eigen::MatrixX<T> *mat){
        *mat = mat->array().isNaN().select(replaceVal, mat->array());
    }

    template<typename T>
    static void replaceInf(T replaceVal, Eigen::MatrixX<T> *mat){
        *mat = mat->unaryExpr([replaceVal](double v) { return std::isinf(v) ? replaceVal : v; });
    }

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
    static void fitPolynomialQRDecomposition(
            const Eigen::MatrixXd &points,
            int order,
            QVector<double> *coeff
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
    static void thresholdMatrix(T thresholdValue, Eigen::MatrixX<T> *mat) {
        *mat = (mat->array() < thresholdValue).select(0.0, *mat);
    }

    template <typename T>
    static void thresholdMatrix(T thresholdValue, T fillVal, Eigen::MatrixX<T> *mat) {
        *mat = (mat->array() < thresholdValue).select(fillVal, *mat);
    }

    template <typename T>
    static void thresholdVector(T thresholdValue, Eigen::VectorX<T> *vec) {
        *vec = (vec->array() < thresholdValue).select(0.0, *vec);
    }

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
    static double klDivergence(Eigen::VectorX<T> v1, Eigen::VectorX<T> v2) {

        const double nearZero = 1e-5;
        const double threshold = 0.0;
        thresholdVector(threshold, nearZero, &v1);
        thresholdVector(threshold, nearZero, &v2);

        v1 = v1.array() + nearZero;
        v2 = v2.array() + nearZero;

        v1 = v1.array() / v1.sum();
        v2 = v2.array() / v2.sum();

        Eigen::VectorXd log_ratio = v1.array().log() - v2.array().log();
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

//        topX = std::min(topX, static_cast<int>(vec.size()));
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
    static Err writeMatrixToFile(
            const Eigen::MatrixX<T> &mat,
            const QString &outputFilePath
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(mat.rows() > 0); ree;

        QFile file(outputFilePath);
        if (file.open(QIODevice::ReadWrite)) {

            QTextStream stream(&file);

            for (int i = 0; i < mat.rows(); i++) {

                const Eigen::VectorX<T> v = mat.row(i);

                QStringList rtwQStringList;
                for (int j = 0; j < v.size(); j++) {
                    rtwQStringList.push_back(QString::number(v.coeff(j)));
                }

                stream << rtwQStringList.join(S_GLOBAL_SETTINGS.COMMA);
                stream << S_GLOBAL_SETTINGS.NEWLINE;
            }

        }

        file.close();
        qDebug() << mat.rows() << "Rows written to" << outputFilePath;

        ERR_RETURN
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
    }

    template<typename T>
    static Eigen::MatrixX<double> convertQVectorsToEigenMatrix(const QVector<QVector<T>> &matA) {

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
