//
// Created by anichols on 11/07/2021.
//
#ifndef EIGEN_KERNEL_UTILITIES_H
#define EIGEN_KERNEL_UTILITIES_H

#include "EigenLib_Exports.h"

#include "EigenUtils.h"
#include "EigenSparseUtils.h"
#include "ErrorUtils.h"
#include "MathUtils.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>


class EIGENLIB_EXPORTS EigenKernelUtils
{

public:

    template <typename T>
    static void normalizeVector(Eigen::MatrixX<T> *vec) {
        const double normalizerDemonator = std::max(std::numeric_limits<double>::min(), vec->maxCoeff());
        *vec /= normalizerDemonator;
    }

    template <typename T>
    static Err buildSavitzkyGolayKernel(
            int windowSize,
            int order,
            int derivative,
            int rate,
            Eigen::MatrixX<T> *mat
    ) {

        //https://en.wikipedia.org/wiki/Savitzky–Golay_filter

        ERR_INIT

        e = ErrorUtils::isNotEqual(windowSize % 2 != 1, true, eValueError); ree;

        e = ErrorUtils::isNotEqual(windowSize < order + 2, true, eValueError); ree;

        const int halfWindow = std::floor((windowSize - 1) / 2);
        const int cols = order + 1;
        const int rows = (halfWindow * 2) + 1;

        mat->resize(rows, cols);
        for (int i = -halfWindow; i <= halfWindow; ++i) {
            for (int j = 0; j < cols; j++) {
                mat->coeffRef(i + halfWindow, j) = std::pow(i, j);
            }
        }

        int factorialDerivative = 1;
        for (int i = 1; i <= derivative; ++i) {
            factorialDerivative *= i;
        }

        const Eigen::MatrixX<T> pinv = mat->completeOrthogonalDecomposition().pseudoInverse();
        Eigen::MatrixX<T> kernel = pinv.row(derivative) * std::pow(rate, derivative) * factorialDerivative;

        *mat = kernel.transpose();
        return e;

        ERR_RETURN
    }

    /*!
    * @brief creates a uniform grid spaced equally from start to stop
    * of specified length == num
    */
    template <typename T>
    static Eigen::VectorX<T> lineSpace(
            double start,
            double stop,
            int num,
            bool endPoint
    ) {

        if (num < 0){
            return {};
        }

        const double div = endPoint
                           ? static_cast<double>(num - 1)
                           : static_cast<double>(num);

        const double delta = stop - start;

        QVector<T> r(num);
        std::iota(r.begin(), r.end(), 0.0);

        Eigen::VectorX<T> y = Eigen::Map<Eigen::VectorX<T>, Eigen::Unaligned>(r.data(), r.size());

        if (div > 0 && !MathUtils::tZero(div)){

            double step = delta / div;

            if (std::any_of(r.begin(), r.end(), [](double d){return MathUtils::tZero(d);})){
                y /= div;
                y *= delta;
            }

            else{
                y *= step;
            }
        }

        else{
            y = y * delta;
        }

        y = y.array() + start;

        if (endPoint && num > 1){
            y.coeffRef(y.size() -1) = stop;
        }

        return y;
    }

    template <typename T>
    static Eigen::VectorX<T> buildGaussianFilter1D (
            int filterLength,
            double sigma,
            bool normalize = false
    ) {
        //https://www.geeksforgeeks.org/gaussian-filter-generation-c/
        //code from link was modified to only be 1D

        Eigen::VectorX<T> ls = lineSpace<T>(-sigma, sigma, filterLength, true);

        double r, s = sigma;

        for (int i = 0; i < ls.size(); i++) {

            double &l = ls.coeffRef(i);

            r = sqrt(l * l);
            l = (std::exp( -(r * r) / s)) / (static_cast<float>(M_PI) * s);

        }

        if(normalize && !MathUtils::tZero(ls.maxCoeff())){
            ls /= ls.maxCoeff();
        }

        return ls;
    }

    template <typename T>
    static Eigen::VectorX<T> buildMexicanHatFilter1D(
            int filterLength,
            double width,
            bool normalize = false
    ){

        const double A = 2 / (std::sqrt(3 * width) * std::pow(M_PI,0.25));
        const double wsq = std::pow(width, 2);

        QVector<T> r(filterLength);
        std::iota(r.begin(), r.end(), 0.0);
        //TODO replace w/ EigenUtils::convertQvectorToEigenVector
        Eigen::VectorX<T> vec
                = Eigen::Map<Eigen::VectorX<T>, Eigen::Unaligned>(r.data(), r.size());
        vec = vec.array() - (filterLength - 1.0) / 2;

        const Eigen::VectorX<T> xsq = vec.cwiseProduct(vec);
        const Eigen::VectorX<T> xsqwsq = xsq / wsq;
        const Eigen::VectorX<T> mod = 1 - xsqwsq.array();

        Eigen::VectorX<T> gauss = -xsq / (2 * wsq);
        gauss = Eigen::exp(gauss.array());

        const Eigen::VectorX<T> total = A * mod.cwiseProduct(gauss);

        return normalize ? total / total.maxCoeff() : total;
    }

    /*!
    * \brief Calculates the number of times a filter can be convolved given a vector length,
     * filter length, and a skip value.
    */
    static int calculateNumberOfStrides(
            int vectorLength,
            int strideLength,
            int strideSkip = 1
    );

    template <typename T>
    static Eigen::VectorX<T> convolveVectorWithKernel(
            const Eigen::VectorX<T> &_vec,
            const Eigen::VectorX<T> &_kernel
            ) {

        if (_vec.size() <= _kernel.size()){
            return _vec;
        }

        Eigen::VectorX<T> vec = _vec;
        Eigen::VectorX<T> kernel = _kernel;

        if(kernel.cols() != 1){
            kernel = kernel.transpose();
        }

        const int paddingAmount = kernel.size() - 1;
        Eigen::VectorX<T> vecResized(vec.size() + paddingAmount);

        const int halfWindow = std::floor(paddingAmount / 2.0);
        const T frontPadding = static_cast<T>(EigenUtils::calculateMeanOfVector(vec.head(halfWindow + 1)));
        const T backPadding = static_cast<T>(EigenUtils::calculateMeanOfVector(vec.tail(halfWindow + 1)));
        for(int i = 0; i < halfWindow; i++){
            vecResized.coeffRef(i) = frontPadding;
            vecResized.coeffRef(vecResized.size() - 1 - i) = backPadding;
        }

        vecResized.segment(halfWindow, vec.size()) = vec;

        const int filterLength = static_cast<int>(kernel.size());

        Eigen::MatrixX<T> convolutionMatrix(vec.size(), filterLength);
        convolutionMatrix.setZero();

        for (int i = 0; i < filterLength; ++i) {
            convolutionMatrix.col(i) = vecResized.segment(i, vec.size());
        }

        return convolutionMatrix * kernel;
    }

    template <typename T>
    static Eigen::SparseVector<T> convolveVectorWithKernel(
            const Eigen::SparseVector<T> &vec,
            const Eigen::VectorX<T> &kernel
            ) {

        const int kernelSize = static_cast<int>(kernel.size());
        const int filterCorrection = kernelSize % 2 == 1 ? 1 : 0;
        const int ogVectorSize = static_cast<int>(vec.size()) - kernelSize + filterCorrection;

        Eigen::SparseMatrix<T> mat(ogVectorSize, kernelSize);
        for (int i = 0; i < kernelSize; ++i) {
            const Eigen::SparseVector<T> v = vec.segment(i, ogVectorSize);
            mat.col(i) = v;
        }

        Eigen::VectorX<T> dotProd = mat * kernel;
        return dotProd.sparseView();
    }

    template <typename T>
    static Eigen::SparseVector<T> addPaddingToSparseVector(
            const Eigen::SparseVector<T> &vec,
            int filterLength
            ) {

        const int filterLengthHalved = std::floor(filterLength/2.0);
        const int rows = static_cast<int>(vec.size());
        const int newRows = rows + (filterLengthHalved * 2);

        Eigen::SparseVector<T> paddedVec(newRows);
        paddedVec.setZero();

        for (typename Eigen::SparseVector<T>::InnerIterator it(vec); it; ++it) {
            const T &val = it.value();
            const int &ogIndex = it.index();

            paddedVec.insert(ogIndex + filterLengthHalved) = val;
        }

        return paddedVec;
    }

    template <typename T>
    static Eigen::SparseMatrix<T, Eigen::RowMajor> addPaddingToSparseMatrixRowWise(
            const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat,
            int filterLength
            ) {
        const int filterLengthHalved = std::floor(filterLength/2.0);
        const int rows = static_cast<int>(mat.rows());
        const int newRows = rows + (filterLengthHalved * 2);
        const int cols = static_cast<int>(mat.cols());

        Eigen::SparseMatrix<T, Eigen::RowMajor> paddedMatrix(newRows, cols);
        paddedMatrix.setZero();

        paddedMatrix.middleRows(filterLengthHalved, rows) = mat;

        return paddedMatrix;
    }

    template <typename T>
    static Eigen::MatrixX<T> addPaddingToMatrixRowWise(
            const Eigen::MatrixX<T> &mat,
            int filterLength
    ) {

        const int filterLengthHalved = std::floor(filterLength/2.0);
        const int rows = static_cast<int>(mat.rows());
        const int newRows = rows + (filterLengthHalved * 2);
        const int cols = static_cast<int>(mat.cols());

        Eigen::MatrixX<T> paddedMatrix(newRows, cols);
        paddedMatrix.setZero();

        paddedMatrix.block(filterLengthHalved, 0, rows, cols) = mat;

        return paddedMatrix;
    }

    template <typename T>
    static Eigen::SparseMatrix<T> addPaddingToSparseMatrixColWise(
            const Eigen::SparseMatrix<double> &_mat,
            int filterLength
    ) {

        //This is a hack because assigning middleCols() below is behaving funky unless
        //storage order is changed set to RowMajor.  //TODO Drew Nichols figure this out.
        const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat = _mat;

        const int filterLengthHalved = std::floor(filterLength / 2.0);
        const int rows = static_cast<int>(mat.rows());
        const int cols = static_cast<int>(mat.cols());
        const int newCols = cols + (filterLengthHalved * 2);

        Eigen::SparseMatrix<T> paddedMatrix(rows, newCols);
        paddedMatrix.setZero();

        paddedMatrix.middleCols(filterLengthHalved, cols) = mat;

        return paddedMatrix;
    }

    template <typename T>
    static Eigen::SparseMatrix<T> applyKernelColumnWiseToMatrix(
            const Eigen::SparseMatrix<T> &_mat,
            const Eigen::VectorX<T> &kernel,
            bool matchOriginalMaximum
            ) {

        const int filterLength = static_cast<int>(kernel.size());

        const int minimumAllowableFilterLength = 2;
        if(filterLength < minimumAllowableFilterLength){
            return _mat;
        }

        Eigen::SparseMatrix<T> mat = addPaddingToSparseMatrixColWise<T>(_mat, filterLength);

        const int matRows = static_cast<int>(mat.rows());
        const int matCols = static_cast<int>(mat.cols());

        const int resultMatrixCols = calculateNumberOfStrides(matCols, filterLength);
        Eigen::SparseMatrix<T> resultMatrix(matRows, resultMatrixCols);
        resultMatrix.setZero();

        for (int i = 0; i < filterLength; ++i) {

            const T multiplier = kernel.coeff(i);

            Eigen::SparseMatrix<T> addingMatrix = mat.middleCols(i, resultMatrixCols);

            addingMatrix *= multiplier;
            resultMatrix += addingMatrix;
        }

        if(matchOriginalMaximum){
            const double originalMaximum = EigenSparseUtils::max(_mat);

            double smoothedMaximum = EigenSparseUtils::max(resultMatrix);
            smoothedMaximum = MathUtils::tZero(smoothedMaximum) ? 1 : smoothedMaximum;

            return resultMatrix * (originalMaximum / smoothedMaximum);
        }

        return resultMatrix;
    }

    template <typename T>
    static Eigen::MatrixX<T> applyKernelToEachColumnInMatrix(
            const Eigen::MatrixX<T> &_mat,
            const Eigen::VectorX<T> &kernel
            ) {

        const int filterLength = static_cast<int>(kernel.size());

        const int minimumAllowableFilterLength = 2;
        if(filterLength < minimumAllowableFilterLength){
            return _mat;
        }

        const Eigen::MatrixX<T> matPadded = addPaddingToMatrixRowWise<T>(_mat, filterLength);

        const int matRows = static_cast<int>(matPadded.rows());
        const int matCols = static_cast<int>(matPadded.cols());

        const int resultMatrixRows = calculateNumberOfStrides(matRows, filterLength);
        Eigen::MatrixX<T> resultMatrix(resultMatrixRows, matCols);
        resultMatrix.setZero();

        for (int i = 0; i < filterLength; ++i) {

            const double multiplier = kernel.coeff(i);

            Eigen::MatrixX<T> addingMatrix = matPadded.middleRows(i, resultMatrixRows);

            addingMatrix *= multiplier;
            resultMatrix += addingMatrix;
        }

        return resultMatrix;
    }

    template <typename T>
    static Err savitskyGolaySmooth(
            int windowSize,
            int order = 1,
            int derivative = 0,
            int rate = 1,
            Eigen::VectorX<T> *smoothedVec = Q_NULLPTR
            ) {

        ERR_INIT

        Eigen::MatrixX<T> savitskyGolayKernel;
        e = buildSavitzkyGolayKernel(
                windowSize,
                order,
                derivative,
                rate,
                &savitskyGolayKernel
                ); ree;

        Eigen::VectorX<T> vec = *smoothedVec;

        *smoothedVec = convolveVectorWithKernel<T>(
                vec,
                savitskyGolayKernel
                );

        ERR_RETURN
    }

};


#endif // EIGEN_KERNEL_UTILITIES_H