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

    /*!
    * @brief   Builds a Savitzky-Golay smoothing kernel.
    * @tparam  T: Numeric type used in Eigen Matrix. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   windowSize: The size of the window for the filter, usually a small odd integer.
    * @param   order: The order of the polynomial to be fitted to the points within the window.
    * @param   derivative: The order of the derivative to be computed. Generally, set to 0 for smoothing or 1 to compute the first derivative.
    * @param   rate: Scaling factor usually set to 1.
    * @param   mat: Pointer to Eigen Matrix where the output kernel will be stored.
    * @return  Err: Error object which contains the details about whether the function succeeded or not.
    *
    * Note: This function first validates the arguments and then creates a Vandermonde-like matrix to fit the polynomial.
    * It then calculates the pseudo-inverse of the matrix. By using the row corresponding to the derivative,
    * the kernel for the Savitzky–Golay filter is calculated. Finally, it stores the output kernel in the input matrix pointer.
    */
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
    * @brief   Creates a VectorX of evenly spaced values over a defined interval, similar to linspace MATLAB function.
    * @tparam  T: Numeric type used in Eigen VectorX. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   start: Value of start of sequence.
    * @param   stop: Value of end of sequence.
    * @param   num: Number of values to generate. Must be non-negative.
    * @param   endPoint: If true, 'stop' is the last value in the range, otherwise, 'stop' is not included in the range.
    * @return  Eigen::VectorX<T>: Eigen Vector with equally spaced values.
    *
    * Note: This function first constructs a sequence of increasing numbers starting from 0 to num-1 using std::iota. This sequence is then converted into an Eigen::Vector. Depending on the values of 'div' and 'num', different transformations happen to 'y'. Finally, the start value is added to each number in the sequence and if 'endPoint' is true, the last value of the sequence is set to 'stop'.
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

    /*!
    * @brief   Builds a one-dimensional Gaussian filter.
    * @tparam  T: Numeric type used in Eigen VectorX. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   filterLength: The length of output filter, must be non-negative.
    * @param   sigma: The standard deviation used in the Gaussian function.
    * @param   normalize: If true, scale the filter values so that the maximum value is 1.
    * @return  Eigen::VectorX<T>: Eigen Vector representing the Gaussian filter.
    *
    * Note: This function first constructs a sequence of equally spaced values over the interval [-sigma, sigma]. It then computes a Gaussian filter using this sequence. If 'normalize' is true, all values in filter will be divided by maximum value in the filter to normalize values between 0 and 1.
    */
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

            T &l = ls.coeffRef(i);

            r = sqrt(l * l);
            l = (std::exp( -(r * r) / s)) / (static_cast<float>(M_PI) * s);

        }

        if(normalize && !MathUtils::tZero(ls.maxCoeff())){
            ls /= ls.maxCoeff();
        }

        return ls;
    }

    /*!
    * @brief   Builds a one-dimensional Mexican Hat (Ricker Wavelet) filter.
    * @tparam  T: Numeric type used in Eigen VectorX. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   filterLength: The length of the output filter, must be non-negative.
    * @param   width: The width used in the Mexican Hat function.
    * @param   normalize: If true, scale the filter values so that the maximum value is 1.
    * @return  Eigen::VectorX<T>: Eigen Vector representing the Mexican Hat filter.
    *
    * Note: This function first creates a sequence of numbers from 0 to filterLength-1, and transforms it to center around 0. Then, it calculates the Mexican Hat (Ricker wavelet) function on this sequence. If 'normalize' is true, the result is divided by its maximum value to normalize the values between 0 and 1.
    */
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
    * @brief   Calculates the number of strides that can be taken with a certain stride length and stride skip over a vector of a certain length.
    * @param   vectorLength: The length of the vector that is being covered by strides.
    * @param   strideLength: The length of a single stride.
    * @param   strideSkip: The number of elements to skip over between the start of each stride.
    * @return  int: The total number of strides that will occur. If strideSkip is <= 0, returns -1.
    *
    * Note: This function first checks if strideSkip is > 0. If it's not, the function returns -1 as it's invalid to not skip any elements or go backwards. Otherwise, it calculates the number of strides by using the formula ((vectorLength - strideLength) / strideSkip) + 1 and returns this count.
    */
    static int calculateNumberOfStrides(
            int vectorLength,
            int strideLength,
            int strideSkip = 1
    );

    /*!
    * @brief   Convolves a given vector with a given kernel using the method of padding and flipping.
    * @tparam  T: Numeric type used in Eigen VectorX. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   _vec: The input vector to be convolved.
    * @param   _kernel: The kernel to convolve the vector with.
    * @return  Eigen::VectorX<T>: Eigen Vector representing the convolved output.
    *
    * Note: This function first checks if the length of the vector is smaller than the kernel, if true it returns the original vector. Otherwise, it pads the original vector based on the kernel size.   It then builds a so-called toeplitz (diagonal-constant) matrix for the padded vector. Multiplying this matrix with the kernel provides the convolution result.
    */
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
        const T frontPadding = static_cast<T>((vec.head(halfWindow + 1).mean()));
        const T backPadding = static_cast<T>((vec.tail(halfWindow + 1).mean()));
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

    /*!
    * @brief   Adds padding to a given Sparse Vector which will be used for convolution operation.
    * @tparam  T: Numeric type used in Eigen SparseVector. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   vec: The Sparse Vector to which padding needs to be added.
    * @param   filterLength: The length of the filter which will be convolved with the Sparse Vector.
    * @return  Eigen::SparseVector<T>: Eigen SparseVector with added padding.
    *
    * Note: This function creates a new SparseVector with padding added to the original Sparse Vector.
    * The padding length is calculated from the filter length. The padding is added equally at the beginning and the ending of the Vector.
    * The non-zero elements of the original Sparse Vector are then copied to the right position in the padded Sparse Vector.
    */
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

    /*!
    * @brief   Convolves a given sparse vector with a given kernel.
    * @tparam  T: Numeric type used in Eigen SparseVector. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   vec: The input sparse vector to be convolved.
    * @param   kernel: The kernel to convolve the vector with.
    * @return  Eigen::SparseVector<T>: Eigen SparseVector representing the convolved output.
    *
    * Note: This function constructs a sparse matrix, where each column is a segment of the original sparse vector, offset by the index. This constructs overlapping "windows" of the original vector, matching the size of the kernel. These windows are then each convolved with the kernel by taking the product of the sparse matrix with the kernel, resulting in the convolved sparse vector.
    */
    template <typename T>
    static Eigen::SparseVector<T> convolveVectorWithKernel(
            const Eigen::SparseVector<T> &vec,
            const Eigen::VectorX<T> &kernel
            ) {

        Eigen::SparseVector<T> vecPadded = addPaddingToSparseVector(vec, kernel.size());

        Eigen::SparseMatrix<T> mat(vec.size(), kernel.size());
        for (int i = 0; i < kernel.size(); ++i) {
            const Eigen::SparseVector<T> v = vecPadded.segment(i, vec.size());
            mat.col(i) = v;
        }

        Eigen::VectorX<T> dotProd = mat * kernel;
        return dotProd.sparseView();
    }

    /*!
    * @brief   Adds padding to a given Sparse Matrix row-wise.
    * @tparam  T: Numeric type used in Eigen SparseMatrix. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   mat: The SparseMatrix (with RowMajor storage order) to which padding needs to be added.
    * @param   filterLength: The length of the filter which will be convolved with the SparseMatrix.
    * @return  Eigen::SparseMatrix<T, RowMajor>: Eigen SparseMatrix with RowMajor storage order and added padding.
    *
    * Note: This function creates a new SparseMatrix with added padding row-wise to the original SparseMatrix. The padding length is calculated from the filter length. New rows filled with zeros are added equally at the beginning and the ending of the SparseMatrix, resulting in larger matrix dimensions.
    */
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

    /*!
    * @brief   Adds padding to a given matrix row-wise.
    * @tparam  T: Numeric type used in Eigen Matrix. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   mat: The Matrix to which padding needs to be added.
    * @param   filterLength: The length of the filter which will be convolved with the Matrix.
    * @return  Eigen::MatrixX<T>: Eigen Matrix with added padding.
    *
    * Note: This function creates a new Matrix with padding added to the original Matrix.
    * The padding length is calculated from the filter length. The padding is added equally at the beginning and the ending of each row.
    */
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

    /*!
    * @brief   Adds padding to a given Sparse Matrix column-wise.
    * @tparam  T: Numeric type used in Eigen SparseMatrix. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   _mat: The SparseMatrix to which padding needs to be added.
    * @param   filterLength: The length of the filter which will be convolved with the SparseMatrix.
    * @return  Eigen::SparseMatrix<T>: Eigen SparseMatrix with added padding.
    *
    * Note: This function creates a new SparseMatrix with padding added to the original SparseMatrix.
    * The padding length is calculated from the filter length. The padding is added equally at the beginning and the ending of each column.
    * The non-zero elements of the original SparseMatrix are then copied to the right position in the padded SparseMatrix.
    */
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

    /*!
    * @brief   Applies a given kernel to a Sparse Matrix to each row, i.e. left->right.
    * @tparam  T: Numeric type used in Eigen SparseMatrix. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   _mat: The SparseMatrix to which the kernel needs to be applied.
    * @param   kernel: The kernel that will be applied.
    * @param   matchOriginalMaximum: If true, the result matrix will be scaled so that its maximum equals the original matrix's maximum.
    * @return  Eigen::SparseMatrix<T>: Eigen SparseMatrix which is the result of applying the kernel to the input SparseMatrix.
    *
    * Note: This function creates a new SparseMatrix from the input SparseMatrix by padding it column-wise. The kernel, which is a Vector, is then applied to each column. If 'matchOriginalMaximum' is true, the resulting values are scaled such that the maximum value of the result matrix matches the maximum value of the original matrix.
    */
    template <typename T>
    static Eigen::SparseMatrix<T> applyKernelToEachRowInMatrix(
            const Eigen::SparseMatrix<T> &_mat,
            const Eigen::VectorX<T> &kernel,
            bool matchOriginalMaximum = false
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

    /*!
    * @brief   Applies a given kernel to each column of a Matrix.
    * @tparam  T: Numeric type used in Eigen Matrix. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   _mat: The Matrix to which the kernel needs to be applied.
    * @param   kernel: The kernel that will be applied.
    * @return  Eigen::MatrixX<T>: Eigen Matrix which is the result of applying the kernel to the input Matrix.
    *
    * Note: This function creates a new Matrix from the input Matrix by padding it row-wise. The kernel, which is a Vector, is then applied to each column of the matrix by adding padded Matrix rows weighted by kernel coefficients. The result is effectively the convolution of each column of the input matrix with the kernel.
    */
    template <typename T>
    static Eigen::MatrixX<T> applyKernelToEachColumnInMatrix(
            const Eigen::MatrixX<T> &_mat,
            const Eigen::VectorX<T> &kernel
            ) {

        const int filterLength = static_cast<int>(kernel.size());

        constexpr int minimumAllowableFilterLength = 2;
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

    /*!
    * @brief   Applies Savitzky-Golay smoothing to an input vector.
    * @tparam  T: Numeric type used in Eigen Vector. It should support mathematical operations such as multiplication, division, addition, subtraction, and power.
    * @param   windowSize: The size of the window used for smoothing (must be odd and greater than order).
    * @param   order: The order of the polynomial used in smoothing (must be less than windowSize).
    * @param   derivative: The order of the derivative to compute.
    * @param   rate: The rate of the lowess window increment (typically it should be equal to 1).
    * @param   smoothedVec: The input Vector that needs to be smoothed. After the function, it contains the smoothed data.
    * @return  Err: Error code, indicating the success or failure of the function.
    *
    * Note: This function first builds a Savitzky-Golay filter using the supplied parameters (window size, order of the polynomial, derivative order, rate). The function then applies this filter to the input vector for smoothing.
    */
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