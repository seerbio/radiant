//
// Created by anichols on 10/21/22.
//

#include "EigenKernelUtils.h"

#include "EigenUtils.h"
#include "EigenSparseUtils.h"

#include <QtCore>

Err EigenKernelUtils::buildSavitzkyGolayKernel(
        int windowSize,
        int order,
        int derivative,
        int rate,
        Eigen::MatrixX<double> *mat
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

    const Eigen::MatrixX<double> pinv = mat->completeOrthogonalDecomposition().pseudoInverse();
    Eigen::MatrixX<double> kernel = pinv.row(derivative) * std::pow(rate, derivative) * factorialDerivative;

    *mat = kernel.transpose();
    return e;

    ERR_RETURN
}


Eigen::VectorXd EigenKernelUtils::lineSpace(
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

    QVector<double> r(num);
    std::iota(r.begin(), r.end(), 0.0);
    //TODO replace w/ EigenUtils::convertQvectorToEigenVector
    Eigen::VectorXd y = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());

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


Eigen::VectorXd EigenKernelUtils::buildGaussianFilter1D (
        int filterLength,
        double sigma,
        bool normalize
) {
    //https://www.geeksforgeeks.org/gaussian-filter-generation-c/
    //code from link was modified to only be 1D

    Eigen::VectorXd ls = lineSpace(-sigma, sigma, filterLength);

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


Eigen::VectorXd EigenKernelUtils::buildMexicanHatFilter1D(
        int filterLength,
        double width,
        bool normalize
){

    const double A = 2 / (std::sqrt(3 * width) * std::pow(M_PI,0.25));
    const double wsq = std::pow(width, 2);

    QVector<double> r(filterLength);
    std::iota(r.begin(), r.end(), 0.0);
    //TODO replace w/ EigenUtils::convertQvectorToEigenVector
    Eigen::VectorXd vec
            = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(r.data(), r.size());
    vec = vec.array() - (filterLength - 1.0) / 2;

    const Eigen::VectorXd xsq = vec.cwiseProduct(vec);
    const Eigen::VectorXd xsqwsq = xsq / wsq;
    const Eigen::VectorXd mod = 1 - xsqwsq.array();

    Eigen::VectorXd gauss = -xsq / (2 * wsq);
    gauss = Eigen::exp(gauss.array());

    const Eigen::VectorXd total = A * mod.cwiseProduct(gauss);

    return normalize ? total / total.maxCoeff() : total;
}


int EigenKernelUtils::calculateNumberOfStrides(
        int vectorLength,
        int strideLength,
        int strideSkip
) {
    if (strideSkip <= 0) {
        return -1;
    }

    return static_cast<int>((vectorLength - strideLength) / strideSkip) + 1;
}


Eigen::VectorXd EigenKernelUtils::convolveVectorWithKernel(
        const Eigen::VectorXd &_vec,
        const Eigen::VectorXd &_kernel
) {

    if (_vec.size() <= _kernel.size()){
        return _vec;
    }

    Eigen::VectorXd vec = _vec;
    Eigen::VectorXd kernel = _kernel;

    if(kernel.cols() != 1){
        kernel = kernel.transpose();
    }

    const int paddingAmount = kernel.size() - 1;
    Eigen::VectorXd vecResized(vec.size() + paddingAmount);

    const int halfWindow = std::floor(paddingAmount / 2.0);
    const float frontPadding = static_cast<float>(EigenUtils::calculateMeanOfVector(vec.head(halfWindow + 1)));
    const float backPadding = static_cast<float>(EigenUtils::calculateMeanOfVector(vec.tail(halfWindow + 1)));
    for(int i = 0; i < halfWindow; i++){
        vecResized.coeffRef(i) = frontPadding;
        vecResized.coeffRef(vecResized.size() - 1 - i) = backPadding;
    }

    vecResized.segment(halfWindow, vec.size()) = vec;

    const int filterLength = static_cast<int>(kernel.size());

    Eigen::MatrixXd convolutionMatrix(vec.size(), filterLength);
    convolutionMatrix.setZero();

    for (int i = 0; i < filterLength; ++i) {
        convolutionMatrix.col(i) = vecResized.segment(i, vec.size());
    }

    return convolutionMatrix * kernel;
}


Eigen::VectorXd EigenKernelUtils::convolveVectorWithKernel(
        const QVector<double> &_vec,
        const Eigen::VectorXd &_kernel
) {

    std::vector<double> vec = _vec.toStdVector();
    //TODO replace w/ EigenUtils::convertQvectorToEigenVector
    Eigen::VectorX<double> ev
        = Eigen::Map<Eigen::VectorX<double>, Eigen::Unaligned>(vec.data(), vec.size());


    return convolveVectorWithKernel(ev, _kernel);
}


Eigen::SparseVector<double> EigenKernelUtils::convolveVectorWithKernel(
        const Eigen::SparseVector<double> &vec,
        const Eigen::VectorXd &kernel
) {

    const int kernelSize = static_cast<int>(kernel.size());
    const int filterCorrection = kernelSize % 2 == 1 ? 1 : 0;
    const int ogVectorSize = static_cast<int>(vec.size()) - kernelSize + filterCorrection;

    Eigen::SparseMatrix<double> mat(ogVectorSize, kernelSize);
    for (int i = 0; i < kernelSize; ++i) {
        const Eigen::SparseVector<double> v = vec.segment(i, ogVectorSize);
        mat.col(i) = v;
    }

    Eigen::VectorX<double> dotProd = mat * kernel;
    return dotProd.sparseView();
}


Eigen::SparseVector<double> EigenKernelUtils::addPaddingToSparseVector(
        const Eigen::SparseVector<double> &vec,
        int filterLength
) {

    const int filterLengthHalved = std::floor(filterLength/2.0);
    const int rows = static_cast<int>(vec.size());
    const int newRows = rows + (filterLengthHalved * 2);

    Eigen::SparseVector<double> paddedVec(newRows);
    paddedVec.setZero();

    for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it) {
        const double &val = it.value();
        const int &ogIndex = it.index();

        paddedVec.insert(ogIndex + filterLengthHalved) = val;
    }

    return paddedVec;
}


Eigen::SparseMatrix<double, Eigen::RowMajor> EigenKernelUtils::addPaddingToSparseMatrixRowWise(
        const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat,
        int filterLength
) {
    const int filterLengthHalved = std::floor(filterLength/2.0);
    const int rows = static_cast<int>(mat.rows());
    const int newRows = rows + (filterLengthHalved * 2);
    const int cols = static_cast<int>(mat.cols());

    Eigen::SparseMatrix<double, Eigen::RowMajor> paddedMatrix(newRows, cols);
    paddedMatrix.setZero();

    paddedMatrix.middleRows(filterLengthHalved, rows) = mat;

    return paddedMatrix;
}


Eigen::SparseMatrix<double> EigenKernelUtils::applyKernelRowWiseToMatrix(
        const Eigen::SparseMatrix<double, Eigen::RowMajor> &_mat,
        const Eigen::VectorXd &kernel,
        bool matchOriginalMaximum
) {

    const int filterLength = static_cast<int>(kernel.size());

    const int minimumAllowableFilterLength = 2;
    if(filterLength < minimumAllowableFilterLength){
        return _mat;
    }

    Eigen::SparseMatrix<double, Eigen::RowMajor> mat
            = addPaddingToSparseMatrixRowWise(_mat, filterLength);

    const int matRows = static_cast<int>(mat.rows());
    const int matCols = static_cast<int>(mat.cols());

    const int resultMatrixRows = calculateNumberOfStrides(matRows, filterLength);
    Eigen::SparseMatrix<double> resultMatrix(resultMatrixRows, matCols);
    resultMatrix.setZero();

    for (int i = 0; i < filterLength; ++i) {

        const double multiplier = kernel.coeff(i);

        Eigen::SparseMatrix<double> addingMatrix = mat.middleRows(i, resultMatrixRows);

        addingMatrix *= multiplier;
        resultMatrix += addingMatrix;

    }

    if(matchOriginalMaximum){

        const double originalMaximum = EigenSparseUtils::max(_mat);
        double smoothedMaximum = EigenSparseUtils::max(resultMatrix);
        smoothedMaximum = MathUtils::tZero(smoothedMaximum) ? 1 : smoothedMaximum;

        return resultMatrix * (originalMaximum/smoothedMaximum);
    }

    return resultMatrix;
}


Eigen::SparseMatrix<double> EigenKernelUtils::addPaddingToSparseMatrixColWise(
        const Eigen::SparseMatrix<double> &_mat,
        int filterLength
) {

    //This is a hack because assigning middleCols() below is behaving funky unless
    //storage order is changed set to RowMajor.  //TODO Drew Nichols figure this out.
    const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat = _mat;

    const int filterLengthHalved = std::floor(filterLength / 2.0);
    const int rows = static_cast<int>(mat.rows());
    const int cols = static_cast<int>(mat.cols());
    const int newCols = cols + (filterLengthHalved * 2);

    Eigen::SparseMatrix<double> paddedMatrix(rows, newCols);
    paddedMatrix.setZero();

    paddedMatrix.middleCols(filterLengthHalved, cols) = mat;

    return paddedMatrix;
}


Eigen::SparseMatrix<double> EigenKernelUtils::applyKernelColumnWiseToMatrix(
        const Eigen::SparseMatrix<double> &_mat,
        const Eigen::VectorXd &kernel,
        bool matchOriginalMaximum
) {

    const int filterLength = static_cast<int>(kernel.size());

    const int minimumAllowableFilterLength = 2;
    if(filterLength < minimumAllowableFilterLength){
        return _mat;
    }

    Eigen::SparseMatrix<double> mat = addPaddingToSparseMatrixColWise(_mat, filterLength);

    const int matRows = static_cast<int>(mat.rows());
    const int matCols = static_cast<int>(mat.cols());

    const int resultMatrixCols = calculateNumberOfStrides(matCols, filterLength);
    Eigen::SparseMatrix<double> resultMatrix(matRows, resultMatrixCols);
    resultMatrix.setZero();

    for (int i = 0; i < filterLength; ++i) {

        const double multiplier = kernel.coeff(i);

        Eigen::SparseMatrix<double> addingMatrix = mat.middleCols(i, resultMatrixCols);

        addingMatrix *= multiplier;
        resultMatrix += addingMatrix;
    }

    if(matchOriginalMaximum){
        const double originalMaximum = EigenSparseUtils::max(_mat);

        double smoothedMaximum = EigenSparseUtils::max(resultMatrix);
        smoothedMaximum = MathUtils::tZero(smoothedMaximum) ? 1 : smoothedMaximum;

        return resultMatrix * (originalMaximum/smoothedMaximum);
    }

    return resultMatrix;
}

Eigen::MatrixX<double> EigenKernelUtils::applyKernelRowWiseToMatrix(
        const Eigen::MatrixX<double> &_mat,
        const Eigen::VectorXd &kernel
        ) {

    const int filterLength = static_cast<int>(kernel.size());

    const int minimumAllowableFilterLength = 2;
    if(filterLength < minimumAllowableFilterLength){
        return _mat;
    }

    Eigen::MatrixX<double> matPadded = addPaddingToMatrixRowWise(_mat, filterLength);

    const int matRows = static_cast<int>(matPadded.rows());
    const int matCols = static_cast<int>(matPadded.cols());

    const int resultMatrixRows = calculateNumberOfStrides(matRows, filterLength);
    Eigen::MatrixX<double> resultMatrix(resultMatrixRows, matCols);
    resultMatrix.setZero();

    for (int i = 0; i < filterLength; ++i) {

        const double multiplier = kernel.coeff(i);

        Eigen::MatrixX<double> addingMatrix = matPadded.middleRows(i, resultMatrixRows);

        addingMatrix *= multiplier;
        resultMatrix += addingMatrix;
    }

    return resultMatrix;
}

Eigen::MatrixX<double> EigenKernelUtils::addPaddingToMatrixRowWise(
        const Eigen::MatrixX<double> &mat,
        int filterLength
        ) {

    const int filterLengthHalved = std::floor(filterLength/2.0);
    const int rows = static_cast<int>(mat.rows());
    const int newRows = rows + (filterLengthHalved * 2);
    const int cols = static_cast<int>(mat.cols());

    Eigen::MatrixX<double> paddedMatrix(newRows, cols);
    paddedMatrix.setZero();

    paddedMatrix.block(filterLengthHalved, 0, rows, cols) = mat;

    return paddedMatrix;
}
