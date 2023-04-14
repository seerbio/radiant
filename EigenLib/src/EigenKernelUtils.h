//
// Created by anichols on 11/07/2021.
//
#ifndef EIGEN_KERNEL_UTILITIES_H
#define EIGEN_KERNEL_UTILITIES_H

#include "EigenLib_Exports.h"
#include "ErrorUtils.h"
#include "MathUtils.h"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>


class EIGENLIB_EXPORTS EigenKernelUtils
{

public:

    template <typename EigenMatrix>
    static void normalizeVector(EigenMatrix *vec) {
        const double normalizerDemonator = std::max(std::numeric_limits<double>::min(), vec->maxCoeff());
        *vec /= normalizerDemonator;
    }


    static Err buildSavitzkyGolayKernel(
            int windowSize,
            int order,
            int derivative,
            int rate,
            Eigen::VectorXd *mat
    );


    /*!
    * @brief creates a uniform grid spaced equally from start to stop
    * of specified length == num
    */
    static Eigen::VectorXd lineSpace(
            double start,
            double stop,
            int num,
            bool endPoint = true
    );


    static Eigen::VectorXd buildGaussianFilter1D (
            int filterLength,
            double sigma,
            bool normalize = true
    );


    static Eigen::VectorXd buildMexicanHatFilter1D (
            int filterLength,
            double width,
            bool normalize = true
    );


    /*!
    * \brief Calculates the number of times a filter can be convolved given a vector length,
     * filter length, and a skip value.
    */
    static int calculateNumberOfStrides(
            int vectorLength,
            int strideLength,
            int strideSkip = 1
    );


    static Eigen::VectorXd convolveVectorWithKernel(
            const QVector<double> &vec,
            const Eigen::VectorXd &_kernel
    );


    static Eigen::VectorXd convolveVectorWithKernel(
            const Eigen::VectorXd &_vec,
            const Eigen::VectorXd &_kernel
    );


    static Eigen::SparseVector<double> convolveVectorWithKernel(
            const Eigen::SparseVector<double> &vec,
            const Eigen::VectorXd &kernel
    );


    static Eigen::SparseVector<double> addPaddingToSparseVector(
            const Eigen::SparseVector<double> &vec,
            int filterLength
    );


    static Eigen::SparseMatrix<double, Eigen::RowMajor> addPaddingToSparseMatrixRowWise(
            const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat,
            int filterLength
    );


    static Eigen::SparseMatrix<double> applyKernelRowWiseToMatrix(
            const Eigen::SparseMatrix<double, Eigen::RowMajor> &_mat,
            const Eigen::VectorXd &kernel,
            bool matchOriginalMaximum
    );


    static Eigen::SparseMatrix<double> addPaddingToSparseMatrixColWise(
            const Eigen::SparseMatrix<double> &_mat,
            int filterLength
    );


    static Eigen::SparseMatrix<double> applyKernelColumnWiseToMatrix(
            const Eigen::SparseMatrix<double> &_mat,
            const Eigen::VectorXd &kernel,
            bool matchOriginalMaximum
    );

    static Err savitskyGolaySmooth(
            int windowSize,
            int order = 1,
            int derivative = 0,
            int rate = 1,
            Eigen::VectorX<double> *smoothedVec = Q_NULLPTR
            ) {

        ERR_INIT

        Eigen::VectorX<double> savitskyGolayKernel;
        e = buildSavitzkyGolayKernel(
                windowSize,
                order,
                derivative,
                rate,
                &savitskyGolayKernel
                ); ree;

        *smoothedVec = convolveVectorWithKernel(
                *smoothedVec,
                savitskyGolayKernel
                );

        ERR_RETURN
    }

};


#endif // EIGEN_KERNEL_UTILITIES_H