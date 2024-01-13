#ifndef EIGENMATRIXDISKUTILS_H
#define EIGENMATRIXDISKUTILS_H

#include "Error.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace Error;

class EigenMatrixDiskUtils
{
public:

    /*!
    * @brief   Saves an Eigen::SparseMatrix<double, Eigen::ColMajor> to a disk file.
    * @param   matrix: Matrix object to save.
    * @param   filePath: Path of the file to save the matrix to.
    * @return  Err: Error object which contains the details about whether error has occurred during the IO operations.
    *
    * Note: This function first serializes the input matrix using the QDataStream class and saves the serialized data to a file. If the file doesn't exist, it will be created. If it already exists, its contents will be replaced. Error checking is performed at each step using the ErrorUtils::isTrue() and ErrorUtils::isNotEqual() methods, and if errors occur during the process, the function will exit early and return those errors.
    */
    static Err saveEigenMatrix(
            const Eigen::SparseMatrix<double, Eigen::ColMajor> &matrix,
            const QString &filePath
                    );

    /*!
    * @brief   Loads an Eigen::SparseMatrix<double, Eigen::ColMajor> from a disk file.
    * @param   filePath: Path of the file to load the matrix from.
    * @param   matrix: Pointer to Matrix object where the data will be loaded.
    * @return  Err: Error object which contains the details about whether error has occurred during the IO operations.
    *
    * Note: This function first opens the file in read-only mode and reads all the data from it. Then it creates a QDataStream object and sets its version. The data from the QDataStream object is deserialized into the provided matrix object. Error checking is performed at each step using the ErrorUtils::isTrue() and ErrorUtils::isNotEqual() methods, and if errors occur during the process, the function will exit early and return those errors.
    */
    static Err loadEigenMatrix(
            const QString &filePath,
            Eigen::SparseMatrix<double, Eigen::ColMajor> *matrix
            );

    /*!
    * @brief   Saves an Eigen::MatrixXd to a disk file.
    * @param   matrix: Dense MatrixXd object to save.
    * @param   filePath: Path of the file to save the matrix to.
    * @return  Err: Error object which contains the details about whether error has occurred during the IO operations.
    *
    * Note: This function first converts the input dense matrix to a sparse matrix view using Eigen's .sparseView() function. Then, it calls the saveEigenMatrix() function to save the sparse matrix to the disk. If error occurs during the process, it is returned to the caller.
    */
    static Err saveEigenMatrix(
            const Eigen::MatrixXd &matrix,
            const QString &filePath
            );

    /*!
    * @brief   Loads an Eigen::MatrixXd from a disk file.
    * @param   filePath: Path of the file to load the matrix from.
    * @param   matrix: Pointer to MatrixXd object where the data will be loaded.
    * @return  Err: Error object which contains the details about whether error has occurred during the IO operations.
    *
    * Note: This function first creates an empty SparseMatrix. It then calls the loadEigenMatrix() function to load the sparse matrix from the disk. The loaded sparse matrix is then converted to a dense matrix and assigned to the input matrix pointer. If error occurs during the process, it is returned to the caller.
    */
    static Err loadEigenMatrix(
            const QString &filePath,
            Eigen::MatrixXd *matrix
            );

};

#endif // EIGENMATRIXDISKUTILS_H


