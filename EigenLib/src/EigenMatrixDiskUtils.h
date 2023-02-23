#ifndef EIGENMATRIXDISKUTILS_H
#define EIGENMATRIXDISKUTILS_H

#include "Error.h"

#include <Eigen/Core>
#include <Eigen/Sparse>

using namespace Error;

class EigenMatrixDiskUtils
{
public:

    static Err saveEigenMatrix(const Eigen::SparseMatrix<double, Eigen::ColMajor> &matrix, const QString &filePath);


    static Err loadEigenMatrix(const QString &filePath, Eigen::SparseMatrix<double, Eigen::ColMajor> *matrix);


    static Err saveEigenMatrix(const Eigen::MatrixXd &matrix, const QString &filePath);


    static Err loadEigenMatrix(const QString &filePath, Eigen::MatrixXd *matrix);

};

#endif // EIGENMATRIXDISKUTILS_H


