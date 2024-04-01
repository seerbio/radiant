//
// Created by anichols on 11/07/2021.
//

#ifndef EIGENSPARSEUTILS_H
#define EIGENSPARSEUTILS_H

#include "EigenLib_Exports.h"
#include "GlobalSettings.h"
#include "MathUtils.h"

#include <QMap>
#include <QPointF>
#include <QVector>

#include <Eigen/Sparse>

#include <iostream>
#include <vector>


class EIGENLIB_EXPORTS EigenSparseUtils {

public:

    /*!
    * @brief   Returns the maximum coefficient in the input SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support comparison operations.
    * @param   mat: The Eigen::SparseMatrix<T> for which the maximum coefficient will be determined.
    * @return  T: The maximum coefficient in the SparseMatrix. If the matrix is empty, it will return a default-constructed value of type T.
    *
    * Notes:
    * 1. This function first evaluates if the input SparseMatrix is populated. If mat.size() == 0, it immediately returns a default-constructed value of type T.
    * 2. If the input SparseMatrix is populated, the function will cast it to a DenseMatrix using Eigen's ternary matrix conversion mechanism. After the matrix is cast, it will proceed to find the maximum coefficient.
    * 3. The performance will rely on the number of elements in the matrix, as the maxCoeff function from Eigen's library performs a full pass through the matrix coefficients.
    */
    template <typename T>
    static T max(const Eigen::SparseMatrix<T> &mat){

        if(mat.size() == 0){
            return T();
        }

        const Eigen::MatrixX<T> castMat = mat;
        return castMat.maxCoeff();
    }

    /*!
    * @brief   Returns the maximum coefficient of the input SparseMatrix stored in RowMajor order.
    * @tparam  T: Numeric type used in SparseMatrix. It should support comparison operations.
    * @param   vec: The RowMajor Eigen::SparseMatrix<T> for which the maximum coefficient will be determined.
    * @return  T: The maximum coefficient in the SparseMatrix. If the matrix is empty, it will return a default-constructed value of type T.
    *
    * Note: This function first re-interprets the RowMajor SparseMatrix as a ColumnMajor SparseMatrix, which is the default storage order in Eigen. It then calls the max function on this reinterpretation to find the maximum coefficient. The performance will depend on the number of non-zero coefficients in the matrix.
    */
    template <typename T>
    static T max(const Eigen::SparseMatrix<T, Eigen::RowMajor> &vec){

        const Eigen::SparseMatrix<T> &matColMajor = vec;
        return max(matColMajor);
    }

    /*!
    * @brief   Returns the maximum coefficient of the input SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support comparison operations.
    * @param   vec: The SparseVector from the Eigen linear algebra library, for which the maximum coefficient will be determined.
    * @return  T: The maximum coefficient in the SparseVector. If the vector is empty, it will return a default-constructed value of type T.
    *
    * Note: This function converts the given SparseVector to a SparseMatrix (which is the default storage order in Eigen) before calling the max function to find the maximum coefficient. The performance will depend on the number of non-zero coefficients in the vector.
    */
    template <typename T>
    static T max(const Eigen::SparseVector<T> &vec){

        const Eigen::SparseMatrix<T> &matrixColMajor = vec;
        return max(matrixColMajor);
    }

    /*!
    * @brief   Returns the maximum coefficient of the input RowMajor SparseVector.
    * @tparam  T: Numeric type used in SparseMatrix. It should support comparison operations.
    * @param   vec: The RowMajor Eigen::SparseVector<T> for which the maximum coefficient will be determined.
    * @return  T: The maximum coefficient in the SparseVector. If the vector is empty, it will return a default-constructed value of type T.
    *
    * Note: This function first transposes the given RowMajor SparseVector to form a SparseMatrix. It then calls the max function to find the maximum coefficient. The performance will depend on the number of non-zero coefficients in the vector.
    */
    template <typename T>
    static T max(const Eigen::SparseVector<T, Eigen::RowMajor> &vec){

        const Eigen::SparseMatrix<T, Eigen::RowMajor> &matRowMajor = vec.transpose();
        return max(matRowMajor);
    }


    /*!
    * @brief   Returns the minimum coefficient in the input SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support comparison operations.
    * @param   mat: The Eigen::SparseMatrix<T> for which the minimum coefficient will be determined.
    * @return  T: The minimum coefficient in the SparseMatrix. If the matrix is empty, it will return a default-constructed value of type T.
    *
    * Note: This function first checks if the input SparseMatrix is populated. If mat.size() == 0, it immediately returns a default-constructed value of type T. If the SparseMatrix is populated, the function will cast it to a DenseMatrix using Eigen's ternary matrix conversion mechanism. After the matrix is cast, it will proceed to find the minimum coefficient. The performance will depend on the matrix's number of elements as the minCoeff function performs a full pass through the matrix coefficients.
    */
    template <typename T>
    static T min(const Eigen::SparseMatrix<T> &mat){

        if(mat.size() == 0){
            return T();
        }

        const Eigen::MatrixX<T> castMat = mat;
        return castMat.minCoeff();
    }

    /*!
    * @brief   Returns the minimum coefficient of the input RowMajor SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support comparison operations.
    * @param   vec: The RowMajor Eigen::SparseMatrix<T> for which the minimum coefficient will be determined.
    * @return  T: The minimum coefficient in the SparseMatrix. If the matrix is empty, it will return a default-constructed value of type T.
    *
    * Note: This function first re-interprets the RowMajor SparseMatrix as a ColumnMajor SparseMatrix, which is the default storage order in Eigen. It then calls the min function on this reinterpretation to find the minimum coefficient. The performance will depend on the number of non-zero coefficients in the matrix.
    */
    template <typename T>
    static T min(const Eigen::SparseMatrix<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T> &matColMajor = vec;
        return min(matColMajor);
    }

    /*!
    * @brief   Returns the minimum coefficient of the input SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support comparison operations.
    * @param   vec: The Eigen::SparseVector<T> for which the minimum coefficient will be determined.
    * @return  T: The minimum coefficient in the SparseVector. If the vector is empty, it will return a default-constructed value of type T.
    *
    * Note: This function first converts the given SparseVector to a SparseMatrix (which is the default storage order in Eigen) before calling the min function to find the minimum coefficient. The performance will depend on the number of non-zero coefficients in the vector.
    */
    template <typename T>
    static T min(const Eigen::SparseVector<T> &vec){
        const Eigen::SparseMatrix<T> &matrixColMajor = vec;
        return min(matrixColMajor);
    }

    /*!
    * @brief   Returns the minimum coefficient of the input RowMajor SparseVector.
    * @tparam  T: Numeric type used in SparseMatrix. It should support comparison operations.
    * @param   vec: The RowMajor Eigen::SparseVector<T> for which the minimum coefficient will be determined.
    * @return  T: The minimum coefficient in the SparseVector. If the vector is empty, it will return a default-constructed value of type T.
    *
    * Note: This function first transposes the given RowMajor SparseVector to form a SparseMatrix. It then calls the min function to find the minimum coefficient. The performance will depend on the number of non-zero coefficients in the vector.
    */
    template <typename T>
    static T min(const Eigen::SparseVector<T, Eigen::RowMajor> &vec) {
        const Eigen::SparseMatrix<T, Eigen::RowMajor> &matRowMajor = vec.transpose();
        return min(matRowMajor);
    }


    /*!
    * @brief   Returns the mean of the non-zero coefficients in the input SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support operations of addition and division.
    * @param   mat: The Eigen::SparseMatrix<T> for which the mean will be computed.
    * @return  double: The mean value of the non-zero coefficients in the SparseMatrix. If the matrix is empty (i.e., has zero size, or all of its coefficients are zero), it will return a default-constructed value of double.
    *
    * Note: This function first checks if the input matrix has any non-zero coefficients. If it does, it compresses the input SparseMatrix to ensure that it is in compressed form (which allows the function to skip over zero entries during computation), then casts the sum of the non-zero coefficients to double, and divides it by the count of non-zero coefficients. The performance will depend on the number of non-zero coefficients in the matrix.
    */
    template<typename T>
    static double mean(const Eigen::SparseMatrix<T> &mat){
        const int nonZeros = mat.nonZeros();

        if(mat.size() == 0 || nonZeros == 0){
            return double();
        }

        Eigen::SparseMatrix<T> copy = mat;
        if (!copy.isCompressed()) {
            copy.makeCompressed();
        }

        const auto coeffsSum = static_cast<double>(copy.coeffs().sum());
        return coeffsSum / nonZeros;
    }

    /*!
    * @brief   Returns the mean of the non-zero coefficients of the input RowMajor SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support operations of addition and division.
    * @param   vec: The RowMajor Eigen::SparseMatrix<T> for which the mean will be computed.
    * @return  double: The mean value of the non-zero coefficients in the SparseMatrix. If the matrix is empty (i.e., has zero size, or all of its coefficients are zero), it will return a default-constructed value of double.
    *
    * Note: This function first re-interprets the RowMajor SparseMatrix as a ColumnMajor SparseMatrix, which is the default storage order in Eigen. It then calls the mean function on this reinterpretation to compute the mean coefficient. The performance will depend on the number of non-zero coefficients in the matrix.
    */
    template <typename T>
    static double mean(const Eigen::SparseMatrix<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T> &matColMajor = vec;
        return mean(matColMajor);
    }

    /*!
    * @brief   Returns the mean of the non-zero coefficients of the input SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support operations of addition and division.
    * @param   vec: The Eigen::SparseVector<T> for which the mean will be computed.
    * @return  double: The mean value of the non-zero coefficients in the SparseVector. If the vector is empty (i.e., has zero size, or all of its coefficients are zero), it will return a default-constructed value of double.
    *
    * Note: This function first converts the given SparseVector to a SparseMatrix (which is the default storage order in Eigen) before calling the mean function to compute the mean coefficient. The performance will depend on the number of non-zero coefficients in the vector.
    */
    template <typename T>
    static double mean(const Eigen::SparseVector<T> &vec){
        const Eigen::SparseMatrix<T> &matrixColMajor = vec;
        return mean(matrixColMajor);
    }

    /*!
    * @brief   Returns the mean of the non-zero coefficients of the input SparseVector represented in RowMajor format.
    * @tparam  T: Numeric type used in SparseVector. It should support operations of addition and division.
    * @param   vec: The SparseVector from the Eigen library, for which the mean will be computed.
    * @return  double: The mean value of the non-zero coefficients in the SparseVector. If the vector is empty (i.e., has zero size, or all of its coefficients are zero), it will return a default-constructed value of double.
    *
    * Note: This function first transposes the input RowMajor SparseVector to form a RowMajor SparseMatrix. It then calls the mean function to compute the mean of non-zero coefficients. The performance will depend on the number of non-zero coefficients in the vector.
    */
    template <typename T>
    static double mean(const Eigen::SparseVector<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T, Eigen::RowMajor> &matRowMajor = vec.transpose();
        return mean(matRowMajor);
    }

    /*!
    * @brief   Returns the standard deviation of the non-zero coefficients in the input SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support operations of addition, subtraction, multiplication, and division.
    * @param   mat: The Eigen::SparseMatrix<T> for which the standard deviation will be computed.
    * @return  double: The standard deviation of the non-zero coefficients in the SparseMatrix. If the matrix is empty (i.e., has zero size, or all of its coefficients are zero), it will return a default-constructed value of double.
    *
    * Note: This function first checks if the input matrix has any non-zero coefficients. If it does, it compresses the input SparseMatrix to ensure that it is in compressed form (which allows the function to skip over zero entries during computation), calculates the mean of the non-zero coefficients, computes the sum of the squared differences from the mean, and finally takes the square root of this quantity divided by the count of non-zero coefficients to compute the standard deviation. The performance will depend on the number of non-zero coefficients in the matrix.
    */
    template <typename T>
    static double stDev(const Eigen::SparseMatrix<T> &mat){
        const int nonZeros = mat.nonZeros();

        if(mat.size() == 0 || nonZeros == 0){
            return double();
        }

        Eigen::SparseMatrix<T> copy = mat;
        if (!copy.isCompressed()) {
            copy.makeCompressed();
        }

        const double meanOfVec = mean(copy);
        const Eigen::VectorXd diffVec = copy.coeffs().array().template cast<double>() - meanOfVec;
        return std::sqrt(diffVec.cwiseProduct(diffVec).sum() / nonZeros);
    }

    /*!
    * @brief   Returns the standard deviation of the non-zero coefficients of the input RowMajor SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support operations of addition, subtraction, multiplication and division.
    * @param   vec: The RowMajor Eigen::SparseMatrix<T> for which the standard deviation will be computed.
    * @return  double: The standard deviation of the non-zero coefficients in the SparseMatrix. If the matrix is empty (i.e., has zero size, or all of its coefficients are zero), it will return a default-constructed value of double.
    *
    * Note: This function first re-interprets the RowMajor SparseMatrix as a ColumnMajor SparseMatrix, which is the default storage order in Eigen. It then calls the stDev function on this reinterpretation to compute the standard deviation. The performance will depend on the number of non-zero coefficients in the matrix.
    */
    template <typename T>
    static double stDev(const Eigen::SparseMatrix<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T> &matColMajor = vec;
        return stDev(matColMajor);
    }

    /*!
    * @brief   Returns the standard deviation of the non-zero coefficients of the input SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support operations of addition, subtraction, multiplication, and division.
    * @param   vec: The Eigen::SparseVector<T> for which the standard deviation will be computed.
    * @return  double: The standard deviation of the non-zero coefficients in the SparseVector. If the vector is empty (i.e., has zero size, or all of its coefficients are zero), it will return a default-constructed value of double.
    *
    * Note: This function first converts the given SparseVector to a SparseMatrix (which is the default storage order in Eigen) before calling the stDev function to compute the standard deviation. The performance will depend on the number of non-zero coefficients in the vector.
    */
    template <typename T>
    static double stDev(const Eigen::SparseVector<T> &vec){
        const Eigen::SparseMatrix<T> &matrixColMajor = vec;
        return stDev(matrixColMajor);
    }

    /*!
    * @brief   Returns the standard deviation of the non-zero coefficients of the input RowMajor SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support operations of addition, subtraction, multiplication and division.
    * @param   vec: The RowMajor Eigen::SparseVector<T> for which the standard deviation will be computed.
    * @return  double: The standard deviation of the non-zero coefficients in the SparseVector. If the vector is empty (i.e., has zero size, or all of its coefficients are zero), it will return a default-constructed value of double.
    *
    * Note: This function first transposes the input RowMajor SparseVector to form a RowMajor SparseMatrix. It then calls the stDev function to compute the standard deviation. The performance will depend on the number of non-zero coefficients in the SparseVector.
    */
    template <typename T>
    static double stDev(const Eigen::SparseVector<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T, Eigen::RowMajor> &matRowMajor = vec.transpose();
        return stDev(matRowMajor);
    }

    /*!
    * @brief   Checks if given index is within the bounds of the input SparseVector.
    * @tparam  T: Numeric type used in SparseVector.
    * @param   vec: The Eigen::SparseVector<T> whose bounds will be checked against the provided index.
    * @param   index: The index that will be checked for its validity.
    * @return  bool: Returns true if the index lies within the bounds of the SparseVector, false otherwise.
    *
    * Note: This function evaluates the validity of the index by checking if it's non-negative and is less than the size of the input SparseVector. The function returns the logical result of this comparison.
    */
    template <typename T>
    static bool validIndex(
            const Eigen::SparseVector<T> &vec,
            int index
            ){
        return 0 <= index && index < vec.size();
    }

    /*!
    * @brief   Checks if the given index is within the bounds of the input RowMajor SparseVector.
    * @tparam  T: Numeric type used in SparseVector.
    * @param   vec: The RowMajor Eigen::SparseVector<T> whose bounds will be checked against the provided index.
    * @param   index: The index that will be checked for its validity.
    * @return  bool: Returns true if the index lies within the bounds of the SparseVector, false otherwise.
    *
    * Note: This function first reinterprets the RowMajor SparseVector as a SparseVector (default is column major in Eigen). It then calls the validIndex function to check if the provided index lies within the bounds of the reinterpretation. The performance will depend on the size of the SparseVector.
    */
    template <typename T>
    static bool validIndex(
            const Eigen::SparseVector<T, Eigen::RowMajor> &vec,
            int index
            ){
        const Eigen::SparseVector<T> &vecRowMajor = vec;
        return validIndex(vecRowMajor, index);
    }

    /*!
    * @brief   Checks if given row and column indices are within the bounds of the input SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix.
    * @param   mat: The Eigen::SparseMatrix<T> whose bounds will be checked against the provided row and column indices.
    * @param   row: The row index that will be checked for its validity.
    * @param   col: The column index that will be checked for its validity.
    * @return  bool: Returns true if both indices lie within the bounds of the SparseMatrix, false otherwise.
    *
    * Note: This function first checks if the row index is non-negative and is less than the number of rows in the SparseMatrix. It also checks if the column index is non-negative and is less than the number of columns in the SparseMatrix. The function returns true if both these conditions are met, otherwise it returns false.
    */
    template <typename T>
    static bool validIndex(
            const Eigen::SparseMatrix<T> &mat,
            int row,
            int col
            ){
        return (0 <= row && row < mat.rows())
            && (0 <= col && col < mat.cols());
    }

    /*!
    * @brief   Checks if given row and column indices are within the bounds of the input RowMajor SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix.
    * @param   mat: The RowMajor Eigen::SparseMatrix<T> whose bounds will be checked against the provided row and column indices.
    * @param   row: The row index that will be checked for its validity.
    * @param   col: The column index that will be checked for its validity.
    * @return  bool: Returns true if both indices lie within the bounds of the SparseMatrix, false otherwise.
    *
    * Note: This function first re-interprets the RowMajor SparseMatrix as a ColumnMajor SparseMatrix (which is the default storage order in Eigen). It then calls the validIndex function to check if the provided row and column indices lie within the bounds of the reinterpretation.
    */
    template <typename T>
    static bool validIndex(
            const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat,
            int row,
            int col
            ) {
        const Eigen::SparseMatrix<T> &matRowMajor = mat;
        return validIndex(matRowMajor,
                          row,
                          col);
    }

    /*!
    * @brief   Returns the median of the non-zero coefficients of the input SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support operations of addition, subtraction, multiplication, and division.
    * @param   vec: The Eigen::SparseVector<T> for which the median will be computed.
    * @return  double: The median value of the non-zero coefficients in the SparseVector. For an even number of non-zero coefficients, the average of the two middle elements is returned. For an odd number, the middle element is returned.
    *
    * Note: This function first converts the SparseVector coefficients to an Eigen::VectorXd and cast them to double. It then moves the data of the VectorXd to a std::vector of double and computes the median using the MathUtils::median() function. The performance will depend on the number of non-zero coefficients in the SparseVector.
    */
    template <typename T>
    static double median(const Eigen::SparseVector<T> &vec){
        const Eigen::VectorXd v = vec.coeffs().template cast<double>();
        std::vector<double> qvec(v.data(), v.data() + v.size());
        return MathUtils::median(qvec);
    }

    /*!
    * @brief   Returns the median of the non-zero coefficients of the input RowMajor SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support operations of addition, subtraction, multiplication, and division.
    * @param   vec: The RowMajor Eigen::SparseVector<T> for which the median will be computed.
    * @return  double: The median value of the non-zero coefficients in the SparseVector. For an even number of non-zero coefficients, the average of the two middle elements is returned. For an odd number, the middle element is returned.
    *
    * Note: This function first transposes the RowMajor SparseVector to form a SparseVector (default storage order in Eigen is column-major). It then calls the median function to compute the median coefficient. The performance will depend on the number of non-zero coefficients in the SparseVector.
    */
    template <typename T>
    static double median(const Eigen::SparseVector<T, Eigen::RowMajor> &vec){
        const Eigen::SparseVector<T> &matRowMajor = vec.transpose();
        return median(matRowMajor);
    }

    /*!
    * @brief   Returns the median of the non-zero coefficients in the input SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support operations of addition, subtraction, multiplication, and division.
    * @param   mat: The Eigen::SparseMatrix<T> for which the median will be computed.
    * @return  double: The median value of the non-zero coefficients in the SparseMatrix. For an even number of non-zero coefficients, the average of the two middle elements is returned. For an odd number, the middle element is returned.
    *
    * Note: This function first checks if the SparseMatrix is in compressed form, if it's not, it compresses it. It then converts its coefficients to an Eigen::VectorXd, casts them to double, and transfers the data of the VectorXd to a std::vector. The median is then calculated and returned using the MathUtils::median() function. The performance will depend on the number of non-zero coefficients in the SparseMatrix.
    */
    template <typename T>
    static double median(const Eigen::SparseMatrix<T> &mat){
        Eigen::SparseMatrix<T> copy = mat;
        if (!copy.isCompressed()) {
            copy.makeCompressed();
        }

        const Eigen::VectorXd v = copy.coeffs().template cast<double>();
        std::vector<double> qvec(v.data(), v.data() + v.size());
        return MathUtils::median(qvec);
    }

    /*!
    * @brief   Returns the median of the non-zero coefficients of the input RowMajor SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support  operations of addition, subtraction, multiplication, and division.
    * @param   mat: The RowMajor Eigen::SparseMatrix<T> for which the median will be computed.
    * @return  double: The median value of the non-zero coefficients in the SparseMatrix. For an even number of non-zero coefficients, the average of the two middle elements is returned. For an odd number, the middle element is returned.
    *
    * Note: This function first re-interprets the RowMajor SparseMatrix as a ColumnMajor SparseMatrix, which is the default storage order in Eigen. It then calls the median function on this reinterpretation to compute the median coefficient. The performance will depend on the number of non-zero coefficients in the SparseMatrix.
    */
    template <typename T>
    static double median(const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat){
        const Eigen::SparseMatrix<T> &matColMajor = mat;
        return median(matColMajor);
    }

    /*!
    * @brief   Removes elements below a certain threshold value from the input SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support operations of comparison (<).
    * @param   thresholdValue: The value below which (exclusive) elements in the SparseVector will be removed.
    * @param   vec: Pointer to Eigen::SparseVector<T>, where elements below a certain threshold value will be removed.
    *
    * Note: This function creates a new SparseVector of the same size as vec. It then iterates over each element in the original vector and only adds it to the new vector if its value is above the thresholdValue. Finally, it assigns the new SparseVector to the original. Performance will depend on the number of elements found above the thresholdValue.
    */
    template <typename T>
    static void removeElementsBelowThreshold(
            double thresholdValue,
            Eigen::SparseVector<T> *vec
            ){

        Eigen::SparseVector<T> newVec(vec->size());
        for (typename Eigen::SparseVector<T>::InnerIterator it(*vec); it; ++it) {

            if (it.value() < thresholdValue) {
                continue;
            }

            newVec.coeffRef(it.index()) = it.value();
        }

        *vec = newVec;
    }

    /*!
    * @brief   Removes elements below a certain threshold value from the input RowMajor SparseVector.
    * @tparam  T: Numeric type used in SparseVector. It should support operations of comparison (<).
    * @param   thresholdValue: The value below which (exclusive) elements in the SparseVector will be removed.
    * @param   vec: Pointer to RowMajor Eigen::SparseVector<T>, where elements below a certain threshold value will be removed.
    *
    * Note: This function creates a new SparseVector (default column-major storage order in Eigen) of the same size as vec. It then iterates using InnerIterator over each element in the original vector (in the order of inner dimension) and only adds it to the new vector if its value is above the thresholdValue. Finally, it assigns the new vector to the original, converting it to the default storage order (column-major). Performance will depend on the number of elements found above the thresholdValue.
    */
    template <typename T>
    static void removeElementsBelowThreshold(
            double thresholdValue,
            Eigen::SparseVector<T, Eigen::RowMajor> *vec
            ){

        Eigen::SparseVector<T> newVec(vec->size());
        for (typename Eigen::SparseVector<T, Eigen::RowMajor>::InnerIterator it(*vec); it; ++it) {

            if (it.value() < thresholdValue) {
                continue;
            }

            newVec.coeffRef(it.index()) = it.value();
        }

        *vec = newVec;
    }

    /*!
    * @brief   Removes elements below a certain threshold value from the input SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support operations of comparison (<).
    * @param   thresholdValue: The value below which (exclusive) elements in the SparseMatrix will be removed.
    * @param   mat: Pointer to Eigen::SparseMatrix<T>, where elements below a certain threshold value will be removed.
    *
    * Note: This function creates a new SparseMatrix of the same size as mat. It then iterates over each element in the original matrix and only adds it to the new matrix if its value is above the thresholdValue. Finally, it assigns the new SparseMatrix to the original. Performance will depend on the number of elements found above the thresholdValue.
    */
    template <typename T>
    static void removeElementsBelowThreshold(
            double thresholdValue,
            Eigen::SparseMatrix<T> *mat
            ){

        Eigen::SparseMatrix<T> newMat(mat->rows(), mat->cols());

        for (int i = 0; i < mat->outerSize(); ++i) {

            for (typename Eigen::SparseMatrix<T>::InnerIterator it(*mat, i); it; ++it){

                if (it.value() < thresholdValue) {
                    continue;
                }

                newMat.coeffRef(it.row(), it.col()) = it.value();
            }
        }

        *mat = newMat;
    }

    /*!
    * @brief   Removes elements below a certain threshold value from the RowMajor SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix. It should support operations of comparison (<).
    * @param   thresholdValue: The value less than which (exclusive) elements in the SparseMatrix will be removed.
    * @param   mat: Pointer to Eigen::SparseMatrix<T, Eigen::RowMajor>, elements less than the threshold will be removed from this matrix.
    *
    * Note: This function creates a new SparseMatrix of the same size as mat. The function then iterates over each element in the original matrix using InnerIterator and adds it to the new matrix only if the value is greater than the thresholdValue. In the end, it replaces the original matrix with the new matrix having removed elements below the set threshold. Performance depends on the number of elements found above the thresholdValue.
    */
    template <typename T>
    static void removeElementsBelowThreshold(
            T thresholdValue,
            Eigen::SparseMatrix<T, Eigen::RowMajor> *mat
            ){

        Eigen::SparseMatrix<T, Eigen::RowMajor> newMat(mat->rows(), mat->cols());

        for (int i = 0; i < mat->outerSize(); ++i) {

            for (typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(*mat, i); it; ++it){

                if (it.value() < thresholdValue) {
                    continue;
                }

                newMat.coeffRef(it.row(), it.col()) = it.value();
            }
        }

        *mat = newMat;
    }

    /*!
    * @brief   Rolls elements of the input SparseVector by the provided roll distance.
    * @tparam  T: Numeric type used in SparseVector.
    * @param   vec: The Eigen::SparseVector<T> where elements will be rolled.
    * @param   rollDistance: Integer roll distance, positive values will roll elements to the right, negative values will roll them to the left.
    * @return  Eigen::SparseVector<T>: New SparseVector where elements have been rotated by rollDistance.
    *
    * Note: This function first creates a new SparseVector of the same size as vec. It then iterates over every element in the input SparseVector. For each element, it calculates a "rolled" index which is the sum of the current index and the roll distance. It then verifies with the validIndex function if the rolled index lies within the bounds of SparseVector. If true, it assigns the value at the current index to the rolled index in the new SparseVector. This function returns the new rolled SparseVector.
    */
    template <typename T>
    static Eigen::SparseVector<T> roll(
            const Eigen::SparseVector<T> &vec,
            int rollDistance
            ){

        Eigen::SparseVector<T> returnVec(vec.size());

        for (typename Eigen::SparseVector<T>::InnerIterator it(vec); it; ++it) {

            const int rolledIndex = it.index() + rollDistance;

            if (validIndex(returnVec, rolledIndex)) {

                returnVec.coeffRef(rolledIndex) = it.value();
            }
        }

        return returnVec;
    }

    /*!
    * @brief   Rolls elements of the input RowMajor SparseVector by the provided roll distance.
    * @tparam  T: Numeric type used in SparseVector.
    * @param   vec: The RowMajor Eigen::SparseVector<T> where elements will be rolled.
    * @param   rollDistance: Integer roll distance, positive values will roll elements to the right, negative values will roll them to the left.
    * @return  Eigen::SparseVector<T>: New SparseVector (default is column-major in Eigen) where elements have been rotated by rollDistance.
    *
    * Note: This function first creates a new SparseVector (default is column-major in Eigen) of the same size as vec. It then iterates over every element in the input SparseVector. For each element, it calculates a "rolled" index which is the sum of the current index and the roll distance. It then verifies with the validIndex function if the rolled index lies within the bounds of the SparseVector. If true, it assigns the value at the current index to the rolled index in the new SparseVector. This function returns the new rolled SparseVector.
    */
    template <typename T>
    static Eigen::SparseVector<T> roll(
            const Eigen::SparseVector<T, Eigen::RowMajor> &vec,
            int rollDistance
            )  {

        Eigen::SparseVector<T> returnVec(vec.size());

        for (typename Eigen::SparseVector<T, Eigen::RowMajor>::InnerIterator it(vec); it; ++it) {

            const int rolledIndex = it.index() + rollDistance;

            if (validIndex(returnVec, rolledIndex)) {

                returnVec.coeffRef(rolledIndex) = it.value();
            }
        }

        return returnVec;
    }

    class SparseMatrixPoint{
    public:
        int row = -1;
        int col = -1;
        double value = 0;
    };

    /*!
    * @brief   Obtains the apexes of values in the input SparseVector.
    * @tparam  T: Numeric type used in SparseVector.
    * @param   vec: The Eigen::SparseVector<T> from which the apex values will be obtained.
    * @param   precision: The precision level for checking the equality of values, default value is 1e4.
    * @return  QMap<int, T>: QMap containing indices where the apex values are located and their corresponding values.
    *
    * Note: This function first gets the size of SparseVector, then it iterates over each element in the input SparseVector. For each element, it checks if it is a local apex (greater than or equal to its surrounding values when multiplied by the precision value). The local apex indices and their respective values are collected in a QMap and returned.
    */
    template <typename T>
    static QMap<int, T> apexes(
            const Eigen::SparseVector<T> &vec,
            int precision = 1e4
                    ){

        QMap<int, T> apexIndicies;

        const int vecSize = static_cast<int>(vec.size());
        for (typename Eigen::SparseVector<T>::InnerIterator it(vec); it; ++it) {

            const int index = it.index();
            if(index < 1 || index >= vecSize - 1){
                continue;
            }

            const auto leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));
            const auto rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));
            const auto centerPointValue = static_cast<int>(std::round(it.value() * precision));

            if(centerPointValue > leftPointValue && centerPointValue >= rightPointValue){
                apexIndicies.insert(index,  it.value());
            }
        }

        return apexIndicies;
    }

    // helper method apexes(const Eigen::SparseMatrix<T> &mat)
    static void sortApexesHiLoValue(QVector<SparseMatrixPoint> *apxs);

    // helper method apexes(const Eigen::SparseMatrix<T> &mat)
    template <typename T>
    static QMap<int, QMap<int, T> >  matrixApexesByColumn(const Eigen::SparseMatrix<T> &mat){
        QMap<int, QMap<int, T> > columnApexes;
        QSet<int> visitedColumns;
        for (int i = 0; i < mat.outerSize(); ++i) {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it) {

                const int column = it.col();

                if(visitedColumns.contains(column)){
                    continue;
                }
                visitedColumns.insert(column);

                const Eigen::SparseVector<T> colVec = mat.col(it.col());
                const QMap<int, T> colApexes = apexes(colVec);

                columnApexes.insert(column, colApexes);
            }
        }

        return columnApexes;
    }

    /*!
    * @brief   Obtains the apexes of values in the input SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix.
    * @param   mat: The Eigen::SparseMatrix<T> from which the apex values will be obtained.
    * @param   precision: The precision level for checking the equality of values. The default value is 1e4.
    * @return  QVector<SparseMatrixPoint>: Vector containing points where the apex values are located and their corresponding values.
    *
    * Note: This function first gets the apexes by column from the SparseMatrix. It then iterates over each column apex and every row coefficient in that column. For each coefficient, it checks if it is a local apex (greater than or equal to its surrounding values when multiplied by the precision value). The local apex points (row and column indices along with their respective values) are collected in a QVector and returned. The function also adjusts the edge buffer of the matrix and sorts the apexes in a high-to-low value order.
    */
    template <typename T>
    static QVector<SparseMatrixPoint> apexes(
            const Eigen::SparseMatrix<T> &mat,
            int precision = 1e4
                    ){

        const int matrixEdgeBuffer = 2;

        const QMap<int, QMap<int, T> > columnApexes = matrixApexesByColumn(mat);

        QVector<SparseMatrixPoint> returnApexes;
        for (auto it = columnApexes.begin(); it != columnApexes.end(); it++){

            int column = it.key();
            column = std::max(1, column);
            column = std::min(column, static_cast<int>(mat.cols()) - matrixEdgeBuffer);

            const  QMap<int, double> colApexes = it.value();
            for(auto it2 = colApexes.begin(); it2!= colApexes.end(); it2++){

                int row = it2.key();
                row = std::max(1, row);
                row = std::min(row, static_cast<int>(mat.rows()) - matrixEdgeBuffer);

                const int colApex = static_cast<int>(std::round(it2.value() * precision));

                const int left = static_cast<int>(std::round(mat.coeff(row , column - 1) * precision));
                const int right = static_cast<int>(std::round(mat.coeff(row, column + 1) * precision));
                const int leftUp = static_cast<int>(std::round(mat.coeff(row - 1, column - 1) * precision));
                const int rightUp = static_cast<int>(std::round(mat.coeff(row - 1, column + 1) * precision));
                const int leftDown = static_cast<int>(std::round(mat.coeff(row + 1, column - 1) * precision));
                const int rightDown = static_cast<int>(std::round(mat.coeff(row + 1, column + 1) * precision));

                if(colApex > left
                    && colApex >= right
                    && colApex > leftUp
                    && colApex >= rightUp
                    && colApex > leftDown
                    && colApex >= rightDown
                ){
                    SparseMatrixPoint a;
                    a.row = row;
                    a.col = column;
                    a.value = static_cast<double>(it2.value());
                    returnApexes.push_back(a);
                }
            }
        }

        sortApexesHiLoValue(&returnApexes);
        return returnApexes;
    }

    /*!
    * @brief   Obtains the apexes of values in the input RowMajor SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix.
    * @param   mat: The RowMajor Eigen::SparseMatrix<T> from which the apex values will be obtained.
    * @param   precision: The precision level for checking the equality of values. The default value is 1e4.
    * @return  QVector<SparseMatrixPoint>: Vector containing points where the apex values are located and their corresponding values.
    *
    * Note: This function first re-interprets the RowMajor SparseMatrix as a ColumnMajor SparseMatrix (default storage order in Eigen). It then uses apexes() function to calculate apexes from the re-interpreted SparseMatrix. The performance will depend on the number of apexes found in the SparseMatrix.
    */
    template <typename T>
    static QVector<SparseMatrixPoint> apexes(
            const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat,
            int precision = 1e4
                    ){
        const Eigen::SparseMatrix<T>  matConverted = mat.template cast<T>();
        return apexes(matConverted, precision);
    }

    /*!
    * @brief   Builds a comb filter in the form of a RowMajor SparseMatrix.
    * @tparam  T: Numeric type used in SparseMatrix and QVector. It should support operations of addition, subtraction, comparison, and casting to T2.
    * @tparam  T2: Numeric type used in precision.
    * @param   teethValues: QVector<T> that defines the value of teeth in comb filter.
    * @param   tolerance: The tolerance level for creating a range around each tooth in comb filter.
    * @param   maxValue: The maximum possible value after hashing. It sets the column size of output SparseMatrix.
    * @param   precision: The precision level for hashing.
    * @return  Eigen::SparseMatrix<T, Eigen::RowMajor>: RowMajor SparseMatrix representing the output comb filter.
    *
    * Note: This function first calculates hashed tolerance and max value using MathUtils::hashDecimal() method. It then creates a SparseMatrix of row size equal to teethValues.size() and column size equal to hashedMaxValue. Next, it iterates over teethValues and for each value, it calculates a range (hashedVal - hashedTolerance, hashedVal + hashedTolerance) and inserts `1.0` in output SparseMatrix for that row in range. Finally it returns the SparseMatrix.
    */
    template <typename T, typename T2>
    static Eigen::SparseMatrix<T, Eigen::RowMajor> buildCombFilter(
            const QVector<T> &teethValues,
            T tolerance,
            T maxValue,
            T2 precision
            ) {

        const int hashedTolerance = MathUtils::hashDecimal(tolerance, precision);
        const int hashedMaxValue = MathUtils::hashDecimal(maxValue, precision);

        Eigen::SparseMatrix<T, Eigen::RowMajor> combFilter(teethValues.size(), hashedMaxValue);

        combFilter.setZero();

        for (int row = 0; row < teethValues.size(); row++) {

            const T val = teethValues.at(row);
            const int hashedVal = MathUtils::hashDecimal(val, precision);

            const int startCol = std::max(hashedVal - hashedTolerance, 0);
            const int endCol = std::min(hashedVal + hashedTolerance + 1, hashedMaxValue - 1);

            for(int col = startCol; col < endCol; col++) {

                combFilter.insert(row, col) = static_cast<T>(1.0);
            }
        }

        return combFilter;
    }

    template <typename Point, typename T2>
    static Err pointsToSparseVector(
            const QVector<Point> &points,
            T2 maxValue,
            int precision,
            Eigen::SparseVector<T2> *vec
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(points); ree;
        e = ErrorUtils::isAboveThreshold(
                maxValue,
                static_cast<T2>(0.0),
                ErrorUtilsParam::ExcludeThreshold
                ); ree

        e = ErrorUtils::isAboveThreshold(
                precision,
                0,
                ErrorUtilsParam::ExcludeThreshold
                ); ree

        const int vecSize = MathUtils::hashDecimal(maxValue, precision);
        vec->resize(vecSize);
        vec->setZero();

        for (const Point &p : points) {

            const int xHashed = MathUtils::hashDecimal(p.x(), precision);
            if (xHashed < 0 || xHashed >= vecSize) {
                continue;
            }

            vec->coeffRef(xHashed) += p.y();
        }

        ERR_RETURN
    }

    template <typename Point, typename T>
    static Err sparseVectorToPoints(
            const Eigen::SparseVector<T> &vec,
            int precision,
            QVector<Point> *points

    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(vec.size() > 0); ree;
        e = ErrorUtils::isAboveThreshold(
                precision,
                0,
                ErrorUtilsParam::ExcludeThreshold
        ); ree

        for (typename Eigen::SparseVector<T>::InnerIterator it(vec); it; ++it) {
            const T xUnHashed = MathUtils::unHashDecimal<T>(it.index(), precision);
            points->push_back({xUnHashed, it.value()});
        }

        ERR_RETURN
    }


};


#endif //EIGENSPARSEUTILS_H
