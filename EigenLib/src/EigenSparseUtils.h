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
    * @brief return maximum of sparse data structure
    */
    template <typename T>
    static T max(const Eigen::SparseMatrix<T> &mat){

        if(mat.size() == 0){
            return T();
        }

        const Eigen:MatrixX<T> castVec = mat;
        return castVec.maxCoeff();
    }


    template <typename T>
    static T max(const Eigen::SparseMatrix<T, Eigen::RowMajor> &vec){

        const Eigen::SparseMatrix<T> &matColMajor = vec;
        return max(matColMajor);
    }


    template <typename T>
    static T max(const Eigen::SparseVector<T> &vec){

        const Eigen::SparseMatrix<T> &matrixColMajor = vec;
        return max(matrixColMajor);
    }


    template <typename T>
    static T max(const Eigen::SparseVector<T, Eigen::RowMajor> &vec){

        const Eigen::SparseMatrix<T, Eigen::RowMajor> &matRowMajor = vec.transpose();
        return max(matRowMajor);
    }


    /*!
    * @brief return minimum value of sparse data structure
    */
    template <typename T>
    static T min(const Eigen::SparseMatrix<T> &mat){

        if(mat.size() == 0){
            return T();
        }

        const Eigen::MatrixX<T> castVec = mat;
        return castVec.minCoeff();
    }


    template <typename T>
    static T min(const Eigen::SparseMatrix<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T> &matColMajor = vec;
        return min(matColMajor);
    }


    template <typename T>
    static T min(const Eigen::SparseVector<T> &vec){
        const Eigen::SparseMatrix<T> &matrixColMajor = vec;
        return min(matrixColMajor);
    }


    template <typename T>
    static T min(const Eigen::SparseVector<T, Eigen::RowMajor> &vec) {
        const Eigen::SparseMatrix<T, Eigen::RowMajor> &matRowMajor = vec.transpose();
        return min(matRowMajor);
    }


    /*!
    * @brief Calculates mean.
    */
    template<typename T>
    static double mean(const Eigen::SparseMatrix<T> &mat){
        if(mat.size() == 0){
            return double();
        }

        const T coeffsSum = mat.coeffs().sum();
        return coeffsSum / static_cast<double>(mat.nonZeros());
    }


    template <typename T>
    static double mean(const Eigen::SparseMatrix<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T> &matColMajor = vec;
        return mean(matColMajor);
    }


    template <typename T>
    static double mean(const Eigen::SparseVector<T> &vec){
        const Eigen::SparseMatrix<T> &matrixColMajor = vec;
        return mean(matrixColMajor);
    }


    template <typename T>
    static double mean(const Eigen::SparseVector<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T, Eigen::RowMajor> &matRowMajor = vec.transpose();
        return mean(matRowMajor);
    }


    /*!
    * @brief Calculates standard deviation.
    */
    template <typename T>
    static double stDev(const Eigen::SparseMatrix<T> &mat){

        if(mat.size() == 0){
            return double();
        }

        const double meanOfVec = mean(mat);
        const Eigen::VectorXd diffVec = mat.coeffs().array().template cast<double>() - meanOfVec;
        return std::sqrt(diffVec.cwiseProduct(diffVec).sum() / static_cast<double>(mat.nonZeros()));
    }


    template <typename T>
    static double stDev(const Eigen::SparseMatrix<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T> &matColMajor = vec;
        return stDev(matColMajor);
    }


    template <typename T>
    static double stDev(const Eigen::SparseVector<T> &vec){
        const Eigen::SparseMatrix<T> &matrixColMajor = vec;
        return stDev(matrixColMajor);
    }


    template <typename T>
    static double stDev(const Eigen::SparseVector<T, Eigen::RowMajor> &vec){
        const Eigen::SparseMatrix<T, Eigen::RowMajor> &matRowMajor = vec.transpose();
        return stDev(matRowMajor);
    }


    /*!
    * @brief Checks to see if an index is in range of sparse vector.
    */
    template <typename T>
    static bool validIndex(const Eigen::SparseVector<T> &vec,
                           int index){
        return 0 <= index && index < vec.size();
    }


    template <typename T>
    static bool validIndex(const Eigen::SparseVector<T, Eigen::RowMajor> &vec,
                           int index){
        const Eigen::SparseVector<T> &vecRowMajor = vec;
        return validIndex(vecRowMajor, index);
    }


    template <typename T>
    static bool validIndex(const Eigen::SparseMatrix<T> &mat,
                           int row,
                           int col){
        return (0 <= row && row < mat.rows())
            && (0 <= col && col < mat.cols());
    }


    template <typename T>
    static bool validIndex(const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat,
                           int row,
                           int col) {
        const Eigen::SparseMatrix<T> &matRowMajor = mat;
        return validIndex(matRowMajor,
                          row,
                          col);
    }


    /*!
    * @brief Calculates median of sparse data structures.
    */
    template <typename T>
    static double median(const Eigen::SparseVector<T> &vec){
        const Eigen::VectorXd v = vec.coeffs().template cast<double>();
        std::vector<double> qvec(v.data(), v.data() + v.size());
        return MathUtils::median(qvec);
    }


    template <typename T>
    static double median(const Eigen::SparseVector<T, Eigen::RowMajor> &vec){
        const Eigen::SparseVector<T> &matRowMajor = vec.transpose();
        return median(matRowMajor);
    }


    template <typename T>
    static double median(const Eigen::SparseMatrix<T> &mat){
        const Eigen::VectorXd v = mat.coeffs().template cast<double>();
        std::vector<double> qvec(v.data(), v.data() + v.size());
        return MathUtils::median(qvec);
    }


    template <typename T>
    static double median(const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat){
        const Eigen::SparseMatrix<T> &matColMajor = mat;
        return median(matColMajor);
    }


    /*!
     * @brief Returns a sparse vector w/ all values below thresholdValue removed.
     */
    template <typename T>
    static void removeElementsBelowThreshold(double thresholdValue,
                                             Eigen::SparseVector<T> *vec){
        Eigen::SparseVector<T> newVec(vec->size());
        for (typename Eigen::SparseVector<T>::InnerIterator it(*vec); it; ++it) {

            if (it.value() < thresholdValue) {
                continue;
            }

            newVec.coeffRef(it.index()) = it.value();
        }

        *vec = newVec;
    }


    template <typename T>
    static void removeElementsBelowThreshold(double thresholdValue,
                                             Eigen::SparseVector<T, Eigen::RowMajor> *vec){
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
    * @brief Returns a sparse matrix w/ all values below thresholdValue removed.
    */
    template <typename T>
    static void removeElementsBelowThreshold(double thresholdValue,
                                             Eigen::SparseMatrix<T> *mat){
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


    template <typename T>
    static void removeElementsBelowThreshold(T thresholdValue,
                                             Eigen::SparseMatrix<T, Eigen::RowMajor> *mat){
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
     * @brief Translates all values of indexes by the roll distance.  This function does not wrap
     * values over the size of the data structure if roll distance plus an index is greater.
     *
     */
    template <typename T>
    static Eigen::SparseVector<T> roll(const Eigen::SparseVector<T> &vec,
                                                                        int rollDistance){
        Eigen::SparseVector<T> returnVec(vec.size());

        for (typename Eigen::SparseVector<T>::InnerIterator it(vec); it; ++it) {

            const int rolledIndex = it.index() + rollDistance;

            if (validIndex(returnVec, rolledIndex)) {

                returnVec.coeffRef(rolledIndex) = it.value();
            }
        }

        return returnVec;
    }


    template <typename T>
    static Eigen::SparseVector<T> roll(const Eigen::SparseVector<T, Eigen::RowMajor> &vec,
                                       int rollDistance)  {

        Eigen::SparseVector<T> returnVec(vec.size());

        for (typename Eigen::SparseVector<T, Eigen::RowMajor>::InnerIterator it(vec); it; ++it) {

            const int rolledIndex = it.index() + rollDistance;

            if (validIndex(returnVec, rolledIndex)) {

                returnVec.coeffRef(rolledIndex) = it.value();
            }
        }

        return returnVec;
    }


    template <typename T>
    static Eigen::SparseVector<int> markVectorApexes(
            const Eigen::SparseVector<T> &vec,
            int precision = 1e4
    ) {

        const int vecSize = static_cast<int>(vec.size());

        Eigen::SparseVector<int> resultVec(vecSize);
        resultVec.setZero();

        for (typename Eigen::SparseVector<T>::InnerIterator it(vec); it; ++it) {

            const int index = it.index();

            int centerPointValue = -1;
            int leftPointValue = -1;
            int rightPointValue = -1;

            if (index < 1) {

                leftPointValue = static_cast<int>(std::round(it.value() * precision));
                rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));

                if (leftPointValue > rightPointValue) {
                    resultVec.insert(index) = 1;
                }

                continue;
            }

            else if (index >= vecSize - 1) {

                leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1)  * precision));
                rightPointValue = static_cast<int>(std::round(it.value() * precision));

                if (leftPointValue < rightPointValue) {
                    resultVec.insert(index) = 1;
                }

                continue;
            }

            centerPointValue = static_cast<int>(std::round(it.value() * precision));
            leftPointValue = static_cast<int>(std::round(vec.coeff(index - 1) * precision));
            rightPointValue = static_cast<int>(std::round(vec.coeff(index + 1)  * precision));


            if (centerPointValue > leftPointValue && centerPointValue >= rightPointValue) {
                resultVec.insert(index) = 1;
            }
        }

        return resultVec;
    }


    class SparseMatrixPoint{
    public:
        int row = -1;
        int col = -1;
        double value = 0;
    };


    template<typename T>
    static Eigen::SparseMatrix<int> markMatrixApexes(
            const Eigen::SparseMatrix<T> &mat,
            int precision = 1e4
    ) {

        //NOTE: This method is slower than the apexes() method below.  Figure out why.

        Eigen::SparseMatrix<int> resultMat(mat.rows(), mat.cols());

        QSet<int> visitedColumns;
        for (int i = 0; i < mat.outerSize(); ++i) {
            for (typename Eigen::SparseMatrix<T>::InnerIterator it(mat, i); it; ++it) {

                const int column = it.col();

                if(visitedColumns.contains(column)){
                    continue;
                }

                visitedColumns.insert(column);

                const Eigen::SparseVector<T> colVec = mat.col(column);
                const Eigen::SparseVector<int> markedColVec = markVectorApexes(colVec);

                for (Eigen::SparseVector<int>::InnerIterator itt(markedColVec); itt; ++itt) {
                    resultMat.coeffRef(itt.index(), column) = itt.value();
                }
            }
        }

        const Eigen::SparseMatrix<T, Eigen::RowMajor> _mat = mat;

        Eigen::SparseMatrix<int, Eigen::RowMajor> resultMatRows(mat.rows(), mat.cols());

        QSet<int> visitedRows;
        for (int i = 0; i < _mat.outerSize(); ++i) {
            for (typename Eigen::SparseMatrix<T, Eigen::RowMajor>::InnerIterator it(_mat, i); it; ++it) {

                const int row = it.row();

                if(visitedRows.contains(row)){
                    continue;
                }

                visitedRows.insert(row);

                const Eigen::SparseVector<T> rowVec = _mat.row(row);
                const Eigen::SparseVector<int> markedRowVec = markVectorApexes(rowVec);

                for (Eigen::SparseVector<int>::InnerIterator itt(markedRowVec); itt; ++itt) {
                    resultMatRows.coeffRef(row, itt.index()) = itt.value();
                }
            }
        }

        resultMat += Eigen::SparseMatrix<int>(resultMatRows);

        return resultMat;
    }


    template<typename T>
    static QVector<SparseMatrixPoint> findMatrixApexes(const Eigen::SparseMatrix<T> &mat) {

        const Eigen::SparseMatrix<int> matMarkedApexes = markMatrixApexes<T>(mat);

        QVector<SparseMatrixPoint> resultVec;

        for (int i = 0; i < matMarkedApexes.outerSize(); ++i) {
            for (typename Eigen::SparseMatrix<int>::InnerIterator it(matMarkedApexes, i); it; ++it) {

                const int val = it.value();

                if (val < 2) {
                    continue;
                }

                SparseMatrixPoint smp;
                smp.row = static_cast<int>(it.row());
                smp.col = static_cast<int>(it.col());
                smp.value = static_cast<double>(mat.coeff(smp.row, smp.col));

                resultVec.push_back(smp);
            }
        }

        return resultVec;
    }

//    template<typename T>
//    static QVector<QPair<Index, T>> findSparseVectorApexes(const Eigen::SparseVector<T> &vec) {
//
//        QVector<QPair<Index, T>> apexPoints;
//
//        const int vecSize = static_cast<int>(vec.size());
//
//        for (typename Eigen::SparseVector<T>::InnerIterator it(vec); it; ++it) {
//
//            const int index = it.index();
//            const T val = it.value();
//
//            if(index < 1 || index >= vecSize - 1){
//                continue;
//            }
//
//            const T leftPointVal = vec.coeff(index - 1);
//            const T rightPointVal = vec.coeff(index + 1);
//
//            if(val > leftPointVal && val > rightPointVal){
//                apexPoints.push_back({static_cast<int>(index), val});
//            }
//        }
//
//        return apexPoints;
//    }

    /*!
    * \brief Finds all apexes in a sparse vector.
     *
     *  Returns a QMap where the key is the index of the apex and the value is the value of
     *  the apex.
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


//    template <typename T>
//    static QMap<int, T> apexes(const Eigen::SparseVector<T, Eigen::RowMajor> &vec,
//                               int precision = 1e4){
//        Eigen::SparseVector<double> convertTemplateVec = vec.template cast<double>();
//        return apexes(convertTemplateVec, precision);
//    }


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


    template <typename T>
    static QVector<SparseMatrixPoint> apexes(const Eigen::SparseMatrix<T> &mat,
                                             int precision = 1e4){

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


    template <typename T>
    static QVector<SparseMatrixPoint> apexes(
            const Eigen::SparseMatrix<T, Eigen::RowMajor> &mat,
            int precision = 1e4
                    ){
        const Eigen::SparseMatrix<T>  matConverted = mat.template cast<T>();
        return apexes(matConverted, precision);
    }


    template <typename T, typename T2>
    static Eigen::SparseMatrix<T, Eigen::RowMajor> buildCombFilter(
            const QVector<T> &teethValues,
            T tolerance,
            T maxValue,
            T2 precision) {

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


    template <typename T1, typename  T2, typename T3>
    static Eigen::SparseVector<T2> vectorize(
            const QVector<T1> x,
            const QVector<T2> y,
            T1 maxXVal,
            T3 precision
            ) {

        if (x.size() != y.size()) {
            return {};
        }

        const int buffer = 1;
        const int hashedMaxXVal = MathUtils::hashDecimal(maxXVal + buffer, precision);

        Eigen::SparseVector<T2> vec(hashedMaxXVal);
        vec.setZero();

        for (int i = 0; i < x.size(); i++) {

            const T1 xVal = x.at(i);
            const T2 yVal = y.at(i);
            const int insertIndex = MathUtils::hashDecimal<T1>(xVal, precision);

            if (insertIndex < 0 || insertIndex >= hashedMaxXVal) {
                continue;
            }

            vec.coeffRef(insertIndex) = yVal;
        }

        return vec;
    }


    static Eigen::SparseVector<double> vectorizeFromScan(
            const QVector<QPointF> &scan,
            double maxXVal,
            int precision
            ) {

        const int buffer = 10;
        const int hashedMaxXVal = MathUtils::hashDecimal(maxXVal + buffer, precision);

        Eigen::SparseVector<double> vec(hashedMaxXVal);
        vec.setZero();

        for (const QPointF &p : scan) {

            const double xVal = p.x();
            const double yVal = p.y();
            const int insertIndex = MathUtils::hashDecimal(xVal, precision);

            if (insertIndex < 0 || insertIndex >= hashedMaxXVal) {
                continue;
            }

            vec.coeffRef(insertIndex) = yVal;
        }

        return vec;
    }

    static Eigen::SparseMatrix<double, Eigen::RowMajor> loadFrameToSparseMatrixRowMajor(
            const QMap<int, QVector<QPointF>> &frame,
            int precision,
            double maxRowValue
            );

    static Eigen::SparseMatrix<double, Eigen::ColMajor> loadFrameToSparseMatrixColMajor(
            const QMap<int, QVector<QPointF>> &frame,
            int precision,
            double maxRowValue
    );

    static QMap<int, QVector<QPointF>> loadSparseMatrixToFrame(
            const Eigen::SparseMatrix<double> &mat,
            int precision
            );

};


#endif //EIGENSPARSEUTILS_H
