//
// Created by Drucifer on 12/31/2021.
//

#ifndef PARALLELUTILS_H
#define PARALLELUTILS_H

#include "ErrorUtils.h"
#include "UtilsLib_Exports.h"

#include <QtConcurrent/QtConcurrent>
#include <QVector>

class UTILSLIB_EXPORTS ParallelUtils {

public:

    /*!
    * @brief  Gets the number of processors available in the current system.
    * @return int: The number of processors available in the system as reported by the Qt framework.
    *
    * Note: This method wraps around QThread::idealThreadCount() provided by the Qt framework to get the ideal number of threads that can run concurrently for the current system.
    */
    static int numberOfAvailableSystemProcessors() {
        return QThread::idealThreadCount();
    }

    /*!
    * @brief  Distributes elements of an input vector into many sub-vectors ("tranches"), typically for parallel processing.
    * @tparam T: The datatype of the elements in the input and the output vectors.
    * @param input: The input vector to be divided into tranche vectors.
    * @param desiredTrancheSize: The desired number of sub-vectors (tranches). If set to -1, it will default to the number of available system processors.
    * @param output: A pointer to the container that will store the resulting tranche vectors.
    * @return Err: Error status after the checks and operations. Error code will be returned if the input vector is empty or desiredTrancheSize is zero.
    *
    * Note: If there are 'N' desired tranche, each tranche will have approximately input.size()/N elements. Elements are distributed in round-robin fashion, and 'output' vector will have 'N' vectors.
    */
    template <typename T>
    static Err trancheVectorForParallelization(
            const QVector<T> &input,
            int desiredTrancheSize,
            QVector<QVector<T>> *output
    ) {

        ERR_INIT

        output->clear();
        output->reserve(desiredTrancheSize);

        e = ErrorUtils::isNotEmpty(input); ree;
        e = ErrorUtils::isNotEqual(desiredTrancheSize, 0); ree;

        if (desiredTrancheSize == -1) {
            desiredTrancheSize = numberOfAvailableSystemProcessors();
        }

        output->resize(desiredTrancheSize);

        int currentVecIndex = 0;
        for (const T &inp : input) {

            (*output)[currentVecIndex++].push_back(inp) ;

            if (currentVecIndex == desiredTrancheSize){
                currentVecIndex = 0;
            }
        }

        ERR_RETURN;
    }

    /*!
    * @brief  Distributes elements of an input vector into many sub-vectors ("tranches") in their original order, typically for parallel processing.
    * @tparam T: The datatype of the elements in the input and the output vectors.
    * @param input: The input vector to be divided into tranches.
    * @param desiredTrancheSegments: The desired number of sub-vectors (tranches). If set to -1, it will default to the number of available system processors.
    * @param trancheBuffer: The overlapping-yield buffer from the current batch of inputs.
    * @param output: A pointer to the container that will store the resulting tranches.
    * @return Err: Error status after the checks and operations. Error code will be returned if the input vector is empty or desiredTrancheSegments is zero.
    *
    * Note: This method differs from "trancheVectorForParallelization" in preservation of the input order in tranches
    * and provides an optional 'buffer' between tranches to smoothly transition between different processing zones. Depending on the division,
    * there are times when an extra tranche will be creader, i.e. desiredTrancheSegments + 1
    */
    template <typename T>
    static Err trancheVectorForParallelizationInOrder(
            const QVector<T> &input,
            int desiredTrancheSegments,
            int trancheBuffer,
            QVector<QVector<T>> *output
    ) {

        ERR_INIT

        output->clear();
        output->reserve(desiredTrancheSegments);

        e = ErrorUtils::isNotEmpty(input); ree;
        e = ErrorUtils::isNotEqual(desiredTrancheSegments, 0); ree;

        if (desiredTrancheSegments == -1) {
            desiredTrancheSegments = numberOfAvailableSystemProcessors();
        }

        const int minTrancheSize = 1;
        const int trancheSize = std::max(
                    static_cast<int>(std::round(input.size() / static_cast<double>(desiredTrancheSegments))),
                    minTrancheSize
        );

        QVector<T> currentTranchVec;
        for (int i = 0; i < input.size(); i++) {

            const T &inp = input.at(i);

            if (i % trancheSize == 0 && i > 0) {

                for (int j = i; j < i + trancheBuffer; j++) {

                    if (j >= input.size()) {
                        break;
                    }

                    const T &inpBuffer = input.at(j);
                    currentTranchVec.push_back(inpBuffer);
                }

                output->push_back(currentTranchVec);
                currentTranchVec.clear();
            }

            currentTranchVec.push_back(inp);
        }

        if(!currentTranchVec.isEmpty()) {
            output->push_back(currentTranchVec);
        }

        ERR_RETURN;
    }

    /*!
    * @brief  Converts a QMap to a QVector.
    * @tparam T: The datatype of the elements in the QMap and the QVector. This can be a numeric datatype or any other datatype that supports assignment.
    * @param m: The QMap to be converted. The keys of the QMap are assumed to be indices for the QVector.
    * @param vecSize: The desired size of the output QVector.
    * @return QVector<T>: The converted QVector. For any Teuchos::Range1D where key < vecSize, map value is assigned to the corresponding vector index. All other indices are assigned '0.0'.
    *
    * Note: Be aware that if the QMap contains keys that are greater than or equal to vecSize, those keys (and their associated values) will be ignored during the conversion process.
    */
    template <typename T>
    static QVector<T> convertMapToVector(const QMap<int, T> &m, int vecSize){

        QVector<T> vec(vecSize, 0.0);
        vec.reserve(vecSize);

        for (auto it = m.begin(); it != m.end(); it++) {

            const int key = it.key();
            const T val = it.value();

            if (key >= vecSize) {
                continue;
            }

            vec[key] = val;
        }

        return vec;
    }

    /*!
    * @brief  Divides a QMap into multiple segments ("tranches") and stores them in a QVector for easier parallel processing.
    * @tparam T: The datatype of the map keys.
    * @tparam U: The datatype of the map values.
    * @param map: The QMap to be divided into segments.
    * @param numberOfProcesses: The number of desired tranche segments.
    * @param tranchedMaps: A pointer to the QVector that will store the tranche segments.
    * @return Err: Error status after the checks and operations. Error code will be returned if there are errors in the operations, and 'eNoError' otherwise.
    *
    * Note: This method is a specialized wrapper for 'trancheVectorForParallelization' that allows for QMaps to be divided into segments.
    */
    template <typename T, typename U>
    static Err trancheMapForParallelization(
            const QMap<T, U> &map,
            int numberOfProcesses,
            QVector<QVector<QPair<T, U>>> *tranchedMaps
            ) {

        ERR_INIT

        QVector<QPair<T, U>> pairs = ParallelUtils::convertMapToVectorPairs(map);

        e = trancheVectorForParallelization(
                pairs,
                numberOfProcesses,
                tranchedMaps
                ); ree;

        ERR_RETURN
    }

    /*!
    * @brief  Divides values of a QMap (QVector elements) into multiple segments ("tranches") and stores them in a QVector of QMaps for easier parallel processing.
    * @tparam T: The datatype of the map keys.
    * @tparam U: The datatype of the map values.
    * @param map: The QMap with QVector values to be divided into segments.
    * @param numberOfTranches: The number of desired tranche segments.
    * @param tranchedMaps: A pointer to the QVector of QMap that will store the tranche segments.
    * @return Err: Error status after the checks and operations. Error code will be returned if there are errors in the operations, and 'eNoError' otherwise.
    *
    * Note: This method is a specialised version of 'trancheMapForParallelization' that segments QMap values (which are QVectors) based on the numberOfTranches.
    */
    template <typename T, typename U>
    static Err trancheMapValueVectorsByKeyForParallelization(
            const QMap<T, QVector<U>> &map,
            int numberOfTranches,
            QVector<QMap<T, QVector<U>>> *tranchedMaps
    ) {

        ERR_INIT

        tranchedMaps->resize(numberOfTranches);
        tranchedMaps->reserve(numberOfTranches);

        for (auto it = map.begin(); it != map.end(); it++) {
            const T key = it.key();
            const QVector<U> &val = it.value();

            QVector<QVector<U>> valTranched;
            e = trancheVectorForParallelization(
                    val,
                    numberOfTranches,
                    &valTranched
            ); ree;

            for (int i = 0; i < valTranched.size(); i++) {
                e = ErrorUtils::isTrue(i < tranchedMaps->size());
                (*tranchedMaps)[i][key] = valTranched[i];
            }
        }

        ERR_RETURN
    }

    /*!
    * @brief  Converts a QMap into a QVector of QPairs.
    * @tparam T: The datatype of the map keys.
    * @tparam U: The datatype of the map values.
    * @param map: The QMap to be converted to a QVector of QPairs.
    * @return QVector<QPair<T, U>>: A QVector of QPairs, created from the input QMap. Each pair in the vector represents an original key-value pair from the map.
    *
    * Note: This conversion method doesn't change the order of the elements, it just changes the container of the elements from QMap to QVector.
    */
    template <typename T, typename U>
    static QVector<QPair<T, U>> convertMapToVectorPairs(const QMap<T, U> &map) {

        QVector<QPair<T, U>> output;
        for (auto it = map.begin(); it != map.end(); it++) {
            const T &t = it.key();
            const U &u = it.value();

            output.push_back({t, u});
        }

        return output;
    }

    /*!
    * @brief  Combines two vectors into a single vector of pairs.
    * @tparam T: The datatype of the elements in the first vector.
    * @tparam U: The datatype of the elements in the second vector.
    * @param z1: The first vector to zip.
    * @param z2: The second vector to zip.
    * @param zipResult: A pointer to a QVector that will store the result of zipping.
    * @return Err: Error status after the checks and operations. Error code will be 'eSizeError' if two vectors are not the same size, and 'eNoError' otherwise.
    *
    * Note: The two input vectors should be of the same size. Each individual element from the two vectors is combined into a QPair and stored in zipResult QVector.
    */
    template <typename T, typename U>
    static Err zip(
            const QVector<T> &z1,
            const QVector<U> &z2,
            QVector<QPair<T, U>> *zipResult
            ) {

        ERR_INIT

        e = ErrorUtils::isEqual(z1.size(), z2.size()); ree;

        for (int i = 0; i < z1.size(); i++) {
            zipResult->push_back({z1.at(i), z2.at(i)});
        }

        ERR_RETURN
    }

    /*!
    * @brief  Combines two vectors into a single vector of QPointF objects.
    * @tparam T: The datatype of the elements in the first vector. This should be a numeric type that can be converted to qreal, the type Qt uses for floating point values.
    * @tparam U: The datatype of the elements in the second vector. This should be a numeric type that can be converted to qreal.
    * @param z1: The first vector to zip.
    * @param z2: The second vector to zip.
    * @param zipResult: A pointer to a QVector that will store the result of zipping as QPointF.
    * @return Err: Error status after the checks and operations. Error code will be 'eSizeError' if two vectors are not the same size, and 'eNoError' otherwise.
    *
    * Note: The two input vectors should be of the same size. Each individual element from the two vectors is combined into a QPointF and stored in zipResult QVector.
    */
    template <typename T, typename U>
    static Err zip(
            const QVector<T> &z1,
            const QVector<U> &z2,
            QVector<QPointF> *zipResult
    ) {

        ERR_INIT

        e = ErrorUtils::isEqual(z1.size(), z2.size()); ree;

        for (int i = 0; i < z1.size(); i++) {
            zipResult->push_back({z1.at(i), z2.at(i)});
        }

        ERR_RETURN
    }

    /*!
    * @brief  Combines two vectors of the same type into a single vector of pairs.
    * @tparam T: The datatype of the elements in the two vectors.
    * @param z1: The first vector to zip.
    * @param z2: The second vector to zip.
    * @param zipResult: A pointer to a QVector that will store the result of zipping.
    * @return Err: Error status after the checks and operations. Error code will be 'eSizeError' if two vectors are not the same size, and 'eNoError' otherwise.
    *
    * Note: This function requires the two input vectors to be of the same type and size. Each individual element from the two vectors is combined into a QPair and stored in zipResult QVector.
    */
    template <typename T>
    static Err zip(
            const QVector<T> &z1,
            const QVector<T> &z2,
            QVector<QPair<T, T>> *zipResult
            ) {

        ERR_INIT

        e = ErrorUtils::isEqual(z1.size(), z2.size()); ree;

        for (int i = 0; i < z1.size(); i++) {
            zipResult->push_back({z1.at(i), z2.at(i)});
        }

        ERR_RETURN
    }

    /*!
    * @brief  Splits a vector of QPointF into a pair of double vectors.
    * @param points: The QVector of QPointF to be split.
    * @return QPair<QVector<double>, QVector<double>>: A pair of vectors. The first vector in the pair contains the x-coordinates of the points, and the second vector contains the y-coordinates of the points.
    *
    * Note: This operation is the opposite of zipping operation where two vectors are combined into one. Here, one vector is split into two vectors. The QPointF's x and y values are used to form two distinct vectors.
    */
    template <typename PointXY>
    static QPair<QVector<double>, QVector<double>> unZip(const QVector<PointXY> &points) {

        QVector<double> v1;
        v1.reserve(points.size());

        QVector<double> v2;
        v1.reserve(points.size());

        for (const PointXY &p : points) {
            v1.push_back(p.x());
            v2.push_back(p.y());
        }

        return {v1, v2};
    }

    /*!
    * @brief  Splits a vector of QPairs into two separate vectors.
    * @tparam T: The datatype of the first elements in the pairs.
    * @tparam U: The datatype of the second elements in the pairs.
    * @param vecOfPairs: The QVector of QPairs to be split.
    * @return QPair<QVector<T>, QVector<U>>: A pair of vectors. The first vector in the pair contains the first
    *   elements from the pairs, and the second vector contains the second elements from the pairs.
    *
    * Note: This operation is the opposite of a zip operation where two vectors are combined into one as a
    *   QVector of QPairs. Here, one vector of QPairs is split into two vectors. The QPair's first and second values are used to form two distinct vectors.
    */
    template <typename T, typename U>
    static QPair<QVector<T>, QVector<U>> unZip(const QVector<QPair<T, U>> &vecOfPairs) {

        QVector<T> v1;
        v1.reserve(vecOfPairs.size());

        QVector<U> v2;
        v2.reserve(vecOfPairs.size());

        for (const QPair<T, U> &pr : vecOfPairs) {
            v1.push_back(pr.first);
            v2.push_back(pr.second);
        }

        return {v1, v2};
    }

    /*!
    * @brief  Splits a vector of QPairs into two separate vectors.
    * @tparam T: The datatype of the first elements in the pairs.
    * @tparam U: The datatype of the second elements in the pairs.
    * @param vecOfPairs: The QVector of QPairs to be split.
    * @return QPair<QVector<T>, QVector<U>>: A pair of vectors. The first vector in the pair contains the first elements from the pairs, and the second vector contains the second elements from the pairs.
    *
    * Note: This operation is the opposite of zip operation where two vectors are combined into one as QVector of QPairs. Here, one vector is split into two vectors. The QPair's first and second values are used to form two distinct vectors.
    */
    template <typename T>
    static QMap<int, T> convertVectorToMap(const QVector<T> &vec) {

        QMap<int, T> vMap;

        for (const T &v : vec) {
            vMap.insert(vMap.size(), v);
        }

        return vMap;
    }

    /*!
    * @brief  Converts a QMap to a QVector of QPointF objects.
    * @tparam T: The datatype of the map keys. It should be a numeric type that can be converted to qreal, the type Qt uses for floating point values.
    * @tparam U: The datatype of the map values. It should be a numeric type that can be converted to qreal.
    * @param vec: The QMap to be converted into a QVector of QPointF. The keys of QMap will be used as x coordinates and the values as y coordinates for QPointF.
    * @return QVector<QPointF>: A QVector of QPointF objects.
    *
    * Note: This method transforms key-value pairs in QMap to x-y coordinates in QPointF. Ideal when QPointF's x and y have different significance and are stored in QMap accordingly.
    */
    template <typename T, typename U>
    static QVector<QPointF> convertMapToPoints(const QMap<T, U> &vec) {

        QVector<QPointF> points;

        for (auto it = vec.begin(); it != vec.end(); it++) {

            const U x = static_cast<U>(it.key());
            const U y = it.value();

            points.push_back({x, y});
        }

        return points;
    }

};


#endif //PARALLELUTILS_H
