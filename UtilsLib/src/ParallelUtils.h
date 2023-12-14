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

    static int numberOfAvailableSystemProcessors() {
        return QThread::idealThreadCount();
    }

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

    static QPair<QVector<double>, QVector<double>> unZip(const QVector<QPointF> &points);

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

    template <typename T>
    static QMap<int, T> convertVectorToMap(const QVector<T> &vec) {

        QMap<int, T> vMap;

        for (const T &v : vec) {
            vMap.insert(vMap.size(), v);
        }

        return vMap;
    }

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
