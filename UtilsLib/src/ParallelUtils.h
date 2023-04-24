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
    static Err tranchVectorForParallelization(
            const QVector<T> &input,
            int desiredTrancheSize,
            QVector<QVector<T>> *output
    ) {

        ERR_INIT

        output->clear();
        output->reserve(desiredTrancheSize);

        e = ErrorUtils::isNotEmpty(input); ree;

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
    static Err tranchVectorForParallelizationInOrder(
            const QVector<T> &input,
            int desiredTrancheSegments,
            int trancheBuffer,
            QVector<QVector<T>> *output
    ) {

        ERR_INIT

        output->clear();
        output->reserve(desiredTrancheSegments);

        e = ErrorUtils::isNotEmpty(input); ree;
        if (desiredTrancheSegments == -1) {
            desiredTrancheSegments = numberOfAvailableSystemProcessors();
        }

        const int minTrancheSize = 1;
        const int tranchSize = std::max(
                    static_cast<int>(std::round(input.size() / static_cast<double>(desiredTrancheSegments))),
                    minTrancheSize
        );

        QVector<T> currentTranchVec;
        for (int i = 0; i < input.size(); i++) {

            const T &inp = input.at(i);

            if (i % tranchSize == 0 && i > 0) {

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

    template <typename T, typename U>
    static Err tranchMapForParallelizationInOrder(
            const QMap<T, U> &map,
            int numberOfProcesses,
            QVector<QMap<T, U>> *tranchedMaps
            ) {

        ERR_INIT

        QVector<QPair<T, U>> pairs
                = ParallelUtils::convertMapToVectorPairs(map);

        const int itemsPerTranche = static_cast<int>(std::round(map.size() / static_cast<double>(numberOfProcesses)));

        QMap<T,U> trancheMap;
        for (int i = 0; i < pairs.size(); i++) {

            if (i % itemsPerTranche == 0 && i > 0 && tranchedMaps->size() < numberOfProcesses - 1) {
                tranchedMaps->push_back(trancheMap);
                trancheMap.clear();
            }

            const QPair<T, U> &pr = pairs.at(i);
            trancheMap.insert(pr.first, pr.second);
        }

        if (!trancheMap.isEmpty()) {
            tranchedMaps->push_back(trancheMap);
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
            const T &z1,
            const U &z2,
            QVector<QPair<T, U>> *zipResult
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(z1); ree;
        e = ErrorUtils::isEqual(z1.size(), z2.size());

        for (int i = 0; i < z1.size(); i++) {
            zipResult->push_back({z1.at(i), z2.at(i)});
        }

        ERR_RETURN
    }

    template <typename T, typename U>
    static Err zip(
            const T &z1,
            const U &z2,
            QVector<QPointF> *zipResult
    ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(z1); ree;
        e = ErrorUtils::isEqual(z1.size(), z2.size());

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

};


#endif //PARALLELUTILS_H
