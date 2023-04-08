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
    static QVector<QPair<T, U>> convertMapToVectorPairs(const QMap<T, U> &map) {

        QVector<QPair<T, U>> output;
        for (auto it = map.begin(); it != map.end(); it++) {
            const T &t = it.key();
            const U &u = it.value();

            output.push_back({t, u});
        }

        return output;
    }

};


#endif //PARALLELUTILS_H
