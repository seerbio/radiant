//
// Created by anichols on 10/21/22.
//

#include "EigenKernelUtils.h"

#include <QtCore>


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
