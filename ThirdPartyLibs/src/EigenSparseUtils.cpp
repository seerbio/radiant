//
// Created by anichols on 11/07/2021.
//

#include "EigenSparseUtils.h"


void  EigenSparseUtils::sortApexesHiLoValue(std::vector<SparseMatrixPoint> *apexes)
{
    const auto sortLogic = [](const SparseMatrixPoint &l, const SparseMatrixPoint &r){
        return l.value < r.value;
    };
    std::sort(apexes->rbegin(), apexes->rend(), sortLogic);
}