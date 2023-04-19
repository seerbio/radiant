//
// Created by anichols on 4/18/23.
//

#ifndef PYTHIADIACPP_NEARESTNEIGBHORSSEARCH_H
#define PYTHIADIACPP_NEARESTNEIGBHORSSEARCH_H

#include "Error.h"
#include "MachineLrnLib_Exports.h"

#include "Eigen/Dense"

using namespace Error;


struct NNSearchResult {
    std::vector<long> indexes;
    std::vector<double> distancesSquared;
    std::vector<double> values;

    NNSearchResult() = default;
    ~NNSearchResult() = default;

    NNSearchResult(
            const std::vector<long> &indexes,
            const std::vector<double> &distancesSquared,
            const std::vector<double> &values
            )
            : indexes(indexes)
            , distancesSquared(distancesSquared)
            , values(values) {}

};

class MACHINELRNLIB_EXPORTS NearestNeighborsSearch {

public:
    NearestNeighborsSearch();
    ~NearestNeighborsSearch();

    Err init(
            const QVector<QPair<double, QVector<double>>> &valuesVsTreePoints,
            int maxTreeLeafSize = 30
            );

    Err kNearestNeighborsSearch(
            const QVector<QVector<double>> &searchPointCoors,
            int k,
            QVector<NNSearchResult> *result
            );

    Err radiusSearch(
            const QVector<QVector<double>> &searchPointCoors,
            double searchRadiusSquared,
            QVector<NNSearchResult> *result
    );


private:

    Q_DISABLE_COPY(NearestNeighborsSearch) class Private;
    const QScopedPointer<Private> d_ptr;



};


#endif //PYTHIADIACPP_NEARESTNEIGBHORSSEARCH_H
