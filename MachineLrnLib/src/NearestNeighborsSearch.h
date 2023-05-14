//
// Created by anichols on 4/18/23.
//

#ifndef PYTHIADIACPP_NEARESTNEIGBHORSSEARCH_H
#define PYTHIADIACPP_NEARESTNEIGBHORSSEARCH_H

#include "Error.h"
#include "GlobalSettings.h"
#include "MachineLrnLib_Exports.h"

using namespace Error;

struct NNSearchResult {
    QVector<double> searchedCoor;
    double distanceSquared;
    double values;

    NNSearchResult() = default;
    ~NNSearchResult() = default;

    NNSearchResult(
            const QVector<double> &searchedCoor,
            double distanceSquared,
            double value
            )
            : searchedCoor(searchedCoor)
            , distanceSquared(distanceSquared)
            , values(value) {}
};

class MACHINELRNLIB_EXPORTS NearestNeighborsSearch {

public:
    NearestNeighborsSearch();
    ~NearestNeighborsSearch();

    Err init(
            const QVector<QPair<double, Coors>> &valuesVsTreePoints,
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

    [[nodiscard]] int kdTreeSize() const;

    Err writeNearestNeighbors(
            const QString &msDataFilePath,
            QString *calibrationMatFilePath,
            QString *calibarationCalFilePath
            );

    Err readNearestNeighbors(
            const QString &calFilePath,
            const QString &matFilePath
            );


private:

    Q_DISABLE_COPY(NearestNeighborsSearch) class Private;
    const QScopedPointer<Private> d_ptr;



};


#endif //PYTHIADIACPP_NEARESTNEIGBHORSSEARCH_H
