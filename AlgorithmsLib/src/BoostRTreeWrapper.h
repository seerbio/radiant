//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_BOOSTRTREEWRAPPER_H
#define PYTHIADIACPP_BOOSTRTREEWRAPPER_H

#include "AlgorithmsLib_Exports.h"

#include "BoostRTreeAbstractFactory.h"
#include "Error.h"
#include "GlobalSettings.h"

class RTreePointData2D;
class RTreeBoxData2D;

class ALGORITHMSLIB_EXPORTS BoostRTreeWrapper {

public:

    BoostRTreeWrapper() = default;
    ~BoostRTreeWrapper();

    Err init(const QVector<RTreePointData2D> &rTreeData);

    Err init(const QVector<RTreeBoxData2D> &rTreeData);

    Err getPoints(
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            QVector<RTreePointData2D> *vals
    );

    Err getBoxes(
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            QVector<RTreeBoxData2D> *vals
    );


private:

    QSharedPointer<BoostRTreeAbstractFactory> m_rTree;


};


#endif //PYTHIADIACPP_BOOSTRTREEWRAPPER_H
