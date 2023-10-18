//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_BOOSTRTREEABSTRACTFACTORY_H
#define PYTHIADIACPP_BOOSTRTREEABSTRACTFACTORY_H


#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"

using namespace Error;

struct RTreePointData2D {
    double x = -1.0;
    double y = -1.0;
    double val = -1.0;
};

struct RTreeBoxData2D {
    double xMin = -1.0;
    double xMax = -1.0;
    double yMin = -1.0;
    double yMax = -1.0;
    double val = -1.0;
};

class ALGORITHMSLIB_EXPORTS BoostRTreeAbstractFactory {

public:

    virtual ~BoostRTreeAbstractFactory() {}

    virtual Err init(const QVector<RTreePointData2D> &rTreeDataPoint2D) = 0;

    virtual Err init(const QVector<RTreeBoxData2D> &rTreeDataPoint2D) = 0;

};


#endif //PYTHIADIACPP_BOOSTRTREEABSTRACTFACTORY_H
