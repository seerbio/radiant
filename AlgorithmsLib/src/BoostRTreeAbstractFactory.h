//
// Created by anichols on 10/18/23.
//

#ifndef PYTHIADIACPP_BOOSTRTREEABSTRACTFACTORY_H
#define PYTHIADIACPP_BOOSTRTREEABSTRACTFACTORY_H


#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"

using namespace Error;


struct ALGORITHMSLIB_EXPORTS RTreePointData2D {
    double x = -1.0;
    double y = -1.0;
    double val = -1.0;

    RTreePointData2D() = default;

    RTreePointData2D(
            double x,
            double y,
            double val
            )
            : x(x)
            , y(y)
            , val(val) {}



};


struct ALGORITHMSLIB_EXPORTS RTreeBoxData2D {
    double xMin = -1.0;
    double xMax = -1.0;
    double yMin = -1.0;
    double yMax = -1.0;
    double val = -1.0;

    RTreeBoxData2D() = default;

    RTreeBoxData2D(
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            double val
            )
            : xMin(xMin)
            , xMax(xMax)
            , yMin(xMin)
            , yMax(xMax)
            , val(val) {}
};


class ALGORITHMSLIB_EXPORTS BoostRTreeAbstractFactory {

public:

    virtual ~BoostRTreeAbstractFactory() = default;

    virtual Err init(const QVector<RTreePointData2D> &rTreeDataPoint2D) = 0;

    virtual Err init(const QVector<RTreeBoxData2D> &rTreeDataPoint2D) = 0;

    virtual Err getPoints(
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            QVector<RTreePointData2D> *vals
            ) = 0;

    virtual Err getBoxes(
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            QVector<RTreeBoxData2D> *vals
    ) = 0;

};


#endif //PYTHIADIACPP_BOOSTRTREEABSTRACTFACTORY_H
