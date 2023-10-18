//
// Created by anichols on 10/18/23.
//

#include "BoostRTreePoints.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


class Q_DECL_HIDDEN BoostRTreePoints::Private
{
    using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, double> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QVector<RTreePointData2D> &rTreeDataPoint2D);


private:

    RTree *m_rTree;


};

BoostRTreePoints::Private::Private() {}

BoostRTreePoints::Private::~Private() {}

Err BoostRTreePoints::Private::init(const QVector<RTreePointData2D> &rTreeDataPoint2D) {
    return eNetworkError;
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


Err BoostRTreePoints::init(const QVector<RTreePointData2D> &rTreeDataPoint2D) {
    return eNetworkError;
}

Err BoostRTreePoints::init(const QVector<RTreeBoxData2D> &rTreeDataPoint2D) {
    qDebug() << QStringLiteral("This is a points rTree. Trying to use box rTree data");
    return eFunctionNotImplemented;
}
