//
// Created by anichols on 10/18/23.
//

#include "BoostRTreeWrapper.h"

#include "BoostRTreeAbstractFactory.h"


#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;



class Q_DECL_HIDDEN BoostRTreeWrapper::Private
{
    using rTreeCoor2D = bg::model::point<double, 2, bg::cs::cartesian>;
    using rTreeBox2D = bg::model::box<rTreeCoor2D>;
    using rTreeEntryPoint2DInt = std::pair<rTreeCoor2D, int> ;
    using RTreePonitsInt = bgi::rtree<rTreeEntryPoint2DInt, bgi::dynamic_quadratic>;
    using rTreeEntryPoint2DDouble = std::pair<rTreeCoor2D, double> ;
    using RTreePointsDouble = bgi::rtree<rTreeEntryPoint2DDouble, bgi::dynamic_quadratic>;


public:

    Private();
    ~Private();

    Err init(const QVector<RTreePointData2D> &rTreeData);

    Err init(const QVector<RTreeBoxData2D> &rTreeData);



private:

    QSharedPointer<BoostRTreeAbstractFactory> *m_rTree;


};

BoostRTreeWrapper::Private::Private() {}

BoostRTreeWrapper::Private::~Private() {}

Err BoostRTreeWrapper::Private::init(const QVector<RTreePointData2D> &rTreeData) {
    ERR_INIT

    ERR_RETURN
}

Err BoostRTreeWrapper::Private::init(const QVector<RTreeBoxData2D> &rTreeData) {
    ERR_INIT

    ERR_RETURN
}



///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

BoostRTreeWrapper::~BoostRTreeWrapper() {}

Err BoostRTreeWrapper::init(const QVector<RTreePointData2D> &rTreeData) {
    ERR_INIT

    ERR_RETURN
}

Err BoostRTreeWrapper::init(const QVector<RTreeBoxData2D> &rTreeData) {
    ERR_INIT

    ERR_RETURN
}
