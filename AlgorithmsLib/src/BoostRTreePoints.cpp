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

    Err getPoints(
            double xMin,
            double xMax,
            double yMin,
            double yMax,
            QVector<RTreePointData2D> *vals
    );


private:

    bool m_isInit;
    RTree *m_rTree;

};

BoostRTreePoints::Private::Private()
: m_rTree(Q_NULLPTR)
, m_isInit(false)
{}

BoostRTreePoints::Private::~Private() {
    delete m_rTree;
}

Err BoostRTreePoints::Private::init(const QVector<RTreePointData2D> &rTreeDataPoint2D) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(rTreeDataPoint2D); ree;

    std::vector<rTreePoint> cloudLoader;

    std::transform(
            rTreeDataPoint2D.begin(),
            rTreeDataPoint2D.end(),
            std::back_inserter(cloudLoader),
            [](const RTreePointData2D &p){

                const rTreeCoor coor(p.x, p.y);
                return rTreePoint(coor, p.val);
            }
    );

    const int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    m_isInit = true;
    ERR_RETURN
}

Err BoostRTreePoints::Private::getPoints(
        double xMin,
        double xMax,
        double yMin,
        double yMax,
        QVector<RTreePointData2D> *vals
) {

    ERR_INIT

    e = ErrorUtils::isTrue(xMin < xMax); ree;
    e = ErrorUtils::isTrue(yMin < yMax); ree;
    e = ErrorUtils::isTrue(m_isInit); ree;

    const rTreeSearchBox queryBox(
            rTreeCoor(xMin, yMin),
            rTreeCoor(xMax, yMax)
            );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    for (const rTreePoint &rtp : result) {

        RTreePointData2D dp;
        dp.x = rtp.first.get<0>();
        dp.y = rtp.first.get<1>();
        dp.val = rtp.second;

        vals->push_back(dp);
    }

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

BoostRTreePoints::BoostRTreePoints() : d_ptr(new Private()) {}

BoostRTreePoints::~BoostRTreePoints() {}

Err BoostRTreePoints::init(const QVector<RTreePointData2D> &rTreeDataPoint2Ds) {
    ERR_INIT
    e = d_ptr->init(rTreeDataPoint2Ds); ree;
    ERR_RETURN
}

Err BoostRTreePoints::init(const QVector<RTreeBoxData2D> &rTreeBoxPoint2D) {
    qDebug() << QStringLiteral("This is a points rTree. Trying to use box rTree data");
    return eFunctionNotImplemented;
}

Err BoostRTreePoints::getPoints(
        double xMin,
        double xMax,
        double yMin,
        double yMax,
        QVector<RTreePointData2D> *vals
        ) {

    ERR_INIT
    e = d_ptr->getPoints(
            xMin,
            xMax,
            yMin,
            yMax,
            vals
            ); ree;
    ERR_RETURN
}

Err BoostRTreePoints::getBoxes(
        double xMin,
        double xMax,
        double yMin,
        double yMax,
        QVector<RTreeBoxData2D> *vals
        ) {
    qDebug() << QStringLiteral("This is a points rTree. Trying to use box rTree data");
    return eFunctionNotImplemented;
}
