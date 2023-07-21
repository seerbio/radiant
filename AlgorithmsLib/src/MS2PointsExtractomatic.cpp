//
// Created by anichols on 7/7/23.
//

#include "MS2PointsExtractomatic.h"

#include "ErrorUtils.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN MS2PointsExtractomatic::Private
{
    using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, int> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QVector<ScanPoints> &scanPoints);

    Err findScanPoints(
            double mzMin,
            double mzMax,
            QMap<Id, ScanPoint> *foundIdVsScanPoints
    );

    Err removeScanPoints(const QMap<Id, ScanPoint> &scanPoints);

    Err updateScanPoints(const QMap<Id, ScanPoint> &scanPoints);

    QMap<Id, ScanPoint> idVsScanPoints();

private:

    RTree *m_rTree;
    QMap<Id, ScanPoint> m_idVsScanPoints;

    double m_defaultYVal;

};

MS2PointsExtractomatic::Private::Private()
: m_rTree(nullptr)
, m_defaultYVal(0.0)
{}

MS2PointsExtractomatic::Private::~Private() {delete m_rTree;}

Err MS2PointsExtractomatic::Private::init(const QVector<ScanPoints> &scanPoints) {
    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;

    std::vector<rTreePoint> cloudLoader;

    int idCurrent = 0;
    for (const ScanPoints &sps : scanPoints) {

        for (const ScanPoint &sp : sps) {
            rTreeCoor coor(sp.x(), m_defaultYVal);
            cloudLoader.emplace_back(coor, idCurrent);
            m_idVsScanPoints.insert(idCurrent, sp);
            idCurrent++;
        }

    }

    const int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    e = ErrorUtils::isFalse(m_rTree->empty()); ree

    ERR_RETURN
}

Err MS2PointsExtractomatic::Private::findScanPoints(
        double mzMin,
        double mzMax,
        QMap<Id, ScanPoint> *foundIdVsScanPoints
        ) {
    ERR_INIT

    const rTreeSearchBox queryBox(
            rTreeCoor(mzMin, m_defaultYVal),
            rTreeCoor(mzMax, m_defaultYVal)
    );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    for (const rTreePoint &rtp : result) {
        foundIdVsScanPoints->insert(rtp.second, m_idVsScanPoints.value(rtp.second));
    }

    ERR_RETURN
}

Err MS2PointsExtractomatic::Private::removeScanPoints(const QMap<Id, ScanPoint> &scanPoints) {
    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree

    for (auto it = scanPoints.begin(); it != scanPoints.end(); it++) {

        const Id id = it.key();
        const ScanPoint &sp = it.value();

        rTreeCoor coor(sp.x(), m_defaultYVal);

        const bool isDeleted = m_rTree->remove(rTreePoint(coor, id));
        e = ErrorUtils::isTrue(isDeleted); ree

        m_idVsScanPoints.remove(id);
        e = ErrorUtils::isFalse(m_idVsScanPoints.contains(id)); ree
    }

    ERR_RETURN
}

Err MS2PointsExtractomatic::Private::updateScanPoints(const QMap<Id, ScanPoint> &scanPoints) {
    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree

    for (auto it = scanPoints.begin(); it != scanPoints.end(); it++) {

        const Id id = it.key();
        const ScanPoint &sp = it.value();

        e = removeScanPoints(scanPoints); ree;

        rTreeCoor coor(sp.x(), m_defaultYVal);
        m_rTree->insert(rTreePoint(coor, id));
        m_idVsScanPoints[id] = sp;
    }

    ERR_RETURN
}

QMap<Id, ScanPoint> MS2PointsExtractomatic::Private::idVsScanPoints() {
    return m_idVsScanPoints;
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


MS2PointsExtractomatic::MS2PointsExtractomatic()
: d_ptr(new Private()) {}

MS2PointsExtractomatic::~MS2PointsExtractomatic() {}

Err MS2PointsExtractomatic::init(const QVector<ScanPoints> &scanPoints) {
    ERR_INIT
    e = d_ptr->init(scanPoints); ree
    ERR_RETURN
}

Err MS2PointsExtractomatic::findScanPoints(
        double mzMin,
        double mzMax,
        QMap<Id, ScanPoint> *foundIdVsScanPoints
        ) {
    ERR_INIT
    e = d_ptr->findScanPoints(
            mzMin,
            mzMax,
            foundIdVsScanPoints
            ); ree
    ERR_RETURN
}

Err MS2PointsExtractomatic::removeScanPoints(const QMap<Id, ScanPoint> &scanPoints) {
    ERR_INIT
    e = d_ptr->removeScanPoints(scanPoints); ree
    ERR_RETURN
}

Err MS2PointsExtractomatic::updateScanPoints(const QMap<Id, ScanPoint> &scanPoints) {
    ERR_INIT
    e = d_ptr->updateScanPoints(scanPoints); ree
    ERR_RETURN
}

QMap<Id, ScanPoint> MS2PointsExtractomatic::idVsScanPoints() {
    return d_ptr->idVsScanPoints();
}
