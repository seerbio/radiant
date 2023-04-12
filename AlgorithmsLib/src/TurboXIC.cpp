//
// Created by anichols on 12/18/22.
//

#include "TurboXIC.h"

#include "ErrorUtils.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN TurboXIC::Private
{
    using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, double> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber);

    XICPoints extractPoints(
            double mzMin,
            double mzMax,
            ScanNumber scanNumberMin,
            ScanNumber scanNumberMax
    );

private:

    RTree *m_rTree;

};


TurboXIC::Private::Private() {}


TurboXIC::Private::~Private() {
    delete m_rTree;
}


Err TurboXIC::Private::init(const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPointsByScanNumber); ree;

    std::vector<rTreePoint> cloudLouder;

    for (auto it = scanPointsByScanNumber.begin(); it != scanPointsByScanNumber.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const ScanPoints &scanPoints = it.value();

        for (const ScanPoint &sp : scanPoints) {
            rTreeCoor coor(static_cast<double>(scanNumber), sp.x());
            cloudLouder.emplace_back(coor, sp.y());
        }
    }

    const int maxElements = 16;
    m_rTree = new RTree(cloudLouder, bgi::dynamic_quadratic(maxElements));

    ERR_RETURN
}


XICPoints TurboXIC::Private::extractPoints(
        double mzMin,
        double mzMax,
        ScanNumber scanNumberMin,
        ScanNumber scanNumberMax
) {

    XICPoints xicPoints;

    const rTreeSearchBox queryBox(
            rTreeCoor(scanNumberMin, mzMin),
            rTreeCoor(scanNumberMax, mzMax)
    );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    for (const rTreePoint &rtp : result) {
        xicPoints[static_cast<int>(rtp.first.get<0>())] += rtp.second;
    }

    return xicPoints;
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


TurboXIC::TurboXIC() : d_ptr(new Private()) {}


TurboXIC::~TurboXIC() {}


Err TurboXIC::init(const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber) {

    ERR_INIT
    e = d_ptr->init(scanPointsByScanNumber); ree;
    ERR_RETURN
}


XICPoints TurboXIC::extractPoints(
        double mzMin,
        double mzMax,
        ScanNumber scanNumberMin,
        ScanNumber scanNumberMax
) {
    return d_ptr->extractPoints(
            mzMin,
            mzMax,
            scanNumberMin,
            scanNumberMax
    );
}
