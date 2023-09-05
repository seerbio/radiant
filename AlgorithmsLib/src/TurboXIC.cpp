//
// Created by anichols on 12/18/22.
//

#include "TurboXIC.h"

#include "EigenSparseUtils.h"
#include "ErrorUtils.h"
#include "MsUtils.h"

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

    XICPoints extractPointsXIC(
            double mzMin,
            double mzMax,
            ScanNumber scanNumberMin,
            ScanNumber scanNumberMax
    );

    ScanPoints extractSpectrum(
            double mzMin,
            double mzMax,
            ScanNumber scanNumberMin,
            ScanNumber scanNumberMax
    );

    Err getRTreeLimits(
            double *scanNumberMin,
            double *scanNumberMax,
            double *mzMin,
            double *mzMax
            );

    bool isInit();

private:

    RTree *m_rTree;
    int m_defaultPrecision;

};


TurboXIC::Private::Private()
: m_defaultPrecision(4)
, m_rTree(Q_NULLPTR) {}


TurboXIC::Private::~Private() {
    delete m_rTree;
}


Err TurboXIC::Private::init(const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPointsByScanNumber); ree;

    std::vector<rTreePoint> cloudLoader;

    for (auto it = scanPointsByScanNumber.begin(); it != scanPointsByScanNumber.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const ScanPoints &scanPoints = it.value();

        for (const ScanPoint &sp : scanPoints) {
            rTreeCoor coor(static_cast<double>(scanNumber), sp.x());
            cloudLoader.emplace_back(coor, sp.y());
        }
    }

    const int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    ERR_RETURN
}


XICPoints TurboXIC::Private::extractPointsXIC(
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
        xicPoints.scanNumbersVsIntensityVals[static_cast<int>(rtp.first.get<0>())] += rtp.second;
        xicPoints.scanNumberVsMzVals[static_cast<int>(rtp.first.get<0>())].push_back(rtp.first.get<1>());
    }

    return xicPoints;
}

Err TurboXIC::Private::getRTreeLimits(
        double *scanNumberMin,
        double *scanNumberMax,
        double *mzMin,
        double *mzMax
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(m_rTree->empty()); ree

    *scanNumberMin = m_rTree->bounds().min_corner().get<0>();
    *scanNumberMax = m_rTree->bounds().max_corner().get<0>();
    *mzMin = m_rTree->bounds().min_corner().get<1>();
    *mzMax = m_rTree->bounds().max_corner().get<1>();

    ERR_RETURN
}

ScanPoints TurboXIC::Private::extractSpectrum(
        double mzMin,
        double mzMax,
        ScanNumber scanNumberMin,
        ScanNumber scanNumberMax
        ) {

    const rTreeSearchBox queryBox(
            rTreeCoor(scanNumberMin, mzMin),
            rTreeCoor(scanNumberMax, mzMax)
    );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    const int buffer = 1;
    const int vecSize = MathUtils::hashDecimal(
            mzMax + buffer,
            m_defaultPrecision
            );

    Eigen::SparseVector<double> vec(vecSize);
    vec.setZero();

    for (const rTreePoint &rtp : result) {

        const double mz = rtp.first.get<1>();
        const int mzHashed = MathUtils::hashDecimal(mz, m_defaultPrecision);

        vec.coeffRef(mzHashed) += rtp.second;
    }

    ScanPoints scanPoints;

    for (Eigen::SparseVector<double>::InnerIterator it(vec); it; ++it) {

        const auto mz = MathUtils::unHashDecimal<double>(it.index(), m_defaultPrecision);
        const double intensity = it.value();

        scanPoints.push_back({mz, intensity});
    }

    std::sort(scanPoints.begin(), scanPoints.end(), MsUtilsNamespace::sortAscMz);

    return scanPoints;
}

bool TurboXIC::Private::isInit() {
    return !m_rTree->empty();
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


XICPoints TurboXIC::extractPointsXIC(
        double mzMin,
        double mzMax,
        ScanNumber scanNumberMin,
        ScanNumber scanNumberMax
) {
    return d_ptr->extractPointsXIC(
            mzMin,
            mzMax,
            scanNumberMin,
            scanNumberMax
    );
}

Err TurboXIC::getRTreeLimits(
        double *scanNumberMin,
        double *scanNumberMax,
        double *mzMin,
        double *mzMax
        ) const {
    ERR_INIT

    e = d_ptr->getRTreeLimits(
            scanNumberMin,
            scanNumberMax,
            mzMin,
            mzMax
            ); ree

    ERR_RETURN
}

ScanPoints TurboXIC::extractSpectrum(
        double mzMin,
        double mzMax,
        ScanNumber scanNumberMin,
        ScanNumber scanNumberMax
        ) const {

    return d_ptr->extractSpectrum(mzMin, mzMax, scanNumberMin, scanNumberMax);
}

bool TurboXIC::isInit() {
    return false;
}
