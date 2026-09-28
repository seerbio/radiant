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
#include <boost/iterator/function_output_iterator.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN TurboXIC::Private
{
    using rTreeScanNumber = float;
    using rTreeIntensity = float;
    using rTreeCoor = bg::model::point<float, 1, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, std::pair<rTreeScanNumber , rTreeIntensity>> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);

    Err init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints);

    XICPoints extractPointsXIC(
            float mzMin,
            float mzMax,
            bool restrictScans = false,
            ScanNumber scanNumberMin = 0,
            ScanNumber scanNumberMax = 0,
            bool includeEndpoints = true
    ) const;

    Err getRTreeLimits(
            float *mzMin,
            float *mzMax
            ) const;

    bool isInit() const;

private:

    RTree *m_rTree;

};


TurboXIC::Private::Private()
: m_rTree(Q_NULLPTR) {}

TurboXIC::Private::~Private() {
    delete m_rTree;
}


Err TurboXIC::Private::init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanNumberVsScanPoints); ree;

    const int scanPointsCount = std::accumulate(
            scanNumberVsScanPoints.begin(),
            scanNumberVsScanPoints.end(),
            0,
            [](int sum, ScanPoints *sp){return sum + sp->size();}
            );

    QElapsedTimer et;
    et.start();

    std::vector<rTreePoint> cloudLoader;
    cloudLoader.reserve(scanPointsCount);
    for (auto it = scanNumberVsScanPoints.begin(); it != scanNumberVsScanPoints.end(); ++it) {

        const ScanNumber scanNumber = it.key();
        ScanPoints *scanPoints = it.value();

        for (ScanPoint &sp : *scanPoints) {
            rTreeCoor coor(sp.x());
            const std::pair<rTreeScanNumber, rTreeIntensity> pr(static_cast<float>(scanNumber), sp.y());
            cloudLoader.emplace_back(coor, pr);
        }
    }

    std::sort(cloudLoader.begin(), cloudLoader.end(), [](const rTreePoint &l, const rTreePoint &r){
        return l.first.get<0>() < r.first.get<0>();
    });

    delete m_rTree;

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    ERR_RETURN
}

Err TurboXIC::Private::init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(*scanNumberVsScanPoints); ree;

    const int scanPointsCount = std::accumulate(
            scanNumberVsScanPoints->begin(),
            scanNumberVsScanPoints->end(),
            0,
            [](int sum, ScanPoints *sp){return sum + sp->size();}
    );

    std::vector<rTreePoint> cloudLoader;
    cloudLoader.reserve(scanPointsCount);
    for (auto it = scanNumberVsScanPoints->begin(); it != scanNumberVsScanPoints->end(); it++) {

        const ScanNumber scanNumber = it.key();
        ScanPoints *scanPoints = it.value();

        for (const ScanPoint &sp : *scanPoints) {
            rTreeCoor coor(sp.x());
            std::pair<rTreeScanNumber, rTreeIntensity> pr(static_cast<float>(scanNumber), sp.y());
            cloudLoader.emplace_back(coor, pr);
        }
    }

    std::sort(cloudLoader.begin(), cloudLoader.end(), [](const rTreePoint &l, const rTreePoint &r){
        return l.first.get<0>() < r.first.get<0>();
    });

    delete m_rTree;

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    ERR_RETURN
}

XICPoints TurboXIC::Private::extractPointsXIC(
        float mzMin,
        float mzMax,
        bool restrictScans,
        ScanNumber scanNumberMin,
        ScanNumber scanNumberMax,
        bool includeEndpoints
) const {

    const rTreeSearchBox queryBox(
            (rTreeCoor(mzMin)),
            rTreeCoor(mzMax)
    );

    XICPoints xicPoints;
    const auto append = [&xicPoints](const rTreePoint &rtp) {
        const auto &pr = rtp.second;
        XICPoint xp;
        xp.mz = rtp.first.get<0>();
        xp.intensity = pr.second;
        xp.scanNumber = static_cast<ScanNumber>(pr.first);
        xicPoints.push_back(xp);
    };
    const auto output = boost::make_function_output_iterator(append);
    if (restrictScans) {
        const auto selectedScan = [=](const rTreePoint &point) {
            // Match the old filter after its float-to-ScanNumber conversion.
            const auto scan = static_cast<ScanNumber>(point.second.first);
            return includeEndpoints
                ? scanNumberMin <= scan && scan <= scanNumberMax
                : scanNumberMin < scan && scan < scanNumberMax;
        };
        m_rTree->query(bgi::intersects(queryBox) && bgi::satisfies(selectedScan), output);
    } else {
        m_rTree->query(bgi::intersects(queryBox), output);
    }

    return xicPoints;
}

Err TurboXIC::Private::getRTreeLimits(
        float *mzMin,
        float *mzMax
        ) const {

    ERR_INIT

    e = ErrorUtils::isFalse(m_rTree->empty()); ree

    *mzMin = m_rTree->bounds().min_corner().get<0>();
    *mzMax = m_rTree->bounds().max_corner().get<0>();

    ERR_RETURN
}

bool TurboXIC::Private::isInit() const {

    if (!m_rTree) {
        return false;
    }

    return !m_rTree->empty();
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


TurboXIC::TurboXIC() : d_ptr(QScopedPointer<Private>(new Private)) {}

TurboXIC::~TurboXIC() {
}


Err TurboXIC::init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints) const {
    ERR_INIT
    e = d_ptr->init(scanNumberVsScanPoints); ree;
    ERR_RETURN
}

Err TurboXIC::init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints) const {
    ERR_INIT
    e = d_ptr->init(scanNumberVsScanPoints); ree;
    ERR_RETURN
}

XICPoints TurboXIC::extractPointsXIC(
        float mzMin,
        float mzMax
) const {
    return d_ptr->extractPointsXIC(
            mzMin,
            mzMax
    );
}

XICPoints TurboXIC::extractPointsXIC(
        float mzMin, float mzMax, ScanNumber scanNumberMin,
        ScanNumber scanNumberMax, bool includeEndpoints) const {
    return d_ptr->extractPointsXIC(
        mzMin, mzMax, true, scanNumberMin, scanNumberMax, includeEndpoints);
}

Err TurboXIC::getRTreeLimits(
        float *mzMin,
        float *mzMax
        ) const {
    ERR_INIT

    e = d_ptr->getRTreeLimits(
            mzMin,
            mzMax
            ); ree

    ERR_RETURN
}

bool TurboXIC::isInit() {
    return d_ptr->isInit();
}

void TurboXIC::filterXICPointsByScanNumber(
    ScanNumber scanNumberMin,
    ScanNumber scanNumberMax,
    XICPoints* xicPoints
    ) {

    const auto terminatorLogic = [scanNumberMin, scanNumberMax](const XICPoint &p) {
        return !(scanNumberMin <= p.scanNumber && p.scanNumber <= scanNumberMax);
    };

    const auto terminator = std::remove_if(xicPoints->begin(), xicPoints->end(), terminatorLogic);

    xicPoints->erase(terminator, xicPoints->end());
}
