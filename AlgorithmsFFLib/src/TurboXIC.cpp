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
    using rTreeScanNumber = float;
    using rTreeIntensity = float;
    using rTreeCoor = bg::model::point<float, 1, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, std::pair<rTreeScanNumber , rTreeIntensity>*> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;
    
    struct Ms2Point {
        float mz = -1.0;
        float intensity = -1.0;
        FrameIndex frameIndex = -1;
    };

public:

    Private();
    ~Private();

    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);

    Err init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints);

    XICPoints extractPointsXIC(
            float mzMin,
            float mzMax
    );

    Err getRTreeLimits(
            float *mzMin,
            float *mzMax
            );

    bool isInit();

private:

    RTree *m_rTree;
    QVector<std::pair<rTreeScanNumber , rTreeIntensity>> m_scanNumberVsIntensityPtrs;

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

    m_scanNumberVsIntensityPtrs.resize(scanPointsCount);
    m_scanNumberVsIntensityPtrs.reserve(scanPointsCount);
    int pointCounter = 0;

    std::vector<rTreePoint> cloudLoader;
    for (auto it = scanNumberVsScanPoints.begin(); it != scanNumberVsScanPoints.end(); it++) {

        const ScanNumber scanNumber = it.key();
        ScanPoints *scanPoints = it.value();

        for (ScanPoint &sp : *scanPoints) {
            rTreeCoor coor(sp.x());
            std::pair<rTreeScanNumber, rTreeIntensity> pr(static_cast<float>(scanNumber), sp.y());
            m_scanNumberVsIntensityPtrs[pointCounter] = pr;
            cloudLoader.push_back(std::make_pair(coor, &m_scanNumberVsIntensityPtrs[pointCounter++]));
        }
    }

    std::sort(cloudLoader.begin(), cloudLoader.end(), [](const rTreePoint &l, const rTreePoint &r){
        return l.first.get<0>() < r.first.get<0>();
    });

    const int maxElements = 16;
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

    m_scanNumberVsIntensityPtrs.resize(scanPointsCount);
    m_scanNumberVsIntensityPtrs.reserve(scanPointsCount);
    int pointCounter = 0;

    std::vector<rTreePoint> cloudLoader;
    for (auto it = scanNumberVsScanPoints->begin(); it != scanNumberVsScanPoints->end(); it++) {

        const ScanNumber scanNumber = it.key();
        ScanPoints *scanPoints = it.value();

        for (const ScanPoint &sp : *scanPoints) {
            rTreeCoor coor(sp.x());
            std::pair<rTreeScanNumber, rTreeIntensity> pr(static_cast<float>(scanNumber), sp.y());
            m_scanNumberVsIntensityPtrs[pointCounter] = pr;
            cloudLoader.push_back(std::make_pair(coor, &m_scanNumberVsIntensityPtrs[pointCounter++]));
        }
    }

    std::sort(cloudLoader.begin(), cloudLoader.end(), [](const rTreePoint &l, const rTreePoint &r){
        return l.first.get<0>() < r.first.get<0>();
    });

    const int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    ERR_RETURN
}

XICPoints TurboXIC::Private::extractPointsXIC(
        float mzMin,
        float mzMax
) {

    XICPoints xicPoints;

    const rTreeSearchBox queryBox(
            (rTreeCoor(mzMin)),
            rTreeCoor(mzMax)
    );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    for (const rTreePoint &rtp : result) {
        std::pair<rTreeScanNumber , rTreeIntensity> *pr = rtp.second;
        xicPoints.scanNumbersVsScanPoints[static_cast<int>(pr->first)].push_back({rtp.first.get<0>(), pr->second});
    }

    return xicPoints;
}

Err TurboXIC::Private::getRTreeLimits(
        float *mzMin,
        float *mzMax
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(m_rTree->empty()); ree

    *mzMin = m_rTree->bounds().min_corner().get<0>();
    *mzMax = m_rTree->bounds().max_corner().get<0>();

    ERR_RETURN
}

bool TurboXIC::Private::isInit() {
    return !m_rTree->empty();
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


TurboXIC::TurboXIC() : d_ptr(new Private()) {}


TurboXIC::~TurboXIC() {}


Err TurboXIC::init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints) {

    ERR_INIT
    e = d_ptr->init(scanNumberVsScanPoints); ree;
    ERR_RETURN
}

Err TurboXIC::init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints) {
    ERR_INIT
    e = d_ptr->init(scanNumberVsScanPoints); ree;
    ERR_RETURN
}


XICPoints TurboXIC::extractPointsXIC(
        float mzMin,
        float mzMax
) {
    return d_ptr->extractPointsXIC(
            mzMin,
            mzMax
    );
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
