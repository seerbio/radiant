//
// Created by andrewnichols on 8/15/24.
//

#include "TurboXIC2D.h"
#include "TurboXIC.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN TurboXIC2D::Private
{
    using rTreeScanNumber = float;
    using rTreeIntensity = float;
    using rTreeCoor = bg::model::point<float, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, rTreeIntensity> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);

    XICPoints extractPointsXIC(
        float dim1Min,
        float dim1Max,
        float dim2Min,
        float dim2Max
    ) const;

    Err getRTreeLimits(
            float *dim1Min,
            float *dim1Max,
            float *dim2Min,
            float *dim2Max
            ) const;

    bool isInit() const;

private:

    RTree *m_rTree;

};

TurboXIC2D::Private::Private() : m_rTree(Q_NULLPTR) {}

TurboXIC2D::Private::~Private() {
    delete m_rTree;
}

Err TurboXIC2D::Private::init(const QMap<ScanNumber, ScanPoints*>& scanNumberVsScanPoints) {

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
            rTreeCoor coor(sp.x(), static_cast<float>(scanNumber));
            cloudLoader.emplace_back(coor, sp.y());
        }
    }

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    ERR_RETURN;
}

XICPoints TurboXIC2D::Private::extractPointsXIC(
    float dim1Min,
    float dim1Max,
    float dim2Min,
    float dim2Max
    ) const {
    ERR_INIT

    const rTreeSearchBox queryBox(
            rTreeCoor(dim1Min, dim2Min),
            rTreeCoor(dim1Max, dim2Max)
    );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    const ulong loaderSize = result.size();
    XICPoint xicPointsLoader[loaderSize];

    for (int i = 0; i <  loaderSize; i++) {

        const rTreePoint &rtp = result[i];

        XICPoint xp;
        xp.mz = rtp.first.get<0>();
        xp.intensity = rtp.second;
        xp.scanNumber = static_cast<ScanNumber>(rtp.first.get<1>());

        xicPointsLoader[i] = xp;
    }

    XICPoints xicPoints = {xicPointsLoader, xicPointsLoader + loaderSize};
    std::sort(
        xicPoints.begin(),
        xicPoints.end(),
        [](const XICPoint &l, const XICPoint &r){return l.scanNumber < r.scanNumber;}
        );

    return xicPoints;
}

Err TurboXIC2D::Private::getRTreeLimits(
    float* dim1Min,
    float* dim1Max,
    float* dim2Min,
    float* dim2Max
    ) const {

    ERR_INIT

    e = ErrorUtils::isFalse(m_rTree->empty()); ree

    *dim1Min = m_rTree->bounds().min_corner().get<0>();
    *dim1Max = m_rTree->bounds().max_corner().get<0>();
    *dim2Min = m_rTree->bounds().max_corner().get<1>();
    *dim2Max = m_rTree->bounds().max_corner().get<1>();

    ERR_RETURN
}

bool TurboXIC2D::Private::isInit() const {
    if (!m_rTree) {
        return false;
    }

    return !m_rTree->empty();
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

TurboXIC2D::TurboXIC2D() : d_ptr(QScopedPointer<Private>(new Private)) {
}

TurboXIC2D::~TurboXIC2D() {}

Err TurboXIC2D::init(const QMap<ScanNumber, ScanPoints*>& scanNumberVsScanPoints) const {
    ERR_INIT
    e = d_ptr->init(scanNumberVsScanPoints); ree;
    ERR_RETURN;
}

XICPoints TurboXIC2D::extractPointsXIC(
    float dim1Min,
    float dim1Max,
    float dim2Min,
    float dim2Max
    ) {
    return d_ptr->extractPointsXIC(
        dim1Min,
        dim1Max,
        dim2Min,
        dim2Max
        );
}

Err TurboXIC2D::getRTreeLimits(
    float* dim1Min,
    float* dim1Max,
    float* dim2Min,
    float* dim2Max
    ) const {

    ERR_INIT
    e = d_ptr->getRTreeLimits(dim1Min, dim1Max, dim2Min, dim2Max); ree;
    ERR_RETURN
}

bool TurboXIC2D::isInit() const {
    return d_ptr->isInit();
}
