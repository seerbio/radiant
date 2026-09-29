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
#include <boost/container/small_vector.hpp>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

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

    struct IndexedPoint {
        float mz;
        float scan;
        float intensity;
        std::uint32_t traversalOrder;
    };

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

    void finishInit(std::vector<rTreePoint> &cloudLoader);
    void prepareMassBins();

    RTree *m_rTree;
    std::vector<IndexedPoint> m_pointsByMz;
    static constexpr double massBinScale = 8.0;
    double m_massBinBase = 0.0;
    std::vector<std::uint32_t> m_massBinStarts;
    std::vector<std::uint32_t> m_massBinTraversal;

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

    finishInit(cloudLoader);

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

    finishInit(cloudLoader);

    ERR_RETURN
}

void TurboXIC::Private::finishInit(std::vector<rTreePoint> &cloudLoader) {
    std::sort(cloudLoader.begin(), cloudLoader.end(), [](const rTreePoint &l, const rTreePoint &r){
        return l.first.get<0>() < r.first.get<0>();
    });

    delete m_rTree;

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    m_pointsByMz.clear();
    m_massBinStarts.clear();
    m_massBinTraversal.clear();
    if (cloudLoader.empty()
        || cloudLoader.size() > std::numeric_limits<std::uint32_t>::max()
        || !std::all_of(cloudLoader.begin(), cloudLoader.end(), [](const rTreePoint &point) {
            return std::isfinite(point.first.get<0>());
        })) {
        return;
    }

    // Every spatial query walks surviving tree nodes in the same order.
    // Record that order before sorting by mass so lookup can restore the
    // exact original sequence, including equal-mass observations.
    m_pointsByMz.reserve(cloudLoader.size());
    m_rTree->query(bgi::intersects(m_rTree->bounds()),
        boost::make_function_output_iterator([this](const rTreePoint &point) {
            m_pointsByMz.push_back({
                point.first.get<0>(), point.second.first, point.second.second,
                static_cast<std::uint32_t>(m_pointsByMz.size())});
        }));
    std::sort(m_pointsByMz.begin(), m_pointsByMz.end(),
        [](const IndexedPoint &left, const IndexedPoint &right) {
            return left.mz < right.mz;
        });
    prepareMassBins();
}

void TurboXIC::Private::prepareMassBins() {
    // Power-of-two scaling is exact for stored float masses. The directory
    // only bounds the existing binary searches; it never rounds a query or
    // changes which observations are selected.
    m_massBinBase = std::floor(static_cast<double>(m_pointsByMz.front().mz) * massBinScale);
    const double lastBin = std::floor(static_cast<double>(m_pointsByMz.back().mz) * massBinScale);
    const double span = lastBin - m_massBinBase;
    constexpr std::size_t maxBins = 65536;
    if (!(span >= 0.0 && span < maxBins)) return;
    const std::size_t binCount = static_cast<std::size_t>(span) + 1;
    m_massBinStarts.resize(binCount + 1);
    std::size_t point = 0;
    for (std::size_t bin = 0; bin < binCount; ++bin) {
        while (point < m_pointsByMz.size()
               && std::floor(static_cast<double>(m_pointsByMz[point].mz) * massBinScale)
                      - m_massBinBase < static_cast<double>(bin)) ++point;
        m_massBinStarts[bin] = static_cast<std::uint32_t>(point);
    }
    m_massBinStarts[binCount] = static_cast<std::uint32_t>(m_pointsByMz.size());

    // Stable counting scatter groups points by mass bin while preserving the
    // original tree traversal order within every bin. Building this index is
    // linear in the point/bin counts; it stores one extra 32-bit index per point.
    std::vector<std::uint32_t> traversal(m_pointsByMz.size());
    for (std::size_t i = 0; i < m_pointsByMz.size(); ++i)
        traversal[m_pointsByMz[i].traversalOrder] = static_cast<std::uint32_t>(i);
    auto next = m_massBinStarts;
    m_massBinTraversal.resize(m_pointsByMz.size());
    for (const std::uint32_t index : traversal) {
        const auto bin = static_cast<std::size_t>(
            std::floor(static_cast<double>(m_pointsByMz[index].mz) * massBinScale)
            - m_massBinBase);
        m_massBinTraversal[next[bin]++] = index;
    }
}

XICPoints TurboXIC::Private::extractPointsXIC(
        float mzMin,
        float mzMax,
        bool restrictScans,
        ScanNumber scanNumberMin,
        ScanNumber scanNumberMax,
        bool includeEndpoints
) const {

    if (!m_pointsByMz.empty() && std::isfinite(mzMin) && std::isfinite(mzMax)
        && mzMin <= mzMax) {
        auto rangeBegin = m_pointsByMz.cbegin();
        auto rangeEnd = m_pointsByMz.cend();
        std::size_t firstBin = 0, pastLastBin = 0;
        if (!m_massBinStarts.empty()) {
            const std::size_t binCount = m_massBinStarts.size() - 1;
            const auto clampBin = [binCount](double bin) -> std::size_t {
                if (bin <= 0.0) return 0;
                if (bin >= static_cast<double>(binCount)) return binCount;
                return static_cast<std::size_t>(bin);
            };
            const double lowerBin = std::floor(static_cast<double>(mzMin) * massBinScale) - m_massBinBase;
            const double upperBin = std::floor(static_cast<double>(mzMax) * massBinScale) - m_massBinBase + 1.0;
            firstBin = clampBin(lowerBin);
            pastLastBin = clampBin(upperBin);
            rangeBegin += m_massBinStarts[firstBin];
            rangeEnd = m_pointsByMz.cbegin() + m_massBinStarts[pastLastBin];
        }
        const auto first = std::lower_bound(
            rangeBegin, rangeEnd, mzMin,
            [](const IndexedPoint &point, float mass) { return point.mz < mass; });
        const auto last = std::upper_bound(
            first, rangeEnd, mzMax,
            [](float mass, const IndexedPoint &point) { return mass < point.mz; });

        // Broad queries are cheaper in tree traversal order. Narrow mass
        // windows avoid the tree's recursive visitors and restore that order
        // only among the selected points.
        constexpr std::ptrdiff_t maxIndexedRange = 4096;
        if (last - first <= maxIndexedRange) {
            boost::container::small_vector<const IndexedPoint*, 64> selected;
            const auto appendSelected = [&](const IndexedPoint &point) {
                const auto scan = static_cast<ScanNumber>(point.scan);
                if (restrictScans && !(includeEndpoints
                    ? scanNumberMin <= scan && scan <= scanNumberMax
                    : scanNumberMin < scan && scan < scanNumberMax)) {
                    return;
                }
                selected.push_back(&point);
            };
            // Only scan a whole bin when it contains at most twice as many
            // points as the exact mass range. This bounds the work by O(k)
            // and avoids sorting k matches on every dense, narrow query.
            const auto massMatches = last - first;
            const bool scanOrderedBin = pastLastBin == firstBin + 1
                && massMatches >= 16 && rangeEnd - rangeBegin <= 2 * massMatches;
            if (scanOrderedBin) {
                const auto firstIndex = static_cast<std::uint32_t>(first - m_pointsByMz.cbegin());
                const auto matchCount = static_cast<std::uint32_t>(massMatches);
                for (auto offset = m_massBinStarts[firstBin];
                     offset < m_massBinStarts[pastLastBin]; ++offset) {
                    const auto index = m_massBinTraversal[offset];
                    // The sorted index interval already encodes the exact
                    // inclusive mass bounds, including equal-mass points.
                    if (index - firstIndex < matchCount) appendSelected(m_pointsByMz[index]);
                }
            } else {
                for (auto it = first; it != last; ++it) appendSelected(*it);
                std::sort(selected.begin(), selected.end(),
                    [](const IndexedPoint *left, const IndexedPoint *right) {
                        return left->traversalOrder < right->traversalOrder;
                    });
            }
            XICPoints xicPoints;
            xicPoints.reserve(selected.size());
            for (const IndexedPoint *point : selected) {
                XICPoint xp;
                xp.mz = point->mz;
                xp.intensity = point->intensity;
                xp.scanNumber = static_cast<ScanNumber>(point->scan);
                xicPoints.push_back(xp);
            }
            return xicPoints;
        }
    }

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
