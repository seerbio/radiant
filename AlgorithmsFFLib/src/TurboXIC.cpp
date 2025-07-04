//
// Created by anichols on 12/18/22.
//

#include "TurboXIC.h"

#include "EigenSparseUtils.h"
#include "ErrorUtils.h"
#include "MsUtils.h"

#include <immintrin.h>

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
    using rTreePoint = std::pair<rTreeCoor, std::pair<rTreeScanNumber , rTreeIntensity>>;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);

    Err init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints);

    XICPoints extractPointsXIC(
            float mzMin,
            float mzMax
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
        float mzMax
) const {

    const rTreeSearchBox queryBox(
            (rTreeCoor(mzMin)),
            rTreeCoor(mzMax)
    );

    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    const ulong loaderSize = result.size();
    XICPoint xicPointsLoader[loaderSize];

    for (int i = 0; i <  result.size(); i++) {

        const rTreePoint &rtp = result[i];
        const std::pair<rTreeScanNumber , rTreeIntensity> &pr = rtp.second;

        XICPoint xp;
        xp.mz = rtp.first.get<0>();
        xp.intensity = pr.second;
        xp.scanNumber = static_cast<ScanNumber>(pr.first);

        xicPointsLoader[i] = xp;
    }

    XICPoints xicPoints = {xicPointsLoader, xicPointsLoader + loaderSize};
    // std::sort(
    //     xicPoints.begin(),
    //     xicPoints.end(),
    //     [](const XICPoint &l, const XICPoint &r){return l.scanNumber < r.scanNumber;}
    //     );

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


TurboXIC::TurboXIC()
: d_ptr(QScopedPointer<Private>(new Private))
, m_mzVals(nullptr)
, m_indexesMz(nullptr)
, m_indexesI(nullptr)
, m_alignasPiecesOfEight(0)
// , m_indexes(nullptr)
, m_scanPointsCountAlignas(0)
{}

TurboXIC::~TurboXIC() {
	delete m_mzVals;
}


Err TurboXIC::init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints) {
    ERR_INIT
    // e = d_ptr->init(scanNumberVsScanPoints); ree;
    e = initAVX(scanNumberVsScanPoints); ree;
    ERR_RETURN
}

// Err TurboXIC::init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints) const {
//     ERR_INIT
//     e = d_ptr->init(scanNumberVsScanPoints); ree;
//     ERR_RETURN
// }

namespace {

	ulong calculateNextAlignasSize(ulong vecSize) {
		while (vecSize % 32 != 0) {
			vecSize++;
		}
		return vecSize;
	}

}//namespace
Err TurboXIC::initAVX(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints) {
	ERR_INIT

	e = ErrorUtils::isNotEmpty(scanNumberVsScanPoints); ree;
	e = ErrorUtils::isEqual(
		static_cast<short>(sizeof(float)),
		static_cast<short>(sizeof(uint32_t))
		); ree;

	const ulong scanPointsCount = std::accumulate(
		scanNumberVsScanPoints.begin(),
		scanNumberVsScanPoints.end(),
		0,
		[](ulong sum, ScanPoints *sp){return sum + sp->size();}
		);

	e = ErrorUtils::isBelowThreshold(
		scanPointsCount,
		static_cast<ulong>(std::numeric_limits<uint32_t>::max()),
		ErrorUtilsParam::IncludeThreshold
		); ree;

	m_scanPointsCountAlignas = calculateNextAlignasSize(scanPointsCount);
	delete m_mzVals;
	auto* mzVals = static_cast<float*>(std::aligned_alloc(32, m_scanPointsCountAlignas * sizeof(float)));
	m_mzVals = mzVals;

	m_alignasPiecesOfEight = m_scanPointsCountAlignas / 8;

	delete m_indexesMz;
	auto* indexesMz = static_cast<float*>(std::aligned_alloc(32, m_alignasPiecesOfEight * sizeof(float)));
	m_indexesMz = indexesMz;

	delete m_indexesI;
	auto* indexesI = static_cast<uint32_t*>(std::aligned_alloc(32, m_alignasPiecesOfEight * sizeof(uint32_t)));
	m_indexesI = indexesI;

	m_xicPoints.reserve(scanPointsCount);
	for (auto it = scanNumberVsScanPoints.begin(); it != scanNumberVsScanPoints.end(); ++it) {
		const ScanNumber scanNumber = it.key();
		const ScanPoints *scanPoints = it.value();
		for(const ScanPoint &sp : *scanPoints) {
			XICPoint xicPoint;
			xicPoint.mz = sp.x();
			xicPoint.intensity = sp.y();
			xicPoint.scanNumber = scanNumber;
			xicPoint.ionMobilityIndex = -1;
			m_xicPoints.push_back(xicPoint);
		}
	}
	std::sort(
		m_xicPoints.begin(),
		m_xicPoints.end(),
		[](const XICPoint &l, const XICPoint &r){return l.mz < r.mz;}
		);

	int indexCounter = 0;
	for (int i = 0; i <  m_xicPoints.size(); i++) {
		m_mzVals[i] = m_xicPoints[i].mz;
		if (i % 8 == 0) {
			m_indexesMz[indexCounter] = m_xicPoints[i].mz;
			m_indexesI[indexCounter++] = i;
		}
	}

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

XICPointsPntrs TurboXIC::extractPointsXICAVX(float mzMin, float mzMax) {

	XICPointsPntrs results;
	results.reserve(1000);

	const __m256 mz_min_v = _mm256_set1_ps(mzMin);
	const __m256 mz_max_v = _mm256_set1_ps(mzMax);

	size_t left = 0, right = m_alignasPiecesOfEight;
	while (left < right) {
		size_t mid = left + (right - left) / 2;
		if (m_indexesMz[mid] < mzMin)
			left = mid + 1;
		else
			right = mid;
	}

	for (size_t i = m_indexesI[left & ~0x7]; i + 8 <= m_scanPointsCountAlignas; i += 8) {
		__m256 mz_v = _mm256_load_ps(&m_mzVals[i]);

		__m256 cmp_ge = _mm256_cmp_ps(mz_v, mz_min_v, _CMP_GE_OS);
		__m256 cmp_le = _mm256_cmp_ps(mz_v, mz_max_v, _CMP_LE_OS);
		__m256 mask = _mm256_and_ps(cmp_ge, cmp_le);

		int movemask = _mm256_movemask_ps(mask);

		if (movemask == 0) {
			__m256 cmp_gt_max = _mm256_cmp_ps(mz_v, mz_max_v, _CMP_GT_OS);
			int gt_mask = _mm256_movemask_ps(cmp_gt_max);
			if (gt_mask == 0xFF) {
				break;
			}
			continue;
		}

		for (int bit = 0; bit < 8; ++bit) {
			if (movemask & (1 << bit)) {
				results.push_back(&m_xicPoints[i + bit]);
			}
		}
	}

	return results;
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
	return !m_xicPoints.empty();
    // return d_ptr->isInit();
}

void TurboXIC::filterXICPointsByScanNumber(
    ScanNumber scanNumberMin,
    ScanNumber scanNumberMax,
    XICPointsPntrs* xicPoints
    ) {

    const auto terminatorLogic = [scanNumberMin, scanNumberMax](XICPoint *p) {
        return !(scanNumberMin <= p->scanNumber && p->scanNumber <= scanNumberMax);
    };

    const auto terminator = std::remove_if(xicPoints->begin(), xicPoints->end(), terminatorLogic);

    xicPoints->erase(terminator, xicPoints->end());
}
