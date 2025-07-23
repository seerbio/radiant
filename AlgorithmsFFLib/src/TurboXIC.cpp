//
// Created by anichols on 12/18/22.
//

#include "TurboXIC.h"

#include "AVXUtils.h"
#include "ErrorUtils.h"
#include "MsUtils.h"

#include <immintrin.h>


TurboXIC::TurboXIC()
: m_mzVals(nullptr)
, m_indexesMz(nullptr)
, m_indexesI(nullptr)
, m_alignasPiecesOfEight(0)
, m_scanPointsCountAlignas(0)
{}

TurboXIC::~TurboXIC() {
	delete m_mzVals;
	delete m_indexesMz;
	delete m_indexesI;
}

Err TurboXIC::init(const QVector<MsScan> &msScans) {
	ERR_INIT

	e = ErrorUtils::isNotEmpty(msScans); ree;
	e = ErrorUtils::isEqual(
		static_cast<short>(sizeof(float)),
		static_cast<short>(sizeof(uint32_t))
		); ree;

	const ulong scanPointsCount = std::accumulate(
		msScans.begin(),
		msScans.end(),
		0,
		[](ulong sum, const MsScan &s){return sum + s.msScanInfoPntr->pointCount;}
		);

	e = ErrorUtils::isBelowThreshold(
		scanPointsCount,
		static_cast<ulong>(std::numeric_limits<uint32_t>::max()),
		ErrorUtilsParam::IncludeThreshold
		); ree;

	m_scanPointsCountAlignas = AVXUtils::calculateNextAlignedBlockSize(
		scanPointsCount,
		AVXUtils::AVX2_ALIGNAS_SIZE
		);

	delete m_mzVals;
	auto* mzVals = static_cast<float*>(std::aligned_alloc(AVXUtils::AVX2_ALIGNAS_SIZE, m_scanPointsCountAlignas * sizeof(float)));
	m_mzVals = mzVals;

	m_alignasPiecesOfEight = m_scanPointsCountAlignas / 8;

	delete m_indexesMz;
	auto* indexesMz = static_cast<float*>(std::aligned_alloc(AVXUtils::AVX2_ALIGNAS_SIZE, m_alignasPiecesOfEight * sizeof(float)));
	m_indexesMz = indexesMz;

	delete m_indexesI;
	auto* indexesI = static_cast<uint32_t*>(std::aligned_alloc(AVXUtils::AVX2_ALIGNAS_SIZE, m_alignasPiecesOfEight * sizeof(uint32_t)));
	m_indexesI = indexesI;

	m_xicPoints.reserve(scanPointsCount);
	for (const MsScan &s : msScans) {

		for(int i = 0; i < s.msScanInfoPntr->pointCount; i++) {
			XICPoint xicPoint;
			xicPoint.mz = s.mzVals[i];
			xicPoint.intensity = s.intensityVals[i];
			xicPoint.scanNumber = s.msScanInfoPntr->scanNumber;
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

Err TurboXIC::init(const QString &filePath) {

	ERR_INIT

	QFile file(filePath);
	if (!file.open(QIODevice::ReadOnly)) {
		qWarning() << "Failed to open file for reading:" << filePath;
		rrr(eFileError);
	}

	QDataStream in(&file);
	in.setVersion(QDataStream::Qt_5_15);  // Ensure compatibility with the saved version

	in >> m_scanPointsCountAlignas;
	in >> m_alignasPiecesOfEight;

	delete m_mzVals;
	auto* mzVals
		= static_cast<float*>(std::aligned_alloc(AVXUtils::AVX2_ALIGNAS_SIZE, m_scanPointsCountAlignas * sizeof(float)));
	m_mzVals = mzVals;

	delete m_indexesMz;
	auto* indexesMz
		= static_cast<float*>(std::aligned_alloc(AVXUtils::AVX2_ALIGNAS_SIZE, m_alignasPiecesOfEight * sizeof(float)));
	m_indexesMz = indexesMz;

	delete m_indexesI;
	auto* indexesI
		= static_cast<uint32_t*>(std::aligned_alloc(AVXUtils::AVX2_ALIGNAS_SIZE, m_alignasPiecesOfEight * sizeof(uint32_t)));
	m_indexesI = indexesI;

	QVector<float> mzValsLoad;
	in >> mzValsLoad;
	AVXUtils::copyAVXFloatToAligned(mzValsLoad.data(), m_mzVals, m_scanPointsCountAlignas);

	QVector<float> indexesMzLoad;
	in >> indexesMzLoad;
	AVXUtils::copyAVXFloatToAligned(indexesMzLoad.data(), m_indexesMz, m_alignasPiecesOfEight);

	QVector<uint32_t> indexesILoad;
	in >> indexesILoad;
	AVXUtils::copyAVXIntToAligned(indexesILoad.data(), m_indexesI, m_alignasPiecesOfEight);

	in >> m_xicPoints;

	file.close();
	// qDebug() << "TurboXIC data loaded successfully from" << filePath;

	ERR_RETURN
}

XICPointsPntrs TurboXIC::extractPointsXIC(float mzMin, float mzMax) {

	XICPointsPntrs results;

	constexpr int reserveEstimate = 1000;
	results.reserve(reserveEstimate);

	const __m256 mz_min_v = _mm256_set1_ps(mzMin);
	const __m256 mz_max_v = _mm256_set1_ps(mzMax);

	size_t left = 0;
	size_t right = m_alignasPiecesOfEight;
	while (left < right) {
		size_t mid = left + (right - left) / 2;
		if (m_indexesMz[mid] < mzMin) {
			left = mid + 1;
		}
		else {
			right = mid;
		}
	}

	for (
		size_t i = m_indexesI[left];
		i + AVXUtils::AVX2_FLOAT_REGISTER_SIZE <= m_scanPointsCountAlignas;
		i += AVXUtils::AVX2_FLOAT_REGISTER_SIZE
		) {

		__m256 mz_v = _mm256_load_ps(&m_mzVals[i]);

		__m256 cmp_ge = _mm256_cmp_ps(mz_v, mz_min_v, _CMP_GE_OS);
		__m256 cmp_le = _mm256_cmp_ps(mz_v, mz_max_v, _CMP_LE_OS);
		__m256 mask = _mm256_and_ps(cmp_ge, cmp_le);

		const int movemask = _mm256_movemask_ps(mask);

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

bool TurboXIC::isInit() const {
	return !m_xicPoints.empty();
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


Err TurboXIC::write(const QString &filePath) const {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(filePath); ree;
	e = ErrorUtils::isNotEmpty(m_xicPoints); ree;

	QFile file(filePath);
	if (!file.open(QIODevice::WriteOnly)) {
		qWarning() << "Failed to open file for writing:" << filePath;
		rrr(eFileError);
	}

	QDataStream out(&file);
	out.setVersion(QDataStream::Qt_5_15);

	QVector<float> mzValsVector(m_mzVals, m_mzVals + m_scanPointsCountAlignas);
	QVector<float> indexesMz(m_indexesMz, m_indexesMz + m_alignasPiecesOfEight);
	QVector<int> indexesI(m_indexesI, m_indexesI + m_alignasPiecesOfEight);

	out << m_scanPointsCountAlignas;
	out << m_alignasPiecesOfEight;
	out << mzValsVector;
	out << indexesMz;
	out << indexesI;
	out << m_xicPoints;

	file.close();
	// qDebug() << "Turbo Data saved successfully to" << filePath;

	ERR_RETURN
}
