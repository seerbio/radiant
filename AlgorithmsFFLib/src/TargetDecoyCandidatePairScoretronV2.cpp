//
// Created by andrewnichols on 7/26/25.
//

#include "TargetDecoyCandidatePairScoretronV2.h"

#include "AVXUtils.h"
#include "MS2Ion.h"
#include "TargetDecoyCandidatePair.h"

TargetDecoyCandidatePairScoretronV2::TargetDecoyCandidatePairScoretronV2()
: m_mzTargetKeyCurrent("null")
, m_turboXICCurrent(nullptr)
, m_xicSizeMaxAlignas(-1)
, m_ms2IonsCount(-1)
, m_frameIndexTargetMax(-1)
{}

TargetDecoyCandidatePairScoretronV2::~TargetDecoyCandidatePairScoretronV2() {
	for (float* v : m_xicsAlignasIntensity) {
		delete v;
	}

	for (float* v : m_xicsAlignasIntensityShadows) {
		delete v;
	}

	for (float* v : m_xicsAlignasMz) {
		delete v;
	}

	for (float* v : m_xicsAlignasMzShadows) {
		delete v;
	}
}

Err TargetDecoyCandidatePairScoretronV2::init(
	const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXICs,
	const PythiaParameters &pythiaParameters,
	int meanFrameScanCountMS2,
	int ms2IonsCount
	) {
	ERR_INIT

	e = ErrorUtils::isNotEmpty(mzTargetKeyVsTurboXICs); ree;
	e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
	e = ErrorUtils::isGreaterThanZero(meanFrameScanCountMS2); ree;
	e = ErrorUtils::isGreaterThanZero(ms2IonsCount); ree;

	m_mzTargetKeyVsTurboXICs = mzTargetKeyVsTurboXICs;
	m_pythiaParameters = pythiaParameters;
	m_xicSizeMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		meanFrameScanCountMS2,
		S_GLOBAL_SETTINGS.AVX2_BYTE_SIZE_ALIGNAS
		);
	m_ms2IonsCount = ms2IonsCount;

	m_xicsAlignasIntensity.resize(m_ms2IonsCount);
	m_xicsAlignasIntensityShadows.resize(m_ms2IonsCount);
	m_xicsAlignasMz.resize(m_ms2IonsCount);
	m_xicsAlignasMzShadows.resize(m_ms2IonsCount);
	for (int i = 0; i < m_ms2IonsCount; i++) {

		auto* alignedIntensityVals = static_cast<float*>(std::aligned_alloc(
			S_GLOBAL_SETTINGS.AVX2_BYTE_SIZE_ALIGNAS,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasIntensity[i] = alignedIntensityVals;

		auto* alignedMzVals = static_cast<float*>(std::aligned_alloc(
			S_GLOBAL_SETTINGS.AVX2_BYTE_SIZE_ALIGNAS,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasMz[i] = alignedMzVals;

		auto* alignedIntensityValsShadows = static_cast<float*>(std::aligned_alloc(
			S_GLOBAL_SETTINGS.AVX2_BYTE_SIZE_ALIGNAS,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasIntensityShadows[i] = alignedIntensityValsShadows;

		auto* alignedMzValsShadows = static_cast<float*>(std::aligned_alloc(
			S_GLOBAL_SETTINGS.AVX2_BYTE_SIZE_ALIGNAS,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasMzShadows[i] = alignedMzValsShadows;
	}

	ERR_RETURN
}

namespace {



}//namespace
Err TargetDecoyCandidatePairScoretronV2::scoreTargetDecoyCandidatePairPntr(
	const QPair<MzTargetKey, TargetDecoyCandidatePair*> &mzTargetKeyVsTdcpPntr
	) {

	ERR_INIT

	const MzTargetKey &mzTargetKey = mzTargetKeyVsTdcpPntr.first;
	if (mzTargetKey != m_mzTargetKeyCurrent) {
		m_turboXICCurrent = m_mzTargetKeyVsTurboXICs.value(mzTargetKey);
		e = ErrorUtils::isTrue(m_turboXICCurrent->isInit()); ree;
		m_mzTargetKeyCurrent = mzTargetKey;
	}

	TargetDecoyCandidatePair *tdcpPntr = mzTargetKeyVsTdcpPntr.second;

	QVector<MS2Ion> ms2IonsTarget = tdcpPntr->ms2IonsTarget();

	ms2IonsTarget.resize(std::min(ms2IonsTarget.size(), m_ms2IonsCount));

	const QVector<MS2Ion> &ms2IonsDecoy = tdcpPntr->ms2IonsDecoy(ms2IonsTarget);

	e = scoreMS2Ions(ms2IonsTarget); ree;
	e = scoreMS2Ions(ms2IonsDecoy); ree;

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::scoreMS2Ions(const QVector<MS2Ion> &ms2Ions) {

	ERR_INIT

	e = loadMS2IonArrays(ms2Ions); ree;



	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::loadMS2IonArrays(const QVector<MS2Ion> &ms2Ions) {

	ERR_INIT

	m_frameIndexTargetMax = -1;
	zeroOutArrays();

	int arrayIndex = 0;
	for (const MS2Ion &ms2Ion : ms2Ions) {

		const float mzTol = MathUtils::calculatePPM(
			ms2Ion.mz,
			static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
			);

		const float mzMin = ms2Ion.mz - mzTol;
		const float mzMax = ms2Ion.mz + mzTol;
		const XICPointsPntrs xicPointsPntrs = m_turboXICCurrent->extractPointsXIC(
			mzMin,
			mzMax
			);

		const auto isotopeDistanceThomsons = static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge);
		const float mzShadow = ms2Ion.mz - isotopeDistanceThomsons;
		const float mzShadowTol = MathUtils::calculatePPM(
			mzShadow,
			static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
			);

		const float mzShadowMin = mzShadow - mzShadowTol;
		const float mzShadowMax = mzShadow + mzShadowTol;
		const XICPointsPntrs xicPointsPntrsShadows = m_turboXICCurrent->extractPointsXIC(
			mzShadowMin,
			mzShadowMax
			);

		float* arrIntensity = m_xicsAlignasIntensity[arrayIndex];
		float* arrMz = m_xicsAlignasMz[arrayIndex];
		for (const XICPoint *xicPoint : xicPointsPntrs) {

			m_frameIndexTargetMax = std::max(m_frameIndexTargetMax, xicPoint->frameIndex);

			arrIntensity[xicPoint->frameIndex] += xicPoint->intensity;

			if (arrMz[xicPoint->frameIndex] > 1) {
				arrMz[xicPoint->frameIndex] += xicPoint->mz;
				arrMz[xicPoint->frameIndex] * 0.5;
				continue;
			}

			arrMz[xicPoint->frameIndex] = xicPoint->mz;
		}

		float* arrIntensityShadows = m_xicsAlignasIntensityShadows[arrayIndex];
		float* arrMzShadows = m_xicsAlignasMzShadows[arrayIndex];
		for (const XICPoint *xicPoint : xicPointsPntrsShadows) {

			arrIntensityShadows[xicPoint->frameIndex] += xicPoint->intensity;

			if (arrMzShadows[xicPoint->frameIndex] > 1) {
				arrMzShadows[xicPoint->frameIndex] += xicPoint->mz;
				arrMzShadows[xicPoint->frameIndex] * 0.5;
				continue;
			}

			arrMzShadows[xicPoint->frameIndex] = xicPoint->mz;
		}

		arrayIndex++;
	}

	ERR_RETURN
}

void TargetDecoyCandidatePairScoretronV2::zeroOutArrays() {

	for (float *f : m_xicsAlignasIntensity) {
		std::memset(f, 0, m_xicSizeMaxAlignas * sizeof(float));
	}

	for (float *f : m_xicsAlignasIntensityShadows) {
		std::memset(f, 0, m_xicSizeMaxAlignas * sizeof(float));
	}

	for (float *f : m_xicsAlignasMz) {
		std::memset(f, 0, m_xicSizeMaxAlignas * sizeof(float));
	}

	for (float *f : m_xicsAlignasMzShadows) {
		std::memset(f, 0, m_xicSizeMaxAlignas * sizeof(float));
	}
}

Err TargetDecoyCandidatePairScoretronV2::smoothMS2IonArrays() {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(m_frameIndexTargetMax);

	int m_xicSizeTargetMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		m_frameIndexTargetMax,
		S_GLOBAL_SETTINGS.AVX2_BYTE_SIZE_ALIGNAS
		);


	ERR_RETURN
}
