//
// Created by andrewnichols on 7/26/25.
//

#include "TargetDecoyCandidatePairScoretronV2.h"

#include "AVXUtils.h"
#include "EigenUtils.h"
#include "EigenKernelUtils.h"
#include "MS2Ion.h"
#include "TargetDecoyCandidatePair.h"

TargetDecoyCandidatePairScoretronV2::TargetDecoyCandidatePairScoretronV2()
: m_mzTargetKeyCurrent("null")
, m_turboXICCurrent(nullptr)
, m_xicSizeMaxAlignas(-1)
, m_ms2IonsCount(-1)
, m_frameIndexTargetMax(-1)
, m_intensityVec(nullptr)
, m_ionCountVec(nullptr)
, m_integrationVecCosineSim(nullptr)
, m_productVec(nullptr)
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

	delete m_intensityVec;
	delete m_ionCountVec;
	delete m_integrationVecCosineSim;
	delete m_productVec;
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
		AVXUtils::AVX2_ALIGNAS_SIZE
		);
	m_ms2IonsCount = ms2IonsCount;

	m_xicsAlignasIntensity.resize(m_ms2IonsCount);
	m_xicsAlignasIntensityShadows.resize(m_ms2IonsCount);
	m_xicsAlignasMz.resize(m_ms2IonsCount);
	m_xicsAlignasMzShadows.resize(m_ms2IonsCount);
	for (int i = 0; i < m_ms2IonsCount; i++) {

		auto* alignedIntensityVals = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasIntensity[i] = alignedIntensityVals;

		auto* alignedMzVals = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasMz[i] = alignedMzVals;

		auto* alignedIntensityValsShadows = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasIntensityShadows[i] = alignedIntensityValsShadows;

		auto* alignedMzValsShadows = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasMzShadows[i] = alignedMzValsShadows;
	}

	auto* alignedIntensityVec = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_intensityVec = alignedIntensityVec;

	auto* alignedIonCountVec = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_ionCountVec = alignedIonCountVec;

	auto* alignedIntegrationVecCosineSim = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_integrationVecCosineSim = alignedIntegrationVecCosineSim;

	auto* alignedProductVec = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_productVec = alignedProductVec;

	constexpr int order = 1;
	constexpr int derivative = 0;
	constexpr int rate = 1;
	Eigen::MatrixX<float> kernel;
	e = EigenKernelUtils::buildSavitzkyGolayKernel(
		m_pythiaParameters.filterLengthMS2,
		order,
		derivative,
		rate,
		&kernel
		); ree;
	const Eigen::VectorX<float> kernelVec(kernel);
	m_savitzkyGolayKernel = EigenUtils::convertEigenVectorToQVector(kernelVec);

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

	if (m_frameIndexTargetMax < 0) {
		//TODO figure out what to return here. Most likely blank candidate score
		ERR_RETURN
	}

	if (m_pythiaParameters.subtractShadows) {
		e = subtractShadowsArrays(); ree;
	}

	e = smoothMS2IonArrays(); ree;



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

Err TargetDecoyCandidatePairScoretronV2::subtractShadowsArrays() {

	ERR_INIT

	const size_t xicSizeTargetMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		m_frameIndexTargetMax,
		AVXUtils::AVX2_ALIGNAS_SIZE
		);

	for (int i = 0; i < m_ms2IonsCount; i++) {
		e = AVXUtils::subtractArraysAVX2(
			m_xicsAlignasIntensity[i],
			m_xicsAlignasIntensityShadows[i],
			xicSizeTargetMaxAlignas
			); ree;
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

	std::memset(m_intensityVec, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_ionCountVec, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_integrationVecCosineSim, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_productVec, 0, m_xicSizeMaxAlignas * sizeof(float));
}

Err TargetDecoyCandidatePairScoretronV2::smoothMS2IonArrays() {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(m_frameIndexTargetMax); ree;
	e = ErrorUtils::isNotEmpty(m_savitzkyGolayKernel); ree;

	const size_t xicSizeTargetMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		m_frameIndexTargetMax,
		AVXUtils::AVX2_ALIGNAS_SIZE
		);

	for (int i = 0; i < m_pythiaParameters.smoothCountMS2; i++) {

		e = AVXUtils::convolveWithKernelAVXFloat(
			m_savitzkyGolayKernel,
			xicSizeTargetMaxAlignas,
			m_xicsAlignasIntensity[0],
			m_xicsAlignasIntensity[1],
			m_xicsAlignasIntensity[2],
			m_xicsAlignasIntensity[3],
			m_xicsAlignasIntensity[4],
			m_xicsAlignasIntensity[5],
			m_xicsAlignasIntensity[6],
			m_xicsAlignasIntensity[7]
			); ree;

		if (m_ms2IonsCount == S_GLOBAL_SETTINGS.MAX_MS2_IONS) {
			e = AVXUtils::convolveWithKernelAVXFloat(
					m_savitzkyGolayKernel,
					xicSizeTargetMaxAlignas,
					m_xicsAlignasIntensity[8],
					m_xicsAlignasIntensity[9],
					m_xicsAlignasIntensity[10],
					m_xicsAlignasIntensity[11],
					m_xicsAlignasIntensity[12],
					m_xicsAlignasIntensity[13],
					m_xicsAlignasIntensity[14],
					m_xicsAlignasIntensity[15]
					); ree;
		}
	}

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::buildLocationVectors() {
	ERR_INIT


	ERR_RETURN
}
