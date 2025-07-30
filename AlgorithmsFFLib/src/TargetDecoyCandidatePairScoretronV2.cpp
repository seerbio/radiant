//
// Created by andrewnichols on 7/26/25.
//

#include "TargetDecoyCandidatePairScoretronV2.h"

#include "AVXUtils.h"
#include "EigenUtils.h"
#include "EigenKernelUtils.h"
#include "MS2Ion.h"
#include "MsReaderMzMLLazyLoad.h"
#include "TargetDecoyCandidatePair.h"

TargetDecoyCandidatePairScoretronV2::TargetDecoyCandidatePairScoretronV2()
: m_msReaderPointerAcc(nullptr)
, m_mzTargetKeyCurrent("null")
, m_msFrameV2Current(nullptr)
, m_xicSizeMaxAlignas(-1)
, m_ms2IonsCount(-1)
, m_intensityVec(nullptr)
, m_ionCountVec(nullptr)
, m_integrationVecCosineSim(nullptr)
, m_productVec(nullptr)
, m_targetFrameIndexMax(-1)
, m_mzMs2Min(-1)
, m_mzMs2Max(-1)
, m_intensityVecMax(0)
, m_minMs2IonsFoundCount(-1)
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

	delete m_msFrameV2Current;
}

Err TargetDecoyCandidatePairScoretronV2::init(
	const QMap<MzTargetKey, QVector<MsScanInfo*>> &mzTargetKeyVsMsScanInfos,
	const PythiaParameters &pythiaParameters,
	int ms2IonsCount,
	float minMs2IonsFoundCount,
	MsReaderPointerAcc *msReaderPointerAcc
	) {
	ERR_INIT

	e = ErrorUtils::isNotEmpty(mzTargetKeyVsMsScanInfos); ree;
	e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); ree;
	e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
	e = ErrorUtils::isGreaterThanZero(ms2IonsCount); ree;
	e = ErrorUtils::isGreaterThanZero(minMs2IonsFoundCount); ree;

	m_msReaderPointerAcc = msReaderPointerAcc;
	m_mzTargetKeyVsMsScanInfos = mzTargetKeyVsMsScanInfos;
	m_pythiaParameters = pythiaParameters;
	m_xicSizeMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		m_msReaderPointerAcc->ptr->meanFrameScanCountMS2(),
		AVXUtils::AVX2_ALIGNAS_SIZE
		);
	m_ms2IonsCount = ms2IonsCount;
	m_mzMs2Min = m_msReaderPointerAcc->ptr->mzMs2Min();
	m_mzMs2Max = m_msReaderPointerAcc->ptr->mzMs2Max();
	m_minMs2IonsFoundCount = minMs2IonsFoundCount;

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

Err TargetDecoyCandidatePairScoretronV2::scoreTargetDecoyCandidatePairPntr(
	const QPair<MzTargetKey, TargetDecoyCandidatePair*> &mzTargetKeyVsTdcpPntr
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(m_msReaderPointerAcc->isInit()); ree;

	const MzTargetKey &mzTargetKey = mzTargetKeyVsTdcpPntr.first;
	if (mzTargetKey != m_mzTargetKeyCurrent) {

		delete m_msFrameV2Current;

		const QVector<MsScanInfo*> &msScanInfosCopy = m_mzTargetKeyVsMsScanInfos.value(mzTargetKey);
		QVector<MsScan> msScans;
		e = m_msReaderPointerAcc->ptr->extractScanPoints(
			msScanInfosCopy,
			&msScans
			); rtee;

		e = ErrorUtils::isTrue(
			msScans.front().msScanInfoPntr->targetKey() ==
			msScans.back().msScanInfoPntr->targetKey()
			); rtee;

		auto* msFrame = new MsFrameV2();
		e = msFrame->init(msScans); rtee;
		m_msFrameV2Current = msFrame;

		e = ErrorUtils::isTrue(m_msFrameV2Current->isInit()); ree;
		m_mzTargetKeyCurrent = mzTargetKey;

	}

	TargetDecoyCandidatePair *tdcpPntr = mzTargetKeyVsTdcpPntr.second;

	constexpr float mzMs2MinMin = 200.0f;
	QVector<MS2Ion> ms2IonsTarget = tdcpPntr->ms2IonsTarget(
		std::max(m_mzMs2Min, mzMs2MinMin),
		m_mzMs2Max
		);
	const int ms2IonsTargetSizeOG = ms2IonsTarget.size(); //TODO, add this to feature vec in CandScores
	ms2IonsTarget.resize(std::min(ms2IonsTarget.size(), m_ms2IonsCount));

	const QVector<MS2Ion> &ms2IonsDecoy = tdcpPntr->ms2IonsDecoy(ms2IonsTarget);

	e = scoreMS2Ions(ms2IonsTarget); ree;
	e = scoreMS2Ions(ms2IonsDecoy); ree;

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::scoreMS2Ions(const QVector<MS2Ion> &ms2Ions) {

	ERR_INIT

	e = loadMS2IonArrays(ms2Ions); ree;

	if (m_targetFrameIndexMax < 0) {
		//TODO figure out what to return here. Most likely blank candidate score
		ERR_RETURN
	}

	if (m_pythiaParameters.subtractShadows) {
		e = subtractShadowsArrays(); ree;
	}

	e = smoothMS2IonArrays(); ree;

	e = buildLocationVectors(ms2Ions); ree;

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::loadMS2IonArrays(const QVector<MS2Ion> &ms2Ions) {

	ERR_INIT

	m_targetFrameIndexMax = -1;
	m_intensityVecMax = 0;
	zeroOutArrays();

	int arrayIndex = 0;
	for (const MS2Ion &ms2Ion : ms2Ions) {

		const float mzTol = MathUtils::calculatePPM(
			ms2Ion.mz,
			static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
			);

		const float mzMin = ms2Ion.mz - mzTol;
		const float mzMax = ms2Ion.mz + mzTol;
		XICPointsPntrs xicPointsPntrs;
		e = m_msFrameV2Current->getTurboXICPntr()->extractPointsXIC(
			mzMin,
			mzMax,
			&xicPointsPntrs
			); ree;

		const auto isotopeDistanceThomsons = static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge);
		const float mzShadow = ms2Ion.mz - isotopeDistanceThomsons;
		const float mzShadowTol = MathUtils::calculatePPM(
			mzShadow,
			static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
			);

		const float mzShadowMin = mzShadow - mzShadowTol;
		const float mzShadowMax = mzShadow + mzShadowTol;
		XICPointsPntrs xicPointsPntrsShadows;
		e = m_msFrameV2Current->getTurboXICPntr()->extractPointsXIC(
			mzShadowMin,
			mzShadowMax,
			&xicPointsPntrsShadows
			); ree;

		float* arrIntensity = m_xicsAlignasIntensity[arrayIndex];
		float* arrMz = m_xicsAlignasMz[arrayIndex];
		for (const XICPoint *xicPoint : xicPointsPntrs) {

			m_targetFrameIndexMax = std::max(m_targetFrameIndexMax, xicPoint->frameIndex);

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
		m_targetFrameIndexMax,
		AVXUtils::AVX2_ALIGNAS_SIZE
		);

	for (int i = 0; i < m_ms2IonsCount; i++) {
		constexpr bool zeroNegatives = true;
		e = AVXUtils::subtractArraysAVX2(
			m_xicsAlignasIntensity[i],
			m_xicsAlignasIntensityShadows[i],
			xicSizeTargetMaxAlignas,
			zeroNegatives
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

	e = ErrorUtils::isGreaterThanZero(m_targetFrameIndexMax); ree;
	e = ErrorUtils::isNotEmpty(m_savitzkyGolayKernel); ree;

	const size_t xicSizeTargetMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		m_targetFrameIndexMax,
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

Err TargetDecoyCandidatePairScoretronV2::buildLocationVectors(const QVector<MS2Ion> &ms2Ions) {

	ERR_INIT

	for (int i = 0; i < m_targetFrameIndexMax; i += AVXUtils::AVX2_FLOAT_REGISTER_SIZE) {

		__m256 v0 = _mm256_load_ps(m_xicsAlignasIntensity[0] + i);
		__m256 v1 = _mm256_load_ps(m_xicsAlignasIntensity[1] + i);
		__m256 v2 = _mm256_load_ps(m_xicsAlignasIntensity[2] + i);
		__m256 v3 = _mm256_load_ps(m_xicsAlignasIntensity[3] + i);
		__m256 v4 = _mm256_load_ps(m_xicsAlignasIntensity[4] + i);
		__m256 v5 = _mm256_load_ps(m_xicsAlignasIntensity[5] + i);
		__m256 v6 = _mm256_load_ps(m_xicsAlignasIntensity[6] + i);
		__m256 v7 = _mm256_load_ps(m_xicsAlignasIntensity[7] + i);

		__m256 sum1 = _mm256_add_ps(v0, v1);
		__m256 sum2 = _mm256_add_ps(v2, v3);
		__m256 sum3 = _mm256_add_ps(v4, v5);
		__m256 sum4 = _mm256_add_ps(v6, v7);
		__m256 intermediate1 = _mm256_add_ps(sum1, sum2);
		__m256 intermediate2 = _mm256_add_ps(sum3, sum4);
		__m256 finalSum1 = _mm256_add_ps(intermediate1, intermediate2);

		constexpr float thresholdValue = 0.1;
		constexpr float replaceValue = 1.0;

		AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v0);
		AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v1);
		AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v2);
		AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v3);
		AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v4);
		AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v5);
		AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v6);
		AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v7);

		__m256 sumCount1 = _mm256_add_ps(v0, v1);
		__m256 sumCount2 = _mm256_add_ps(v2, v3);
		__m256 sumCount3 = _mm256_add_ps(v4, v5);
		__m256 sumCount4 = _mm256_add_ps(v6, v7);
		__m256 intermediateCount1 = _mm256_add_ps(sumCount1, sumCount2);
		__m256 intermediateCount2 = _mm256_add_ps(sumCount3, sumCount4);
		__m256 finalSumCount1 = _mm256_add_ps(intermediateCount1, intermediateCount2);

		if (m_ms2IonsCount == S_GLOBAL_SETTINGS.MAX_MS2_IONS) {
			__m256 v8 = _mm256_load_ps(m_xicsAlignasIntensity[8] + i);
			__m256 v9 = _mm256_load_ps(m_xicsAlignasIntensity[9] + i);
			__m256 v10 = _mm256_load_ps(m_xicsAlignasIntensity[10] + i);
			__m256 v11 = _mm256_load_ps(m_xicsAlignasIntensity[11] + i);
			__m256 v12 = _mm256_load_ps(m_xicsAlignasIntensity[12] + i);
			__m256 v13 = _mm256_load_ps(m_xicsAlignasIntensity[13] + i);
			__m256 v14 = _mm256_load_ps(m_xicsAlignasIntensity[14] + i);
			__m256 v15 = _mm256_load_ps(m_xicsAlignasIntensity[15] + i);

			__m256 sum5 = _mm256_add_ps(v8, v9);
			__m256 sum6 = _mm256_add_ps(v10, v11);
			__m256 sum7 = _mm256_add_ps(v12, v13);
			__m256 sum8 = _mm256_add_ps(v14, v15);
			__m256 intermediate3 = _mm256_add_ps(sum5, sum6);
			__m256 intermediate4 = _mm256_add_ps(sum7, sum8);
			__m256 finalSum2 = _mm256_add_ps(intermediate3, intermediate4);
			__m256 finalSumForReal = _mm256_add_ps(finalSum1, finalSum2);

			AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v8);
			AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v9);
			AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v10);
			AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v11);
			AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v12);
			AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v13);
			AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v14);
			AVXUtils::replaceArrayValuesAVXGreaterThan(thresholdValue, replaceValue, v15);

			__m256 sumCount5 = _mm256_add_ps(v8, v9);
			__m256 sumCount6 = _mm256_add_ps(v10, v11);
			__m256 sumCount7 = _mm256_add_ps(v12, v13);
			__m256 sumCount8 = _mm256_add_ps(v14, v15);
			__m256 intermediateCount3 = _mm256_add_ps(sumCount5, sumCount6);
			__m256 intermediateCount4 = _mm256_add_ps(sumCount7, sumCount8);
			__m256 finalSumCount2 = _mm256_add_ps(intermediateCount3, intermediateCount4);
			__m256 finalSumForRealCount = _mm256_add_ps(finalSumCount1, finalSumCount2);

			m_intensityVecMax = std::max(m_intensityVecMax, AVXUtils::maxFloat(finalSumForReal));
			_mm256_store_ps(m_intensityVec + i, finalSumForReal);
			_mm256_store_ps(m_ionCountVec + i, finalSumForRealCount);
		}
		else {
			m_intensityVecMax = std::max(m_intensityVecMax, AVXUtils::maxFloat(finalSum1));
			_mm256_store_ps(m_intensityVec + i, finalSum1);
			_mm256_store_ps(m_ionCountVec + i, finalSumCount1);
		}

	}

	e = buildIntegrationVecCosineSim(ms2Ions); ree;
	e = buildProductVec(); ree;

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::buildIntegrationVecCosineSim(const QVector<MS2Ion> &ms2Ions) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(ms2Ions); ree;
	e = ErrorUtils::isGreaterThanZero(m_targetFrameIndexMax); ree;

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float scanIntensities[m_ms2IonsCount] = {0};
	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float theoIntensities[m_ms2IonsCount] = {0};
	for (int i = 0; i < std::min(m_ms2IonsCount, ms2Ions.size()); i++) {
		theoIntensities[i] = ms2Ions[i].intensity;
	}

	for (int i = 0; i < m_targetFrameIndexMax; i++) {
		for (int j = 0; j < std::min(m_ms2IonsCount, ms2Ions.size()); j++) {
			scanIntensities[j] = m_xicsAlignasIntensity[j][i];
		}

		const float cosineSim = AVXUtils::cosineSimilarityAVX(scanIntensities, theoIntensities, m_ms2IonsCount);
		m_integrationVecCosineSim[i] = cosineSim;

	}

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::buildProductVec() const {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(m_targetFrameIndexMax); ree;

	const size_t xicSizeTargetMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		m_targetFrameIndexMax,
		AVXUtils::AVX2_ALIGNAS_SIZE
		);

	for (int i = 0; i < xicSizeTargetMaxAlignas; i += AVXUtils::AVX2_FLOAT_REGISTER_SIZE) {

		__m256 ionCount = _mm256_load_ps(m_ionCountVec + i);
		const __m256 thresholds = _mm256_set1_ps(m_minMs2IonsFoundCount);
		const __m256 mask = _mm256_cmp_ps(ionCount, thresholds, _CMP_LT_OQ);

		if (AVXUtils::isAllOnes(mask)) {
			continue;
		}

		ionCount = _mm256_blendv_ps(ionCount, _mm256_setzero_ps(), mask);
		const __m256 cosineSim = _mm256_load_ps(m_integrationVecCosineSim + i);
		const __m256 cosineSimSqrt = _mm256_sqrt_ps(cosineSim);

		const __m256 product = _mm256_mul_ps(ionCount, cosineSimSqrt);
		_mm256_store_ps(m_productVec + i, product);
	}

	ERR_RETURN
}
