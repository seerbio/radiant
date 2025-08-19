//
// Created by andrewnichols on 8/14/25.
//

#include "CandidateScorertronV2.h"

#include "AVXUtils.h"
#include "PeakIntegratomatic.h"


CandidateScorertronV2::CandidateScorertronV2()
: m_ms2IonsCount(-1)
, m_integrationArraySizeMax(-1)
, m_productVecIntegration(nullptr)
, m_peakLength(-1)
{}

CandidateScorertronV2::~CandidateScorertronV2() {
	for (float* f : m_xicsAlignasIntensityIntegration) {
		delete f;
	}
	for (float* f : m_xicsAlignasIntensityIntegrationTight1) {
		delete f;
	}
	delete m_productVecIntegration;
}

Err CandidateScorertronV2::init(
	const QVector<CandidateScoresFeatureManager::Features> &featuresCalibration,
	int ms2IonCount
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(featuresCalibration); ree;
	e = ErrorUtils::isGreaterThanZero(ms2IonCount); ree;

	m_features = featuresCalibration;
	m_ms2IonsCount = ms2IonCount;
	m_xicsAlignasIntensityIntegration.resize(m_ms2IonsCount);
	m_xicsAlignasIntensityIntegrationTight1.resize(m_ms2IonsCount);

	m_integrationArraySizeMax
		= AVXUtils::calculateNextAlignedBlockSize(100, AVXUtils::AVX2_ALIGNAS_SIZE);;

	for (const CandidateScoresFeatureManager::Features &f : m_features) {
		m_featuresVsFeaturesArrayIndex[f] = m_featuresVsFeaturesArrayIndex.size();
	}

	for (int i = 0; i < m_ms2IonsCount; i++) {
		auto* alignedIntensityIntegrationVals = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_integrationArraySizeMax * sizeof(float))
			);
		m_xicsAlignasIntensityIntegration[i] = alignedIntensityIntegrationVals;

		auto* alignedIntensityIntegrationValsTight1 = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_integrationArraySizeMax * sizeof(float))
			);
		m_xicsAlignasIntensityIntegrationTight1[i] = alignedIntensityIntegrationValsTight1;
	}

	auto* productIntegrationVec = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_integrationArraySizeMax * sizeof(float))
		);
	m_productVecIntegration = productIntegrationVec;

	ERR_RETURN
}

bool CandidateScorertronV2::isInit() const {
	return !m_featuresVsFeaturesArrayIndex.isEmpty();
}

Err CandidateScorertronV2::scoreCandidate(
	const PeakIntegrationIndexes &pii,
	const QVector<MS2Ion> &ms2IonsFull,
	const QVector<float*> &xicsAlignasIntensity,
	const QVector<float*> &xicsAlignasIntensityTight1,
	float* productVec,
	CandidateScoresV2 *candidateScores
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(ms2IonsFull); ree;
	e = ErrorUtils::isNotEmpty(m_features); ree;
	e = ErrorUtils::isNotEmpty(m_featuresVsFeaturesArrayIndex); ree;

	QVector<MS2Ion> ms2IonsTrunc = ms2IonsFull;
	ms2IonsTrunc.resize(m_ms2IonsCount);

	e = copyToPeakVecs(
		pii,
		xicsAlignasIntensity,
		xicsAlignasIntensityTight1,
		productVec
		); ree;

	e = calculateScores(candidateScores); ree;

	ERR_RETURN
}

void CandidateScorertronV2::zeroOutArrays() {

	m_peakLength = -1;

	for (float *f : m_xicsAlignasIntensityIntegration) {
		std::memset(f, 0, m_integrationArraySizeMax * sizeof(float));
	}
	for (float *f : m_xicsAlignasIntensityIntegrationTight1) {
		std::memset(f, 0, m_integrationArraySizeMax * sizeof(float));
	}
	std::memset(m_productVecIntegration, 0, m_integrationArraySizeMax * sizeof(float));
}

Err CandidateScorertronV2::copyToPeakVecs(
	const PeakIntegrationIndexes &pii,
	const QVector<float*> &xicsAlignasIntensity,
	const QVector<float*> &xicsAlignasIntensityTight1,
	float* productVec
	) {
	ERR_INIT

	e = ErrorUtils::isTrue(pii.second > pii.first); ree;

	zeroOutArrays();

	m_peakLength = pii.second - pii.first + 1;
	e = ErrorUtils::isBelowThreshold(
		m_peakLength,
		m_integrationArraySizeMax,
		ErrorUtilsParam::ExcludeThreshold
		); ree;

	for (int i = 0; i < m_ms2IonsCount; i++) {
		std::memcpy(
			m_xicsAlignasIntensityIntegration[i],
			xicsAlignasIntensity[i] + pii.first,
			m_peakLength * sizeof(float)
			);

		std::memcpy(
			m_xicsAlignasIntensityIntegrationTight1[i],
			xicsAlignasIntensityTight1[i] + pii.first,
			m_peakLength * sizeof(float)
			);
	}

	std::memcpy(
		m_productVecIntegration,
		productVec + pii.first,
		m_peakLength * sizeof(float)
		);

	ERR_RETURN
}

Err CandidateScorertronV2::calculateScores(CandidateScoresV2 *candidateScoresV2) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_featuresVsFeaturesArrayIndex); ree;
	e = ErrorUtils::isGreaterThanZero(m_peakLength); ree;

	e = calculateCorrelationScoresMS2(candidateScoresV2); ree;
	e = calculateCorrelationScoresMS2Tight1(candidateScoresV2); ree;


	ERR_RETURN
}

Err CandidateScorertronV2::calculateCorrelationScoresMS2(CandidateScoresV2 *candScores) {

	ERR_INIT

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArr[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArr, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArrUpper[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArrUpper, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	float cosineSimSum = 0.0;
	float CosineSimSumTop16 = 0.0;
	float cosineSimSumGreaterThan80 = 0.0;

	if (m_ms2IonsCount > S_GLOBAL_SETTINGS.MIN_MS2_IONS) {
		AVXUtils::cosineSimilarityAVXParallel(
			m_productVecIntegration,
			m_xicsAlignasIntensityIntegration[8],
			m_xicsAlignasIntensityIntegration[9],
			m_xicsAlignasIntensityIntegration[10],
			m_xicsAlignasIntensityIntegration[11],
			m_xicsAlignasIntensityIntegration[12],
			m_xicsAlignasIntensityIntegration[13],
			m_xicsAlignasIntensityIntegration[14],
			m_xicsAlignasIntensityIntegration[15],
			m_peakLength,
			resultArrUpper
			);

		for (int i = 0; i < AVXUtils::AVX2_FLOAT_REGISTER_SIZE; i++) {
			const float cosineSimToAnchorI = resultArrUpper[i];
			candScores->featuresArray[
				CandidateScoresFeatureManager::CosineSimToAnchor1
				+ i
				+ AVXUtils::AVX2_FLOAT_REGISTER_SIZE
				] = cosineSimToAnchorI;
			CosineSimSumTop16 += cosineSimToAnchorI;
			if (cosineSimToAnchorI > 0.8f) {
				cosineSimSumGreaterThan80 += cosineSimToAnchorI;
			}
		}
	}

	AVXUtils::cosineSimilarityAVXParallel(
		m_productVecIntegration,
		m_xicsAlignasIntensityIntegration[0],
		m_xicsAlignasIntensityIntegration[1],
		m_xicsAlignasIntensityIntegration[2],
		m_xicsAlignasIntensityIntegration[3],
		m_xicsAlignasIntensityIntegration[4],
		m_xicsAlignasIntensityIntegration[5],
		m_xicsAlignasIntensityIntegration[6],
		m_xicsAlignasIntensityIntegration[7],
		m_peakLength,
		resultArr
		);

	for (int i = 0; i < AVXUtils::AVX2_FLOAT_REGISTER_SIZE; i++) {
		const float cosineSimToAnchorI = resultArr[i];
		candScores->featuresArray[CandidateScoresFeatureManager::CosineSimToAnchor1 + i] = cosineSimToAnchorI;
		cosineSimSum += cosineSimToAnchorI;
		if (cosineSimToAnchorI > 0.8f) {
			cosineSimSumGreaterThan80 += cosineSimToAnchorI;
		}
	}

	candScores->featuresArray[CandidateScoresFeatureManager::CosineSimSumTop8] = cosineSimSum;
	candScores->featuresArray[CandidateScoresFeatureManager::CosineSimSumGreaterThan80] = cosineSimSumGreaterThan80;

	if (m_ms2IonsCount > S_GLOBAL_SETTINGS.MIN_MS2_IONS) {
		candScores->featuresArray[CandidateScoresFeatureManager::CosineSimSumTop16] = cosineSimSum + CosineSimSumTop16;
	}

	ERR_RETURN
}

Err CandidateScorertronV2::calculateCorrelationScoresMS2Tight1(CandidateScoresV2 *candScores) {
	ERR_INIT

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArr[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArr, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArrUpper[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArrUpper, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	float cosineSimSum = 0.0;
	float CosineSimSumTop16 = 0.0;
	float cosineSimSumGreaterThan80 = 0.0;

	if (m_ms2IonsCount > S_GLOBAL_SETTINGS.MIN_MS2_IONS) {
		AVXUtils::cosineSimilarityAVXParallel(
			m_productVecIntegration,
			m_xicsAlignasIntensityIntegrationTight1[8],
			m_xicsAlignasIntensityIntegrationTight1[9],
			m_xicsAlignasIntensityIntegrationTight1[10],
			m_xicsAlignasIntensityIntegrationTight1[11],
			m_xicsAlignasIntensityIntegrationTight1[12],
			m_xicsAlignasIntensityIntegrationTight1[13],
			m_xicsAlignasIntensityIntegrationTight1[14],
			m_xicsAlignasIntensityIntegrationTight1[15],
			m_peakLength,
			resultArrUpper
			);

		for (int i = 0; i < AVXUtils::AVX2_FLOAT_REGISTER_SIZE; i++) {
			const float cosineSimToAnchorI = resultArrUpper[i];
			candScores->featuresArray[
				CandidateScoresFeatureManager::CosineSimToAnchorTight1_1
				+ i
				+ AVXUtils::AVX2_FLOAT_REGISTER_SIZE
				] = cosineSimToAnchorI;
			CosineSimSumTop16 += cosineSimToAnchorI;
			if (cosineSimToAnchorI > 0.8f) {
				cosineSimSumGreaterThan80 += cosineSimToAnchorI;
			}
		}
	}

	AVXUtils::cosineSimilarityAVXParallel(
		m_productVecIntegration,
		m_xicsAlignasIntensityIntegrationTight1[0],
		m_xicsAlignasIntensityIntegrationTight1[1],
		m_xicsAlignasIntensityIntegrationTight1[2],
		m_xicsAlignasIntensityIntegrationTight1[3],
		m_xicsAlignasIntensityIntegrationTight1[4],
		m_xicsAlignasIntensityIntegrationTight1[5],
		m_xicsAlignasIntensityIntegrationTight1[6],
		m_xicsAlignasIntensityIntegrationTight1[7],
		m_peakLength,
		resultArr
		);

	for (int i = 0; i < AVXUtils::AVX2_FLOAT_REGISTER_SIZE; i++) {
		const float cosineSimToAnchorI = resultArr[i];
		candScores->featuresArray[CandidateScoresFeatureManager::CosineSimToAnchorTight1_1 + i] = cosineSimToAnchorI;
		cosineSimSum += cosineSimToAnchorI;
		if (cosineSimToAnchorI > 0.8f) {
			cosineSimSumGreaterThan80 += cosineSimToAnchorI;
		}
	}

	candScores->featuresArray[CandidateScoresFeatureManager::CosineSimSumTop8Tight1] = cosineSimSum;
	candScores->featuresArray[CandidateScoresFeatureManager::CosineSimSumGreaterThan80Tight1] = cosineSimSumGreaterThan80;

	if (m_ms2IonsCount > S_GLOBAL_SETTINGS.MIN_MS2_IONS) {
		candScores->featuresArray[CandidateScoresFeatureManager::CosineSimSumTop16Tight1] = cosineSimSum + CosineSimSumTop16;
	}

	ERR_RETURN
}
