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
	const QVector<FTR> &featuresCalibration,
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

	for (const FTR &f : m_features) {
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
	const CandidateScorertronV2Input &input,
	CandidateScoresV2 *candidateScores
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(input.ms2IonsFull); ree;
	e = ErrorUtils::isNotEmpty(m_features); ree;
	e = ErrorUtils::isNotEmpty(m_featuresVsFeaturesArrayIndex); ree;

	QVector<MS2Ion> ms2IonsTrunc = input.ms2IonsFull;
	ms2IonsTrunc.resize(m_ms2IonsCount);

	e = copyToPeakVecs(
		input.pii,
		input.xicsAlignasIntensity,
		input.xicsAlignasIntensityTight1,
		input.productVec
		); ree;

	e = calculateScores(
		input,
		candidateScores
		); ree;

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

	// std::cout << m_peakLength << " fldsjf " << m_integrationArraySizeMax << std::endl;
	e = ErrorUtils::isTrue(m_peakLength < m_integrationArraySizeMax); ree;

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

Err CandidateScorertronV2::calculateScores(
	const CandidateScorertronV2Input &input,
	CandidateScoresV2 *candidateScoresV2
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_featuresVsFeaturesArrayIndex); ree;
	e = ErrorUtils::isGreaterThanZero(m_peakLength); ree;

	e = calculateRTCorrelationScoresMS2(candidateScoresV2); ree;
	e = calculateRTCorrelationScoresMS2Tight1(candidateScoresV2); ree;
	e = calculateFragmentCorrelationScoresMS2(
		input.pii,
		input.ms2IonsFull,
		candidateScoresV2
		); ree;

	e = calculateRTCorrelationScoresMS1(
		input,
		candidateScoresV2
		); ree;


	ERR_RETURN
}

Err CandidateScorertronV2::calculateRTCorrelationScoresMS2(CandidateScoresV2 *candScores) {

	ERR_INIT

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArr[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArr, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArrUpper[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArrUpper, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	float cosineSimSum = 0.0;
	float CosineSimSumTop16 = 0.0;
	float cosineSimSumGreaterThan80 = 0.0;
	QVector<float> cosineSims(8, 0.0f);

	if (m_ms2IonsCount > S_GLOBAL_SETTINGS.MIN_MS2_IONS) {
		AVXUtils::cosineSimilarityIntraAVXParallel(
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
				FTR::CosineSimToAnchor1
				+ i
				+ AVXUtils::AVX2_FLOAT_REGISTER_SIZE
				] = cosineSimToAnchorI;
			CosineSimSumTop16 += cosineSimToAnchorI;
			if (cosineSimToAnchorI > 0.8f) {
				cosineSimSumGreaterThan80 += cosineSimToAnchorI;
			}
		}
	}

	AVXUtils::cosineSimilarityIntraAVXParallel(
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

	std::memcpy(cosineSims.data(), resultArr, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	for (int i = 0; i < AVXUtils::AVX2_FLOAT_REGISTER_SIZE; i++) {
		const float cosineSimToAnchorI = resultArr[i];
		cosineSims.push_back(cosineSimToAnchorI);
		candScores->featuresArray[FTR::CosineSimToAnchor1 + i] = cosineSimToAnchorI;
		cosineSimSum += cosineSimToAnchorI;
		if (cosineSimToAnchorI > 0.8f) {
			cosineSimSumGreaterThan80 += cosineSimToAnchorI;
		}
	}

	candScores->featuresArray[FTR::CosineSimSumTop8] = cosineSimSum;
	candScores->featuresArray[FTR::CosineSimSumMeanCorrelation] = MathUtils::mean(cosineSims);
	candScores->featuresArray[FTR::CosineSimSumStDevCorrelation] = MathUtils::stDev(cosineSims);
	candScores->featuresArray[FTR::CosineSimSumGreaterThan80] = cosineSimSumGreaterThan80;

	if (m_ms2IonsCount > S_GLOBAL_SETTINGS.MIN_MS2_IONS) {
		candScores->featuresArray[FTR::CosineSimSumTop16] = cosineSimSum + CosineSimSumTop16;
	}

	ERR_RETURN
}

Err CandidateScorertronV2::calculateRTCorrelationScoresMS2Tight1(CandidateScoresV2 *candScores) {
	ERR_INIT

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArr[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArr, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArrUpper[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArrUpper, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	float cosineSimSum = 0.0;
	float CosineSimSumTop16 = 0.0;
	float cosineSimSumGreaterThan80 = 0.0;

	if (m_ms2IonsCount > S_GLOBAL_SETTINGS.MIN_MS2_IONS) {
		AVXUtils::cosineSimilarityIntraAVXParallel(
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
				FTR::CosineSimToAnchorTight1_1
				+ i
				+ AVXUtils::AVX2_FLOAT_REGISTER_SIZE
				] = cosineSimToAnchorI;
			CosineSimSumTop16 += cosineSimToAnchorI;
			if (cosineSimToAnchorI > 0.8f) {
				cosineSimSumGreaterThan80 += cosineSimToAnchorI;
			}
		}
	}

	AVXUtils::cosineSimilarityIntraAVXParallel(
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
		candScores->featuresArray[FTR::CosineSimToAnchorTight1_1 + i] = cosineSimToAnchorI;
		cosineSimSum += cosineSimToAnchorI;
		if (cosineSimToAnchorI > 0.8f) {
			cosineSimSumGreaterThan80 += cosineSimToAnchorI;
		}
	}

	candScores->featuresArray[FTR::CosineSimSumTop8Tight1] = cosineSimSum;
	candScores->featuresArray[FTR::CosineSimSumGreaterThan80Tight1] = cosineSimSumGreaterThan80;

	if (m_ms2IonsCount > S_GLOBAL_SETTINGS.MIN_MS2_IONS) {
		candScores->featuresArray[FTR::CosineSimSumTop16Tight1] = cosineSimSum + CosineSimSumTop16;
	}

	ERR_RETURN
}

Err CandidateScorertronV2::calculateFragmentCorrelationScoresMS2(
	const PeakIntegrationIndexes &pii,
	const QVector<MS2Ion> &ms2IonsFull,
	CandidateScoresV2 *candScores
	) {
	ERR_INIT

	QVector<float> theoIntensities(AVXUtils::AVX2_FLOAT_REGISTER_SIZE, 0.0f);
	for (int i = 0; i < std::min(m_ms2IonsCount, ms2IonsFull.size()); i++) {
		theoIntensities[i] = ms2IonsFull[i].intensity;
	}

	const __m256 theoFragmentIntensities = _mm256_set_ps(
		theoIntensities[7],
		theoIntensities[6],
		theoIntensities[5],
		theoIntensities[4],
		theoIntensities[3],
		theoIntensities[2],
		theoIntensities[1],
		theoIntensities[0]
		);

	QVector<float> cosineCumulative(m_peakLength, 0.0f);
	for (int i = 0; i < m_peakLength; i++) {

		const __m256 empericalIntensitiesFrameIndex = _mm256_set_ps(
			m_xicsAlignasIntensityIntegration[7][i],
			m_xicsAlignasIntensityIntegration[6][i],
			m_xicsAlignasIntensityIntegration[5][i],
			m_xicsAlignasIntensityIntegration[4][i],
			m_xicsAlignasIntensityIntegration[3][i],
			m_xicsAlignasIntensityIntegration[2][i],
			m_xicsAlignasIntensityIntegration[1][i],
			m_xicsAlignasIntensityIntegration[0][i]
			);

		const float cosineSimFrameIndex = AVXUtils::cosineSimilarityAVX(
			empericalIntensitiesFrameIndex,
			theoFragmentIntensities
			);
		cosineCumulative[i] = cosineSimFrameIndex;
	}

	const float cosineCumulativeMean = MathUtils::mean(cosineCumulative);
	const float cosineCumulativeStDev = MathUtils::stDev(cosineCumulative);


	candScores->featuresArray[FTR::CosineSimSpectrumOverTime] = cosineCumulativeMean;
	candScores->featuresArray[FTR::CosineSimSpectrumOverTimeCubed] = static_cast<float>(std::pow(cosineCumulativeMean, 3));
	candScores->featuresArray[FTR::CosineSimSpectrumOverTimeStDev] = cosineCumulativeStDev;


	ERR_RETURN
}

Err CandidateScorertronV2::calculateRTCorrelationScoresMS1(
	const CandidateScorertronV2Input &input,
	CandidateScoresV2 *candScores
	) {

	ERR_INIT

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float resultArr[AVXUtils::AVX2_FLOAT_REGISTER_SIZE];
	std::memset(resultArr, 0, AVXUtils::AVX2_FLOAT_REGISTER_SIZE * sizeof(float));

	float* ms1MonoIntegrate = m_xicsAlignasIntensityIntegration[0];
	std::memcpy(
		ms1MonoIntegrate,
		input.ms1MonoIsotopeVec + input.pii.first,
		m_peakLength * sizeof(float)
		);

	float* ms1C13Integrate = m_xicsAlignasIntensityIntegration[1];
	std::memcpy(
		ms1C13Integrate,
		input.ms1C13Vec + input.pii.first,
		m_peakLength * sizeof(float)
		);

	float* ms1C132Integrate = m_xicsAlignasIntensityIntegration[2];
	std::memcpy(
		ms1C132Integrate,
		input.ms1C132Vec + input.pii.first,
		m_peakLength * sizeof(float)
		);

	float* ms1PreMonoShadowIntegrate = m_xicsAlignasIntensityIntegration[3];
	std::memcpy(
		ms1PreMonoShadowIntegrate,
		input.ms1PreMonoShadowVec + input.pii.first,
		m_peakLength * sizeof(float)
		);

	AVXUtils::cosineSimilarityIntraAVXParallel(
		m_productVecIntegration,
		ms1MonoIntegrate,
		ms1C13Integrate,
		ms1C132Integrate,
		ms1PreMonoShadowIntegrate,
		m_xicsAlignasIntensityIntegration[4],
		m_xicsAlignasIntensityIntegration[5],
		m_xicsAlignasIntensityIntegration[6],
		m_xicsAlignasIntensityIntegration[7],
		m_peakLength,
		resultArr
		);

	candScores->featuresArray[FTR::CosineSimToAnchorMS1MonoIsotope] = resultArr[0];
	candScores->featuresArray[FTR::CosineSimToAnchorMS1C13] = resultArr[1];
	candScores->featuresArray[FTR::CosineSimToAnchorMS1C132] = resultArr[2];
	candScores->featuresArray[FTR::CosineSimToAnchorMS1PreMonoShadow] = resultArr[3];
	candScores->featuresArray[FTR::CosineSimSumDiffMonoVsPreMonoShadow] = resultArr[0] - resultArr[3];
	candScores->featuresArray[FTR::CosineSimSumDiffMonoVsPreMonoShadowAbs]
					= std::abs(candScores->featuresArray[FTR::CosineSimSumDiffMonoVsPreMonoShadow]);

	ERR_RETURN
}
