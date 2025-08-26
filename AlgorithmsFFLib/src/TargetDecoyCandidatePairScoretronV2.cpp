//
// Created by andrewnichols on 7/26/25.
//

#include "TargetDecoyCandidatePairScoretronV2.h"

#include "AVXUtils.h"
#include "CandidateScorertronV2.h"
#include "EigenUtils.h"
#include "EigenKernelUtils.h"
#include "MS2Ion.h"
#include "MsReaderMzMLLazyLoad.h"
#include "ObjectCSVWriters.h"
#include "PeakIntegratomatic.h"
#include "TargetDecoyCandidatePair.h"

TargetDecoyCandidatePairScoretronV2::TargetDecoyCandidatePairScoretronV2()
: m_mzTargetKeyCurrent("null")
, m_msFrameV2Current(nullptr)
, m_msFrameV2MS1(nullptr)
, m_xicSizeMaxAlignas(-1)
, m_ms2IonsCount(-1)
, m_intensityVec(nullptr)
, m_ionCountVec(nullptr)
, m_productVec(nullptr)
, m_mzMs1MonoIsotopeVecIntensity(nullptr)
, m_mzMs1C13VecIntensity(nullptr)
, m_mzMs1C132VecIntensity(nullptr)
, m_mzMs1MonoIsotopeShadowVecIntensity(nullptr)
, m_mzMs1MonoIsotopeVecMz(nullptr)
, m_mzMs1C13VecMz(nullptr)
, m_mzMs1C132VecMz(nullptr)
, m_mzMs1MonoIsotopeShadowVecMz(nullptr)
, m_targetFrameIndexMax(-1)
, m_xicSizeTargetMaxAlignas(-1)
, m_mzMs2Min(-1)
, m_mzMs2Max(-1)
, m_intensityVecMax(0)
, m_minMs2IonsFoundCount(-1)
, m_msCalibratomatic(nullptr)
, m_smoothingKernel({0.25, 0.5, 0.25})
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

	for (float* v : m_xicsAlignasIntensityTight1) {
		delete v;
	}

	for (float* v : m_xicsAlignasIntensityShadowsTight1) {
		delete v;
	}

	for (float* v : m_xicsAlignasMzTight1) {
		delete v;
	}

	for (float* v : m_xicsAlignasMzShadowsTight1) {
		delete v;
	}

	delete m_intensityVec;
	delete m_ionCountVec;
	delete m_productVec;

	delete m_mzMs1MonoIsotopeVecIntensity;
	delete m_mzMs1C13VecIntensity;
	delete m_mzMs1C132VecIntensity;
	delete m_mzMs1MonoIsotopeShadowVecIntensity;
	delete m_mzMs1MonoIsotopeVecMz;
	delete m_mzMs1C13VecMz;
	delete m_mzMs1C132VecMz;
	delete m_mzMs1MonoIsotopeShadowVecMz;
}

Err TargetDecoyCandidatePairScoretronV2::init(
	const QMap<MzTargetKey, MsFrameV2*> &mzTargetKeyVsMsFramesMS2Pntrs,
	const QVector<FTR> &featuresCalibration,
	const PythiaParameters &pythiaParameters,
	int ms2IonsCount,
	float minMs2IonsFoundCount,
	MsFrameV2 *msFrameV2MS1
	) {
	ERR_INIT

	e = ErrorUtils::isNotEmpty(mzTargetKeyVsMsFramesMS2Pntrs); ree;
	e = ErrorUtils::isNotEmpty(featuresCalibration); ree;
	e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
	e = ErrorUtils::isGreaterThanZero(ms2IonsCount); ree;
	e = ErrorUtils::isGreaterThanZero(minMs2IonsFoundCount); ree;
	// e = ErrorUtils::isTrue(msFrameV2MS1->isInit()); ree;

	m_mzTargetKeyVsMsFramesMS2Pntrs = mzTargetKeyVsMsFramesMS2Pntrs;
	m_pythiaParameters = pythiaParameters;
	m_features = featuresCalibration;

	constexpr int buffer = 10;
	m_xicSizeMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		mzTargetKeyVsMsFramesMS2Pntrs.first()->frameIndexSize() + buffer,
		AVXUtils::AVX2_ALIGNAS_SIZE
		);
	m_mzMs2Min = mzTargetKeyVsMsFramesMS2Pntrs.first()->mzMin();
	m_mzMs2Max = mzTargetKeyVsMsFramesMS2Pntrs.first()->mzMax();
	m_ms2IonsCount = ms2IonsCount;
	m_minMs2IonsFoundCount = minMs2IonsFoundCount;
	m_msFrameV2MS1 = msFrameV2MS1;

	m_xicsAlignasIntensity.resize(m_ms2IonsCount);
	m_xicsAlignasIntensityShadows.resize(m_ms2IonsCount);
	m_xicsAlignasMz.resize(m_ms2IonsCount);
	m_xicsAlignasMzShadows.resize(m_ms2IonsCount);
	m_xicsAlignasIntensityTight1.resize(m_ms2IonsCount);
	m_xicsAlignasIntensityShadowsTight1.resize(m_ms2IonsCount);
	m_xicsAlignasMzTight1.resize(m_ms2IonsCount);
	m_xicsAlignasMzShadowsTight1.resize(m_ms2IonsCount);
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

		auto* alignedIntensityValsTight1 = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
		m_xicsAlignasIntensityTight1[i] = alignedIntensityValsTight1;

		auto* alignedMzValsTight1 = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasMzTight1[i] = alignedMzValsTight1;

		auto* alignedIntensityValsShadowsTight1 = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasIntensityShadowsTight1[i] = alignedIntensityValsShadowsTight1;

		auto* alignedMzValsShadowsTight1 = static_cast<float*>(std::aligned_alloc(
			AVXUtils::AVX2_ALIGNAS_SIZE,
			m_xicSizeMaxAlignas * sizeof(float))
			);
		m_xicsAlignasMzShadowsTight1[i] = alignedMzValsShadowsTight1;
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

	auto* alignedProductVec = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_productVec = alignedProductVec;

	auto* mzMs1MonoIsotopeVecIntensity = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_mzMs1MonoIsotopeVecIntensity = mzMs1MonoIsotopeVecIntensity;

	auto* mzMs1C13VecIntensity = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_mzMs1C13VecIntensity = mzMs1C13VecIntensity;

	auto* mzMs1C132VecIntensity = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_mzMs1C132VecIntensity = mzMs1C132VecIntensity;

	auto* mzMs1MonoIsotopeShadowVecIntensity = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_mzMs1MonoIsotopeShadowVecIntensity = mzMs1MonoIsotopeShadowVecIntensity;

	auto* mzMs1MonoIsotopeVecMz = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_mzMs1MonoIsotopeVecMz = mzMs1MonoIsotopeVecMz;

	auto* mzMs1C13VecMz = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_mzMs1C13VecMz = mzMs1C13VecMz;

	auto* mzMs1C132VecMz = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_mzMs1C132VecMz = mzMs1C132VecMz;

	auto* mzMs1MonoIsotopeShadowVecMz = static_cast<float*>(std::aligned_alloc(
		AVXUtils::AVX2_ALIGNAS_SIZE,
		m_xicSizeMaxAlignas * sizeof(float))
		);
	m_mzMs1MonoIsotopeShadowVecMz = mzMs1MonoIsotopeShadowVecMz;

	e = m_candidateScoretronV2.init(
		featuresCalibration,
		m_ms2IonsCount
		); ree;

	// constexpr int order = 1;
	// constexpr int derivative = 0;
	// constexpr int rate = 1;
	// Eigen::MatrixX<float> kernel;
	// e = EigenKernelUtils::buildSavitzkyGolayKernel(
	// 	m_pythiaParameters.filterLengthMS2,
	// 	order,
	// 	derivative,
	// 	rate,
	// 	&kernel
	// 	); ree;
	// const Eigen::VectorX<float> kernelVec(kernel);
	// m_smoothingKernel = EigenUtils::convertEigenVectorToQVector(kernelVec);



	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::scoreTargetDecoyCandidatePairPntr(
	const QPair<MzTargetKey, TargetDecoyCandidatePair*> &mzTargetKeyVsTdcpPntr,
	QPair<CandidateScoresV2Target, CandidateScoresV2Decoy> *scoresTargetVsDecoyPair
	) {

	ERR_INIT

	const MzTargetKey &mzTargetKey = mzTargetKeyVsTdcpPntr.first;
	if (mzTargetKey != m_mzTargetKeyCurrent) {
		e = ErrorUtils::contains(mzTargetKey, m_mzTargetKeyVsMsFramesMS2Pntrs); ree;
		m_msFrameV2Current = m_mzTargetKeyVsMsFramesMS2Pntrs.value(mzTargetKey);
		e = ErrorUtils::isTrue(m_msFrameV2Current->isInit()); ree;
		m_mzTargetKeyCurrent = mzTargetKey;
		m_mzMs2Min = m_msFrameV2Current->mzMin();
		m_mzMs2Max = m_msFrameV2Current->mzMax();
	}

	TargetDecoyCandidatePair *tdcpPntr = mzTargetKeyVsTdcpPntr.second;

	constexpr float mzMs2MinMin = 200.0f;
	const QVector<MS2Ion> ms2IonsTargetFullLength = tdcpPntr->ms2IonsTarget(
		std::max(m_mzMs2Min, mzMs2MinMin),
		m_mzMs2Max
		);
	const QVector<MS2Ion> &ms2IonsDecoyFullLength = tdcpPntr->ms2IonsDecoy(ms2IonsTargetFullLength);

	CandidateScoresV2 candidateScoresV2Target;
	e = scoreMS2Ions(
		ms2IonsTargetFullLength,
		false,
		tdcpPntr,
		&candidateScoresV2Target
		); ree;

// #define CHECK_VEC
#ifdef CHECK_VEC
	QVector<float> vectorFromPointer(m_xicSizeMaxAlignas, 0);

	std::copy_n(m_ionCountVec, m_xicSizeMaxAlignas, vectorFromPointer.data());
	const float max = *std::max_element(vectorFromPointer.begin(), vectorFromPointer.end());

	if (max > 7.9) {

		qDebug() << max << "Maximulmmml";

		e = ObjectCSVWriters::writeVectorToFile(vectorFromPointer, "testes.csv"); ree;
		e = ObjectCSVWriters::writeRawPointerToFile(m_intensityVec, m_xicSizeMaxAlignas, "testesIntz.csv"); ree;
		e = ObjectCSVWriters::writeRawPointerToFile(m_ionCountVec, m_xicSizeMaxAlignas, "testesCnt.csv"); ree;
		e = ObjectCSVWriters::writeRawPointerToFile(m_productVec, m_xicSizeMaxAlignas, "testesProd.csv"); ree;

		for (int i = 0; i < m_xicsAlignasIntensity.size(); i++) {
			e = ObjectCSVWriters::writeRawPointerToFile(
				m_xicsAlignasIntensity[i],
				m_xicSizeMaxAlignas,
				"intz_" + QString::number(i) + ".csv"
				); ree;
		}

		// for (const QPair<int, float> &ap : m_productVecApexes) {
		//
		// 	if (m_ionCountVec[ap.first] < 3) {
		// 		continue;
		// 	}
		//
		// 	std::cout << ap.first << "," << std::endl;
		// }

		throw std::runtime_error("Error in TargetDecoyCandidatePairScoretronV2");
	}
#endif

	CandidateScoresV2 candidateScoresV2Decoy;
	e = scoreMS2Ions(
		ms2IonsDecoyFullLength,
		true,
		tdcpPntr,
		&candidateScoresV2Decoy
		); ree;

	*scoresTargetVsDecoyPair = {candidateScoresV2Target, candidateScoresV2Decoy};

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::scoreMS2Ions(
	const QVector<MS2Ion> &ms2IonsFull,
	bool isDecoy,
	TargetDecoyCandidatePair *tdcp,
	CandidateScoresV2 *candidateScoresV2Best
	) {

	ERR_INIT

	QVector<MS2Ion> ms2IonsTrunc = ms2IonsFull;
	ms2IonsTrunc.resize(m_ms2IonsCount);

	e = loadMS2IonArrays(ms2IonsTrunc); ree;

	if (m_targetFrameIndexMax < 0) {
		//TODO figure out what to return here. Most likely blank candidate score
		ERR_RETURN
	}

	// e = smoothMS2IonArrays(); ree;
	if (m_pythiaParameters.subtractShadows) {
		e = subtractShadowsArrays(); ree;
	}

	e = buildLocationVectors(); ree;
	if (m_msFrameV2MS1->isInit()) {
		e = buildMs1Vec(isDecoy, tdcp); ree;
	}
	e = buildProductVec(); ree;
	e = smoothMS1Arrays(); ree;
	e = scoreProductVecApexes(
		ms2IonsFull,
		isDecoy,
		tdcp,
		candidateScoresV2Best
		);ree;

	{
	// if (m_pythiaParameters.subtractShadows) {
	// 	constexpr bool zeroNegatives = true;
	// 	e = AVXUtils::subtractArraysAVX2(
	// 		m_mzMs1MonoIsotopeVecIntensity,
	// 		m_mzMs1MonoIsotopeShadowVecIntensity,
	// 		m_xicSizeTargetMaxAlignas,
	// 		zeroNegatives
	// 		); ree;
	// }
	// e = scoreProductVecApexes(); ree;
	}

	// e = connectApexes(); ree;

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::loadMS2IonArrays(const QVector<MS2Ion> &ms2Ions) {

	ERR_INIT

	zeroOutArrays();

	int arrayIndex = 0;
	for (const MS2Ion &ms2Ion : ms2Ions) {

		const float mzTol = MathUtils::calculatePPM(
			ms2Ion.mz,
			static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
			);

		const float mzMin = ms2Ion.mz - mzTol;
		const float mzMax = ms2Ion.mz + mzTol;
		const float mzTolTight1 = mzTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION;
		const float mzMinTight1 = ms2Ion.mz - mzTolTight1;
		const float mzMaxTight1 = ms2Ion.mz + mzTolTight1;

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
		const float mzShadowTolTight1 = mzShadowTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION;
		const float mzShadowMinTight1 = mzShadow - mzShadowTolTight1;
		const float mzShadowMaxTight1 = mzShadow + mzShadowTolTight1;

		XICPointsPntrs xicPointsPntrsShadows;
		e = m_msFrameV2Current->getTurboXICPntr()->extractPointsXIC(
			mzShadowMin,
			mzShadowMax,
			&xicPointsPntrsShadows
			); ree;

		float* arrIntensity = m_xicsAlignasIntensity[arrayIndex];
		float* arrMz = m_xicsAlignasMz[arrayIndex];
		float* arrIntensityTight1 = m_xicsAlignasIntensityTight1[arrayIndex];
		float* arrMzTight1 = m_xicsAlignasMzTight1[arrayIndex];
		for (const XICPoint *xicPoint : xicPointsPntrs) {

			m_targetFrameIndexMax = std::max(m_targetFrameIndexMax, xicPoint->frameIndex);

			const float mz = xicPoint->mz;
			const FrameIndex frameIndex = xicPoint->frameIndex;

			arrIntensity[frameIndex] += xicPoint->intensity;
			if (arrMz[frameIndex] > 1) {
				arrMz[frameIndex] += mz;
				arrMz[frameIndex] *= 0.5;
				continue;
			}
			arrMz[xicPoint->frameIndex] = mz;

			if (mz < mzMinTight1 || mz > mzMaxTight1) {
				continue;
			}

			arrIntensityTight1[frameIndex] += xicPoint->intensity;
			if (arrMzTight1[frameIndex] > 1) {
				arrMzTight1[frameIndex] += mz;
				arrMzTight1[frameIndex] *= 0.5;
				continue;
			}
			arrMzTight1[xicPoint->frameIndex] = mz;
		}

		float* arrIntensityShadows = m_xicsAlignasIntensityShadows[arrayIndex];
		float* arrMzShadows = m_xicsAlignasMzShadows[arrayIndex];
		float* arrIntensityShadowsTight1 = m_xicsAlignasIntensityShadowsTight1[arrayIndex];
		float* arrMzShadowsTight1 = m_xicsAlignasMzShadowsTight1[arrayIndex];

		for (const XICPoint *xicPoint : xicPointsPntrsShadows) {

			const float mz = xicPoint->mz;
			const FrameIndex frameIndex = xicPoint->frameIndex;

			arrIntensityShadows[frameIndex] += xicPoint->intensity;
			if (arrMzShadows[frameIndex] > 1) {
				arrMzShadows[frameIndex] += mz;
				arrMzShadows[frameIndex] *= 0.5;
				continue;
			}
			arrMzShadows[xicPoint->frameIndex] = mz;

			if (mz < mzShadowMinTight1 || mz > mzShadowMaxTight1) {
				continue;
			}

			arrIntensityShadowsTight1[frameIndex] += xicPoint->intensity;
			if (arrMzShadowsTight1[frameIndex] > 1) {
				arrMzShadowsTight1[frameIndex] += mz;
				arrMzShadowsTight1[frameIndex] *= 0.5;
				continue;
			}
			arrMzShadowsTight1[xicPoint->frameIndex] = mz;
		}

		arrayIndex++;
	}

	m_xicSizeTargetMaxAlignas = AVXUtils::calculateNextAlignedBlockSize(
		m_targetFrameIndexMax,
		AVXUtils::AVX2_ALIGNAS_SIZE
		);

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::subtractShadowsArrays() {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(m_xicSizeTargetMaxAlignas); ree;

	for (int i = 0; i < m_ms2IonsCount; i++) {
		constexpr bool zeroNegatives = true;
		e = AVXUtils::subtractArraysAVX2(
			m_xicsAlignasIntensity[i],
			m_xicsAlignasIntensityShadows[i],
			m_xicSizeTargetMaxAlignas,
			zeroNegatives
			); ree;

		e = AVXUtils::subtractArraysAVX2(
			m_xicsAlignasIntensityTight1[i],
			m_xicsAlignasIntensityShadowsTight1[i],
			m_xicSizeTargetMaxAlignas,
			zeroNegatives
			); ree;
	}

	ERR_RETURN
}

void TargetDecoyCandidatePairScoretronV2::zeroOutArrays() {

	m_targetFrameIndexMax = -1;
	m_xicSizeTargetMaxAlignas = -1;
	m_intensityVecMax = 0;

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

	for (float *f : m_xicsAlignasIntensityTight1) {
		std::memset(f, 0, m_xicSizeMaxAlignas * sizeof(float));
	}

	for (float *f : m_xicsAlignasIntensityShadowsTight1) {
		std::memset(f, 0, m_xicSizeMaxAlignas * sizeof(float));
	}

	for (float *f : m_xicsAlignasMzTight1) {
		std::memset(f, 0, m_xicSizeMaxAlignas * sizeof(float));
	}

	for (float *f : m_xicsAlignasMzShadowsTight1) {
		std::memset(f, 0, m_xicSizeMaxAlignas * sizeof(float));
	}

	std::memset(m_intensityVec, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_ionCountVec, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_productVec, 0, m_xicSizeMaxAlignas * sizeof(float));

	std::memset(m_mzMs1MonoIsotopeVecIntensity, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_mzMs1C13VecIntensity, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_mzMs1C132VecIntensity, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_mzMs1MonoIsotopeShadowVecIntensity, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_mzMs1MonoIsotopeVecMz, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_mzMs1C13VecMz, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_mzMs1C132VecMz, 0, m_xicSizeMaxAlignas * sizeof(float));
	std::memset(m_mzMs1MonoIsotopeShadowVecMz, 0, m_xicSizeMaxAlignas * sizeof(float));
}

Err TargetDecoyCandidatePairScoretronV2::smoothMS2IonArrays() {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(m_xicSizeTargetMaxAlignas); ree;
	e = ErrorUtils::isNotEmpty(m_smoothingKernel); ree;

	e = AVXUtils::convolveEightVecsWithKernelAVXFloat(
		m_smoothingKernel,
		m_xicSizeTargetMaxAlignas,
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
		e = AVXUtils::convolveEightVecsWithKernelAVXFloat(
				m_smoothingKernel,
				m_xicSizeTargetMaxAlignas,
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

	e = AVXUtils::convolveEightVecsWithKernelAVXFloat(
		m_smoothingKernel,
		m_xicSizeTargetMaxAlignas,
		m_xicsAlignasIntensityTight1[0],
		m_xicsAlignasIntensityTight1[1],
		m_xicsAlignasIntensityTight1[2],
		m_xicsAlignasIntensityTight1[3],
		m_xicsAlignasIntensityTight1[4],
		m_xicsAlignasIntensityTight1[5],
		m_xicsAlignasIntensityTight1[6],
		m_xicsAlignasIntensityTight1[7]
		); ree;

	if (m_ms2IonsCount == S_GLOBAL_SETTINGS.MAX_MS2_IONS) {
		e = AVXUtils::convolveEightVecsWithKernelAVXFloat(
				m_smoothingKernel,
				m_xicSizeTargetMaxAlignas,
				m_xicsAlignasIntensityTight1[8],
				m_xicsAlignasIntensityTight1[9],
				m_xicsAlignasIntensityTight1[10],
				m_xicsAlignasIntensityTight1[11],
				m_xicsAlignasIntensityTight1[12],
				m_xicsAlignasIntensityTight1[13],
				m_xicsAlignasIntensityTight1[14],
				m_xicsAlignasIntensityTight1[15]
				); ree;
	}

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::buildLocationVectors() {

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

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::buildProductVec() const {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(m_targetFrameIndexMax); ree;
	e = ErrorUtils::isGreaterThanZero(m_intensityVecMax); ree;

	for (int i = 0; i < m_xicSizeTargetMaxAlignas; i += AVXUtils::AVX2_FLOAT_REGISTER_SIZE) {

		__m256 ionCount = _mm256_load_ps(m_ionCountVec + i);
		const __m256 thresholds = _mm256_set1_ps(m_minMs2IonsFoundCount);
		const __m256 mask = _mm256_cmp_ps(ionCount, thresholds, _CMP_LT_OQ);

		if (AVXUtils::isAllOnes(mask)) {
			continue;
		}
		ionCount = _mm256_blendv_ps(ionCount, _mm256_setzero_ps(), mask);

		const __m256 intensity = _mm256_load_ps(m_intensityVec + i);
		const __m256 product = _mm256_mul_ps(ionCount, intensity);

		_mm256_store_ps(m_productVec + i, product);
	}

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::buildMs1Vec(bool isDecoy, TargetDecoyCandidatePair *tdcp) const {

	ERR_INIT

	e = ErrorUtils::isTrue(m_msFrameV2MS1->isInit()); ree;
	e = ErrorUtils::isGreaterThanZero(m_xicSizeTargetMaxAlignas); ree;

	const auto isotopeDistanceThomsons = static_cast<float>(S_GLOBAL_SETTINGS.ISO_DIFF / tdcp->charge());

	const float mzMs1MonoIsotope = tdcp->mz(isDecoy);
	const float mzMs1C13 = mzMs1MonoIsotope + isotopeDistanceThomsons;
	const float mzMs1C132 = mzMs1C13 + isotopeDistanceThomsons;
	const float mzMs1MonoIsotopeShadow = mzMs1MonoIsotope - isotopeDistanceThomsons;
	const auto ms2ExtractionWidthPPM = static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM);

	const float massTolMonoIsotope = MathUtils::calculatePPM(mzMs1MonoIsotope, ms2ExtractionWidthPPM);
	const float massTolMs1C13 = MathUtils::calculatePPM(mzMs1C13, ms2ExtractionWidthPPM);
	const float massTolMs1C132 = MathUtils::calculatePPM(mzMs1C132, ms2ExtractionWidthPPM);
	const float massTolMs1MonoIsotopeShadow = MathUtils::calculatePPM(mzMs1C13 , ms2ExtractionWidthPPM);

	TurboXIC *turboXICMS1 = m_msFrameV2MS1->getTurboXICPntr();

	XICPointsPntrs xicPointsPntrsMono;
	XICPointsPntrs xicPointsPntrsC13;
	XICPointsPntrs xicPointsPntrsC132;
	XICPointsPntrs xicPointsPntrsMonoShadow;

	e = turboXICMS1->extractPointsXIC(
		mzMs1MonoIsotope - massTolMonoIsotope,
		mzMs1MonoIsotope + massTolMonoIsotope,
		&xicPointsPntrsMono
		); ree;
	e = fitMS1XICToVecs(
		xicPointsPntrsMono,
		m_mzMs1MonoIsotopeVecIntensity,
		m_mzMs1MonoIsotopeVecMz
		); ree;

	e = turboXICMS1->extractPointsXIC(
		mzMs1C13 - massTolMs1C13,
		mzMs1C13 + massTolMs1C13,
		&xicPointsPntrsC13
		); ree;
	e = fitMS1XICToVecs(
		xicPointsPntrsC13,
		m_mzMs1C13VecIntensity,
		m_mzMs1C13VecMz
		); ree;

	e = turboXICMS1->extractPointsXIC(
		mzMs1C132 - massTolMs1C132,
		mzMs1C132 + massTolMs1C132,
		&xicPointsPntrsC132
		); ree;
	e = fitMS1XICToVecs(
		xicPointsPntrsC132,
		m_mzMs1C132VecIntensity,
		m_mzMs1C132VecMz
		); ree;

	e = turboXICMS1->extractPointsXIC(
		mzMs1MonoIsotopeShadow - massTolMs1MonoIsotopeShadow,
		mzMs1MonoIsotopeShadow + massTolMs1MonoIsotopeShadow,
		&xicPointsPntrsMonoShadow
		); ree;
	e = fitMS1XICToVecs(
		xicPointsPntrsMonoShadow,
		m_mzMs1MonoIsotopeShadowVecIntensity,
		m_mzMs1MonoIsotopeShadowVecMz
		); ree;

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::fitMS1XICToVecs(
	const XICPointsPntrs &xicPointsPntrs,
	float* vecIntensity,
	float *vecMz
	) const {

	ERR_INIT

	e = ErrorUtils::isTrue(m_msFrameV2MS1->isInit()); ree;
	e = ErrorUtils::isTrue(m_msFrameV2Current->isInit()); ree;

	for (XICPoint *xicPoint : xicPointsPntrs) {
		const float ms1ScanTime = m_msFrameV2MS1->getScanTimeFromFrameIndex(xicPoint->frameIndex);
		const FrameIndex frameIndexMS2 = m_msFrameV2Current->getFrameIndex(ms1ScanTime);

		if (frameIndexMS2 > m_targetFrameIndexMax  || frameIndexMS2 < 0) {
			continue;
		}

		vecIntensity[frameIndexMS2] += xicPoint->intensity;
		vecMz[frameIndexMS2] = xicPoint->mz;
	}

	ERR_RETURN
}

Err TargetDecoyCandidatePairScoretronV2::smoothMS1Arrays() const {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(m_xicSizeTargetMaxAlignas); ree;
	e = ErrorUtils::isNotEmpty(m_smoothingKernel); ree;

	alignas(AVXUtils::AVX2_ALIGNAS_SIZE) float dummy[m_xicSizeTargetMaxAlignas];
	e = AVXUtils::convolveEightVecsWithKernelAVXFloat(
		m_smoothingKernel,
		m_xicSizeTargetMaxAlignas,
		m_ionCountVec,
		m_intensityVec,
		m_productVec,
		m_mzMs1MonoIsotopeVecIntensity,
		m_mzMs1C13VecIntensity,
		m_mzMs1C132VecIntensity,
		m_mzMs1MonoIsotopeShadowVecIntensity,
		dummy
		); ree;

	ERR_RETURN
}

namespace {
	void filterApexesByIonCount(
		const float* ionCountVec,
		float thresholdValue,
		QVector<QPair<int, float>> *productVecApexes
		) {

		const auto terminatorLogic = [thresholdValue, ionCountVec](const QPair<int, float> &pr) {
			return ionCountVec[pr.first] < thresholdValue;
		};

		const auto terminator = std::remove_if(
			productVecApexes->begin(),
			productVecApexes->end(),
			terminatorLogic
			);

		productVecApexes->erase(terminator, productVecApexes->end());
	}

	QPair<float, float> calcMedianAndStdDev (
		float *ionCountVec,
		int size
		) {

		QVector<float> ionCountsNoZeros;
		ionCountsNoZeros.reserve(size);

		for (int i = 0; i < size; i++) {
			const float cnt = ionCountVec[i];
			if (MathUtils::tZero(cnt)) {
				continue;
			}
			ionCountsNoZeros.push_back(cnt);
		}

		return {
			static_cast<float>(MathUtils::median(ionCountsNoZeros)),
			static_cast<float>(MathUtils::stDev(ionCountsNoZeros))
		};
	}

}//namespace
Err TargetDecoyCandidatePairScoretronV2::scoreProductVecApexes(
	const QVector<MS2Ion> &ms2IonsFull,
	bool isDecoy,
	TargetDecoyCandidatePair *tdcp,
	CandidateScoresV2 *candidateScoresBest
	) {

	ERR_INIT

	e = ErrorUtils::isGreaterThanZero(m_xicSizeMaxAlignas); ree;
	e = ErrorUtils::isTrue(m_candidateScoretronV2.isInit()); ree;

	m_productVecApexes = AVXUtils::findApexesAVX2(
		m_productVec,
		m_xicSizeTargetMaxAlignas
		);

	QPair<float, float> medianVsStDevIonCount = calcMedianAndStdDev(
		m_ionCountVec,
		m_xicSizeTargetMaxAlignas
		);

	const float ionCountThreshold = medianVsStDevIonCount.first * (2 * medianVsStDevIonCount.second);
	filterApexesByIonCount(
		m_ionCountVec,
		ionCountThreshold,
		&m_productVecApexes
		);

	PeakIntegratomaticParameters params;
	params.stopThresholdFraction = 0.1;

	PeakIntegratomatic peakIntegratomatic;
	e = peakIntegratomatic.init(params); ree;

	QVector<int> apexes;
	std::transform(
		m_productVecApexes.begin(),
		m_productVecApexes.end(),
		std::back_inserter(apexes),
		[](const QPair<int, float> &pr) { return pr.first;}
		);

	QVector<std::tuple<PeakIntegrationIndexes, Intensity, FrameIndex>> peakIntegrationIndexesVsIntensityVsApex;
	e = peakIntegratomatic.simpleIntegrator(
		apexes,
		m_productVec,
		m_xicSizeTargetMaxAlignas,
		&peakIntegrationIndexesVsIntensityVsApex
		); ree;

	// using T = std::tuple<PeakIntegrationIndexes, Intensity, FrameIndex>;
	// std::sort(
	// 	peakIntegrationIndexesVsIntensityVsApex.rbegin(),
	// 	peakIntegrationIndexesVsIntensityVsApex.rend(),
	// 	[](const T &l, const T &r){return std::get<1>(l) < std::get<1>(r);}
	// 	);
	// constexpr int top10 = 10;
	// peakIntegrationIndexesVsIntensityVsApex.resize(std::min(top10, peakIntegrationIndexesVsIntensityVsApex.size()));
	// std::sort(
	// 	peakIntegrationIndexesVsIntensityVsApex.begin(),
	// 	peakIntegrationIndexesVsIntensityVsApex.end(),
	// 	[](const T &l, const T &r){return std::get<2>(l) < std::get<2>(r);}
	// 	);

	CandidateScorertronV2Input input;
	input.ms2IonsFull = ms2IonsFull;
	input.xicsAlignasIntensity = m_xicsAlignasIntensity;
	input.xicsAlignasIntensityTight1 = m_xicsAlignasIntensityTight1;
	input.productVec = m_productVec;
	input.ms1MonoIsotopeVec = m_mzMs1MonoIsotopeVecIntensity;
	input.ms1C13Vec = m_mzMs1C13VecIntensity;
	input.ms1C132Vec = m_mzMs1C132VecIntensity;
	input.ms1PreMonoShadowVec = m_mzMs1MonoIsotopeShadowVecIntensity;

	for (const std::tuple<PeakIntegrationIndexes, Intensity, FrameIndex> &pii : peakIntegrationIndexesVsIntensityVsApex) {

		input.pii = std::get<0>(pii);

		if(std::get<0>(pii).second - std::get<0>(pii).first + 1 > 100) {
			qDebug() << "ACHTUNG!!!";
			const QString filename = "prod_vec_" + QString::number(std::get<0>(pii).first) + "_" + QString::number(std::get<0>(pii).second) + ".csv";
			e = ObjectCSVWriters::writeRawPointerToFile(m_productVec, m_xicSizeMaxAlignas, filename); ree;
			continue;
		}

		CandidateScoresV2 candidateScores;
		candidateScores.isDecoy = isDecoy;
		candidateScores.targetDecoyCandidatePair = tdcp;
		candidateScores.initFeaturesArray();
		candidateScores.frameIndex = std::get<2>(pii);
		candidateScores.scanNumber = m_msFrameV2Current->getScanNumber(candidateScores.frameIndex);
		candidateScores.scanTime = m_msFrameV2Current->getScanTimeFromFrameIndex(candidateScores.frameIndex);
		candidateScores.ionCountVectorMean = medianVsStDevIonCount.first / static_cast<float>(m_ms2IonsCount);
		candidateScores.ionCountVectorStDev = medianVsStDevIonCount.second  / static_cast<float>(m_ms2IonsCount);

		e = m_candidateScoretronV2.scoreCandidate(
			input,
			&candidateScores
			); ree;

		if (
			candidateScores.featuresArray[FTR::CosineSimSumTop8]
			> candidateScoresBest->featuresArray[FTR::CosineSimSumTop8]
			) {
			*candidateScoresBest = candidateScores;
		}

		// if (candidateScores.featuresArray[FTR::CosineSimSumTop8] < 7.5) {
		// 	continue;
		// }
		//
		// qDebug()
		// << candidateScores.isDecoy
		// << candidateScores.featuresArray[FTR::CosineSimSumTop8]
		// << candidateScores.featuresArray[FTR::CosineSimSumTop16]
		// << candidateScores.featuresArray[FTR::CosineSimSumGreaterThan80]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor1]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor2]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor3]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor4]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor5]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor6]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor7]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor8]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor9]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor10]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor11]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor12]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor13]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor14]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchor15]
		// << candidateScores.featuresArray[FTR::CosineSimSumTop8Tight1]
		// << candidateScores.featuresArray[FTR::CosineSimSumTop16Tight1]
		// << candidateScores.featuresArray[FTR::CosineSimSumGreaterThan80Tight1]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_1]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_2]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_3]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_4]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_5]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_6]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_7]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_8]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_9]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_10]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_11]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_12]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_13]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_14]
		// << candidateScores.featuresArray[FTR::CosineSimToAnchorTight1_15]
		// ;

	}

	ERR_RETURN
}
