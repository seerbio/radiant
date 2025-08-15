//
// Created by andrewnichols on 8/14/25.
//

#ifndef CANDIDATESCORERTRONV2_H
#define CANDIDATESCORERTRONV2_H

#include "AlgorithmsFFLib_Exports.h"
#include "CandidateScoresFeatureManager.h"
#include "CandidateScoresV2.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsFrameV2.h"

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS CandidateScorertronV2 {

public:

	CandidateScorertronV2();
	~CandidateScorertronV2();

	Err init(
		const QVector<CandidateScoresFeatureManager::Features> &featuresCalibration,
		float ms2ExtractionWidthPPM,
		float minMs2IonsFoundCountThreshold,
		MsFrameV2* msFrameV2MS2,
		MsFrameV2* msFrameV2MS1
		);

	Err scoreCandidate(
		const QVector<MS2Ion> &ms2IonsTrunc,
		bool isDecoy,
		TargetDecoyCandidatePair *targetDecoyCandidatePair,
		CandidateScoresV2 *candidateScores
		);

private:

	Err initArrays();

	Err loadMS2IonArrays(
		const QVector<MS2Ion> &ms2Ions,
		float ms2ExtractionWidthPPM,
		bool subtractShadows,
		MsFrameV2* msFrameV2MS2
		);

	void zeroOutArrays();

	Err subtractShadowsArrays();

	Err smoothMS2IonArrays();

	Err buildLocationVectors();

	Err buildMs1Vec(
		float ms2ExtractionWidthPPM,
		bool isDecoy,
		TargetDecoyCandidatePair *tdcp
		) const;

	Err fitMS1XICToVecs(
		const XICPointsPntrs &xicPointsPntrs,
		float* vecIntensity,
		float *vecMz
		) const;

	[[nodiscard]] Err buildProductVec() const;

	[[nodiscard]] Err smoothMS1Arrays() const;

	Err scoreProductVecApexes();

	Err calculateScores(CandidateScoresV2 *candidateScoresV2);

	Err calculateCorrelationScores(
		const PeakIntegrationIndexes &pii,
		CandidateScoresV2 *candScores
		);

private:

	size_t m_xicSizeMaxAlignas;
	int m_ms2IonsCount;

	QVector<float*> m_xicsAlignasIntensity;
	QVector<float*> m_xicsAlignasIntensityShadows;
	QVector<float*> m_xicsAlignasMz;
	QVector<float*> m_xicsAlignasMzShadows;

	float* m_intensityVec;
	float* m_ionCountVec;
	float* m_productVec;

	float* m_mzMs1MonoIsotopeVecIntensity;
	float* m_mzMs1C13VecIntensity;
	float* m_mzMs1C132VecIntensity;
	float* m_mzMs1MonoIsotopeShadowVecIntensity;
	float* m_mzMs1MonoIsotopeVecMz;
	float* m_mzMs1C13VecMz;
	float* m_mzMs1C132VecMz;
	float* m_mzMs1MonoIsotopeShadowVecMz;

	int m_targetFrameIndexMax;
	size_t m_xicSizeTargetMaxAlignas;

	float m_intensityVecMax;

	QVector<QPair<int, float>> m_productVecApexes;

	QVector<float> m_smoothingKernel;

	MsFrameV2* m_msFrameV2MS2;
	MsFrameV2* m_msFrameV2MS1;

	float m_minMs2IonsFoundCountThreshold;
	TargetDecoyCandidatePair *m_targetDecoyCandidatePair;
	float m_ms2ExtractionWidthPPM;

	QVector<CandidateScoresFeatureManager::Features> m_features;
	QHash<CandidateScoresFeatureManager::Features, Index> m_featuresVsFeaturesArrayIndex;

	QVector<QPair<PeakIntegrationIndexes, float>> m_peakIntegrationIndexesVsIntensity;
};



#endif //CANDIDATESCORERTRONV2_H
