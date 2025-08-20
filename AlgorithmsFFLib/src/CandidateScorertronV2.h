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

struct CandidateScorertronV2Input {
	PeakIntegrationIndexes pii;
	QVector<MS2Ion> ms2IonsFull;
	QVector<float*> xicsAlignasIntensity;
	QVector<float*> xicsAlignasIntensityTight1;
	float* productVec = nullptr;
	float* ms1MonoIsotopeVec = nullptr;
	float* ms1C13Vec = nullptr;
	float* ms1C132Vec = nullptr;
	float* ms1PreMonoShadowVec = nullptr;

};

class ALGORITHMSFFLIB_EXPORTS CandidateScorertronV2 {

public:

	CandidateScorertronV2();
	~CandidateScorertronV2();

	Err init(
		const QVector<FTR> &featuresCalibration,
		int ms2IonCount
		);

	bool isInit() const;

	Err scoreCandidate(
		const CandidateScorertronV2Input &input,
		CandidateScoresV2 *candidateScores
		);

	void zeroOutArrays();

private:

	Err copyToPeakVecs(
		const PeakIntegrationIndexes &pii,
		const QVector<float*> &xicsAlignasIntensity,
		const QVector<float*> &xicsAlignasIntensityTight1,
		float* productVec
		);

	Err calculateScores(
		const CandidateScorertronV2Input &input,
		CandidateScoresV2 *candidateScoresV2
		);

	Err calculateRTCorrelationScoresMS2(CandidateScoresV2 *candScores);

	Err calculateRTCorrelationScoresMS2Tight1(CandidateScoresV2 *candScores);

	Err calculateFragmentCorrelationScoresMS2(
		const PeakIntegrationIndexes &pii,
		const QVector<MS2Ion> &ms2IonsFull,
		CandidateScoresV2 *candScores
		);

	Err calculateRTCorrelationScoresMS1(
		const CandidateScorertronV2Input &input,
		CandidateScoresV2 *candScores
		);

private:

	int m_ms2IonsCount;
	int m_integrationArraySizeMax;

	int m_peakLength;
	QVector<float*> m_xicsAlignasIntensityIntegration;
	QVector<float*> m_xicsAlignasIntensityIntegrationTight1;
	float* m_productVecIntegration;

	QVector<FTR> m_features;
	QHash<FTR, Index> m_featuresVsFeaturesArrayIndex;
};



#endif //CANDIDATESCORERTRONV2_H
