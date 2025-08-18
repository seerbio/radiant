//
// Created by andrewnichols on 7/26/25.
//

#ifndef TARGETDECOYCANDIDATEPAIRSCORETRONV2_H
#define TARGETDECOYCANDIDATEPAIRSCORETRONV2_H

#include "AlgorithmsFFLib_Exports.h"

#include "CandidateScoresFeatureManager.h"
#include "CandidateScorertronV2.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "MsCalibratomatic.h"
#include "MsFrameV2.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"

using namespace Error;

class MS2Ion;
class TargetDecoyCandidatePair;

class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePairScoretronV2 {

	friend class TargetDecoyCandidatePairScoretronV2Tests;

public:

	TargetDecoyCandidatePairScoretronV2();
	~TargetDecoyCandidatePairScoretronV2();

	Err init(
		const QMap<MzTargetKey, MsFrameV2*> &mzTargetKeyVsMsFramesMS2Pntrs,
		const QVector<CandidateScoresFeatureManager::Features> &featuresCalibration,
		const PythiaParameters &pythiaParameters,
		int ms2IonsCount,
		float minMs2IonsFoundCount,
		MsFrameV2 *msFrameV2MS1
		);

	Err scoreTargetDecoyCandidatePairPntr(
		const QPair<MzTargetKey, TargetDecoyCandidatePair*> &mzTargetKeyVsTdcpPntr
		);

private:

	Err scoreMS2Ions(
		const QVector<MS2Ion> &ms2IonsFull,
		bool isDecoy,
		TargetDecoyCandidatePair *tdcp
		);

	Err loadMS2IonArrays(const QVector<MS2Ion> &ms2Ions);

	Err subtractShadowsArrays();

	void zeroOutArrays();

	Err smoothMS2IonArrays();

	Err buildLocationVectors();

	[[nodiscard]] Err buildProductVec() const;

	[[nodiscard]] Err buildMs1Vec(bool isDecoy, TargetDecoyCandidatePair *tdcp) const;

	Err fitMS1XICToVecs(
		const XICPointsPntrs &xicPointsPntrs,
		float* vecIntensity,
		float *vecMz
		) const;

	[[nodiscard]] Err smoothMS1Arrays() const;

	Err scoreProductVecApexes(
		const QVector<MS2Ion> &ms2IonsFull,
		bool isDecoy,
		TargetDecoyCandidatePair *tdcp
		);


private:

	QMap<MzTargetKey, MsFrameV2*> m_mzTargetKeyVsMsFramesMS2Pntrs;
	PythiaParameters m_pythiaParameters;

	MzTargetKey m_mzTargetKeyCurrent;
	MsFrameV2 *m_msFrameV2Current;
	MsFrameV2 *m_msFrameV2MS1;

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

	QVector<float> m_smoothingKernel;

	float m_mzMs2Min;
	float m_mzMs2Max;
	float m_intensityVecMax;

	float m_minMs2IonsFoundCount;

	QVector<QPair<int, float>> m_productVecApexes;

	MsCalibratomatic *m_msCalibratomatic;

	QVector<CandidateScoresFeatureManager::Features> m_features;

	CandidateScorertronV2 m_candidateScoretronV2;

};



#endif //TARGETDECOYCANDIDATEPAIRSCORETRONV2_H
