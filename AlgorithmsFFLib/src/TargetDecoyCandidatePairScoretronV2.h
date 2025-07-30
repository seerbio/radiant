//
// Created by andrewnichols on 7/26/25.
//

#ifndef TARGETDECOYCANDIDATEPAIRSCORETRONV2_H
#define TARGETDECOYCANDIDATEPAIRSCORETRONV2_H

#include "AlgorithmsFFLib_Exports.h"

#include "GlobalSettings.h"
#include "Error.h"
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
		const QMap<MzTargetKey, QVector<MsScanInfo*>> &mzTargetKeyVsMsScanInfos,
		const PythiaParameters &pythiaParameters,
		int ms2IonsCount,
		float minMs2IonsFoundCount,
		MsReaderPointerAcc *msReaderPointerAcc,
		MsFrameV2 *msFrameV2MS1
		);

	Err scoreTargetDecoyCandidatePairPntr(
		const QPair<MzTargetKey, TargetDecoyCandidatePair*> &mzTargetKeyVsTdcpPntr
		);

private:

	Err scoreMS2Ions(
		const QVector<MS2Ion> &ms2Ions,
		bool isDecoy,
		TargetDecoyCandidatePair *tdcp
		);

	Err loadMS2IonArrays(const QVector<MS2Ion> &ms2Ions);

	Err subtractShadowsArrays();

	void zeroOutArrays();

	Err smoothMS2IonArrays();

	Err buildLocationVectors(const QVector<MS2Ion> &ms2Ions);

	Err buildIntegrationVecCosineSim(const QVector<MS2Ion> &ms2Ions);

	[[nodiscard]] Err buildProductVec() const;

	[[nodiscard]] Err buildMs1Vec(bool isDecoy, TargetDecoyCandidatePair *tdcp) const;

	Err fitMS1XICToVecs(
		const XICPointsPntrs &xicPointsPntrs,
		float* vecIntensity,
		float *vecMz
		) const;

	Err smoothScoreArrays();

private:

	MsReaderPointerAcc *m_msReaderPointerAcc;
	QMap<MzTargetKey, QVector<MsScanInfo*>> m_mzTargetKeyVsMsScanInfos;
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
	float* m_integrationVecCosineSim;
	float* m_productVec;

	float* m_mzMs1MonoIsotopeVecIntensity;
	float* m_mzMs1C13VecIntensity;
	float* m_mzMs1MonoIsotopeShadowVecIntensity;
	float* m_mzMs1MonoIsotopeVecMz;
	float* m_mzMs1C13VecMz;
	float* m_mzMs1MonoIsotopeShadowVecMz;

	int m_targetFrameIndexMax;
	size_t m_xicSizeTargetMaxAlignas;

	QVector<float> m_savitzkyGolayKernel;

	float m_mzMs2Min;
	float m_mzMs2Max;
	float m_intensityVecMax;

	float m_minMs2IonsFoundCount;

};



#endif //TARGETDECOYCANDIDATEPAIRSCORETRONV2_H
