//
// Created by andrewnichols on 7/26/25.
//

#ifndef TARGETDECOYCANDIDATEPAIRSCORETRONV2_H
#define TARGETDECOYCANDIDATEPAIRSCORETRONV2_H

#include "AlgorithmsFFLib_Exports.h"

#include "GlobalSettings.h"
#include "Error.h"
#include "PythiaParameterReader.h"
#include "TurboXIC.h"

using namespace Error;

class MS2Ion;
class TargetDecoyCandidatePair;

class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePairScoretronV2 {

public:

	TargetDecoyCandidatePairScoretronV2();
	~TargetDecoyCandidatePairScoretronV2();

	Err init(
		const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXICs,
		const PythiaParameters &pythiaParameters,
		int meanFrameScanCountMS2,
		int ms2IonsCount
		);

	Err scoreTargetDecoyCandidatePairPntr(
		const QPair<MzTargetKey, TargetDecoyCandidatePair*> &mzTargetKeyVsTdcpPntr
		);

private:

	Err scoreMS2Ions(const QVector<MS2Ion> &ms2Ions);
	Err loadMS2IonArrays(const QVector<MS2Ion> &ms2Ions);
	Err subtractShadowsArrays();
	void zeroOutArrays();
	Err smoothMS2IonArrays();
	Err buildLocationVectors();

private:

	QMap<MzTargetKey, TurboXIC*> m_mzTargetKeyVsTurboXICs;
	PythiaParameters m_pythiaParameters;

	MzTargetKey m_mzTargetKeyCurrent;
	TurboXIC *m_turboXICCurrent;

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

	int m_frameIndexTargetMax;

	QVector<float> m_savitzkyGolayKernel;

};



#endif //TARGETDECOYCANDIDATEPAIRSCORETRONV2_H
