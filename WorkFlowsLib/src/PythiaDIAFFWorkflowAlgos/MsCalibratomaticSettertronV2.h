//
// Created by andrewnichols on 7/24/25.
//

#ifndef MSCALIBRATOMATICSETTERTRONV2_H
#define MSCALIBRATOMATICSETTERTRONV2_H

#include "WorkFlowsLib_Exports.h"

#include "CandidateScoresFeatureManager.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "MsCalibratomatic.h"
#include "MsFrameV2.h"
#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;

class MsReaderPointerAcc;
class TargetDecoyCandidatePairManager;

class WORKFLOWSLIB_EXPORTS MsCalibratomaticSettertronV2 {

public:

	MsCalibratomaticSettertronV2();
	~MsCalibratomaticSettertronV2();

	Err init(
		const QMap<MzTargetKey, MsScanInfo*> &mzTargetKeyVsUniqueMsScanInfoPntrs,
		const QMap<MzTargetKey, MsFrameV2*> &mzTargetKeyVsMsFramesMS2Pntrs,
		const QVector<CandidateScoresFeatureManager::Features> &featuresCalibration,
		TargetDecoyCandidatePairManager *tdcpManager,
		PythiaParameters *pythiaParameters,
		MsFrameV2 *msFrameMS1
		);

	Err buildMsCalibratomatic(MsCalibratomatic *msCalibratomatic);

private:

	Err buildMzTargetKeyVsTargetDecoyCandidatePairPntrs();

private:

	QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> m_mzTargetKeyVsTargetDecoyCandidatePairPntrs;
	TargetDecoyCandidatePairManager *m_tdcpManager;

	PythiaParameters *m_pythiaParameters;
	MsCalibratomatic *m_msCalibratomatic;
	QVector<CandidateScoresFeatureManager::Features> m_featuresCalibration;
	QMap<MzTargetKey, MsScanInfo*> m_mzTargetKeyVsUniqueMsScanInfoPntrs;
	QMap<MzTargetKey, MsFrameV2*> m_mzTargetKeyVsMsFramesMS2Pntrs;
	MsFrameV2 *m_msFrameMS1;

};



#endif //MSCALIBRATOMATICSETTERTRONV2_H
