//
// Created by andrewnichols on 7/24/25.
//

#ifndef MSCALIBRATOMATICSETTERTRONV2_H
#define MSCALIBRATOMATICSETTERTRONV2_H

#include "WorkFlowsLib_Exports.h"

#include "GlobalSettings.h"
#include "Error.h"
#include "MsCalibratomatic.h"
#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TurboXIC.h"

using namespace Error;

class MsReaderPointerAcc;
class TargetDecoyCandidatePairManager;

class WORKFLOWSLIB_EXPORTS MsCalibratomaticSettertronV2 {

public:

	MsCalibratomaticSettertronV2();
	~MsCalibratomaticSettertronV2();

	Err init(
		TargetDecoyCandidatePairManager *tdcpManager,
		MsReaderPointerAcc *msReaderPointerAcc,
		PythiaParameters *pythiaParameters
		);

	Err buildMsCalibratomatic(MsCalibratomatic *msCalibratomatic);

private:

	Err buildMzTargetKeyVsTargetDecoyCandidatePairPntrs();

private:

	QMap<ScanNumber, MsScanInfo> m_msScanInfos;
	QVector<MsScanInfo*> m_msScanInfosPntrs;
	QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> m_mzTargetKeyVsTargetDecoyCandidatePairPntrs;
	TargetDecoyCandidatePairManager *m_tdcpManager;
	MsReaderPointerAcc *m_msReaderPointerAcc;
	PythiaParameters *m_pythiaParameters;
	MsCalibratomatic *m_msCalibratomatic;
	QMap<MzTargetKey, QVector<MsScanInfo*>> m_mzTargetKeyVsMsScanInfos;

};



#endif //MSCALIBRATOMATICSETTERTRONV2_H
