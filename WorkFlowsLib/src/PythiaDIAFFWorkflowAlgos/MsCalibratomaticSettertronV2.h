//
// Created by andrewnichols on 7/24/25.
//

#ifndef MSCALIBRATOMATICSETTERTRONV2_H
#define MSCALIBRATOMATICSETTERTRONV2_H

#include "WorkFlowsLib_Exports.h"

#include "GlobalSettings.h"
#include "Error.h"
#include "MsReaderBase.h"
#include "TurboXIC.h"

using namespace Error;

class MsReaderPointerAcc;

class WORKFLOWSLIB_EXPORTS MsCalibratomaticSettertronV2 {

public:

	MsCalibratomaticSettertronV2() = default;
	~MsCalibratomaticSettertronV2();

	Err init(MsReaderPointerAcc *msReaderPointerAcc);

private:

	QMap<ScanNumber, MsScanInfo> m_msScanInfos;
	QVector<MsScanInfo*> m_msScanInfosPntrs;
	QMap<MzTargetKey, TurboXIC*> m_mzTargetKeyVsTurboXICs;

};



#endif //MSCALIBRATOMATICSETTERTRONV2_H
