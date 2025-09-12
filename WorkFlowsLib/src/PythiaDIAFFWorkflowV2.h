//
// Created by andrewnichols on 7/29/25.
//

#ifndef PYTHIADIAFFWORKFLOWV2_H
#define PYTHIADIAFFWORKFLOWV2_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibReaderRow.h"
#include "MsCalibratomatic.h"
#include "MsFrameV2.h"
#include "PythiaParameterReader.h"
#include "MsCalibratomatic.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;

class MsReaderPointerAcc;

class WORKFLOWSLIB_EXPORTS PythiaDIAFFWorkflowV2 {


  friend class PythiaDIAFFWorkflowTests;

public:

	PythiaDIAFFWorkflowV2();
	~PythiaDIAFFWorkflowV2();

	Err init(
		const PythiaParameters &pythiaParameters,
		const QString &fragLibUri
		);

	Err processFile(
		const QString &msDataFileUri
		);

private:

	Err buildMS1Frame(MsReaderPointerAcc *msReaderPointerAcc);
	Err buildMzTargetKeyVsMsFramesMS2Pntrs(MsReaderPointerAcc *msReaderPointerAcc);


private:

	MsCalibratomatic m_msCalibratomatic;
	QVector<FragLibReaderRow> m_fragLibReaderRows;
	MsFrameV2 m_msFrameMS1;
	QMap<MzTargetKey, MsFrameV2*> m_mzTargetKeyVsMsFramesMS2Pntrs;
	QMap<MzTargetKey, MsScanInfo*> m_mzTargetKeyVsUniqueMsScanInfoPntrs;
	PythiaParameters m_pythiaParameters;
	TargetDecoyCandidatePairManager m_tdcpManager;

};



#endif //PYTHIADIAFFWORKFLOWV2_H
