//
// Created by andrewnichols on 7/29/25.
//

#include "PythiaDIAFFWorkflowV2.h"

#include "FragLibReader.h"
#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertronV2.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderPointerAcc.h"


PythiaDIAFFWorkflowV2::PythiaDIAFFWorkflowV2() {
}

PythiaDIAFFWorkflowV2::~PythiaDIAFFWorkflowV2() {
}

Err PythiaDIAFFWorkflowV2::init(
	const PythiaParameters &pythiaParameters,
	const QString &fragLibUri
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
	e = ErrorUtils::fileExists(fragLibUri); ree;

	m_pythiaParameters = pythiaParameters;

	e = FragLibReader::getFragLibReaderRows(
		fragLibUri,
		&m_fragLibReaderRows
		); ree;

	e = m_tdcpManager.init(
		pythiaParameters,
		&m_fragLibReaderRows
		); ree;


	ERR_RETURN
}

Err PythiaDIAFFWorkflowV2::processFile(const QString &msDataFileUri) {

	ERR_INIT

	e = ErrorUtils::fileExists(msDataFileUri); ree;

	MsReaderPointerAcc msReaderPointerAcc;
	e = msReaderPointerAcc.openFile(msDataFileUri); ree;

	e = buildMS1Frame(&msReaderPointerAcc); ree;

	{
		MsCalibratomaticSettertronV2 setter;
		e = setter.init(
			&m_tdcpManager,
			&msReaderPointerAcc,
			&m_pythiaParameters
			); ree;
	}

	ERR_RETURN
}

Err PythiaDIAFFWorkflowV2::buildMS1Frame(MsReaderPointerAcc *msReaderPointerAcc) {

	ERR_INIT

	constexpr int msLevel = 1;
	const QMap<ScanNumber, MsScanInfo*> ms1ScanInfos = msReaderPointerAcc->ptr->getMsScanInfos(msLevel);
	QVector<MsScan> msScans;
	e = MsReaderMzMLLazyLoad::extractScanPoints(
		msReaderPointerAcc->ptr->filePath(),
		ms1ScanInfos.values().toVector(),
		&msScans
		); ree;

	e = m_msFrameMS1.init(msScans); ree;

	ERR_RETURN
}
