//
// Created by andrewnichols on 7/29/25.
//

#include "PythiaDIAFFWorkflowV2.h"

#include "CandidateScoresFeatureManager.h"
#include "FragLibReader.h"
#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertronV2.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderPointerAcc.h"


PythiaDIAFFWorkflowV2::PythiaDIAFFWorkflowV2() = default;

PythiaDIAFFWorkflowV2::~PythiaDIAFFWorkflowV2() {
	for (MsFrameV2 *f : m_mzTargetKeyVsMsFramesMS2Pntrs) {
		delete f;
	}
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
	e = buildMzTargetKeyVsMsFramesMS2Pntrs(&msReaderPointerAcc); ree;

	{

		const QVector<CandidateScoresFeatureManager::Features> featuresCalibration
			= CandidateScoresFeatureManager::featuresCalibration();

		MsCalibratomaticSettertronV2 msCalibratomaticSettertronV2;
		e = msCalibratomaticSettertronV2.init(
			m_mzTargetKeyVsUniqueMsScanInfoPntrs,
			m_mzTargetKeyVsMsFramesMS2Pntrs,
			featuresCalibration,
			&m_tdcpManager,
			&m_pythiaParameters,
			&m_msFrameMS1
			); ree;
		e = msCalibratomaticSettertronV2.buildMsCalibratomatic(&m_msCalibratomatic); ree;
	}

	ERR_RETURN
}

Err PythiaDIAFFWorkflowV2::buildMS1Frame(MsReaderPointerAcc *msReaderPointerAcc) {

	ERR_INIT

	constexpr int msLevel = 1;
	const QMap<ScanNumber, MsScanInfo*> ms1ScanInfos = msReaderPointerAcc->ptr->getMsScanInfos(msLevel);
	QVector<MsScan> msScans;

	e = msReaderPointerAcc->ptr->extractScanPoints(
		ms1ScanInfos.values().toVector(),
		&msScans
		); ree;

	e = m_msFrameMS1.init(msScans); ree;

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Ms1 Data Frames built";

	ERR_RETURN
}

namespace {

	std::tuple<Err, MsScanInfo*, MsFrameV2*> buildMzTargetKeyVsMsFramesMS2Logic(
		const QVector<MsScanInfo*> &mzTargetKeyMsScanInfoPntrs,
		MsReaderPointerAcc *msReaderPointerAcc
		) {

		ERR_INIT

		e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); rtee;

		QVector<MsScan> msScans;
		e = msReaderPointerAcc->ptr->extractScanPoints(
			mzTargetKeyMsScanInfoPntrs,
			&msScans
			); rtee;

		auto *msFrame = new MsFrameV2;
		e = msFrame->init(msScans);rtee;

		return {e, mzTargetKeyMsScanInfoPntrs.front(), msFrame};
	}

}//namespace
Err PythiaDIAFFWorkflowV2::buildMzTargetKeyVsMsFramesMS2Pntrs(MsReaderPointerAcc *msReaderPointerAcc) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); ree;

	constexpr int msLevel = 2;
	const QMap<ScanNumber, MsScanInfo*> scanNumbersVsMsScanInfoPntrs
									= msReaderPointerAcc->ptr->getMsScanInfos(msLevel);

	QMap<MzTargetKey, QVector<MsScanInfo*>> mzTargetKeysVsMsScanInfoPntrs;
	for (auto it = scanNumbersVsMsScanInfoPntrs.begin(); it != scanNumbersVsMsScanInfoPntrs.end(); it++) {
		mzTargetKeysVsMsScanInfoPntrs[it.value()->targetKey()].push_back(it.value());
	}

#define BUILD_MZ_TARGET_KEY_VS_MS_SCAN_INFO_PNTRS_PARALLEL
#ifdef BUILD_MZ_TARGET_KEY_VS_MS_SCAN_INFO_PNTRS_PARALLEL
	const auto binderLogic = std::bind(
		buildMzTargetKeyVsMsFramesMS2Logic,
		std::placeholders::_1,
		msReaderPointerAcc
		);

	QFuture<std::tuple<Err, MsScanInfo*, MsFrameV2*>> future = QtConcurrent::mapped(
		mzTargetKeysVsMsScanInfoPntrs,
		binderLogic
		);
	future.waitForFinished();

	for (const std::tuple<Err, MsScanInfo*, MsFrameV2*> result : future) {
		e = std::get<0>(result); ree;
		m_mzTargetKeyVsMsFramesMS2Pntrs.insert(std::get<1>(result)->targetKey(), std::get<2>(result));
		m_mzTargetKeyVsUniqueMsScanInfoPntrs.insert(std::get<1>(result)->targetKey(), std::get<1>(result));
	}
#else
	for (const QVector<MsScanInfo*> &msips : mzTargetKeysVsMsScanInfoPntrs) {
		const std::tuple<Err, MzTargetKey, MsFrameV2*> result = buildMzTargetKeyVsMsFramesMS2Logic(msips, msReaderPointerAcc);
		e = std::get<0>(result); ree;
		m_mzTargetKeyVsMsFramesMS2Pntrs.insert(std::get<1>(result), std::get<2>(result));
	}
#endif

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Ms2 Data Frames built";

	ERR_RETURN
}