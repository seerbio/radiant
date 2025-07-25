//
// Created by andrewnichols on 7/24/25.
//

#include "MsCalibratomaticSettertronV2.h"

#include "ErrorUtils.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderPointerAcc.h"


MsCalibratomaticSettertronV2::~MsCalibratomaticSettertronV2() {
	for (const TurboXIC *turboXic : m_mzTargetKeyVsTurboXICs) {
		delete turboXic;
	}
}

namespace {

	Err buildMzTargetKeyVsMsScanInfosTrunc(
		const QVector<MsScanInfo*> &msScanInfosPntrs,
		QMap<MzTargetKey, QVector<MsScanInfo*>> *mzTargetKeyVsMsScanInfos
		) {

		ERR_INIT

		QVector<MsScanInfo*> msScanInfosVec = msScanInfosPntrs;
		const auto msScanInfosThird = static_cast<int>(msScanInfosVec.size() / 3.0);
		msScanInfosVec = msScanInfosVec.mid(msScanInfosThird, msScanInfosThird);

		for (MsScanInfo *msi : msScanInfosVec) {
			(*mzTargetKeyVsMsScanInfos)[msi->targetKey()].push_back(msi);
		}

		ERR_RETURN
	}

	std::tuple<Err, MzTargetKey, TurboXIC*> parallelLoadTurboXIC(
		const QVector<MsScanInfo*> &msScanInfos,
		const QString &filePath
		) {

		ERR_INIT

		QVector<MsScanInfo*> msScanInfosCopy = msScanInfos;
		QVector<MsScan> msScans;
		e = MsReaderMzMLLazyLoad::extractScanPoints(
			filePath,
			msScanInfosCopy,
			&msScans
			); rtee;

		e = ErrorUtils::isTrue(
			msScans.front().msScanInfoPntr->targetKey() ==
			msScans.back().msScanInfoPntr->targetKey()
			); rtee;

		auto *turboXIC = new TurboXIC();
		e = turboXIC->init(msScans); rtee;

		return {e, msScans.back().msScanInfoPntr->targetKey(), turboXIC};
	}

	Err buildTurboXICs(
		const QString &filePath,
		const QVector<MsScanInfo*> &msScanInfosPntrs,
		QMap<MzTargetKey, TurboXIC*> *mzTargetKeyVsTurboXICs
		) {

		ERR_INIT

		QElapsedTimer et;
		et.start();

		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Building calibration TurboXIC";

		QMap<MzTargetKey, QVector<MsScanInfo*>> mzTargetKeyVsMsScanInfos;
		e = buildMzTargetKeyVsMsScanInfosTrunc(
			msScanInfosPntrs,
			&mzTargetKeyVsMsScanInfos
			); ree;

		const auto binderLogic = std::bind(
			parallelLoadTurboXIC,
			std::placeholders::_1,
			filePath
			);

		QFuture<std::tuple<Err, MzTargetKey, TurboXIC*>> futures = QtConcurrent::mapped(
			mzTargetKeyVsMsScanInfos,
			binderLogic
			);
		futures.waitForFinished();

		for (const std::tuple<Err, MzTargetKey, TurboXIC*> &tpl : futures) {
			e = std::get<0>(tpl); ree;
			mzTargetKeyVsTurboXICs->insert(std::get<1>(tpl), std::get<2>(tpl));
		}
		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished building calibration TurboXIC";

		ERR_RETURN
	}

}//namespace
Err MsCalibratomaticSettertronV2::init(MsReaderPointerAcc *msReaderPointerAcc) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); ree;

	constexpr int msLevel = 2;
	m_msScanInfos = msReaderPointerAcc->ptr->getMsScanInfos(msLevel);
	e = ErrorUtils::isNotEmpty(m_msScanInfos); ree;

	m_msScanInfosPntrs.reserve(m_msScanInfos.size());
	for (MsScanInfo &msi : m_msScanInfos) {
		m_msScanInfosPntrs.push_back(&msi);
	}

	e = buildTurboXICs(
		msReaderPointerAcc->ptr->filePath(),
		m_msScanInfosPntrs,
		&m_mzTargetKeyVsTurboXICs
		); ree;







	ERR_RETURN
}
