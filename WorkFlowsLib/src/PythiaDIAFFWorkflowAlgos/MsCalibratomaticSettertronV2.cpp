//
// Created by andrewnichols on 7/24/25.
//

#include "MsCalibratomaticSettertronV2.h"

#include "ErrorUtils.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderPointerAcc.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretronV2.h"


MsCalibratomaticSettertronV2::MsCalibratomaticSettertronV2()
: m_tdcpManager(nullptr)
, m_msReaderPointerAcc(nullptr)
, m_pythiaParameters(nullptr)
, m_msCalibratomatic(nullptr)
{}

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
		const auto msScanInfosThird = static_cast<int>(msScanInfosVec.size() * 0.333333);
		msScanInfosVec = msScanInfosVec.mid(msScanInfosThird, msScanInfosThird);

		QMap<MzTargetKey, bool> temp;
		for (MsScanInfo *msi : msScanInfosVec) {
			temp[msi->targetKey()] = true;
		}

		QVector<MzTargetKey> mzTargetKeys = temp.keys().toVector();
		const auto mzTargetKeysQuartered = static_cast<int>(mzTargetKeys.size() * 0.25);
		mzTargetKeys = mzTargetKeys.mid(mzTargetKeysQuartered, mzTargetKeysQuartered * 2);
		temp.clear();
		for (const MzTargetKey &mzTargetKey : mzTargetKeys) {
			temp[mzTargetKey] = true;
		}

		for (MsScanInfo *msi : msScanInfosVec) {

			const MzTargetKey mzTargetKey = msi->targetKey();
			if (!temp.contains(mzTargetKey)) {
				continue;
			}

			(*mzTargetKeyVsMsScanInfos)[mzTargetKey].push_back(msi);;
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

		mzTargetKeyVsTurboXICs->clear();

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
Err MsCalibratomaticSettertronV2::init(
	TargetDecoyCandidatePairManager *tdcpManager,
	MsReaderPointerAcc *msReaderPointerAcc,
	PythiaParameters *pythiaParameters
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); ree;
	e = ErrorUtils::isTrue(tdcpManager->isInit()); ree;
	e = ErrorUtils::isTrue(pythiaParameters->isValid()); ree;

	m_tdcpManager = tdcpManager;
	m_msReaderPointerAcc = msReaderPointerAcc;
	m_pythiaParameters = pythiaParameters;

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

	e = buildMzTargetKeyVsTargetDecoyCandidatePairPntrs(); ree;

	ERR_RETURN
}

namespace {

	void filterUniqueScanInfosByMzTargetKey(
		const QList<MzTargetKey> &mzTargetKeys,
		QVector<MsScanInfo> *msScanInfos
		) {

		const auto terminatorLogic = [mzTargetKeys](const MsScanInfo &msi) {
			return !mzTargetKeys.contains(msi.targetKey());
		};

		const auto terminator = std::remove_if(
			msScanInfos->begin(),
			msScanInfos->end(),
			terminatorLogic
			);

		msScanInfos->erase(terminator, msScanInfos->end());
	}

	std::tuple<Err, MzTargetKey, QVector<TargetDecoyCandidatePair*>> buildMzTargetKeyVsTargetDecoyCandidatePairPntrsLogic(
		const MsScanInfo &msScanInfo,
		double precursorExtractionWindowThomsons,
		QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsPntrs
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairsPntrs); rtee;
		e = ErrorUtils::isGreaterThanZero(precursorExtractionWindowThomsons); rtee;

		const float mzMin = msScanInfo.precursorTargetMz
						  - (msScanInfo.isoWindowLower + static_cast<float>(precursorExtractionWindowThomsons));

		const float mzMax = msScanInfo.precursorTargetMz
						  + (msScanInfo.isoWindowUpper + static_cast<float>(precursorExtractionWindowThomsons));

		const auto terminatorLogic = [mzMin, mzMax](TargetDecoyCandidatePair *tdcp) {
			const float mzPrecursor = tdcp->mz(false);
			return mzPrecursor < mzMin || mzPrecursor > mzMax;
		};

		const auto terminator = std::remove_if(
			targetDecoyCandidatePairsPntrs.begin(),
			targetDecoyCandidatePairsPntrs.end(),
			terminatorLogic
			);

		targetDecoyCandidatePairsPntrs.erase(terminator, targetDecoyCandidatePairsPntrs.end());

		return {e, msScanInfo.targetKey(), targetDecoyCandidatePairsPntrs};
	}

}//namespace
Err MsCalibratomaticSettertronV2::buildMzTargetKeyVsTargetDecoyCandidatePairPntrs() {

	ERR_INIT

	e = ErrorUtils::isTrue(m_tdcpManager->isInit()); ree;
	e = ErrorUtils::isTrue(m_msReaderPointerAcc->isInit()); ree;
	e = ErrorUtils::isNotEmpty(m_mzTargetKeyVsTurboXICs); ree;

	m_mzTargetKeyVsTargetDecoyCandidatePairPntrs.clear();

	const QList<MzTargetKey> &mzTargetKeys = m_mzTargetKeyVsTurboXICs.keys();
	QVector<MsScanInfo> uniqueMsScanInfosFiltered = m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();
	filterUniqueScanInfosByMzTargetKey(mzTargetKeys, &uniqueMsScanInfosFiltered);

	QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsPntrs;
	e = m_tdcpManager->getTargetDecoyCandidatePairPointers(&targetDecoyCandidatePairsPntrs); rtee;

#define FILTER_PARALLEL_TDCP
#ifdef FILTER_PARALLEL_TDCP
	const auto loadLogicBinder = std::bind(
		buildMzTargetKeyVsTargetDecoyCandidatePairPntrsLogic,
		std::placeholders::_1,
		m_pythiaParameters->precursorExtractionWindowThomsons,
		targetDecoyCandidatePairsPntrs
		);

	QFuture<std::tuple<Err, MzTargetKey, QVector<TargetDecoyCandidatePair*>>> futures = QtConcurrent::mapped(
		uniqueMsScanInfosFiltered,
		loadLogicBinder
		);
	futures.waitForFinished();

	for (const std::tuple<Err, MzTargetKey, QVector<TargetDecoyCandidatePair*>> &tpl : futures) {

		e = std::get<0>(tpl); ree;

		const MzTargetKey mzTargetKey = std::get<1>(tpl);
		const QVector<TargetDecoyCandidatePair*> &tdcpPntrs = std::get<2>(tpl);
		for (TargetDecoyCandidatePair *tdcp : tdcpPntrs) {
			m_mzTargetKeyVsTargetDecoyCandidatePairPntrs.push_back({mzTargetKey, tdcp});
		}
	}
#else
	for (const MsScanInfo &msi : uniqueMsScanInfosFiltered) {
		std::tuple<Err, MzTargetKey, QVector<TargetDecoyCandidatePair*>> tpl = buildMzTargetKeyVsTargetDecoyCandidatePairPntrsLogic(
			msi,
			m_pythiaParameters->precursorExtractionWindowThomsons,
			targetDecoyCandidatePairsPntrs
			);
		e = std::get<0>(tpl); ree;
		m_mzTargetKeyVsTargetDecoyCandidatePairPntrs.insert(std::get<1>(tpl), std::get<2>(tpl));
	}
#endif

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished building TargetPairs";

	ERR_RETURN
}

namespace {

	Err parallelProcessingLogic(
		const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXICs,
		const PythiaParameters &pythiaParameters,
		const QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> &mzTargetKeyVsTargetDecoyCandidatePairPntrs,
		int meanFrameScanCountMS2
		) {

		ERR_INIT


		TargetDecoyCandidatePairScoretronV2 targetDecoyCandidatePairScoretronV2;
		e = targetDecoyCandidatePairScoretronV2.init(
			mzTargetKeyVsTurboXICs,
			pythiaParameters,
			meanFrameScanCountMS2,
			S_GLOBAL_SETTINGS.MIN_MS2_IONS
			); ree;

		for (const QPair<MzTargetKey, TargetDecoyCandidatePair*> &pr : mzTargetKeyVsTargetDecoyCandidatePairPntrs) {


			e = targetDecoyCandidatePairScoretronV2.scoreTargetDecoyCandidatePairPntr(pr); ree;


			// QVector<int> frameIndexVsScanNumber;
			// frameIndexVsScanNumber.reserve(xicPointsPntrs.size());
			// for (const XICPoint *xicPoint : xicPointsPntrs) {
			// 	if (!frameIndexVsScanNumber.isEmpty()) {
			// 		if (frameIndexVsScanNumber.back() == xicPoint->scanNumber) {
			// 			continue;
			// 		}
			// 	}
			// 	frameIndexVsScanNumber.push_back(xicPoint->scanNumber);
			// }



		}


		ERR_RETURN
	}

}//namespace
Err MsCalibratomaticSettertronV2::buildMsCalibratomatic(MsCalibratomatic *msCalibratomatic) {
	ERR_INIT

	m_msCalibratomatic = msCalibratomatic;

	QVector<QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>>> mzTargetDecoyCandidatePairsPntrsTranched;
	e = ParallelUtils::trancheVectorForParallelizationInOrder(
		m_mzTargetKeyVsTargetDecoyCandidatePairPntrs,
		m_pythiaParameters->threadCount,
		0,
		&mzTargetDecoyCandidatePairsPntrsTranched
		); ree;

#define CALIBRATION_PARALLEL_TDCP
#ifdef CALIBRATION_PARALLEL_TDCP
	const auto binderLogic = std::bind(
		parallelProcessingLogic,
		m_mzTargetKeyVsTurboXICs,
		*m_pythiaParameters,
		std::placeholders::_1,
		m_msReaderPointerAcc->ptr->meanFrameScanCountMS2()
		);

	QFuture<Err> futures = QtConcurrent::mapped(
		mzTargetDecoyCandidatePairsPntrsTranched,
		binderLogic
		);
	futures.waitForFinished();
#else
	for (const QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> &pr : mzTargetDecoyCandidatePairsPntrsTranched) {
		e = parallelProcessingLogic(m_mzTargetKeyVsTurboXICs, *m_pythiaParameters, pr); ree;
	}
#endif

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished building Calibration";

	ERR_RETURN
}