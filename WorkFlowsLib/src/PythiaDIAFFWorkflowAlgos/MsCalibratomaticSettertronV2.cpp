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

MsCalibratomaticSettertronV2::~MsCalibratomaticSettertronV2() = default;

namespace {

	Err buildMzTargetKeyVsMsScanInfosTrunc(
		const QVector<MsScanInfo*> &msScanInfosPntrs,
		QMap<MzTargetKey, QVector<MsScanInfo*>> *mzTargetKeyVsMsScanInfos
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(msScanInfosPntrs); ree;

		QVector<MsScanInfo*> msScanInfosVec = msScanInfosPntrs;
		const auto msScanInfosTenthed = static_cast<int>(msScanInfosVec.size() * 0.1);
		msScanInfosVec = msScanInfosVec.mid(msScanInfosTenthed * 1, static_cast<int>(msScanInfosTenthed * 8));

		QMap<MzTargetKey, bool> allTargetKeys;
		for (MsScanInfo *msi : msScanInfosVec) {
			allTargetKeys[msi->targetKey()] = true;
		}

		QVector<MzTargetKey> mzTargetKeys = allTargetKeys.keys().toVector();

		constexpr int skipCount = 5;
		QMap<MzTargetKey, bool> selectTargetKeys;
		for (int i = 0; i < mzTargetKeys.size(); i += skipCount) {
			const MzTargetKey &mzTargetKey = mzTargetKeys[i];
			selectTargetKeys[mzTargetKey] = true;
		}

		for (MsScanInfo *msi : msScanInfosVec) {

			const MzTargetKey mzTargetKey = msi->targetKey();
			if (!selectTargetKeys.contains(mzTargetKey)) {
				continue;
			}

			(*mzTargetKeyVsMsScanInfos)[mzTargetKey].push_back(msi);;
		}

		ERR_RETURN
	}

}//namespace
Err MsCalibratomaticSettertronV2::init(
	TargetDecoyCandidatePairManager *tdcpManager,
	MsReaderPointerAcc *msReaderPointerAcc,
	PythiaParameters *pythiaParameters,
	MsFrameV2 *msFrameMS1
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); ree;
	e = ErrorUtils::isTrue(tdcpManager->isInit()); ree;
	e = ErrorUtils::isTrue(pythiaParameters->isValid()); ree;
	e = ErrorUtils::isTrue(msFrameMS1->isInit()); ree;

	m_tdcpManager = tdcpManager;
	m_msReaderPointerAcc = msReaderPointerAcc;
	m_pythiaParameters = pythiaParameters;
	m_msFrameMS1 = msFrameMS1;

	constexpr int msLevel = 2;
	m_msScanInfosPntrs = m_msReaderPointerAcc->ptr->getMsScanInfos(msLevel).values().toVector();
	e = ErrorUtils::isNotEmpty(m_msScanInfosPntrs); ree;

	e = buildMzTargetKeyVsTargetDecoyCandidatePairPntrs(); ree;

	ERR_RETURN
}

namespace {

	void filterUniqueScanInfosByMzTargetKey(
		const QList<MzTargetKey> &mzTargetKeys,
		QVector<MsScanInfo*> *msScanInfos
		) {

		const auto terminatorLogic = [mzTargetKeys](const MsScanInfo *msi) {
			return !mzTargetKeys.contains(msi->targetKey());
		};

		const auto terminator = std::remove_if(
			msScanInfos->begin(),
			msScanInfos->end(),
			terminatorLogic
			);

		msScanInfos->erase(terminator, msScanInfos->end());
	}

	std::tuple<Err, MzTargetKey, QVector<TargetDecoyCandidatePair*>> buildMzTargetKeyVsTargetDecoyCandidatePairPntrsLogic(
		const MsScanInfo* msScanInfo,
		double precursorExtractionWindowThomsons,
		QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsPntrs
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairsPntrs); rtee;
		e = ErrorUtils::isGreaterThanZero(precursorExtractionWindowThomsons); rtee;

		const float mzMin = msScanInfo->precursorTargetMz
						  - (msScanInfo->isoWindowLower + static_cast<float>(precursorExtractionWindowThomsons));

		const float mzMax = msScanInfo->precursorTargetMz
						  + (msScanInfo->isoWindowUpper + static_cast<float>(precursorExtractionWindowThomsons));

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

		return {e, msScanInfo->targetKey(), targetDecoyCandidatePairsPntrs};
	}

}//namespace
Err MsCalibratomaticSettertronV2::buildMzTargetKeyVsTargetDecoyCandidatePairPntrs() {

	ERR_INIT

	e = ErrorUtils::isTrue(m_tdcpManager->isInit()); ree;
	e = ErrorUtils::isNotEmpty(m_msScanInfosPntrs); ree;

	m_mzTargetKeyVsTargetDecoyCandidatePairPntrs.clear();

	e = buildMzTargetKeyVsMsScanInfosTrunc(
		m_msScanInfosPntrs,
		&m_mzTargetKeyVsMsScanInfos
		); ree;

	const QList<MzTargetKey> &mzTargetKeys = m_mzTargetKeyVsMsScanInfos.keys();
	QVector<MsScanInfo*> uniqueMsScanInfosFiltered = m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();
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
	for (const MsScanInfo *msi : uniqueMsScanInfosFiltered) {
		std::tuple<Err, MzTargetKey, QVector<TargetDecoyCandidatePair*>> tpl = buildMzTargetKeyVsTargetDecoyCandidatePairPntrsLogic(
			msi,
			m_pythiaParameters->precursorExtractionWindowThomsons,
			targetDecoyCandidatePairsPntrs
			);
		const MzTargetKey mzTargetKey = std::get<1>(tpl);
		const QVector<TargetDecoyCandidatePair*> &tdcpPntrs = std::get<2>(tpl);
		for (TargetDecoyCandidatePair *tdcp : tdcpPntrs) {
			m_mzTargetKeyVsTargetDecoyCandidatePairPntrs.push_back({mzTargetKey, tdcp});
		}
	}
#endif

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished building TargetPairs";

	ERR_RETURN
}

namespace {

	Err parallelProcessingLogic(
		const QMap<MzTargetKey, QVector<MsScanInfo*>> &mzTargetKeyVsMsScanInfos,
		const PythiaParameters &pythiaParameters,
		const QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> &mzTargetKeyVsTargetDecoyCandidatePairPntrs,
		MsReaderPointerAcc *msReaderPointerAcc,
		MsFrameV2 *msFrameV2MS1
		) {

		ERR_INIT

		constexpr float minMs2IonsFoundCount = 4.9;

		TargetDecoyCandidatePairScoretronV2 targetDecoyCandidatePairScoretronV2;
		e = targetDecoyCandidatePairScoretronV2.init(
			mzTargetKeyVsMsScanInfos,
			pythiaParameters,
			S_GLOBAL_SETTINGS.MIN_MS2_IONS,
			minMs2IonsFoundCount,
			msReaderPointerAcc,
			msFrameV2MS1
			); ree;

		for (const QPair<MzTargetKey, TargetDecoyCandidatePair*> &pr : mzTargetKeyVsTargetDecoyCandidatePairPntrs) {

			e = targetDecoyCandidatePairScoretronV2.scoreTargetDecoyCandidatePairPntr(pr); ree;

		}


		ERR_RETURN
	}

}//namespace
Err MsCalibratomaticSettertronV2::buildMsCalibratomatic(MsCalibratomatic *msCalibratomatic) {
	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_mzTargetKeyVsMsScanInfos); ree;
	e = ErrorUtils::isTrue(m_pythiaParameters->isValid()); ree;
	e = ErrorUtils::isTrue(m_msReaderPointerAcc->isInit()); ree;

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
		m_mzTargetKeyVsMsScanInfos,
		*m_pythiaParameters,
		std::placeholders::_1,
		m_msReaderPointerAcc,
		m_msFrameMS1
		);

	QFuture<Err> futures = QtConcurrent::mapped(
		mzTargetDecoyCandidatePairsPntrsTranched,
		binderLogic
		);
	futures.waitForFinished();
#else
	for (const QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> &pr : mzTargetDecoyCandidatePairsPntrsTranched) {
		e = parallelProcessingLogic(
			m_mzTargetKeyVsMsScanInfos,
			*m_pythiaParameters,
			pr,
			m_msReaderPointerAcc
			); ree;
	}
#endif

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished building Calibration";

	ERR_RETURN
}