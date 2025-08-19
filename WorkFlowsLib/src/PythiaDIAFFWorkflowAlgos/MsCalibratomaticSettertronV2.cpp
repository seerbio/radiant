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
, m_pythiaParameters(nullptr)
, m_msCalibratomatic(nullptr)
, m_msFrameMS1(nullptr)
{}

MsCalibratomaticSettertronV2::~MsCalibratomaticSettertronV2() = default;

Err MsCalibratomaticSettertronV2::init(
	const QMap<MzTargetKey, MsScanInfo*> &mzTargetKeyVsUniqueMsScanInfoPntrs,
	const QMap<MzTargetKey, MsFrameV2*> &mzTargetKeyVsMsFramesMS2Pntrs,
	const QVector<CandidateScoresFeatureManager::Features> &featuresCalibration,
	TargetDecoyCandidatePairManager *tdcpManager,
	PythiaParameters *pythiaParameters,
	MsFrameV2 *msFrameMS1
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(mzTargetKeyVsUniqueMsScanInfoPntrs); ree;
	e = ErrorUtils::isNotEmpty(mzTargetKeyVsMsFramesMS2Pntrs); ree;
	e = ErrorUtils::isTrue(tdcpManager->isInit()); ree;
	e = ErrorUtils::isTrue(pythiaParameters->isValid()); ree;
	e = ErrorUtils::isTrue(msFrameMS1->isInit()); ree;
	e = ErrorUtils::isNotEmpty(featuresCalibration); ree;

	m_tdcpManager = tdcpManager;
	m_pythiaParameters = pythiaParameters;
	m_msFrameMS1 = msFrameMS1;
	m_mzTargetKeyVsUniqueMsScanInfoPntrs = mzTargetKeyVsUniqueMsScanInfoPntrs;
	m_mzTargetKeyVsMsFramesMS2Pntrs = mzTargetKeyVsMsFramesMS2Pntrs;
	m_featuresCalibration = featuresCalibration;

	e = buildMzTargetKeyVsTargetDecoyCandidatePairPntrs(); ree;

	ERR_RETURN
}

namespace {

	Err buildMsScanInfosPntrsForCalibration(
		const QList<MzTargetKey> &mzTargetKeys,
		const QMap<MzTargetKey, MsScanInfo*> &mzTargetKeyVsUniqueMsScanInfoPntrs,
		QVector<MsScanInfo*> *msScanInfosPntrsForCalibration
		) {

		ERR_INIT

		for (int i = 0; i < mzTargetKeys.size(); ++i) {

			const MzTargetKey &mzTargetKey = mzTargetKeys.at(i);
			// if (constexpr int skipCount = 5; i % skipCount != 0) {
			// 	continue;
			// }

			e = ErrorUtils::contains(mzTargetKey, mzTargetKeyVsUniqueMsScanInfoPntrs);
			msScanInfosPntrsForCalibration->push_back(mzTargetKeyVsUniqueMsScanInfoPntrs.value(mzTargetKey));
		}

		ERR_RETURN
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
	e = ErrorUtils::isNotEmpty(m_mzTargetKeyVsUniqueMsScanInfoPntrs); ree;

	m_mzTargetKeyVsTargetDecoyCandidatePairPntrs.clear();

	const QList<MzTargetKey> &mzTargetKeys = m_mzTargetKeyVsUniqueMsScanInfoPntrs.keys();

	QVector<MsScanInfo*> msScanInfosPntrsForCalibration;
	e = buildMsScanInfosPntrsForCalibration(
		mzTargetKeys,
		m_mzTargetKeyVsUniqueMsScanInfoPntrs,
		&msScanInfosPntrsForCalibration
		); ree;

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
		msScanInfosPntrsForCalibration,
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
	for (const MsScanInfo *msi : msScanInfosPntrsForCalibration) {
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
		const QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> &mzTargetKeyVsTargetDecoyCandidatePairPntrs,
		const QMap<MzTargetKey, MsFrameV2*> &mzTargetKeyVsMsFramesMS2Pntrs,
		const QVector<CandidateScoresFeatureManager::Features> &featuresCalibration,
		const PythiaParameters &pythiaParameters,
		MsFrameV2 *msFrameV2MS1
		) {

		ERR_INIT

		QElapsedTimer et;
		et.start();

		QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePairPntrsSorted
																			= mzTargetKeyVsTargetDecoyCandidatePairPntrs;

		using T = QPair<MzTargetKey, TargetDecoyCandidatePair*>;
		std::sort(
			mzTargetKeyVsTargetDecoyCandidatePairPntrsSorted.begin(),
			mzTargetKeyVsTargetDecoyCandidatePairPntrsSorted.end(),
			[](const T &l, const T &r){return l.first < r.first;}
			);

		constexpr float minMs2IonsFoundCount = 2.22;

		TargetDecoyCandidatePairScoretronV2 targetDecoyCandidatePairScoretronV2;
		e = targetDecoyCandidatePairScoretronV2.init(
			mzTargetKeyVsMsFramesMS2Pntrs,
			featuresCalibration,
			pythiaParameters,
			S_GLOBAL_SETTINGS.MIN_MS2_IONS,
			minMs2IonsFoundCount,
			msFrameV2MS1
			); ree;

		constexpr int skipCountTDCP = 5;
		int counter = 0;
		for (const QPair<MzTargetKey, TargetDecoyCandidatePair*> &pr : mzTargetKeyVsTargetDecoyCandidatePairPntrsSorted) {

			if (counter++ % skipCountTDCP != 0) {
				continue;
			}

			e = targetDecoyCandidatePairScoretronV2.scoreTargetDecoyCandidatePairPntr(pr); ree;
		}

		if (false) { //TODO add verbose variable to parallel processing logic to replace false;
			qDebug()
			<< qPrintable(S_GLOBAL_TIMER.elapsed())
			<< "Finished processing"
			<< mzTargetKeyVsTargetDecoyCandidatePairPntrs.size()
			<< "targets in"
			<< et.elapsed()
			<< "mSec | MzTargetKeys"
			<< mzTargetKeyVsTargetDecoyCandidatePairPntrs.front().first
			<< "to"
			<< mzTargetKeyVsTargetDecoyCandidatePairPntrs.back().first
			;
		}

		ERR_RETURN
	}

}//namespace
Err MsCalibratomaticSettertronV2::buildMsCalibratomatic(MsCalibratomatic *msCalibratomatic) {
	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_mzTargetKeyVsUniqueMsScanInfoPntrs); ree;
	e = ErrorUtils::isNotEmpty(m_featuresCalibration); ree;
	e = ErrorUtils::isTrue(m_pythiaParameters->isValid()); ree;

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
		std::placeholders::_1,
		m_mzTargetKeyVsMsFramesMS2Pntrs,
		m_featuresCalibration,
		*m_pythiaParameters,
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
			pr,
			m_mzTargetKeyVsMsFramesMS2Pntrs,
			m_featuresCalibration,
			*m_pythiaParameters,
			m_msFrameMS1
			); ree;
	}
#endif

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished building Calibration";

	ERR_RETURN
}