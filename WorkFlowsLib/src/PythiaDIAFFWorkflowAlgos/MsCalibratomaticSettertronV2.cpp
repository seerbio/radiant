//
// Created by andrewnichols on 7/24/25.
//

#include "MsCalibratomaticSettertronV2.h"

#include "ErrorUtils.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderPointerAcc.h"
#include "QValueSettertronV2.h"
#include "TargetDecoyCandidatePairManager.h"
#include "TargetDecoyCandidatePairScoretronV2.h"

#include "PythiaDIAFFWorkflowAlgos/PythiaDIAFFWorkflowSharedMethods.h"

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
	const QVector<FTR> &featuresCalibration,
	TargetDecoyCandidatePairManager *tdcpManager,
	PythiaParameters *pythiaParameters,
	MsFrameV2 *msFrameMS1
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(mzTargetKeyVsUniqueMsScanInfoPntrs); ree;
	e = ErrorUtils::isNotEmpty(mzTargetKeyVsMsFramesMS2Pntrs); ree;
	e = ErrorUtils::isTrue(tdcpManager->isInit()); ree;
	e = ErrorUtils::isTrue(pythiaParameters->isValid()); ree;
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

	QPair<Err, QVector<QPair<CandidateScoresV2Target, CandidateScoresV2Decoy>>> parallelProcessingLogic(
		const QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> &mzTargetKeyVsTargetDecoyCandidatePairPntrs,
		const QMap<MzTargetKey, MsFrameV2*> &mzTargetKeyVsMsFramesMS2Pntrs,
		const QVector<FTR> &featuresCalibration,
		const PythiaParameters &pythiaParameters,
		int skipCount,
		MsFrameV2 *msFrameV2MS1
		) {

		ERR_INIT

		QElapsedTimer et;
		et.start();

		QVector<QPair<CandidateScoresV2Target, CandidateScoresV2Decoy>> scoresTargetVsDecoyPairs;

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
			); rree;

		int counter = 0;
		for (const QPair<MzTargetKey, TargetDecoyCandidatePair*> &pr : mzTargetKeyVsTargetDecoyCandidatePairPntrsSorted) {

			if (counter++ % skipCount != 0) {
				continue;
			}
			QPair<CandidateScoresV2Target, CandidateScoresV2Decoy> scoresTargetVsDecoyPair;
			e = targetDecoyCandidatePairScoretronV2.scoreTargetDecoyCandidatePairPntr(
				pr,
				&scoresTargetVsDecoyPair
				); rree;

			if (!scoresTargetVsDecoyPair.first.targetDecoyCandidatePair) {
				// std::cout << "No target found " << scoresTargetVsDecoyPair.first.targetDecoyCandidatePair  << std::endl;
				continue;
			}

			if (!scoresTargetVsDecoyPair.second.targetDecoyCandidatePair) {
				// std::cout << "No decoy found " << scoresTargetVsDecoyPair.second.targetDecoyCandidatePair  << std::endl;
				scoresTargetVsDecoyPair.second.targetDecoyCandidatePair = pr.second;
				scoresTargetVsDecoyPair.second.isDecoy = true;
				scoresTargetVsDecoyPair.second.initFeaturesArray();
			}

			scoresTargetVsDecoyPairs.push_back(scoresTargetVsDecoyPair);

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

		return {e, scoresTargetVsDecoyPairs};
	}

	Err countScoreCandidatesByFDR(
		const QVector<CandidateScoresV2*> &candidateScores,
		double qValueThreshold,
		int *targetCountBelowFDRThreshold
) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(candidateScores); ree;
		e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

		const auto countLogic = [qValueThreshold](const CandidateScoresV2 *cs){
			return cs->qValue < qValueThreshold;
		};

		*targetCountBelowFDRThreshold
				= static_cast<int>(std::count_if(candidateScores.begin(), candidateScores.end(), countLogic));

		ERR_RETURN
	}

	Err outPutFDRCounts(
		const QMap<int, int> &fdrVsCount,
		QString *outputString
		) {

		ERR_INIT
		e = ErrorUtils::isNotEmpty(fdrVsCount); ree;

		QStringList builder;
		for (auto it = fdrVsCount.begin(); it != fdrVsCount.end(); ++it) {
			const QString &k = QString::number(it.key());
			builder.push_back(k + "%: " + QString::number(it.value()));
		}

		*outputString = builder.join(" | ");

		ERR_RETURN
	}

	Err outputFDRResults(
		QVector<CandidateScoresV2*> &candidateScores,
		int verbose,
		QMap<int, int> *fdrVsCount
		) {

		ERR_INIT

		const QVector<double> fdrFractions = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01};
		for (double fdrThresh : fdrFractions) {
			int foundAtThreshold;
			e = countScoreCandidatesByFDR(
					candidateScores,
					fdrThresh,
					&foundAtThreshold
			); ree;
			const int fdrPercent = static_cast<int>(fdrThresh * 100);
			fdrVsCount->insert(fdrPercent, foundAtThreshold);
		}

		QString outputString;
		e = outPutFDRCounts(*fdrVsCount, &outputString); ree;

		if (verbose > -1) {
			qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PSMs found:" << qPrintable(outputString);
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

	QVector<QPair<CandidateScoresV2Target, CandidateScoresV2Decoy>> scoresTargetVsTargetDecoyPairs;
	QVector<CandidateScoresV2*> candidateScores;

	const QVector<int> skipCounts = {4,2,1};
	for (int skipCount : skipCounts) {

		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Calibration skip count:" << skipCount;

#define CALIBRATION_PARALLEL_TDCP
#ifdef CALIBRATION_PARALLEL_TDCP
		const auto binderLogic = std::bind(
			parallelProcessingLogic,
			std::placeholders::_1,
			m_mzTargetKeyVsMsFramesMS2Pntrs,
			m_featuresCalibration,
			*m_pythiaParameters,
			skipCount,
			m_msFrameMS1
			);

		QFuture<QPair<Err, QVector<QPair<CandidateScoresV2Target, CandidateScoresV2Decoy>>>> futures = QtConcurrent::mapped(
			mzTargetDecoyCandidatePairsPntrsTranched,
			binderLogic
			);
		futures.waitForFinished();

		for (const QPair<Err, QVector<QPair<CandidateScoresV2Target, CandidateScoresV2Decoy>>> &pr : futures) {
			e = pr.first; ree;
			scoresTargetVsTargetDecoyPairs.append(pr.second);
		}

		using T = QPair<CandidateScoresV2Target, CandidateScoresV2Decoy>;

		for (T &pr : scoresTargetVsTargetDecoyPairs) {
			candidateScores.push_back(&pr.first);
			candidateScores.push_back(&pr.second);
		};

		std::sort(
			candidateScores.rbegin(),
			candidateScores.rend(),
			[](const CandidateScoresV2 *l, const CandidateScoresV2 *r) {
				return l->featuresArray[FTR::CosineSimSumMeanCorrelation] < r->featuresArray[FTR::CosineSimSumMeanCorrelation];
			});

		e = QValueSettertronV2::setQValueForCandidates(
			QValueSettertronV2::QValueScoreTypeV2::CosineSimSumMeanCorrelation,
			&candidateScores
			); ree;

		QMap<int, int> fdrVsCounts;
		e = outputFDRResults(
		    candidateScores,
		    m_pythiaParameters->verbosity,
		    &fdrVsCounts
		    ); ree;

		constexpr int rtTrainingVolume = 1000;
		if (constexpr int fdrPercent = 1; fdrVsCounts.value(fdrPercent) < rtTrainingVolume) {
			continue;
		}

		e = predictScanTimesForCandidateScores(candidateScores); ree;
		break;

#else
		for (const QVector<QPair<MzTargetKey, TargetDecoyCandidatePair*>> &pr : mzTargetDecoyCandidatePairsPntrsTranched) {
			const QPair<Err, QVector<QPair<CandidateScoresV2Target, CandidateScoresV2Decoy>>> result = parallelProcessingLogic(
				pr,
				m_mzTargetKeyVsMsFramesMS2Pntrs,
				m_featuresCalibration,
				*m_pythiaParameters,
				skipCount,
				m_msFrameMS1
				);
			e = result.first; ree;
		}
#endif
	}//end loop

	for (const CandidateScoresV2 *cs : candidateScores) {
		qDebug() << "SDLFKJDSL" << cs->featuresArray;
	}


	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished building Calibration";

	ERR_RETURN
}

namespace {

	Err predictScanTimesForCandidateScoresParallelLogic(
		MsCalibratomatic *msCalibratomatic,
		CandidateScoresV2 *cs
		) {

		ERR_INIT

		e = msCalibratomatic->predictScanTime(
			cs->targetDecoyCandidatePair->iRt(),
			&cs->scanTimePredicted
			);

		ERR_RETURN
	}

	Err buildRTCurves(
		const QVector<CandidateScoresV2*> &candidateScores,
		const PythiaParameters *pythiaParameters,
		int topN,
		MsCalibratomatic *msCalibratomatic
		) {

		ERR_INIT

		QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
		e = PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
				MSLevelEnum::MS2,
				candidateScores.mid(0, topN),
				pythiaParameters->verbosity,
				&msCalibrationReaderRows
				); ree;

		e = msCalibratomatic->buildRTMapper(msCalibrationReaderRows); ree;
		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "ScanTime StdDev Predicted:"<< msCalibratomatic->scanTimeStDev();

		for (CandidateScoresV2 *cs : candidateScores) {
			if (!cs->targetDecoyCandidatePair) {
				qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Null targetDecoyCandidatePair:" << cs->targetDecoyCandidatePair;
			}
		}

		ERR_RETURN
	}

}//namespace
Err MsCalibratomaticSettertronV2::predictScanTimesForCandidateScores(const QVector<CandidateScoresV2*> &candidateScores) const {
//TODO move this to a place were everyone can use it.
	ERR_INIT

	e = ErrorUtils::isNotEmpty(candidateScores); ree;

	constexpr int topN = 1000;
	e = buildRTCurves(
		candidateScores,
		m_pythiaParameters,
		topN,
		m_msCalibratomatic
		); ree;

	e = ErrorUtils::isTrue(m_msCalibratomatic->isInitRT()); ree;

	const float scanTimeMax = m_mzTargetKeyVsMsFramesMS2Pntrs.last()->scanTimeMax();

	for (CandidateScoresV2 *cs : candidateScores) {
		e = predictScanTimesForCandidateScoresParallelLogic(m_msCalibratomatic, cs); ree;
		cs->featuresArray[FTR::ScanTimeRelative] = cs->scanTime / scanTimeMax;
		cs->featuresArray[FTR::ScanTimeRelativeDelta] = (cs->scanTime - cs->scanTimePredicted) / scanTimeMax;
		cs->featuresArray[FTR::ScanTimeRelativeDeltaAbs] = std::abs(cs->scanTime - cs->scanTimePredicted) / scanTimeMax;
	}

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "ScanTimes prediction finished";

	ERR_RETURN
}








