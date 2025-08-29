//
// Created by andrewnichols on 6/29/25.
//

#include "PythiaDDAWorkflow.h"

#include "DeIsotopotron.h"
#include "Ms2IonFraggertronManager.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderPointerAcc.h"
#include "ObjectCSVWriters.h"
#include "ParallelUtils.h"

Err PythiaDDAWorkflow::init(
	const PythiaParameters &parameters,
	const QString& libraryFilePath
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(parameters.isValid()); ree;
	e = ErrorUtils::fileExists(libraryFilePath); ree;

	m_parameters = parameters;
	m_parameters.print();

	e = FragLibReader::getFragLibReaderRows(
		libraryFilePath,
		&m_fragLibReaderRows
		); ree;

	e = m_tdcpManager.init(
		m_parameters,
		&m_fragLibReaderRows
		); ree;

	int idCounter = 0;
	m_tdcpManager.getTargetDecoyCandidatePairPointers(&m_targetDecoyCandidatePairsPntrs);
	for (TargetDecoyCandidatePair *tdcp : m_targetDecoyCandidatePairsPntrs) {
		tdcp->id = idCounter;
		idCounter += 2;
	}

	qDebug()
	<< qPrintable(S_GLOBAL_TIMER.elapsed())
	<< m_targetDecoyCandidatePairsPntrs.size()
	<< "targets from library of" << m_fragLibReaderRows.size()
	<< "after filter.";

	ERR_RETURN
}

Err PythiaDDAWorkflow::processFile(const QString &msDataFilePath) {

	ERR_INIT

	MsReaderPointerAcc msReaderPtr;
	e = msReaderPtr.openFile(msDataFilePath); ree;

	constexpr int msLevel = 2;
	m_scanNumberVsMsScanInfoMS2 = msReaderPtr.ptr->getMsScanInfos(msLevel); ree;
	e = ErrorUtils::isNotEmpty(m_scanNumberVsMsScanInfoMS2); ree;

	e = extractScansParallel(
		m_scanNumberVsMsScanInfoMS2.values().toVector(),
		&msReaderPtr
		); ree;

	e = performFragging(); ree;

	ERR_RETURN
}

namespace {

	QPair<Err, QVector<MsScanPoint>> extractScansParallelLogic(
		const QVector<MsScanInfo*> &scanInfosPntrs,
		MsReaderPointerAcc *msReaderPtr
		) {
		ERR_INIT

		QVector<MsScan> msScans;
		e = msReaderPtr->ptr->extractScanPoints(
			scanInfosPntrs,
			&msScans
			); rree;

		QVector<MsScanPoint> msScanPoints;
		for (const MsScan &scan : msScans) {
			for (int i = 0; i < scan.msScanInfoPntr->pointCount; i++) {
				MsScanPoint msScanPoint;
				msScanPoint.scanInfoPntr = scan.msScanInfoPntr;
				msScanPoint.mzVal = scan.mzVals[i];
				msScanPoint.intensityVal = scan.intensityVals[i];
				msScanPoints.push_back(msScanPoint);
			}

		}

		constexpr float mzMin = 300; //TODO get values from MsReader for these vals
		constexpr float mzMax = 1600; //TODO get values from MsReader for these vals
		const auto terminatorLogic = [mzMin, mzMax](const MsScanPoint &mssp) {
			return mssp.mzVal < mzMin || mssp.mzVal > mzMax;
		};
		const auto terminator = std::remove_if(msScanPoints.begin(), msScanPoints.end(), terminatorLogic);
		msScanPoints.erase(terminator, msScanPoints.end());

		std::sort(
			msScanPoints.begin(),
			msScanPoints.end(),
			[](const MsScanPoint &l, const MsScanPoint &r) {
				if (MathUtils::tSame(l.scanInfoPntr->precursorTargetMz, r.scanInfoPntr->precursorTargetMz)) {
					return l.mzVal < r.mzVal;
				}
				return l.scanInfoPntr->precursorTargetMz < r.scanInfoPntr->precursorTargetMz;
			});

		return {e, msScanPoints};
	}

}//namespace
Err PythiaDDAWorkflow::extractScansParallel(
	const QVector<MsScanInfo*> &scanInfosPntrs,
	MsReaderPointerAcc *msReaderPtr
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPtr->isInit()); ree;

	QVector<QVector<MsScanInfo*>> scanInfosPntrsTranched;
	e = ParallelUtils::trancheVectorForParallelization(
		scanInfosPntrs,
		m_parameters.threadCount * 4,
		&scanInfosPntrsTranched
		); ree;

	const auto binderLogic = std::bind(
		extractScansParallelLogic,
		std::placeholders::_1,
		msReaderPtr
		);

	QFuture<QPair<Err, QVector<MsScanPoint>>> future = QtConcurrent::mapped(
		scanInfosPntrsTranched,
		binderLogic
		);
	future.waitForFinished();

	for (const QPair<Err, QVector<MsScanPoint>> &res : future) {
		e = res.first; ree;
		m_msScanPointsTranched.push_back(res.second);
	}

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsPoints extracted";

	ERR_RETURN
}

namespace {

	QPair<Err, QVector<QPair<float, float>>> buildMzPrecursorTargetRanges(
		const QVector<QVector<MsScanPoint>> &msScanPointsTranched,
		int threadCount
		) {
		ERR_INIT

		e = ErrorUtils::isNotEmpty(msScanPointsTranched); rree;
		e = ErrorUtils::isNotEmpty(msScanPointsTranched.front()); rree;
		e = ErrorUtils::isNotEmpty(msScanPointsTranched.back()); rree;

		QVector<float> precursorMzValsSample;
		for (const QVector<MsScanPoint> &msScanPoints : msScanPointsTranched) {

			constexpr int skipCount = 50;
			for (int i = 0; i <= msScanPoints.size(); i += skipCount) {

				const MsScanPoint &msScanPoint = msScanPoints[i];
				precursorMzValsSample.append(msScanPoint.scanInfoPntr->precursorTargetMz);
			}
		}
		std::sort(precursorMzValsSample.begin(), precursorMzValsSample.end());
		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Sample mz distribution size" << precursorMzValsSample.size();

		QVector<QVector<float>> precursorMzValsSampleTranched;
		e = ParallelUtils::trancheVectorForParallelizationInOrder(
			precursorMzValsSample,
			threadCount,
			0,
			&precursorMzValsSampleTranched
			); rree;

		QVector<QPair<float, float>> precursorMzValsStartStopTranches;
		for (const QVector<float> &v : precursorMzValsSampleTranched) {
			precursorMzValsStartStopTranches.push_back({v.front(), v.back()});
		}

		return {e, precursorMzValsStartStopTranches};
	}

	QPair<Err, QVector<ProcessingGroup>> rearrangeMsScanPointsByMzPrecursorTargetRanges(
		QVector<QVector<MsScanPoint>> &msScanPointsTranched,
		const QVector<QPair<float, float>> &mzPrecursorTargetRanges
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(msScanPointsTranched); rree;
		e = ErrorUtils::isNotEmpty(mzPrecursorTargetRanges); rree;

		QVector<ProcessingGroup> processingGroups(mzPrecursorTargetRanges.size());
		for (int i = 0; i < mzPrecursorTargetRanges.size(); i++) {
			processingGroups[i].mzPrecursorRangeMinMax = mzPrecursorTargetRanges[i];
		}

		for (QVector<MsScanPoint> &msScanPoints : msScanPointsTranched) {

			const bool isSortedMzAsc = std::is_sorted(
				msScanPoints.begin(),
				msScanPoints.end(),
				[](const MsScanPoint &l, const MsScanPoint &r) {
					return l.scanInfoPntr->precursorTargetMz < r.scanInfoPntr->precursorTargetMz;
				});
			e = ErrorUtils::isTrue(isSortedMzAsc); rree;

			for (ProcessingGroup &pr : processingGroups) {
				const float precursorMzMin = pr.mzPrecursorRangeMinMax.first;
				const float precursorMzMax = pr.mzPrecursorRangeMinMax.second;

				const int precursorMzMinLowerBoundIndex = std::upper_bound(
					msScanPoints.begin(),
					msScanPoints.end(),
					precursorMzMin,
					[](float val, const MsScanPoint &stct) {return stct.scanInfoPntr->precursorTargetMz > val;}
					) - msScanPoints.begin();

				const int precursorMzMaxUpperBoundIndex = std::lower_bound(
					msScanPoints.begin(),
					msScanPoints.end(),
					precursorMzMax,
					[](const MsScanPoint &stct, float val) {return val > stct.scanInfoPntr->precursorTargetMz;}
					) - msScanPoints.begin() - 1;

				const int size = precursorMzMaxUpperBoundIndex - precursorMzMinLowerBoundIndex;
				if (size < 1) {
					continue;
				}

				pr.msScanPoints.reserve(size);
				for(int i = precursorMzMinLowerBoundIndex; i <= precursorMzMaxUpperBoundIndex; i++) {
					pr.msScanPoints.push_back(&msScanPoints[i]);
				}

				e = ErrorUtils::isAboveThreshold(
					pr.msScanPoints.front()->scanInfoPntr->precursorTargetMz,
					precursorMzMin,
					ErrorUtilsParam::IncludeThreshold
					); rree;

				e = ErrorUtils::isBelowThreshold(
					pr.msScanPoints.back()->scanInfoPntr->precursorTargetMz,
					precursorMzMax,
					ErrorUtilsParam::IncludeThreshold
					); rree;

			}
		}

		return {e, processingGroups};
	}

	std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> buildMs2Ions(
		TargetDecoyCandidatePair *tdcp
		) {

		ERR_INIT

		constexpr float mzMin = 300; //TODO get values from MsReader for these vals
		constexpr float mzMax = 1600; //TODO get values from MsReader for these vals
		const QVector<MS2Ion> &ms2IonsTarget = tdcp->ms2IonsTarget(mzMin, mzMax);
		const QVector<MS2Ion> &ms2IonsDecoy = tdcp->ms2IonsDecoy(ms2IonsTarget);

		std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> res(
			tdcp, ms2IonsTarget, ms2IonsDecoy
			);

		return res;
	}


	QPair<Err, QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>>>
										buildLibraryMS2IonsTrancheLogic(const QVector<TargetDecoyCandidatePair*> &tdcps) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(tdcps); rree;

		QFuture<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> futureLib = QtConcurrent::mapped(
			tdcps,
			buildMs2Ions
			);
		futureLib.waitForFinished();

		QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> tdcpMs2IonsTargetDecoyTpl;
		tdcpMs2IonsTargetDecoyTpl.reserve(tdcps.size());
		for (const std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> &r : futureLib) {
			tdcpMs2IonsTargetDecoyTpl.push_back(r);
		}

		return {e, tdcpMs2IonsTargetDecoyTpl};
	}

	Err buildM22IonLibraries(
		QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2IonsLibraryTranche,
		QVector<MS2IonLibrary> *ms2IonLibraries
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(ms2IonsLibraryTranche); ree;
		ms2IonLibraries->clear();

		constexpr int targetDecoyDoubler = 2;
		constexpr int maxLibraryIonCount = 16;
		ms2IonLibraries->reserve(ms2IonsLibraryTranche.size() * targetDecoyDoubler * maxLibraryIonCount);
		for (std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> &tpl : ms2IonsLibraryTranche) {
			TargetDecoyCandidatePair *tdcp = std::get<0>(tpl);

			MS2IonsTarget &ms2IonsTarget = std::get<1>(tpl);
			MS2IonsDecoy &ms2IonsDecoy = std::get<2>(tpl);

			for(int i = 0; i < std::min(maxLibraryIonCount, ms2IonsTarget.size()); ++i) {
				MS2Ion *ms2Ion = &ms2IonsTarget[i];
				MS2IonLibrary ms2IonLibrary;
				ms2IonLibrary.ms2IonPntr = ms2Ion;
				ms2IonLibrary.targeDecoyCandidatePairPntr = tdcp;
				ms2IonLibrary.isDecoy = false;
				ms2IonLibraries->push_back(ms2IonLibrary);
			}

			for(int i = 0; i < std::min(maxLibraryIonCount, ms2IonsDecoy.size()); ++i) {
				MS2Ion *ms2Ion = &ms2IonsDecoy[i];
				MS2IonLibrary ms2IonLibrary;
				ms2IonLibrary.ms2IonPntr = ms2Ion;
				ms2IonLibrary.targeDecoyCandidatePairPntr = tdcp;
				ms2IonLibrary.isDecoy = true;
				ms2IonLibraries->push_back(ms2IonLibrary);
			}
		}

		ERR_RETURN;
	}

	Err checkProcessingGroupRangesAreValid(const QVector<ProcessingGroup> &processingGroups) {
		ERR_INIT

		for (const ProcessingGroup &pg : processingGroups) {

			e = ErrorUtils::isAboveThreshold(
				pg.msScanPoints.front()->scanInfoPntr->precursorTargetMz,
				pg.mzPrecursorRangeMinMax.first,
				ErrorUtilsParam::IncludeThreshold
				); ree;

			e = ErrorUtils::isBelowThreshold(
				pg.msScanPoints.back()->scanInfoPntr->precursorTargetMz,
				pg.mzPrecursorRangeMinMax.second,
				ErrorUtilsParam::IncludeThreshold
				); ree;

			e = ErrorUtils::isAboveThreshold(
				pg.mzPrecursorRangeMinMax.first,
				pg.ms2IonsLibrary.front()->targeDecoyCandidatePairPntr->mz(pg.ms2IonsLibrary.front()->isDecoy),
				ErrorUtilsParam::IncludeThreshold
				); ree;

			e = ErrorUtils::isBelowThreshold(
				pg.mzPrecursorRangeMinMax.second,
				pg.ms2IonsLibrary.back()->targeDecoyCandidatePairPntr->mz(pg.ms2IonsLibrary.back()->isDecoy),
				ErrorUtilsParam::IncludeThreshold
				); ree;

		}

		ERR_RETURN;
	}

	QPair<Err, int> performFraggingLogic(
		const QVector<MsScan*> &msScansPntrs,
		const QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2Ions,
		const PythiaParameters &parameters
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(ms2Ions); rree;
		e = ErrorUtils::isNotEmpty(msScansPntrs); rree;

		Ms2IonFraggertronManager fragger;
		e = fragger.init(ms2Ions); rree;

		for (const MsScan *msScan : msScansPntrs) {

			const MzVals &mzVals = msScan->mzVals ;
			const IntensityVals &intensityVals = msScan->intensityVals ;
			const ScanNumber scanNumber = msScan->msScanInfoPntr->scanNumber;
			const int pointCount = msScan->msScanInfoPntr->pointCount;

			float cutoffIntensity = -1.0;
			constexpr int topNScanPoints = 500;
			if (pointCount > topNScanPoints) {
				QVector<float> intensityValsCopy(intensityVals.data(), intensityVals.data() + pointCount);
				std::sort(intensityValsCopy.rbegin(), intensityValsCopy.rend());
				cutoffIntensity = intensityValsCopy[topNScanPoints - 1];
			}
			for (int i = 0; i < pointCount; ++i) {

				const float intensity = msScan->intensityVals[i];
				if (intensity < cutoffIntensity) {
					continue;
				}

				constexpr float iRTMin = -10000;
				constexpr float iRTMax = 10000;

				const float mzVal = msScan->mzVals[i];
				const float mzTol = MathUtils::calculatePPM(mzVal, static_cast<float>(parameters.ms2ExtractionWidthPPM));
				const float mzMin = mzVal - mzTol;
				const float mzMax = mzVal + mzTol;

				QVector<MS2Frag*> tdPeptideFrags;
				e = fragger.extractMs2Points(
					mzMin,
					mzMax,
					iRTMin,
					iRTMax,
					&tdPeptideFrags
					); rree;
			}
		}

		return {e, -1};
	}

}//namespace
Err PythiaDDAWorkflow::performFragging() {

	ERR_INIT

	QElapsedTimer et;
	et.start();

	QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsPntrsTranched;
	constexpr int libTrancheSize = 4;
	e = ParallelUtils::trancheVectorForParallelization(
		m_targetDecoyCandidatePairsPntrs,
		libTrancheSize,
		&targetDecoyCandidatePairsPntrsTranched
		); ree;


	const QPair<Err, QVector<QPair<float, float>>> mzPrecursorStartVsStopResult = buildMzPrecursorTargetRanges(
		m_msScanPointsTranched,
		m_parameters.threadCount
		);
	e = mzPrecursorStartVsStopResult.first; ree;
	const QVector<QPair<float, float>> &mzPrecursorStartVsStop = mzPrecursorStartVsStopResult.second;

	QPair<Err, QVector<ProcessingGroup>> processingGroupsResult = rearrangeMsScanPointsByMzPrecursorTargetRanges(
		m_msScanPointsTranched,
		mzPrecursorStartVsStop
		);
	e = processingGroupsResult.first; ree;

	QVector<ProcessingGroup> &processingGroups = processingGroupsResult.second;

	for (const QVector<TargetDecoyCandidatePair*> &tdcps : targetDecoyCandidatePairsPntrsTranched) {

		QPair<Err, QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>>> ms2IonsLibraryTrancheResult
																								= buildLibraryMS2IonsTrancheLogic(tdcps);
		e = ms2IonsLibraryTrancheResult.first; ree;
		QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2IonsLibraryTranche
																								= ms2IonsLibraryTrancheResult.second;

		QVector<MS2IonLibrary> ms2IonLibraries;
		e = buildM22IonLibraries(ms2IonsLibraryTranche, &ms2IonLibraries); ree;
		std::sort(
			ms2IonLibraries.begin(),
			ms2IonLibraries.end(),
			[](const MS2IonLibrary &l, const MS2IonLibrary &r) {
				return l.targeDecoyCandidatePairPntr->mz(l.isDecoy) < r.targeDecoyCandidatePairPntr->mz(r.isDecoy);
			});

		for (ProcessingGroup &pg : processingGroups) {

			constexpr float buffer = 2.f;
			const float precursorMzMin = pg.mzPrecursorRangeMinMax.first - m_parameters.precursorExtractionWindowThomsons - buffer;
			const float precursorMzMax = pg.mzPrecursorRangeMinMax.second + m_parameters.precursorExtractionWindowThomsons + buffer;

			const int precursorMzMinLowerBoundIndex = std::upper_bound(
				ms2IonLibraries.begin(),
				ms2IonLibraries.end(),
				precursorMzMin,
				[](float val, const MS2IonLibrary &stct) {return stct.targeDecoyCandidatePairPntr->mz(stct.isDecoy) > val;}
				) - ms2IonLibraries.begin();

			const int precursorMzMaxUpperBoundIndex = std::lower_bound(
				ms2IonLibraries.begin(),
				ms2IonLibraries.end(),
				precursorMzMax,
				[](const MS2IonLibrary &stct, float val) {return val > stct.targeDecoyCandidatePairPntr->mz(stct.isDecoy);}
				) - ms2IonLibraries.begin() - 1;

			const int size = precursorMzMaxUpperBoundIndex - precursorMzMinLowerBoundIndex;
			if (size < 1) {
				continue;
			}

			pg.ms2IonsLibrary.reserve(size);
			for (int i = precursorMzMinLowerBoundIndex; i <= precursorMzMaxUpperBoundIndex; ++i) {
				pg.ms2IonsLibrary.push_back(&ms2IonLibraries[i]);
			}
		}

		e = checkProcessingGroupRangesAreValid(processingGroups); ree;

// #define FRAG_PARALLEL
// #ifdef FRAG_PARALLEL
// 		const auto binderLogic = std::bind(
// 			performFraggingLogic,
// 			std::placeholders::_1,
// 			ms2Ions,
// 			m_parameters
// 			);
//
// 		QFuture<QPair<Err, int>> futureScans = QtConcurrent::mapped(
// 			msScansPntrsTranched,
// 			binderLogic
// 			);
// 		futureScans.waitForFinished();
// #else
// 		for (const QVector<MsScan*> &msScans : msScansPntrsTranched) {
// 			const QPair<Err, int> res = performFraggingLogic(msScans, ms2Ions, m_parameters);
// 			e = res.first; ree;
// 		}
// #endif
//

	}




	ERR_RETURN
}




