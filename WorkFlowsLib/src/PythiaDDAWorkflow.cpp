//
// Created by andrewnichols on 6/29/25.
//

#include "PythiaDDAWorkflow.h"

#include <boost/geometry/index/rtree.hpp>

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
	m_parameters.ms2ExtractionWidthPPM = 7;
	// m_parameters.threadCount = 64; //TODO use higher threadcount to load balans
	m_parameters.print();

	e = FragLibReader::getFragLibReaderRows(
		libraryFilePath,
		&m_fragLibReaderRows
		); ree;

	e = m_tdcpManager.init(
		m_parameters,
		&m_fragLibReaderRows
		); ree;

	m_tdcpManager.getTargetDecoyCandidatePairPointers(&m_targetDecoyCandidatePairsPntrs);

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
	e = buildMsScanPointPntrs(&msReaderPtr); ree;
	e = performFragging(); ree;

	ERR_RETURN
}

Err PythiaDDAWorkflow::buildMsScanPointPntrs(MsReaderPointerAcc *msReaderPtr) {
	ERR_INIT

	constexpr int msLevel = 2;
	m_scanNumberVsMsScanInfoMS2 = msReaderPtr->ptr->getMsScanInfos(msLevel); ree;
	e = ErrorUtils::isNotEmpty(m_scanNumberVsMsScanInfoMS2); ree;

	QVector<MsScanInfo*> msScanInfosMS2 = m_scanNumberVsMsScanInfoMS2.values().toVector();
	std::sort(
		msScanInfosMS2.begin(),
		msScanInfosMS2.end(),
		[](const MsScanInfo *l, const MsScanInfo *r) {
			return l->precursorTargetMz < r->precursorTargetMz;
		});

	e = extractScansParallel(
		msScanInfosMS2,
		msReaderPtr
		); ree;

	const int ms2PointCount = std::accumulate(
		m_msScanPointsTranched.begin(),
		m_msScanPointsTranched.end(),
		0,
		[](int sum, const QVector<MsScanPoint> &v) {return sum + v.size();}
		);

	const int trancheSize = ms2PointCount / m_parameters.threadCount;
	int trancheRemainder = ms2PointCount % m_parameters.threadCount;
	m_msScanPointsPntrsTranched.reserve(m_parameters.threadCount);
	QVector<int> trancheSizes(m_parameters.threadCount,trancheSize);

	int counter = 0;
	while (trancheRemainder > 0) {
		trancheRemainder--;
		trancheSizes[counter]++;
		counter = counter == m_parameters.threadCount ? 0 : counter + 1;
	}

	QVector<MsScanPoint*> currentMsScanPointsPntrsTranche;
	currentMsScanPointsPntrsTranche.reserve(trancheSizes.front());

	for (QVector<MsScanPoint> &v : m_msScanPointsTranched) {

		for (MsScanPoint &p : v) {

			currentMsScanPointsPntrsTranche.push_back(&p);

			if (currentMsScanPointsPntrsTranche.size() >= trancheSizes[m_msScanPointsPntrsTranched.size()]) {

				m_msScanPointsPntrsTranched.push_back(currentMsScanPointsPntrsTranche);

				const auto minMaxPrecursorMz = std::minmax_element(
					currentMsScanPointsPntrsTranche.begin(),
					currentMsScanPointsPntrsTranche.end(),// 	currentMsScanPointsPntrsTranche.end(),
					[](const MsScanPoint *l, const MsScanPoint *r) {
						return l->scanInfoPntr->precursorTargetMz < r->scanInfoPntr->precursorTargetMz;
					});

				const float minPrecursorMz = (*minMaxPrecursorMz.first)->scanInfoPntr->precursorTargetMz;
				const float maxPrecursorMz = (*minMaxPrecursorMz.second)->scanInfoPntr->precursorTargetMz;

				e = ErrorUtils::isTrue(MathUtils::tSame(
						minPrecursorMz,
						currentMsScanPointsPntrsTranche.front()->scanInfoPntr->precursorTargetMz,
						0.5
						));

				e = ErrorUtils::isTrue(MathUtils::tSame(
						maxPrecursorMz,
						currentMsScanPointsPntrsTranche.back()->scanInfoPntr->precursorTargetMz,
						0.5
						));

				m_mzPrecursorStartVsStopResult.push_back({
					minPrecursorMz,
					maxPrecursorMz
				});
				currentMsScanPointsPntrsTranche.clear();
			}
		}
	}

	if (!currentMsScanPointsPntrsTranche.isEmpty()) {
		m_msScanPointsPntrsTranched.back().append(currentMsScanPointsPntrsTranche);
	}
	e = ErrorUtils::isEqual(
		m_parameters.threadCount,
		m_msScanPointsPntrsTranched.size()
		); ree;

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
	e = ParallelUtils::trancheVectorForParallelizationInOrder(
		scanInfosPntrs,
		m_parameters.threadCount * 4,
		0,
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

	Err buildMS2IonLibraries(
		QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2IonsLibraryTranche,
		QVector<MS2IonLibrary> *ms2IonLibraries
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(ms2IonsLibraryTranche); ree;
		ms2IonLibraries->clear();

		constexpr int targetDecoyDoubler = 2;
		constexpr int maxLibraryIonCount = 10;
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

			// e = ErrorUtils::isNotEmpty(pg.ms2IonsLibrary); eee_absorb;
			e = ErrorUtils::isNotEmpty(*pg.msScanPointsPntr); ree;

			e = ErrorUtils::isAboveThreshold(
				pg.msScanPointsPntr->front()->scanInfoPntr->precursorTargetMz,
				pg.mzPrecursorRangeMinMax.first,
				ErrorUtilsParam::IncludeThreshold
				); ree;

			e = ErrorUtils::isBelowThreshold(
				pg.msScanPointsPntr->back()->scanInfoPntr->precursorTargetMz,
				pg.mzPrecursorRangeMinMax.second,
				ErrorUtilsParam::IncludeThreshold
				); ree;

			// e = ErrorUtils::isAboveThreshold(
			// 	pg.mzPrecursorRangeMinMax.first,
			// 	pg.ms2IonsLibrary.front()->targeDecoyCandidatePairPntr->mz(pg.ms2IonsLibrary.front()->isDecoy),
			// 	ErrorUtilsParam::IncludeThreshold
			// 	); ree;
			//
			// e = ErrorUtils::isBelowThreshold(
			// 	pg.mzPrecursorRangeMinMax.second,
			// 	pg.ms2IonsLibrary.back()->targeDecoyCandidatePairPntr->mz(pg.ms2IonsLibrary.back()->isDecoy),
			// 	ErrorUtilsParam::IncludeThreshold
			// 	); ree;

		}

		ERR_RETURN;
	}

	void ms2IonLibrarySortLogic(QVector<MS2IonLibrary> &ms2IonLibraries) {
		std::sort(
			ms2IonLibraries.begin(),
			ms2IonLibraries.end(),
			[](const MS2IonLibrary &l, const MS2IonLibrary &r) {
				return l.targeDecoyCandidatePairPntr->mz(l.isDecoy) < r.targeDecoyCandidatePairPntr->mz(r.isDecoy);
			});
	}

	const auto sortLogicMzAsc = [](const MS2IonLibrary &l, const MS2IonLibrary &r) {
		return l.targeDecoyCandidatePairPntr->mz(l.isDecoy) < r.targeDecoyCandidatePairPntr->mz(r.isDecoy);
	};

	void parallelSortMS2IonLibrariesMzAsc(
		int threadCount,
		QVector<MS2IonLibrary> *ms2IonLibraries
		) {

		auto [minIt, maxIt] = std::minmax_element(
			ms2IonLibraries->begin(),
			ms2IonLibraries->end(),
			sortLogicMzAsc
			);

		const MS2IonLibrary min = *minIt;
		const MS2IonLibrary max = *maxIt;

		const int numTranches = threadCount;
		const int range = max.targeDecoyCandidatePairPntr->mz(max.isDecoy)
						- min.targeDecoyCandidatePairPntr->mz(min.isDecoy) + 1;

		const int trancheSize = range / numTranches;

		QVector<QVector<MS2IonLibrary>> tranches(numTranches);
		for (const MS2IonLibrary &num : *ms2IonLibraries) {
			int trancheIndex = (
				static_cast<int>(num.targeDecoyCandidatePairPntr->mz(num.isDecoy) - min.targeDecoyCandidatePairPntr->mz(min.isDecoy))
				) / trancheSize;
			trancheIndex = std::min(trancheIndex, numTranches - 1);  // Handle edge case
			tranches[trancheIndex].push_back(num);
		}

		QFuture<void> sortedTranchesFuture = QtConcurrent::map(tranches, ms2IonLibrarySortLogic);
		sortedTranchesFuture.waitForFinished();

		ms2IonLibraries->clear();
		for (const QVector<MS2IonLibrary> &l : tranches) {
			ms2IonLibraries->append(l);
		}
	}

	Err addMs2IonsLibraryToProcessingGroups(
		QVector<MS2IonLibrary> &ms2IonLibraries,
		float precursorExtractionWindowThomsons,
		int threadCount,
		QVector<ProcessingGroup> *processingGroups
		) {

		ERR_INIT

		parallelSortMS2IonLibrariesMzAsc(
			threadCount,
			&ms2IonLibraries
			);

		const bool isSorted = std::is_sorted(
			ms2IonLibraries.begin(),
			ms2IonLibraries.end(),
			sortLogicMzAsc
			);
		e = ErrorUtils::isTrue(isSorted); ree;

		for (ProcessingGroup &pg : *processingGroups) {

			e = ErrorUtils::isTrue(pg.ms2IonsLibrary.isEmpty()); ree;

			constexpr float buffer = 2.f;
			const float precursorMzMin = pg.mzPrecursorRangeMinMax.first - precursorExtractionWindowThomsons - buffer;
			const float precursorMzMax = pg.mzPrecursorRangeMinMax.second + precursorExtractionWindowThomsons + buffer;

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

		ERR_RETURN
	}

	QPair<Err, QVector<TallyResultTuple>> collateScanNumberVsOccurrencesTargetDecoyCandidatePairs(
		QHash<TargetDecoyCandidatePair*, QVector<IonSearchResult>> &input
		) {

		ERR_INIT

		struct T {
			Occurrence occurrence = 0;
			QVector<int> ranks;
		};

		QVector<TallyResultTuple> result;

		QHash<ScanNumber, T> occurrencesTarget;
		QHash<ScanNumber, T> occurrencesDecoy;

		for (auto it = input.begin(); it != input.end(); it++) {

			const QVector<IonSearchResult> &isrs = it.value();
			TargetDecoyCandidatePair *targetDecoyCandidatePair = it.key();
			occurrencesTarget.clear();
			occurrencesDecoy.clear();

			for (const IonSearchResult &isr : isrs) {
				const ScanNumber scanNumber = isr.msScanPointPntr->scanInfoPntr->scanNumber;
				if (isr.ms2IonLibraryPntr->isDecoy) {
					T &t = occurrencesDecoy[scanNumber];
					t.occurrence++;
					t.ranks.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->rank);
					continue;
				}

				T &t = occurrencesTarget[scanNumber];
				t.occurrence++;
				t.ranks.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->rank);
			}

			if (occurrencesTarget.isEmpty() && occurrencesDecoy.isEmpty()) {
				continue;
			}

			constexpr int occurenceCountMin = 1;

			QVector<TallyResultTarget> tallyResultsTarget;
			if (!occurrencesTarget.isEmpty()) {
				for (auto itt = occurrencesDecoy.begin(); itt != occurrencesDecoy.end(); itt++) {
					const T &trt =  occurrencesTarget[itt.key()];
					if (trt.occurrence <= occurenceCountMin) {
						continue;
					}

					TallyResultTarget tallyResult;
					tallyResult.scanNumber = itt.key();
					tallyResult.occurrence = trt.occurrence;
					tallyResult.ranks = trt.ranks;

					tallyResultsTarget.push_back(tallyResult);
				}

				if (!tallyResultsTarget.isEmpty()) {
					std::sort(
						tallyResultsTarget.rbegin(),
						tallyResultsTarget.rend(),
						[](const TallyResult &l, const TallyResult &r) {return l.occurrence < r.occurrence;}
					);
				}
			}

			QVector<TallyResultTarget> tallyResultsDecoy;
			if (!occurrencesDecoy.isEmpty()) {

				for (auto itt = occurrencesDecoy.begin(); itt != occurrencesDecoy.end(); itt++) {

					const T &trt =  occurrencesDecoy[itt.key()];
					if (trt.occurrence <= occurenceCountMin) {
						continue;
					}

					TallyResultTarget tallyResult;
					tallyResult.scanNumber = itt.key();
					tallyResult.occurrence = trt.occurrence;
					tallyResult.ranks = trt.ranks;

					tallyResultsDecoy.push_back(tallyResult);
				}

				if (!tallyResultsDecoy.isEmpty()) {
					std::sort(
						tallyResultsDecoy.rbegin(),
						tallyResultsDecoy.rend(),
						[](const TallyResult &l, const TallyResult &r) {return l.occurrence < r.occurrence;}
					);
				}
			}

			if (tallyResultsTarget.isEmpty() && tallyResultsDecoy.isEmpty()) {
				continue;
			}

			result.push_back({targetDecoyCandidatePair, tallyResultsTarget, tallyResultsDecoy});
		}

		return {e, result};
	}

	QPair<Err, QVector<TallyResultTuple>> performFraggingLogic(
		const ProcessingGroup &pgs,
		const PythiaParameters &parameters
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(*pgs.msScanPointsPntr); rree;

		if (pgs.ms2IonsLibrary.isEmpty()) {
			qDebug()
				<< "No ms2IonsLibrary found for precursor range"
				<< pgs.mzPrecursorRangeMinMax.first << "-"
				<< pgs.mzPrecursorRangeMinMax.second
			;
			return {e, {}};
		}

		Ms2IonFraggertronManager fragger;
		e = fragger.init(pgs.ms2IonsLibrary); rree;

		QHash<TargetDecoyCandidatePair*, QVector<IonSearchResult>> ionSearchResults;

		for (MsScanPoint *mssp : *pgs.msScanPointsPntr) {
			const float precursorMzValLower = mssp->scanInfoPntr->precursorTargetMz
											- mssp->scanInfoPntr->isoWindowLower
											- parameters.precursorExtractionWindowThomsons;

			const float precursorMzValUpper = mssp->scanInfoPntr->precursorTargetMz
											+ mssp->scanInfoPntr->isoWindowUpper
											+ parameters.precursorExtractionWindowThomsons;

			const float mzVal = mssp->mzVal;
			const float mzTol = MathUtils::calculatePPM(
				mzVal,
				static_cast<float>(parameters.ms2ExtractionWidthPPM)
				);

			const float mzMin = mzVal - mzTol;
			const float mzMax = mzVal + mzTol;

			QVector<MS2IonLibrary*> tdPeptideFrags;
			e = fragger.extractMs2Points(
				precursorMzValLower,
				precursorMzValUpper,
				mzMin,
				mzMax,
				&tdPeptideFrags
				); rree;

			for (MS2IonLibrary *msil : tdPeptideFrags) {
				IonSearchResult isr;
				isr.ms2IonLibraryPntr = msil;
				isr.msScanPointPntr = mssp;

				ionSearchResults[msil->targeDecoyCandidatePairPntr].push_back(isr);
			}
		}

// #define TR_SHT
#ifdef TR_SHT
		for (auto it = ionSearchResults.begin(); it != ionSearchResults.end(); it++) {
			const TargetDecoyCandidatePair *tcp = it.key();
			QVector<IonSearchResult> &r = it.value();

			if (tcp->peptideStringWithMods() != "LTVTSLQETGLK" || tcp->charge() != 2) {
				continue;
			}

			std::sort(r.begin(), r.end(), [](const IonSearchResult &l, const IonSearchResult &r) {
				return l.msScanPointPntr->scanInfoPntr->scanNumber < r.msScanPointPntr->scanInfoPntr->scanNumber;
			});

			for (auto x : r) {
				qDebug()
					<< x.msScanPointPntr->scanInfoPntr->scanNumber
					<< tcp->peptideStringWithMods()
					// << tcp->mz(x.ms2IonLibraryPntr->isDecoy)
					// << x.msScanPointPntr->scanInfoPntr->precursorTargetMz
					<< x.msScanPointPntr->mzVal
					<< x.ms2IonLibraryPntr->ms2IonPntr->mz
					<< x.ms2IonLibraryPntr->isDecoy
					<< "SLDFJSDL"
				;
			}
		}
#endif

		QPair<Err, QVector<TallyResultTuple>> result
				= collateScanNumberVsOccurrencesTargetDecoyCandidatePairs(ionSearchResults); rree;

		return result;
	}

	int findLongestContinuousSequenceSize(const QSet<int>& set) {
		QSet<int> hashSet = set;

		int maxLength = 0;

		for (int num : set) {
			if (!hashSet.contains(num - 1)) {
				int currentLength = 1;
				int currentNum = num;

				while (hashSet.contains(currentNum + 1)) {
					currentNum++;
					currentLength++;
				}

				maxLength = std::max(maxLength, currentLength);
			}
		}

		return maxLength;
	}

	Err initProcessingGroups(
		const QVector<QPair<float, float>> &mzPrecursorStartVsStopResult,
		QVector<QVector<MsScanPoint*>> *msScanPointsPntrsTranched,
		QVector<ProcessingGroup> *processingGroups
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(*msScanPointsPntrsTranched); ree;
		e = ErrorUtils::isEqual(mzPrecursorStartVsStopResult.size(), msScanPointsPntrsTranched->size()); ree;

		processingGroups->resize(msScanPointsPntrsTranched->size());

		for (int i = 0; i < msScanPointsPntrsTranched->size(); i++) {
			ProcessingGroup processingGroup;
			processingGroup.msScanPointsPntr = &(*msScanPointsPntrsTranched)[i];
			processingGroup.mzPrecursorRangeMinMax = mzPrecursorStartVsStopResult[i];
			(*processingGroups)[i] = processingGroup;
		}

		ERR_RETURN
	}

	Err buildTargetDecoyCandidatePairsPntrsTranched(
		const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairsPntrs,
		QVector<QVector<TargetDecoyCandidatePair*>> *targetDecoyCandidatePairsPntrsTranched
		) {
		ERR_INIT

		e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairsPntrs); ree;
		targetDecoyCandidatePairsPntrsTranched->clear();

		constexpr int batchSize = 1.5e5; //TODO make this settable in params
		const int libTrancheSize = std::max(targetDecoyCandidatePairsPntrs.size() / batchSize, 1);
		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Processing library in" << libTrancheSize << "tranches";

		e = ParallelUtils::trancheVectorForParallelization(
			targetDecoyCandidatePairsPntrs,
			libTrancheSize,
			targetDecoyCandidatePairsPntrsTranched
			); ree;

		ERR_RETURN
	}

}//namespace
Err PythiaDDAWorkflow::performFragging() {

	ERR_INIT

	QElapsedTimer et;
	et.start();

	e = ErrorUtils::isNotEmpty(m_msScanPointsPntrsTranched); ree;
	e = ErrorUtils::isEqual(m_msScanPointsPntrsTranched.size(), m_mzPrecursorStartVsStopResult.size()); ree;

	QVector<ProcessingGroup> processingGroups;
	e = initProcessingGroups(
		m_mzPrecursorStartVsStopResult,
		&m_msScanPointsPntrsTranched,
		&processingGroups
		);

	QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsPntrsTranched;
	e = buildTargetDecoyCandidatePairsPntrsTranched(
		m_targetDecoyCandidatePairsPntrs,
		&targetDecoyCandidatePairsPntrsTranched
		); ree;

	QVector<TallyResultTuple> tallyResultsFinal;
	for (const QVector<TargetDecoyCandidatePair*> &tdcps : targetDecoyCandidatePairsPntrsTranched) {
		const QPair<Err, QVector<TallyResultTuple>> result = processTargetDecoyCandidatePairsPntrsTranch(tdcps, processingGroups); ree;
		e = result.first; ree;
		tallyResultsFinal.append(result.second);
	}

	struct T {
		TargetDecoyCandidatePair *targetDecoyCandidatePair = nullptr;
		bool isDecoy = false;
		Occurrence occurrence;
		float rankMean;
		int rankBest;
	};
	QVector<T> ts;

	for (const TallyResultTuple &tpl : tallyResultsFinal) {
		T t;
		if (std::get<1>(tpl).front().occurrence > 0) {
			t.targetDecoyCandidatePair = std::get<0>(tpl);
			t.occurrence = std::get<1>(tpl).front().occurrence;
			t.rankMean = MathUtils::mean(std::get<1>(tpl).front().ranks);
			t.rankBest = *std::min_element(std::get<1>(tpl).front().ranks.begin(), std::get<1>(tpl).front().ranks.end());
			t.isDecoy = false;
			ts.push_back(t);
		}

		T td;
		if (std::get<2>(tpl).front().occurrence > 0) {
			td.targetDecoyCandidatePair = std::get<0>(tpl);
			td.occurrence = std::get<2>(tpl).front().occurrence;
			td.rankMean = MathUtils::mean(std::get<2>(tpl).front().ranks);
			td.rankBest = *std::min_element(std::get<2>(tpl).front().ranks.begin(), std::get<2>(tpl).front().ranks.end());
			td.isDecoy = true;
			ts.push_back(td);
		}
	}

	std::sort(
		ts.begin(),
		ts.end(),
		[](const T &l, const T &r) {
			if (l.occurrence == r.occurrence) {
				return l.rankMean < r.rankMean;
			}

			return l.occurrence > r.occurrence;
		}
		);

	int counter = 0;
	int decoy = 0;
	for (const T &t : ts) {

		if (t.isDecoy) {
			decoy++;
		}

		const float q = static_cast<float>(decoy) / ++counter;

		std::cout
			<< counter << " fdsafda "
			<< q << " "
			<< t.targetDecoyCandidatePair->peptideStringWithMods().toStdString() << " "
			<< t.occurrence  << " "
			<< t.rankMean  << " "
			<< t.rankBest  << " "
			<< t.isDecoy
			<< std::endl;

		if (q > 0.01 || t.occurrence < 10) break;

		//mean mz thomsons found, higher is more specific than lower
	}

	ERR_RETURN
}

QPair<Err, QVector<TallyResultTuple>> PythiaDDAWorkflow::processTargetDecoyCandidatePairsPntrsTranch(
	const QVector<TargetDecoyCandidatePair*> &tdcps,
	const QVector<ProcessingGroup> &_processingGroups
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(tdcps); rree;
	e = ErrorUtils::isNotEmpty(_processingGroups); rree;

	QVector<ProcessingGroup> processingGroups = _processingGroups;

	QPair<Err, QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>>> ms2IonsLibraryTrancheResult
																		= buildLibraryMS2IonsTrancheLogic(tdcps);
	e = ms2IonsLibraryTrancheResult.first; rree;
	QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2IonsLibraryTranche
																							= ms2IonsLibraryTrancheResult.second;

	QVector<MS2IonLibrary> ms2IonLibraries;
	e = buildMS2IonLibraries(
		ms2IonsLibraryTranche,
		&ms2IonLibraries
		); rree;

	e = addMs2IonsLibraryToProcessingGroups(
		ms2IonLibraries,
		m_parameters.precursorExtractionWindowThomsons,
		m_parameters.threadCount,
		&processingGroups
		); rree;

	e = checkProcessingGroupRangesAreValid(processingGroups); rree;

	QVector<TallyResultTuple> tallyResultsFinal;

#define FRAG_PARALLEL
#ifdef FRAG_PARALLEL
		const auto binderLogic = std::bind(
			performFraggingLogic,
			std::placeholders::_1,
			m_parameters
			);

		QFuture<QPair<Err, QVector<TallyResultTuple>>> future = QtConcurrent::mapped(
			processingGroups,
			binderLogic
			);
		future.waitForFinished();

		for (const QPair<Err, QVector<TallyResultTuple>> &result : future) {
			e = result.first; rree;
			const QVector<TallyResultTuple> &tallyResultTuples = result.second;
			tallyResultsFinal.append(tallyResultTuples);
		}
		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PSMs Found:" << tallyResultsFinal.size();
#else
		for (const ProcessingGroup &pgs : processingGroups) {
			const QPair<Err, QVector<TallyResultTuple>> res = performFraggingLogic(pgs, m_parameters);
			e = res.first; ree;
		}
#endif

	return {e, tallyResultsFinal};
}




