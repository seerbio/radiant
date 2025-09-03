//
// Created by anichols on 6/7/24.
//

#include "Ms2IonFraggertronManager.h"

#include "MsCalibratomatic.h"
#include "TargetDecoyCandidatePair.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>


namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

class Q_DECL_HIDDEN Ms2IonFraggertronManager::Private
{

    struct RTreeMS2IonPoint {
        MS2Ion *ms2Ion = nullptr;
        CandidateScores *candidateScores = nullptr;
    };

    using rTreeScanNumber = float;
    using rTreeIntensity = float;
    using rTreeCoor = bg::model::point<float, 2, bg::cs::cartesian>;
    using rTreeSearchBox = bg::model::box<rTreeCoor>;
    using rTreePoint = std::pair<rTreeCoor, MS2IonLibrary*> ;
    using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

public:

    Private();
    ~Private();

    Err init(const QVector<MS2IonLibrary*> &ms2IonsLibrary);

    Err initTesting(const QVector<CandidateScores*> &candidateScoresPntrs);

    Err buildRTreeInput(QVector<rTreePoint> *cloudLoader);

    Err extractMs2Points(
	    float mzPrecursorMin,
		float mzPrecursorMax,
		float mzMin,
		float mzMax,
		QVector<MS2IonLibrary*> *tdPeptideFrags
        ) const;

private:

	QVector<MS2IonLibrary*> m_ms2IonsLibrary;
    RTree *m_rTree;
    bool m_isInit;

};

Ms2IonFraggertronManager::Private::Private()
: m_rTree(Q_NULLPTR)
, m_isInit(false)
{}

Ms2IonFraggertronManager::Private::~Private() {delete m_rTree;}

Err Ms2IonFraggertronManager::Private::init(
	const QVector<MS2IonLibrary*> &ms2IonsLibrary
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2IonsLibrary); ree;
	m_ms2IonsLibrary = ms2IonsLibrary;

    QVector<rTreePoint> cloudLoader;
    e = buildRTreeInput(&cloudLoader); ree;

    constexpr int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    m_isInit = true;

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::buildRTreeInput(QVector<rTreePoint> *cloudLoader) {

    ERR_INIT

	e = ErrorUtils::isNotEmpty(m_ms2IonsLibrary); ree;

	for (MS2IonLibrary *msil : m_ms2IonsLibrary) {
		rTreeCoor coor(msil->targeDecoyCandidatePairPntr->mz(msil->isDecoy), msil->ms2IonPntr->mz);
		rTreePoint point(coor, msil);
		cloudLoader->push_back(point);
	}

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::extractMs2Points(
	float mzPrecursorMin,
	float mzPrecursorMax,
	float mzMin,
	float mzMax,
	QVector<MS2IonLibrary*> *tdPeptideFrags
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    
    const rTreeSearchBox queryBox(
        rTreeCoor(mzPrecursorMin, mzMin),
        rTreeCoor(mzPrecursorMax, mzMax)
    );
    
    std::vector<rTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    tdPeptideFrags->reserve(static_cast<int>(result.size()));
    for (const auto &[fst, snd] : result) {
        tdPeptideFrags->push_back(snd);
    }

    ERR_RETURN
}

Err Ms2IonFraggertronManager::Private::initTesting(const QVector<CandidateScores*>& candidateScoresPntrs) {

    ERR_INIT

    MsCalibratomatic msCalibratomatic;;
    msCalibratomatic.setScanTimeStDev(0.8);

    e = ErrorUtils::isNotEmpty(candidateScoresPntrs); ree;

    // m_candidateScores = candidateScoresPntrs;
    // m_msCalibratomatic = msCalibratomatic;
    //
    // QVector<rTreePoint> cloudLoader;
    // for (CandidateScores *cs : m_candidateScores) {
    //
    //     float predictedScanTime = cs->targetDecoyCandidatePair->iRt();
    //
    //     const QVector<MS2Ion> &ms2Ions = cs->isDecoy
    //                                    ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
    //                                    : cs->targetDecoyCandidatePair->ms2IonsTarget();
    //
    //     for (const MS2Ion &ms2Ion : ms2Ions) {
    //         rTreeCoor coor(ms2Ion.mz, predictedScanTime);
    //         rTreePoint point(coor, {ms2Ion, cs});
    //         cloudLoader.push_back(point);
    //     }
    // }
    //
    // constexpr int maxElements = 16;
    // m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    //
    // m_isInit = true;

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


Ms2IonFraggertronManager::Ms2IonFraggertronManager()
: d_ptr(QScopedPointer<Private>(new Private()))
, m_msScanPointsTranched(nullptr)
{}

Ms2IonFraggertronManager::~Ms2IonFraggertronManager() = default;

Err Ms2IonFraggertronManager::init(
	const PythiaParameters &parameters,
	const QVector<MS2IonLibrary*> &ms2IonsLibrary,
	MsReaderPointerAcc *msReaderPtr
	) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(ms2IonsLibrary); ree;
	e = ErrorUtils::isTrue(parameters.isValid()); ree;
	e = ErrorUtils::isTrue(msReaderPtr->isInit()); ree;

	e = buildMsScanPointPntrs(msReaderPtr); ree;

	e = d_ptr->init(ms2IonsLibrary); ree;

	ERR_RETURN
}

Err Ms2IonFraggertronManager::buildMsScanPointPntrs(MsReaderPointerAcc *msReaderPtr) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPtr->isInit()); ree;

	constexpr int msLevel = 2;
	QMap<ScanNumber, MsScanInfo*> scanNumberVsMsScanInfoMS2
		= msReaderPtr->ptr->getMsScanInfos(msLevel); ree;
	e = ErrorUtils::isNotEmpty(scanNumberVsMsScanInfoMS2); ree;

	QVector<MsScanInfo*> msScanInfosMS2 = scanNumberVsMsScanInfoMS2.values().toVector();
	std::sort(
		msScanInfosMS2.begin(),
		msScanInfosMS2.end(),
		[](const MsScanInfo *l, const MsScanInfo *r) {
			if (MathUtils::tSame(l->precursorTargetMz, r->precursorTargetMz, 0.01)) {
				return l->scanNumber < r->scanNumber;
			}
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

	for (int i = 0; i < m_msScanPointsPntrsTranched.size() - 1; i++) {
		QVector<MsScanPoint*> &prev = m_msScanPointsPntrsTranched[i];
		QVector<MsScanPoint*> &next = m_msScanPointsPntrsTranched[i+1];

		const MsScanPoint *prevLastPntr = prev.back();
		QVector<MsScanPoint*> frontPntrs;
		for (MsScanPoint *nextFrontPntr : next) {
			if (prevLastPntr->scanInfoPntr->scanNumber != nextFrontPntr->scanInfoPntr->scanNumber) {
				break;
			}
			frontPntrs.push_back(nextFrontPntr);
		}

		prev.append(frontPntrs);
		next.erase(next.begin(), next.begin() + frontPntrs.size());
	}

	for (const QVector<MsScanPoint*> &v : m_msScanPointsPntrsTranched) {

		m_mzPrecursorStartVsStopResult.push_back({
			v.front()->scanInfoPntr->precursorTargetMz,
			v.back()->scanInfoPntr->precursorTargetMz,
		});
	}


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
					return l.scanInfoPntr->scanNumber < r.scanInfoPntr->scanNumber;
				}
				return l.scanInfoPntr->precursorTargetMz < r.scanInfoPntr->precursorTargetMz;
			});

		return {e, msScanPoints};
	}

}//namespace
Err Ms2IonFraggertronManager::extractScansParallel(
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
		constexpr int maxLibraryIonCount = 12;
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

		QVector<TallyResultTuple> result;

		constexpr int expectedHashTableSize = 1000;
		QVector<Tally> tallyResultsTarget;
		tallyResultsTarget.resize(expectedHashTableSize);
		QVector<Tally> tallyResultsDecoy;
		tallyResultsDecoy.resize(expectedHashTableSize);

		Tally tallyTarget;
		Tally tallyDecoy;

		for (auto it = input.begin(); it != input.end(); it++) {

			QVector<IonSearchResult> &isrs = it.value();
			TargetDecoyCandidatePair *targetDecoyCandidatePair = it.key();

			std::sort(
				isrs.begin(),
				isrs.end(),
				[](const IonSearchResult &l, const IonSearchResult &r) {
					return l.msScanPointPntr->scanInfoPntr->scanNumber <
						   r.msScanPointPntr->scanInfoPntr->scanNumber;
				});

			int currentScanNumber = -1;

			for (const IonSearchResult &isr : isrs) {

				const ScanNumber scanNumber = isr.msScanPointPntr->scanInfoPntr->scanNumber;

				if (currentScanNumber != scanNumber) {
					if (currentScanNumber > 0) {
						constexpr int occurenceCountMin = 3;
						if (tallyTarget.occurrence > occurenceCountMin) {
							tallyTarget.scanNumber = currentScanNumber;
							e = EigenUtils::cosineSimilarity(
								EigenUtils::convertQVectorToEigenVector(tallyTarget.intensitiesEmperical),
								EigenUtils::convertQVectorToEigenVector(tallyTarget.intensitiesTheoretical),
								&tallyTarget.cosineSimilarity
								);
							tallyResultsTarget.push_back(tallyTarget);
						}
						if (tallyDecoy.occurrence > occurenceCountMin) {
							tallyDecoy.scanNumber = currentScanNumber;
							e = EigenUtils::cosineSimilarity(
								EigenUtils::convertQVectorToEigenVector(tallyDecoy.intensitiesEmperical),
								EigenUtils::convertQVectorToEigenVector(tallyDecoy.intensitiesTheoretical),
								&tallyDecoy.cosineSimilarity
								);
							tallyResultsDecoy.push_back(tallyDecoy);
						}

						tallyTarget = Tally();
						tallyDecoy = Tally();
					}
					currentScanNumber = scanNumber;
				}

				if (isr.ms2IonLibraryPntr->isDecoy) {
					tallyDecoy.occurrence++;
					tallyDecoy.ranks.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->rank);
					tallyDecoy.intensitiesEmperical.push_back(isr.msScanPointPntr->intensityVal);
					tallyDecoy.intensitiesTheoretical.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->intensity);
					tallyDecoy.intensitiesSum += isr.msScanPointPntr->intensityVal;
					continue;
				}

				tallyTarget.occurrence++;
				tallyTarget.ranks.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->rank);
				tallyTarget.intensitiesEmperical.push_back(isr.msScanPointPntr->intensityVal);
				tallyTarget.intensitiesTheoretical.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->intensity);
				tallyTarget.intensitiesSum += isr.msScanPointPntr->intensityVal;
			}

			if (tallyResultsTarget.isEmpty() && tallyResultsDecoy.isEmpty()) {
				continue;
			}

			const auto sortLogic = [](const Tally &l, const Tally &r) {
				if (l.occurrence == r.occurrence) {
					return l.cosineSimilarity > r.cosineSimilarity;
				}
				return l.occurrence > r.occurrence;
			};
			std::sort(tallyResultsTarget.begin(), tallyResultsTarget.end(), sortLogic);
			std::sort(tallyResultsDecoy.begin(), tallyResultsDecoy.end(), sortLogic);

			QVector<TallyResult> tallyResultsFinalTarget;
			tallyResultsFinalTarget.reserve(tallyResultsTarget.size());
			QVector<TallyResult> tallyResultsFinalDecoy;
			tallyResultsFinalDecoy.reserve(tallyResultsDecoy.size());

			constexpr int topN = 3;
			for (int i = 0; i < std::min(tallyResultsTarget.size(), topN); ++i) {
				TallyResult t;

				const Tally &tally = tallyResultsTarget[i];
				t.scanNumber = tally.scanNumber;
				t.occurrence = tally.occurrence;
				t.isDecoy = false;
				t.cosineSimilarity = tally.cosineSimilarity;
				t.intensitiesSum = tally.intensitiesSum;
				t.ranks = tally.ranks;

				tallyResultsFinalTarget.push_back(t);
			}

			for (int i = 0; i < std::min(tallyResultsDecoy.size(), topN); ++i) {
				TallyResult t;

				const Tally &tally = tallyResultsDecoy[i];
				t.scanNumber = tally.scanNumber;
				t.occurrence = tally.occurrence;
				t.isDecoy = true;
				t.cosineSimilarity = tally.cosineSimilarity;
				t.intensitiesSum = tally.intensitiesSum;
				t.ranks = tally.ranks;

				tallyResultsFinalDecoy.push_back(t);
			}

			result.push_back({targetDecoyCandidatePair, tallyResultsFinalTarget, tallyResultsFinalDecoy});
			tallyResultsTarget.clear();
			tallyResultsDecoy.clear();
		}

		return {e, result};
	}

	QPair<Err, QVector<TallyResultTuple>> performFraggingLogic(
		const ProcessingGroup &pg,
		const PythiaParameters &parameters
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(*pg.msScanPointsPntr); rree;

		if (pg.ms2IonsLibrary.isEmpty()) {
			qDebug()
				<< qPrintable(S_GLOBAL_TIMER.elapsed())
				<< "No ms2IonsLibrary found for precursor range"
				<< pg.mzPrecursorRangeMinMax.first << "-"
				<< pg.mzPrecursorRangeMinMax.second
			;
			return {e, {}};
		}

		Ms2IonFraggertronManager fragger;
		e = fragger.init(pg.ms2IonsLibrary); rree;

		QHash<TargetDecoyCandidatePair*, QVector<IonSearchResult>> ionSearchResults;
		constexpr int batchSize = 2e5; //TODO make this settable in params
		ionSearchResults.reserve(batchSize);

		for (MsScanPoint *mssp : *pg.msScanPointsPntr) {
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

			if (tcp->peptideStringWithMods() != "YGGPNHHLPLPDNWK" || tcp->charge() != 3) {
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
		QElapsedTimer et;
		et.start();
		QPair<Err, QVector<TallyResultTuple>> result
				= collateScanNumberVsOccurrencesTargetDecoyCandidatePairs(ionSearchResults); rree;
		// qDebug() << et.restart() << "SLFJDSL";

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

		constexpr int batchSize = 2e5; //TODO make this settable in params
		const int libTrancheSize = std::max(targetDecoyCandidatePairsPntrs.size() / batchSize, 1);
		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Processing library in" << libTrancheSize << "tranches";

		e = ParallelUtils::trancheVectorForParallelization(
			targetDecoyCandidatePairsPntrs,
			libTrancheSize,
			targetDecoyCandidatePairsPntrsTranched
			); ree;

		ERR_RETURN
	}

	Err combineDuplicateCollates(QVector<TallyResultTuple> *tallyResultsFinal) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(*tallyResultsFinal); ree;

		std::sort(
			tallyResultsFinal->begin(),
			tallyResultsFinal->end(),
			[](const TallyResultTuple &l, const TallyResultTuple &r) {
				if (std::get<1>(l).front().occurrence == std::get<1>(r).front().occurrence) {
					return std::get<2>(l).front().occurrence > std::get<2>(r).front().occurrence;
				}
				return std::get<1>(l).front().occurrence > std::get<1>(r).front().occurrence;
			});

		QHash<TargetDecoyCandidatePair*, TallyResultTuple> tallyResultTuplesCombiner;
		for (const TallyResultTuple &trt : *tallyResultsFinal) {

			TargetDecoyCandidatePair *tdcp = std::get<0>(trt);

			if (tallyResultTuplesCombiner.contains(tdcp)) {

				TallyResultTuple &tallyTupleOld = tallyResultTuplesCombiner[tdcp];
				QVector<TallyResultTarget> &tallyResultTargetOld = std::get<1>(tallyTupleOld);

				const QVector<TallyResultTarget> &tallyResultTargetNew = std::get<1>(trt);

				QHash<ScanNumber, TallyResult> tallyResultTargetNewHash;
				for (const TallyResultTarget &t  : tallyResultTargetNew) {
					tallyResultTargetNewHash[t.scanNumber] = t;
				}
				for (TallyResultTarget &trtOld : tallyResultTargetOld) {
					if (!tallyResultTargetNewHash.contains(trtOld.scanNumber)) {
						continue;
					}

					TallyResult &tn = tallyResultTargetNewHash[trtOld.scanNumber];

					if (tn.occurrence < 1) {
						continue;
					}

					trtOld.occurrence += tn.occurrence;
					trtOld.ranks.append(tn.ranks);

					e = ErrorUtils::isFalse(tallyResultTargetNewHash.contains(trtOld.scanNumber)); eee_absorb;
				}

				std::sort(
					tallyResultTargetOld.rbegin(),
					tallyResultTargetOld.rend(),
					[](const TallyResultTarget &l, const TallyResultTarget &r){return l.occurrence < r.occurrence;}
					);

				const QVector<TallyResultTarget> &tallyResultDecoyNew = std::get<2>(trt);
				QHash<ScanNumber, TallyResult> tallyResultDecoyNewHash;
				for (const TallyResultDecoy &t  : tallyResultDecoyNew) {
					tallyResultDecoyNewHash[t.scanNumber] = t;
				}
				QVector<TallyResultDecoy> &tallyResultDecoyOld = std::get<2>(tallyTupleOld);
				for (TallyResultDecoy &trdOld : tallyResultDecoyOld) {

					if (!tallyResultDecoyNewHash.contains(trdOld.scanNumber)) {
						continue;
					}

					TallyResult &tn = tallyResultDecoyNewHash[trdOld.scanNumber];

					if (tn.occurrence < 1) {
						continue;
					}

					trdOld.occurrence += tn.occurrence;
					trdOld.ranks.append(tn.ranks);

					e = ErrorUtils::isFalse(tallyResultDecoyNewHash.contains(trdOld.scanNumber)); eee_absorb;
				}

				std::sort(
					tallyResultDecoyOld.rbegin(),
					tallyResultDecoyOld.rend(),
					[](const TallyResultTarget &l, const TallyResultTarget &r){return l.occurrence < r.occurrence;}
					);

				continue;
			}

			tallyResultTuplesCombiner.insert(tdcp, trt);
		}

		*tallyResultsFinal = tallyResultTuplesCombiner.values().toVector();

		ERR_RETURN
	}

	void reportResults(const QVector<TallyResultTuple> &tallyResultsFinal) 	{

		struct T {
			TargetDecoyCandidatePair *targetDecoyCandidatePair = nullptr;
			ScanNumber scanNumber = -1;
			bool isDecoy = false;
			Occurrence occurrence = -1;
			float cosineSimilarity = -1;
			float rankMean = -1;
			int rankBest = -1;
		};

		QVector<T> ts;

		for (const TallyResultTuple &tpl : tallyResultsFinal) {

			T t;
			if (std::get<1>(tpl).front().occurrence > 0) {
				t.targetDecoyCandidatePair = std::get<0>(tpl);
				t.occurrence = std::get<1>(tpl).front().occurrence;
				t.scanNumber = std::get<1>(tpl).front().scanNumber;
				t.rankMean = MathUtils::mean(std::get<1>(tpl).front().ranks);
				t.rankBest = *std::min_element(std::get<1>(tpl).front().ranks.begin(), std::get<1>(tpl).front().ranks.end());
				t.isDecoy = false;
				t.cosineSimilarity = std::get<1>(tpl).front().cosineSimilarity;
				ts.push_back(t);
			}

			T td;
			if (std::get<2>(tpl).front().occurrence > 0) {
				td.targetDecoyCandidatePair = std::get<0>(tpl);
				td.occurrence = std::get<2>(tpl).front().occurrence;
				t.scanNumber = std::get<2>(tpl).front().scanNumber;
				td.rankMean = MathUtils::mean(std::get<2>(tpl).front().ranks);
				td.rankBest = *std::min_element(std::get<2>(tpl).front().ranks.begin(), std::get<2>(tpl).front().ranks.end());
				t.cosineSimilarity = std::get<2>(tpl).front().cosineSimilarity;
				td.isDecoy = true;
				ts.push_back(td);
			}
		}

		std::sort(
			ts.begin(),
			ts.end(),
			[](const T &l, const T &r) {
				if (l.occurrence == r.occurrence) {
					return l.cosineSimilarity > r.cosineSimilarity;
				}

				return l.occurrence > r.occurrence;
			}
			);

		int counter = 0;
		int decoy = 0;
		QHash<TargetDecoyCandidatePair*, T> founds;
		for (const T &t : ts) {

			if (t.isDecoy) {
				decoy++;
			}

			const float q = static_cast<float>(decoy) / ++counter;

			// if (founds.contains(t.targetDecoyCandidatePair)) {
			if (true) {
				std::cout
					<< qPrintable(S_GLOBAL_TIMER.elapsed()) << " "
					<< counter << " fdsafda "
					<< q << " "
					<< t.targetDecoyCandidatePair->peptideStringWithMods().toStdString()  << " "
					<< t.targetDecoyCandidatePair->charge()  << " "
					<< t.occurrence << founds.value(t.targetDecoyCandidatePair).occurrence  << " "
					// << t.rankMean  << founds.value(t.targetDecoyCandidatePair).occurrence  << " "
					// << t.rankBest  << " "
					<< t.scanNumber << founds.value(t.targetDecoyCandidatePair).scanNumber  << " "
					<< t.isDecoy
					<< std::endl;
			}

			founds.insert(t.targetDecoyCandidatePair, t);
			if (q > 0.01) break;

			//mean mz thomsons found average, higher is more specific than lower (this is a reminder to add another NN/LDA feature)
		}
	}

}//namespace
Err Ms2IonFraggertronManager::performFragging() {
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

	e = combineDuplicateCollates(&tallyResultsFinal); ree;

	reportResults(tallyResultsFinal);

	ERR_RETURN
}

Err Ms2IonFraggertronManager::extractMs2Points(
	float mzPrecursorMin,
	float mzPrecursorMax,
    float mzMin,
    float mzMax,
	QVector<MS2IonLibrary*> *tdPeptideFrags
    ) const {
    ERR_INIT
    e = d_ptr->extractMs2Points(
    	mzPrecursorMin,
    	mzPrecursorMax,
        mzMin,
        mzMax,
        tdPeptideFrags
        ); ree;
    ERR_RETURN
}

Err Ms2IonFraggertronManager::initTesting(const QVector<CandidateScores*>& candidateScoresPntrs) const {
    ERR_INIT
    e = d_ptr->initTesting(candidateScoresPntrs); ree;
    ERR_RETURN
}
