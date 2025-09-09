
//
// Created by andrewnichols on 9/3/25.
//

#include "MsFraggertron.h"

#include "EigenUtils.h"


#include <iomanip>
namespace {
	constexpr int MAX_LIBRARY_ION_COUNT = 12;
}

MsFraggertron::MsFraggertron()
: m_msReaderPtr(nullptr)
, m_msCalibratomatic(nullptr)
{}

Err MsFraggertron::init(
	const PythiaParameters & params,
	MsReaderPointerAcc *msReaderPtr
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(msReaderPtr->isInit()); ree;
	e = ErrorUtils::isTrue(params.isValid()); ree;

	m_msReaderPtr = msReaderPtr;
	m_parameters = params;

	e = m_precursorMzArrangetron.init(m_msReaderPtr); ree;

	constexpr int msScanIonTranches = 64;
	e = m_precursorMzArrangetron.trancheMsScanPoints(
		msScanIonTranches,
		&m_scanNumberMzIntensityTranched
		); ree;

	ERR_RETURN
}

bool MsFraggertron::isInit() const {
	return !m_scanNumberMzIntensityTranched.isEmpty();
}

namespace {

	constexpr int precursorMzKeyHashingPrecision = 3;

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

	void reportResults(const QVector<TallyResultTuple> &tallyResultsFinal) 	{

		struct T {
			TargetDecoyCandidatePair *targetDecoyCandidatePair = nullptr;
			ScanNumber scanNumber = -1;
			bool isDecoy = false;
			int altCount = -1;
			Occurrence occurrence = -1;
			Occurrence occurrencePrevious = -1;
			float cosineSimilarity = -1;
			float relativeIntensityDifferenceAverage = -1;
			float rankMean = -1;
			int rankBest = -1;
			QVector<int> ranks;
		};

		QVector<T> ts;

		for (const TallyResultTuple &tpl : tallyResultsFinal) {

			T t;
			if (std::get<1>(tpl).front().occurrence > 0) {
				t.targetDecoyCandidatePair = std::get<0>(tpl);
				t.occurrence = std::get<1>(tpl).front().occurrence;
				t.altCount = std::get<1>(tpl).front().totalFound;
				t.occurrencePrevious = t.altCount > 1 ? std::get<1>(tpl)[1].occurrence : 0;
				t.scanNumber = std::get<1>(tpl).front().scanNumber;
				t.rankMean = MathUtils::mean(std::get<1>(tpl).front().ranks);
				t.rankBest = *std::min_element(std::get<1>(tpl).front().ranks.begin(), std::get<1>(tpl).front().ranks.end());
				t.isDecoy = false;
				t.cosineSimilarity = std::get<1>(tpl).front().cosineSimilarity;
				t.relativeIntensityDifferenceAverage = std::get<1>(tpl).front().relativeIntensityDifferenceAverage;
				t.ranks = std::get<1>(tpl).front().ranks;
				ts.push_back(t);
			}

			T td;
			if (std::get<2>(tpl).front().occurrence > 0) {
				td.targetDecoyCandidatePair = std::get<0>(tpl);
				td.occurrence = std::get<2>(tpl).front().occurrence;
				td.altCount = std::get<2>(tpl).front().totalFound;
				td.occurrencePrevious = td.altCount > 1 ? std::get<2>(tpl)[1].occurrence : 0;
				td.scanNumber = std::get<2>(tpl).front().scanNumber;
				td.rankMean = MathUtils::mean(std::get<2>(tpl).front().ranks);
				td.rankBest = *std::min_element(std::get<2>(tpl).front().ranks.begin(), std::get<2>(tpl).front().ranks.end());
				td.cosineSimilarity = std::get<2>(tpl).front().cosineSimilarity;
				td.relativeIntensityDifferenceAverage = std::get<2>(tpl).front().relativeIntensityDifferenceAverage;
				td.isDecoy = true;
				td.ranks = std::get<2>(tpl).front().ranks;
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
			});

		// std::sort(
		// 	ts.begin(),
		// 	ts.end(),
		// 	[](const T &l, const T &r) {
		// 		if (l.targetDecoyCandidatePair->peptideStringWithMods() == r.targetDecoyCandidatePair->peptideStringWithMods()) {
		//
		// 			if (l.isDecoy < r.isDecoy) {
		// 				return l.targetDecoyCandidatePair->charge() == r.targetDecoyCandidatePair->charge();
		// 			}
		//
		// 			return l.targetDecoyCandidatePair->charge() > r.targetDecoyCandidatePair->charge();
		// 		}
		//
		// 		return l.targetDecoyCandidatePair->peptideStringWithMods() > r.targetDecoyCandidatePair->peptideStringWithMods();
		// 	});

		int counter = 0;
		int decoy = 0;
		QHash<TargetDecoyCandidatePair*, T> founds;
		// int dupes = 0;
		// PeptideStringWithMods curr;
		for (const T &t : ts) {

			if (t.isDecoy) {
				decoy++;
			}

			const float q = static_cast<float>(decoy) / ++counter;

			// if (founds.contains(t.targetDecoyCandidatePair)) {
			if (true) {

				// if(counter++ == 10000) {
				// 	break;
				// }
				//
				// if (t.targetDecoyCandidatePair->peptideStringWithMods() != curr) {
				//
				// 	if (dupes > 4) {
				// 		std::cout << "HEREHREHEHE" << std::endl;
				// 	}
				// 	// std::cout << dupes << std::endl;
				// 	dupes = 0;
				// 	curr = t.targetDecoyCandidatePair->peptideStringWithMods();
				// }
				//
				// dupes++;
				std::cout
					<< qPrintable(S_GLOBAL_TIMER.elapsed()) << " "
					<< counter << " "
					<< t.isDecoy << " "
					<< t.altCount << " "
					<< t.scanNumber
					<< " fdsafda "
					<< q << " "
					<< t.targetDecoyCandidatePair->peptideStringWithMods().toStdString()  << " "
					<< t.targetDecoyCandidatePair->charge()  << " "
					<< t.occurrence << " " << t.occurrencePrevious << " "
					// << founds.value(t.targetDecoyCandidatePair).occurrence  << " "
					// << t.rankMean  << founds.value(t.targetDecoyCandidatePair).occurrence  << " "
					// << t.rankBest  << " "
					// << t.scanNumber << " "
					// << founds.value(t.targetDecoyCandidatePair).scanNumber  << " "
					// << std::endl
				;

					for (int r : t.ranks) {
						std::cout << r << ",";
					}
					std::cout << std::endl;


			}

			founds.insert(t.targetDecoyCandidatePair, t);
			if (q > 0.01) break;

			//TODO mean mz thomsons found average, higher is more specific than lower (this is a reminder to add another NN/LDA feature)
		}
	}

}//namespace
Err MsFraggertron::performFragging(const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairsPntrs) {

	ERR_INIT

	QElapsedTimer et;
	et.start();

	e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairsPntrs); ree;
	e = ErrorUtils::isTrue(m_precursorMzArrangetron.isInit()); ree;

	QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePairsPntrsTranched;
	e = buildTargetDecoyCandidatePairsPntrsTranched(
		targetDecoyCandidatePairsPntrs,
		&targetDecoyCandidatePairsPntrsTranched
		); ree;

	QVector<TallyResultTuple> tallyResultsFinal;
	for (const QVector<TargetDecoyCandidatePair*> &tdcps : targetDecoyCandidatePairsPntrsTranched) {
		const QPair<Err, QVector<TallyResultTuple>> result = processTargetDecoyCandidatePairsPntrsTranch(tdcps); ree;
		e = result.first; ree;
		tallyResultsFinal.append(result.second);
	}

	reportResults(tallyResultsFinal);

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

			ms2IonLibraries->reserve(ms2IonsLibraryTranche.size() * targetDecoyDoubler * MAX_LIBRARY_ION_COUNT);
			for (std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy> &tpl : ms2IonsLibraryTranche) {
				TargetDecoyCandidatePair *tdcp = std::get<0>(tpl);

				MS2IonsTarget &ms2IonsTarget = std::get<1>(tpl);
				MS2IonsDecoy &ms2IonsDecoy = std::get<2>(tpl);

				for(int i = 0; i < std::min(MAX_LIBRARY_ION_COUNT, ms2IonsTarget.size()); ++i) {
					MS2Ion *ms2Ion = &ms2IonsTarget[i];
					MS2IonLibrary ms2IonLibrary;
					ms2IonLibrary.ms2IonPntr = ms2Ion;
					ms2IonLibrary.targeDecoyCandidatePairPntr = tdcp;
					ms2IonLibrary.isDecoy = false;
					ms2IonLibraries->push_back(ms2IonLibrary);
				}

				for(int i = 0; i < std::min(MAX_LIBRARY_ION_COUNT, ms2IonsDecoy.size()); ++i) {
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

	void ms2IonLibrarySortLogic(QVector<MS2IonLibrary> &ms2IonLibraries) {
		std::sort(
			ms2IonLibraries.begin(),
			ms2IonLibraries.end(),
			[](const MS2IonLibrary &l, const MS2IonLibrary &r) {

				if (
					constexpr int hashingPrecision = 3;
					MathUtils::hashDecimal(l.targeDecoyCandidatePairPntr->mz(l.isDecoy),hashingPrecision)
					== MathUtils::hashDecimal(r.targeDecoyCandidatePairPntr->mz(r.isDecoy),hashingPrecision)
					) {
					return l.ms2IonPntr->mz < r.ms2IonPntr->mz;
				}

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

	void filterByPPM(Tally *t) {

		const Eigen::VectorX<float> mzEmpericalVec = EigenUtils::convertQVectorToEigenVector(t->mzEmperical);
		Eigen::VectorX<float> intensityEmpericalVecNorm = EigenUtils::convertQVectorToEigenVector(t->intensitiesEmperical);
		intensityEmpericalVecNorm.array() /= intensityEmpericalVecNorm.maxCoeff();

		const Eigen::VectorX<float> mzTheoreticalVec = EigenUtils::convertQVectorToEigenVector(t->mzTheoretical);

		const Eigen::VectorX<float> ppm = 1e6 * ((mzEmpericalVec - mzTheoreticalVec).array() / mzTheoreticalVec.array());
		const Eigen::VectorX<float> ppmDiffsFromMeanAbs = (ppm.array() - ppm.mean()).array().abs();

		using Rank = int;
		using PPMDiff = float;
		QVector<QPair<Rank, PPMDiff>> bestIons(MAX_LIBRARY_ION_COUNT, QPair<Rank, PPMDiff>(-1, std::numeric_limits<float>::max()));
		QVector<int> ogIndexes(MAX_LIBRARY_ION_COUNT, -1);
		for (int i = 0; i < t->occurrence; ++i) {

			ogIndexes[t->ranks[i]] = bestIons[t->ranks[i]].second > ppmDiffsFromMeanAbs.coeff(i)
					   ? i
					   : ogIndexes[t->ranks[i]];

			bestIons[t->ranks[i]] = bestIons[t->ranks[i]].second > ppmDiffsFromMeanAbs.coeff(i)
								  ? QPair<Rank, PPMDiff>(t->ranks[i], ppmDiffsFromMeanAbs.coeff(i))
			                      : bestIons[t->ranks[i]];
		}

		Tally tNew;
		tNew.scanNumber = t->scanNumber;

		for (const QPair<int, float> &pr : bestIons) {

			if (pr.first < 0 || intensityEmpericalVecNorm.coeff(ogIndexes[pr.first]) < 0.01) {
				continue;
			}

			tNew.occurrence++;
			tNew.ranks.push_back(pr.first);
			tNew.intensitiesEmperical.push_back(t->intensitiesEmperical[ogIndexes[pr.first]]);
			tNew.intensitiesTheoretical.push_back(t->intensitiesTheoretical[ogIndexes[pr.first]]);
			tNew.mzEmperical.push_back(t->mzEmperical[ogIndexes[pr.first]]);
			tNew.mzTheoretical.push_back(t->mzTheoretical[ogIndexes[pr.first]]);
		}

		*t = tNew;
	}

	void filterByMedianScanIntensity(
		float medianScanIntensity,
		Tally *t
		) {

		Tally tNew;
		tNew.scanNumber = t->scanNumber;

		for (int i = 0; i < t->mzEmperical.size(); ++i) {

			if (t->intensitiesEmperical[i] < medianScanIntensity) {
				continue;
			}

			tNew.occurrence++;
			tNew.ranks.push_back(t->ranks[i]);
			tNew.intensitiesEmperical.push_back(t->intensitiesEmperical[i]);
			tNew.mzEmperical.push_back(t->mzEmperical[i]);
			tNew.mzTheoretical.push_back(t->mzTheoretical[i]);
			tNew.intensitiesTheoretical.push_back(t->intensitiesTheoretical[i]);
		}

		tNew.intensitiesSum = std::accumulate(tNew.intensitiesEmperical.begin(), tNew.intensitiesEmperical.end(), 0.0f);

		*t = tNew;
	}

	QPair<Err, QVector<TallyResultTuple>> collateScanNumberVsOccurrencesTargetDecoyCandidatePairs(
		QHash<TargetDecoyCandidatePair*, QVector<IonSearchResult2>> &input
		) {

		ERR_INIT

		QVector<TallyResultTuple> result;

		constexpr int expectedHashTableSize = 1000;
		QVector<Tally> tallyResultsTarget;
		tallyResultsTarget.resize(expectedHashTableSize);
		QVector<Tally> tallyResultsDecoy;
		tallyResultsDecoy.resize(expectedHashTableSize);

		for (auto it = input.begin(); it != input.end(); it++) {

			QVector<IonSearchResult2> &isrs = it.value();
			TargetDecoyCandidatePair *targetDecoyCandidatePair = it.key();

			std::sort(
				isrs.begin(),
				isrs.end(),
				[](const IonSearchResult2 &l, const IonSearchResult2 &r) {
					return l.msScanPointPntr->scanInfoPntr->scanNumber <
						   r.msScanPointPntr->scanInfoPntr->scanNumber;
				});

			int currentScanNumber = -1;

			Tally tallyTarget;
			Tally tallyDecoy;

			for (const IonSearchResult2 &isr : isrs) {

				const ScanNumber scanNumber = isr.msScanPointPntr->scanInfoPntr->scanNumber;

				if (currentScanNumber != scanNumber) {
					if (currentScanNumber > 0) {
						constexpr int occurenceCountMin = 3;
						if (tallyTarget.occurrence > occurenceCountMin) {

							tallyTarget.scanNumber = currentScanNumber;

							filterByPPM(&tallyTarget);

							if (isr.msScanPointPntr->scanInfoPntr->pointCount > 1000) {
								filterByMedianScanIntensity(
									isr.msScanPointPntr->scanInfoPntr->medianIntensity,
									&tallyTarget
									);
							}


							if (tallyTarget.occurrence <= occurenceCountMin) {
								continue;
							}

							e = EigenUtils::cosineSimilarity(
								EigenUtils::convertQVectorToEigenVector(tallyTarget.intensitiesEmperical),
								EigenUtils::convertQVectorToEigenVector(tallyTarget.intensitiesTheoretical),
								&tallyTarget.cosineSimilarity
								);

							Eigen::VectorX<float> empIntensityRelative = EigenUtils::convertQVectorToEigenVector(tallyTarget.intensitiesEmperical);
							empIntensityRelative /= empIntensityRelative.maxCoeff();
							Eigen::VectorX<float> theoIntensityRelative = EigenUtils::convertQVectorToEigenVector(tallyTarget.intensitiesTheoretical);
							Eigen::VectorX<float> diffAbs = (empIntensityRelative - theoIntensityRelative).array().abs();
							tallyTarget.relativeIntensityDifferenceAverage = diffAbs.sum() / tallyTarget.occurrence;

							tallyResultsTarget.push_back(tallyTarget);

							// if (
							// 	tallyTarget.occurrence > 11
							// 	// && targetDecoyCandidatePair->charge() == 3
							// 	// && targetDecoyCandidatePair->peptideStringWithMods() == "VEYSAYLDVFSQPEK"
							// 	) {
							// 	Eigen::VectorX<float> intzNorm = EigenUtils::convertQVectorToEigenVector(tallyTarget.intensitiesEmperical);
							// 	intzNorm.array() /= intzNorm.maxCoeff();
							//
							// 	qDebug() << tallyTarget.occurrence << "SDFKJSDL";
							// 	qDebug() << tallyTarget.ranks << "SDFKJSDL";
							// 	qDebug() << tallyTarget.relativeIntensityDifferenceAverage << tallyTarget.cosineSimilarity << "SDFKJSDL";
							// 	qDebug() << tallyTarget.intensitiesEmperical;
							// 	qDebug() << EigenUtils::convertEigenVectorToQVector(intzNorm) << "SDFKJSDL";
							// 	qDebug() << tallyTarget.intensitiesTheoretical << "SDFKJSDL";
							// 	qDebug() << MathUtils::mean(tallyTarget.intensitiesTheoretical) << MathUtils::stDev(tallyTarget.intensitiesTheoretical) <<
							// 	MathUtils::stDev(tallyTarget.intensitiesTheoretical)/ MathUtils::mean(tallyTarget.intensitiesTheoretical) << "mean";
							// 	qDebug() << tallyTarget.mzEmperical << "SDFKJSDL";
							// 	qDebug() << tallyTarget.mzTheoretical << "SDFKJSDL";
							//
							// 	double corr;
							// 	EigenUtils::pearsonCorrelation(intzNorm, EigenUtils::convertQVectorToEigenVector(tallyTarget.intensitiesTheoretical), &corr);
							// 	qDebug() << tallyTarget.cosineSimilarity << corr << "SDFKJSDL";
							// 	const Eigen::VectorX<float> mzEmpericalVec = EigenUtils::convertQVectorToEigenVector(tallyTarget.mzEmperical);
							// 	const Eigen::VectorX<float> mzTheoreticalVec = EigenUtils::convertQVectorToEigenVector(tallyTarget.mzTheoretical);
							// 	const Eigen::VectorX<float> ppm = 1e6 * ((mzEmpericalVec - mzTheoreticalVec).array() / mzTheoreticalVec.array());
							// 	const Eigen::VectorX<float> ppmDiffsFromMeanAbs = (ppm.array() - ppm.mean()).array().abs();
							// 	qDebug() << ppm.mean() << EigenUtils::calculateStDevOfVector(ppm);
							// 	qDebug() << EigenUtils::convertEigenVectorToQVector(ppm);;
							// 	qDebug() << EigenUtils::convertEigenVectorToQVector(ppmDiffsFromMeanAbs);;
							// 	qDebug() << "*********";
							// }
						}
						if (tallyDecoy.occurrence > occurenceCountMin) {

							tallyDecoy.scanNumber = currentScanNumber;

							filterByPPM(&tallyDecoy);

							if (isr.msScanPointPntr->scanInfoPntr->pointCount > 1000) {
								filterByMedianScanIntensity(
									isr.msScanPointPntr->scanInfoPntr->medianIntensity,
									&tallyDecoy
									);
							}

							if (tallyDecoy.occurrence <= occurenceCountMin) {
								continue;
							}

							e = EigenUtils::cosineSimilarity(
								EigenUtils::convertQVectorToEigenVector(tallyDecoy.intensitiesEmperical),
								EigenUtils::convertQVectorToEigenVector(tallyDecoy.intensitiesTheoretical),
								&tallyDecoy.cosineSimilarity
								);

							Eigen::VectorX<float> empIntensityRelative = EigenUtils::convertQVectorToEigenVector(tallyDecoy.intensitiesEmperical);
							empIntensityRelative /= empIntensityRelative.maxCoeff();
							Eigen::VectorX<float> theoIntensityRelative = EigenUtils::convertQVectorToEigenVector(tallyDecoy.intensitiesTheoretical);
							Eigen::VectorX<float> diffAbs = (empIntensityRelative - theoIntensityRelative).array().abs();
							tallyDecoy.relativeIntensityDifferenceAverage = diffAbs.sum() / tallyDecoy.occurrence;

							tallyResultsDecoy.push_back(tallyDecoy);

							// if (
							// 	tallyDecoy.occurrence > 7
							// 	// && targetDecoyCandidatePair->charge() == 3
							// 	// && targetDecoyCandidatePair->peptideStringWithMods() == "VEYSAYLDVFSQPEK"
							// 	) {
							// 		Eigen::VectorX<float> intzNorm = EigenUtils::convertQVectorToEigenVector(tallyDecoy.intensitiesEmperical);
							// 		intzNorm.array() /= intzNorm.maxCoeff();
							//
							// 		qDebug() << tallyDecoy.occurrence << isr.msScanPointPntr->scanInfoPntr->medianIntensity
							// 				<< isr.msScanPointPntr->scanInfoPntr->pointCount << "SDFKJSDL";
							// 		qDebug() << tallyDecoy.ranks << "SDFKJSDL";
							// 		qDebug() << tallyDecoy.relativeIntensityDifferenceAverage << tallyDecoy.cosineSimilarity << "SDFKJSDL";
							// 		qDebug() << tallyDecoy.intensitiesEmperical;
							// 		qDebug() << EigenUtils::convertEigenVectorToQVector(intzNorm) << "SDFKJSDL";
							// 		qDebug() << tallyDecoy.intensitiesTheoretical << "SDFKJSDL";
							// 		qDebug() << MathUtils::mean(tallyDecoy.intensitiesTheoretical) << MathUtils::stDev(tallyDecoy.intensitiesTheoretical) <<
							// 		MathUtils::stDev(tallyDecoy.intensitiesTheoretical)/ MathUtils::mean(tallyDecoy.intensitiesTheoretical) << "mean";
							// 		// qDebug() << tallyDecoy.mzEmperical << "SDFKJSDL";
							// 		// qDebug() << tallyDecoy.mzTheoretical << "SDFKJSDL";
							//
							// 		double corr;
							// 		EigenUtils::pearsonCorrelation(intzNorm, EigenUtils::convertQVectorToEigenVector(tallyDecoy.intensitiesTheoretical), &corr);
							// 		qDebug() << tallyDecoy.cosineSimilarity << corr << "SDFKJSDL";
							// 		// const Eigen::VectorX<float> mzEmpericalVec = EigenUtils::convertQVectorToEigenVector(tallyDecoy.mzEmperical);
							// 		// const Eigen::VectorX<float> mzTheoreticalVec = EigenUtils::convertQVectorToEigenVector(tallyDecoy.mzTheoretical);
							// 		// const Eigen::VectorX<float> ppm = 1e6 * ((mzEmpericalVec - mzTheoreticalVec).array() / mzTheoreticalVec.array());
							// 		// const Eigen::VectorX<float> ppmDiffsFromMeanAbs = (ppm.array() - ppm.mean()).array().abs();
							// 		// qDebug() << ppm.mean() << EigenUtils::calculateStDevOfVector(ppm);
							// 		// qDebug() << EigenUtils::convertEigenVectorToQVector(ppm);;
							// 		// qDebug() << EigenUtils::convertEigenVectorToQVector(ppmDiffsFromMeanAbs);;
							// 		qDebug() << "*********";
							// }
						}

					}

					tallyTarget = Tally();
					tallyDecoy = Tally();
					currentScanNumber = scanNumber;
				}

				if (isr.ms2IonLibraryPntr->isDecoy) {
					tallyDecoy.occurrence++;
					tallyDecoy.ranks.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->rank);
					tallyDecoy.mzEmperical.push_back(isr.msScanPointPntr->mzVal);
					tallyDecoy.mzTheoretical.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->mz);
					tallyDecoy.intensitiesEmperical.push_back(isr.msScanPointPntr->intensityVal);
					tallyDecoy.intensitiesTheoretical.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->intensity);
					tallyDecoy.intensitiesSum += isr.msScanPointPntr->intensityVal;
					continue;
				}

				tallyTarget.occurrence++;
				tallyTarget.ranks.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->rank);
				tallyTarget.mzEmperical.push_back(isr.msScanPointPntr->mzVal);
				tallyTarget.mzTheoretical.push_back(isr.ms2IonLibraryPntr->ms2IonPntr->mz);
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
				t.totalFound = tallyResultsTarget.size();
				t.isDecoy = false;
				t.cosineSimilarity = tally.cosineSimilarity;
				t.relativeIntensityDifferenceAverage = tally.relativeIntensityDifferenceAverage;
				t.intensitiesSum = tally.intensitiesSum;
				t.ranks = tally.ranks;

				tallyResultsFinalTarget.push_back(t);
			}

			for (int i = 0; i < std::min(tallyResultsDecoy.size(), topN); ++i) {
				TallyResult t;

				const Tally &tally = tallyResultsDecoy[i];
				t.scanNumber = tally.scanNumber;
				t.occurrence = tally.occurrence;
				t.totalFound = tallyResultsDecoy.size();
				t.isDecoy = true;
				t.cosineSimilarity = tally.cosineSimilarity;
				t.relativeIntensityDifferenceAverage = tally.relativeIntensityDifferenceAverage;
				t.intensitiesSum = tally.intensitiesSum;
				t.ranks = tally.ranks;

				tallyResultsFinalDecoy.push_back(t);
			}

// #define TS_PSM
#ifdef TS_PSM
			if (
				targetDecoyCandidatePair->peptideStringWithMods() == "MADEAGSEADHEGTHSTK"
				&& targetDecoyCandidatePair->charge() == 2
				) {
					for (const TallyResult &tr : tallyResultsFinalTarget) {
						std::sort(tallyTarget.ranks.begin(), tallyTarget.ranks.end());
						qDebug() << tr.isDecoy << tr.occurrence << tr.ranks << "SDLKJFS";
					}
					for (const TallyResult &tr : tallyResultsFinalDecoy) {
						std::sort(tallyTarget.ranks.begin(), tallyTarget.ranks.end());
						qDebug() << tr.isDecoy << tr.occurrence << tr.ranks << "SDLKJFS";
					}
				}
#endif

			result.push_back({targetDecoyCandidatePair, tallyResultsFinalTarget, tallyResultsFinalDecoy});
			tallyResultsTarget.clear();
			tallyResultsDecoy.clear();
		}

		return {e, result};
	}

	Err indexScanNumberMzIntensities(
		const QVector<ScanNumberMzIntensity*> &scanNumberMzIntensities,
		QVector<QPair<ScanNumberMzIntensity*, QPair<Index, Index>>> *snmiVsIndexes
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(scanNumberMzIntensities); ree;

		QVector<Index> batchIndexes;
		int lastPrecursorMzHashed = -1;
		int lastMs2MzHashed = -1;
		for (int i = 0; i < scanNumberMzIntensities.size(); ++i) {

			ScanNumberMzIntensity *snmi = scanNumberMzIntensities[i];

			const int currentPrecursorMzHashed = MathUtils::hashDecimal(
				snmi->scanInfoPntr->precursorTargetMz,
				precursorMzKeyHashingPrecision
				);
			e = ErrorUtils::isTrue(currentPrecursorMzHashed >= lastPrecursorMzHashed); ree;

			const int currentMs2MzHashed = MathUtils::hashDecimal(
				snmi->mzVal,
				precursorMzKeyHashingPrecision
				);
			if (lastPrecursorMzHashed == currentPrecursorMzHashed) {
				e = ErrorUtils::isTrue(currentMs2MzHashed >= lastMs2MzHashed); ree;
			}

			if (lastMs2MzHashed != currentMs2MzHashed) {
				snmiVsIndexes->push_back({snmi, {batchIndexes.front(), batchIndexes.back()}});
				batchIndexes.clear();
			}

			batchIndexes.push_back(i);

			lastPrecursorMzHashed = currentPrecursorMzHashed;
			lastMs2MzHashed = currentMs2MzHashed;
		}

		ERR_RETURN
	}

	Err indexMS2IonLibraries(
		int libraryPrecursorMzMinLowerBoundIndex,
		int libraryPrecursorMzMaxUpperBoundIndex,
		QVector<MS2IonLibrary> *ms2IonLibraries,
		QHash<MS2IonLibrary*, QVector<MS2IonLibrary*>> *ms2IonLibrariesIndexed
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(*ms2IonLibraries); ree;

		ms2IonLibrariesIndexed->clear();

		int lastHashedPrecursorMz = -1;
		int lastHashedMs2Mz = -1;

		QVector<MS2IonLibrary*> ms2IonLibrariesTemp;
		for (int i = libraryPrecursorMzMinLowerBoundIndex ; i <= libraryPrecursorMzMaxUpperBoundIndex; ++i) {

			MS2IonLibrary &mil = (*ms2IonLibraries)[i];

			constexpr int hashingPrecision = 3;
			const int currentHashedPrecursorMz = MathUtils::hashDecimal(
				mil.targeDecoyCandidatePairPntr->mz(mil.isDecoy),
				hashingPrecision
				);

			const int currentHashedMs2Mz = MathUtils::hashDecimal(
				mil.ms2IonPntr->mz,
				hashingPrecision
				);

			if (
				lastHashedPrecursorMz != currentHashedPrecursorMz
				&& lastHashedMs2Mz != currentHashedMs2Mz
				) {
				ms2IonLibrariesIndexed->insert(ms2IonLibrariesTemp.front(), ms2IonLibrariesTemp);
				ms2IonLibrariesTemp.clear();
			}

			ms2IonLibrariesTemp.push_back(&mil);
			lastHashedPrecursorMz = currentHashedPrecursorMz;
			lastHashedMs2Mz = currentHashedMs2Mz;
		}

		ERR_RETURN
	}

	std::tuple<Err, int, int> getLibraryStartStopIndexes (
		const QVector<ScanNumberMzIntensity*> &scanNumberMzIntensities,
		const PythiaParameters &parameters,
		QVector<MS2IonLibrary> *ms2IonLibraries
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(scanNumberMzIntensities); rtee;
		e = ErrorUtils::isNotEmpty(*ms2IonLibraries); rtee;
		e = ErrorUtils::isTrue(parameters.isValid()); rtee;

		const ScanNumberMzIntensity *snmiFirst = scanNumberMzIntensities.front();
		const ScanNumberMzIntensity *snmiLast = scanNumberMzIntensities.back();

		const float mzTol = MathUtils::calculatePPM(
			snmiFirst->scanInfoPntr->precursorTargetMz,
			static_cast<float>(parameters.ms2ExtractionWidthPPM)
			);

		const float precursorMzValLower = snmiFirst->scanInfoPntr->precursorTargetMz
										- snmiFirst->scanInfoPntr->isoWindowLower
										- mzTol
										- parameters.precursorExtractionWindowThomsons;

		const float precursorMzValUpper = snmiLast->scanInfoPntr->precursorTargetMz
										+ snmiLast->scanInfoPntr->isoWindowUpper
										+ mzTol
										+ parameters.precursorExtractionWindowThomsons;

		const int libraryPrecursorMzMinLowerBoundIndex = std::upper_bound(
			ms2IonLibraries->begin(),
			ms2IonLibraries->end(),
			precursorMzValLower,
			[](float val, const MS2IonLibrary &stct) {return stct.targeDecoyCandidatePairPntr->mz(stct.isDecoy) > val;}
			) - ms2IonLibraries->begin();

		const int libraryPrecursorMzMaxUpperBoundIndex = std::lower_bound(
			ms2IonLibraries->begin(),
			ms2IonLibraries->end(),
			precursorMzValUpper,
			[](const MS2IonLibrary &stct, float val) {return val > stct.targeDecoyCandidatePairPntr->mz(stct.isDecoy);}
			) - ms2IonLibraries->begin() - 1;

		return {e, libraryPrecursorMzMinLowerBoundIndex, libraryPrecursorMzMaxUpperBoundIndex};
	}


	QPair<Err, QVector<TallyResultTuple>> performFraggingLogic(
		const QVector<ScanNumberMzIntensity*> &scanNumberMzIntensities,
		const PythiaParameters &parameters,
		QVector<MS2IonLibrary> *ms2IonLibraries
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(scanNumberMzIntensities); rree;
		e = ErrorUtils::isNotEmpty(*ms2IonLibraries); rree;

		QHash<TargetDecoyCandidatePair*, QVector<IonSearchResult2>> ionSearchResults;
		constexpr int batchSize = 2e5; //TODO make this settable in params
		ionSearchResults.reserve(batchSize);

		const std::tuple<Err, int, int> startStopResult = getLibraryStartStopIndexes(
			scanNumberMzIntensities,
			parameters,
			ms2IonLibraries
			);

		e = std::get<0>(startStopResult); rree;
		const int libraryPrecursorMzMinLowerBoundIndex = std::get<1>(startStopResult);
		const int libraryPrecursorMzMaxUpperBoundIndex = std::get<2>(startStopResult);

		const int size = libraryPrecursorMzMaxUpperBoundIndex - libraryPrecursorMzMinLowerBoundIndex;
		if (size < 1) {
			return {e, {}};
		}

		QHash<MS2IonLibrary*, QVector<MS2IonLibrary*>> ms2IonLibrariesIndexed;
		e = indexMS2IonLibraries(
			libraryPrecursorMzMinLowerBoundIndex,
			libraryPrecursorMzMaxUpperBoundIndex,
			ms2IonLibraries,
			&ms2IonLibrariesIndexed
			); rree;

		QVector<QPair<ScanNumberMzIntensity*, QPair<Index, Index>>> snmiVsIndexes;
		e = indexScanNumberMzIntensities(scanNumberMzIntensities, &snmiVsIndexes); rree;

		QVector<MS2IonLibrary*> ms2IonLibrariesPntrs;
		for (MS2IonLibrary &mil : *ms2IonLibraries) {
			ms2IonLibrariesPntrs.push_back(&mil);
			}

		Ms2IonFraggertronManager fragger;
		e = fragger.init(ms2IonLibrariesPntrs.mid(libraryPrecursorMzMinLowerBoundIndex, size + 1)); rree;

		QVector<QPair<ScanNumberMzIntensity*, QVector<MS2IonLibrary*>>> snmiVsIndexesLibrary;
		for (const QPair<ScanNumberMzIntensity*, QPair<Index, Index>> &pr : snmiVsIndexes) {

			ScanNumberMzIntensity *snmi = pr.first;

			const float mzPrecursorScan = snmi->scanInfoPntr->precursorTargetMz;

			const float mzPrecursorScanTol = MathUtils::calculatePPM(
				mzPrecursorScan,
				static_cast<float>(parameters.ms2ExtractionWidthPPM)
				);

			const float mzPrecursorScanMin = mzPrecursorScan
											- mzPrecursorScanTol
											- snmi->scanInfoPntr->isoWindowLower;

			const float mzPrecursorScanMax = mzPrecursorScan
											+ mzPrecursorScanTol
											+ snmi->scanInfoPntr->isoWindowUpper;

			const float mzMs2Scan = snmi->mzVal;
			const float mzMs2ScanTol = MathUtils::calculatePPM(
				mzMs2Scan,
				static_cast<float>(parameters.ms2ExtractionWidthPPM)
				);
			const float mzMs2ScanMin = mzMs2Scan - mzMs2ScanTol;
			const float mzMs2ScanMax = mzMs2Scan + mzMs2ScanTol;

			QVector<MS2IonLibrary*> tdPeptideFrags;
			e = fragger.extractMs2Points(
				mzPrecursorScanMin,
				mzPrecursorScanMax,
				mzMs2ScanMin,
				mzMs2ScanMax,
				&tdPeptideFrags
				); rree;

			if (!tdPeptideFrags.isEmpty()) {
				snmiVsIndexesLibrary.push_back({snmi, tdPeptideFrags});
			}
		}

		QHash<ScanNumberMzIntensity*, QVector<MS2IonLibrary*>> scanNumberMzIntensityPtrVsMS2IonLibraryPntrs;
		for (const QPair<ScanNumberMzIntensity*, QVector<MS2IonLibrary*>> &pr : snmiVsIndexesLibrary) {
			scanNumberMzIntensityPtrVsMS2IonLibraryPntrs.insert(pr.first, pr.second);
		}

		for (const QPair<ScanNumberMzIntensity*, QPair<Index, Index>> &pr : snmiVsIndexes) {

			ScanNumberMzIntensity *snmi = pr.first;
			const QPair<Index, Index> &snmiVsIndexesRanges = pr.second;

			if (!scanNumberMzIntensityPtrVsMS2IonLibraryPntrs.contains(snmi)) {
				continue;
			}

			const QVector<MS2IonLibrary*> &ms2IonLibrariesSmni = scanNumberMzIntensityPtrVsMS2IonLibraryPntrs[snmi];

			ScanNumberMzIntensity *snmipMin = scanNumberMzIntensities[snmiVsIndexesRanges.first];
			ScanNumberMzIntensity *snmipMax = scanNumberMzIntensities[snmiVsIndexesRanges.second];

			e = ErrorUtils::isTrue(MathUtils::tSame(
				snmipMin->scanInfoPntr->precursorTargetMz,
				snmipMax->scanInfoPntr->precursorTargetMz,
				precursorMzKeyHashingPrecision
				)); rree;

			e = ErrorUtils::isTrue(MathUtils::tSame(
				snmipMin->mzVal,
				snmipMax->mzVal,
				precursorMzKeyHashingPrecision
				)); rree;

			for (int i = snmiVsIndexesRanges.first; i  <= snmiVsIndexesRanges.second; i++) {

				ScanNumberMzIntensity *snmip = scanNumberMzIntensities[i];

				for (MS2IonLibrary *mil : ms2IonLibrariesSmni) {
					TargetDecoyCandidatePair *tdcp = mil->targeDecoyCandidatePairPntr;

					IonSearchResult2 isr2;
					isr2.ms2IonLibraryPntr = mil;
					isr2.msScanPointPntr = snmip;

					ionSearchResults[tdcp].push_back(isr2);
				}
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

		QPair<Err, QVector<TallyResultTuple>> result
				= collateScanNumberVsOccurrencesTargetDecoyCandidatePairs(ionSearchResults); rree;

		return result;
	}

}//namespace
QPair<Err, QVector<TallyResultTuple>> MsFraggertron::processTargetDecoyCandidatePairsPntrsTranch(
	const QVector<TargetDecoyCandidatePair *> &tdcps
	) {
	ERR_INIT

	e = ErrorUtils::isNotEmpty(tdcps); rree;
	e = ErrorUtils::isNotEmpty(m_scanNumberMzIntensityTranched); rree;

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

	parallelSortMS2IonLibrariesMzAsc(
		m_scanNumberMzIntensityTranched.size(),
		&ms2IonLibraries
		);

	QVector<TallyResultTuple> tallyResultsFinal;

#define FRAG_PARALLEL
#ifdef FRAG_PARALLEL
	const auto binderLogic = std::bind(
		performFraggingLogic,
		std::placeholders::_1,
		m_parameters,
		&ms2IonLibraries
		);

	QFuture<QPair<Err, QVector<TallyResultTuple>>> future = QtConcurrent::mapped(
		m_scanNumberMzIntensityTranched,
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
	for (const QVector<ScanNumberMzIntensity*> &pgs : m_scanNumberMzIntensityTranched) {
		auto res = performFraggingLogic(
			pgs,
			m_parameters,
			&ms2IonLibraries
			);
		e = res.first; rree;
	}
#endif

	return {e, tallyResultsFinal};
}
