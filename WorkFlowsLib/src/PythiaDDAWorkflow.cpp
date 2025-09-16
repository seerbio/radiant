//
// Created by andrewnichols on 6/29/25.
//

#include "PythiaDDAWorkflow.h"

#include "DiscriminantScoretron.h"
#include "EigenUtils.h"
#include "FDRCLassifierNeuralNet.h"
#include "Ms2IonFraggertronManager.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderPointerAcc.h"
#include "ObjectCSVWriters.h"
#include "ParallelUtils.h"
#include "QValueSettertron.h"
#include "PythiaDIAFFWorkflowAlgos/PythiaDIAFFWorkflowSharedMethods.h"

Err PythiaDDAWorkflow::init(
	const PythiaParameters &parameters,
	const QString& libraryFilePath
	) {

	ERR_INIT

	e = ErrorUtils::isTrue(parameters.isValid()); ree;
	e = ErrorUtils::fileExists(libraryFilePath); ree;

	m_parameters = parameters;
	m_parameters.ms2ExtractionWidthPPM = 20;
	m_parameters.threadCount = 64; //TODO use higher threadcount to load balans
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

namespace {

	Err serialFilterByValue(
		double threshold,
		QVector<CandidateScoresDDA*> *candidateScoreses
		) {

		ERR_INIT

		e = ErrorUtils::isFalse(candidateScoreses->isEmpty()); ree;

		std::sort(
			candidateScoreses->begin(),
			candidateScoreses->end(),
			[](const CandidateScoresDDA *l, const CandidateScoresDDA *r){
				if (MathUtils::tSame(l->classifierScore, r->classifierScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
					return l->discriminantScore > r->discriminantScore;
				}
				return l->classifierScore < r->classifierScore;
			});

		int counter = 0;
		for (const CandidateScoresDDA *csp : *candidateScoreses) {

			counter++;
			if (csp->qValue >= threshold && !csp->isDecoy) {
				break;
			}
		}

		candidateScoreses->resize(counter);

		ERR_RETURN
	}


	Err convertCandidateScoresDDAToCandidateScoresReaderRow(
		const QVector<CandidateScoresDDA*> &toReport,
		QVector<CandidateScoresReaderRow> *candidateScoresReaderRows
		) {
		ERR_INIT

		e = ErrorUtils::isNotEmpty(toReport); ree;

		candidateScoresReaderRows->clear();
		candidateScoresReaderRows->reserve(toReport.size());

		for (const CandidateScoresDDA *csp : toReport) {
			CandidateScoresReaderRow csrr;
			csrr.peptideStringWithMods = csp->targetDecoyCandidatePair->peptideStringWithMods();
			csrr.isDecoy = csp->isDecoy;
			csrr.discriminantScore = csp->discriminantScore;
			csrr.classifierScore = csp->classifierScore;
			csrr.charge = csp->targetDecoyCandidatePair->charge();
			csrr.cosineSimSum100 = csp->featuresArray[Occurrences];

			const QVector<MS2Ion> ms2Ions = csrr.isDecoy
											? csp->targetDecoyCandidatePair->ms2IonsTarget()
											: csp->targetDecoyCandidatePair->ms2IonsDecoy();

			for (const MS2Ion &ms2Ion : ms2Ions) {
				constexpr int maxRank = 12;

				if (ms2Ion.rank >= maxRank) {
					break;
				}

				switch (ms2Ion.rank) {
				case 0:
					csrr.mzFoundMean1 = ms2Ion.mz;
					csrr.intensityFoundMax1 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 1:
					csrr.mzFoundMean2 = ms2Ion.mz;
					csrr.intensityFoundMax2 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 2:
					csrr.mzFoundMean3 = ms2Ion.mz;
					csrr.intensityFoundMax3 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 3:
					csrr.mzFoundMean4 = ms2Ion.mz;
					csrr.intensityFoundMax4 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 4:
					csrr.mzFoundMean5 = ms2Ion.mz;
					csrr.intensityFoundMax5 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 5:
					csrr.mzFoundMean6 = ms2Ion.mz;
					csrr.intensityFoundMax6 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 6:
					csrr.mzFoundMean7 = ms2Ion.mz;
					csrr.intensityFoundMax7 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 7:
					csrr.mzFoundMean8 = ms2Ion.mz;
					csrr.intensityFoundMax8 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 8:
					csrr.mzFoundMean9 = ms2Ion.mz;
					csrr.intensityFoundMax9 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 9:
					csrr.mzFoundMean10 = ms2Ion.mz;
					csrr.intensityFoundMax10 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 10:
					csrr.mzFoundMean11 = ms2Ion.mz;
					csrr.intensityFoundMax11 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				case 11:
					csrr.mzFoundMean12 = ms2Ion.mz;
					csrr.intensityFoundMax12 = csp->featuresArray[IntensityRank0 + ms2Ion.rank];
					break;
				default:
					break;
				}
			}

			candidateScoresReaderRows->push_back(csrr);
		}

		ERR_RETURN
	}


}//namespace
Err PythiaDDAWorkflow::processFile(const QString &msDataFilePath) {

	ERR_INIT

	MsReaderPointerAcc msReaderPtr;
	e = msReaderPtr.openFile(msDataFilePath); ree;

	e = m_msFraggertron.init(
		m_parameters,
		&msReaderPtr,
		&m_msCalibratomatic
		); ree;

	e = buildMsCalibratomatic(); ree;

	e = optimizeExtractionPPM(); ree;

	QVector<CandidateScoresDDA*> candidateScoresesPntrs;
	e = processAll(&candidateScoresesPntrs); ree;

	QMap<int, int> fdrVsCounts;
	e = FDRCLassifierNeuralNet::outputFDRResults(
		candidateScoresesPntrs,
		1,
		&fdrVsCounts
		); ree;
	const int targetCount1PercentFDRLDA = fdrVsCounts[1];

	QVector<CandidateScoresDDA*> candidateScoreClassifier;
	e = applyNeuralNetClassifier(
		candidateScoresesPntrs,
		S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST,
		&candidateScoreClassifier
		); ree;

	constexpr double fdrQValThreshold = 0.01;
	e = serialFilterByValue(
		fdrQValThreshold,
		&candidateScoreClassifier
		); ree;

	candidateScoresesPntrs.resize(targetCount1PercentFDRLDA);
	const QVector<CandidateScoresDDA*> toReport = candidateScoreClassifier.size() > candidateScoresesPntrs.size()
											? candidateScoreClassifier
											: candidateScoresesPntrs;

	qDebug()
	<< qPrintable(S_GLOBAL_TIMER.elapsed())
	<< "Candidates found by LDA:"
	<< candidateScoresesPntrs.size()
	<< "Candidates found by NN:"
	<< candidateScoreClassifier.size()
	;

	QVector<CandidateScoresReaderRow> candidateScoresReaderRows;
	e = convertCandidateScoresDDAToCandidateScoresReaderRow(
		toReport,
		&candidateScoresReaderRows
		); ree;

	QString resultsFilePath = msReaderPtr.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DDA_FILE_EXTENSION;
	e = ParquetReader::write(candidateScoresReaderRows, resultsFilePath); ree;

	ERR_RETURN
}

namespace {

	constexpr int skipCount = 4;

	QPair<Err, CandidateScoresDDA*> selectBestCandidate(
		const QVector<FeaturesDDA> &featureArray,
		const Eigen::VectorX<float> &weights,
		MsCalibratomatic *msCalibratomatic,
		QVector<CandidateScoresDDA> *csDDAs
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(*csDDAs); rree;

		if (msCalibratomatic->isInitRT()) {
			for (CandidateScoresDDA &cs : *csDDAs) {

				e = msCalibratomatic->predictScanTime(
					cs.targetDecoyCandidatePair->iRt(),
					&cs.scanTimePredicted
					); rree;

				cs.featuresArray[ScanTimeDeltaDDA] = cs.scanTime - cs.scanTimePredicted;
				cs.featuresArray[ScanTimeDeltaDDAAbs] = std::abs(cs.featuresArray[ScanTimeDeltaDDA]);
				cs.featuresArray[ScanTimeDeltaDDAFromMean] = cs.featuresArray[ScanTimeDeltaDDA] - msCalibratomatic->scanTimeMean();
				cs.featuresArray[ScanTimeDeltaDDAFromMeanAbs] = std::abs(cs.featuresArray[ScanTimeDeltaDDAFromMean]);
			}
		}

		CandidateScoresDDA *bestCandidateScoresDDA;
		float bestScore = std::numeric_limits<float>::lowest();
		for (int row = 0; row < csDDAs->size(); row++) {

			const Eigen::VectorX<float> v = EigenUtils::convertQVectorToEigenVector(
				CandidateScoresDDA::selectFeaturesArrayFeatures(
					(*csDDAs)[row].featuresArray,
					featureArray
					)
				);

			const float score = v.dot(weights);
			if (score > bestScore) {
				bestCandidateScoresDDA = &(*csDDAs)[row];
				bestScore = score;
			}
		}

		return {e, bestCandidateScoresDDA};
	}


	Err buildLdaInputData(
		QVector<CandidateScoresDDATuple> &candidateScoresDDATuples,
		const QVector<FeaturesDDA> &features,
		const QVector<float> &weights,
		MsCalibratomatic *msCalibratomatic,
		QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *targetDecoyPairScores,
		QVector<QPair<CandidateScoresDDA*, CandidateScoresDDA*>> *targetDecoyPairCandidateScorePntrs
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(candidateScoresDDATuples); ree;
		targetDecoyPairScores->clear();

		const Eigen::VectorX<float> weightsEigen = EigenUtils::convertQVectorToEigenVector(weights);

		for (CandidateScoresDDATuple &pr : candidateScoresDDATuples) {

			QVector<CandidateScoresDDA> *scoresTarget = &std::get<1>(pr);
			const QPair<Err, CandidateScoresDDA*> bestTargetResult = selectBestCandidate(
				features,
				weightsEigen,
				msCalibratomatic,
				scoresTarget
				); ree;
			e = bestTargetResult.first; ree;
			CandidateScoresDDA *bestTarget = bestTargetResult.second;

			QVector<CandidateScoresDDA> *scoresDecoy = &std::get<2>(pr);
			const QPair<Err, CandidateScoresDDA*> bestDecoyResult = selectBestCandidate(
				features,
				weightsEigen,
				msCalibratomatic,
				scoresDecoy
				); ree;
			e = bestDecoyResult.first; ree;
			CandidateScoresDDA *bestDecoy = bestDecoyResult.second;

			const QVector<float> featuresTarget = bestTarget->selectFeaturesArrayFeatures(features);
			const QVector<float> featuresDecoy = bestDecoy->selectFeaturesArrayFeatures(features);

			targetDecoyPairScores->push_back({featuresTarget, featuresDecoy});
			targetDecoyPairCandidateScorePntrs->push_back({bestTarget, bestDecoy});
		}

		ERR_RETURN
	}

	Err setCandidateScoresDiscriminantScore(
		const QVector<float> &discriminantScores,
		QVector<CandidateScoresDDA*> *candidateScoresesPntrs
		) {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(discriminantScores); ree;
		e = ErrorUtils::isEqual(discriminantScores.size(), candidateScoresesPntrs->size()); ree;

		for (int i = 0; i < discriminantScores.size(); i++) {
			CandidateScoresDDA* cs = (*candidateScoresesPntrs)[i];
			cs->discriminantScore = discriminantScores.at(i);
		}

		ERR_RETURN
	}

}//namespace
Err PythiaDDAWorkflow::buildMsCalibratomatic() {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairsPntrs); ree;

	QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsCalibration;
	for (int i = 0; i < m_targetDecoyCandidatePairsPntrs.size(); i++) {
		if (i % skipCount != 0 ) {
			continue;
		}
		targetDecoyCandidatePairsCalibration.push_back(m_targetDecoyCandidatePairsPntrs[i]);
	}

	QPair<Err, QVector<CandidateScoresDDATuple>> result
		= m_msFraggertron.performFragging(targetDecoyCandidatePairsCalibration);
	e = result.first; ree;

	QVector<CandidateScoresDDATuple> &tallyResults = result.second;

	const QVector<FeaturesDDA> features = DiscriminantScoretron::featuresCalibrationDDA();
	m_weights = DiscriminantScoretron::defaultWeights(features);

	int bestCount = 0;
	for (int i = 0; i < 10; i++) {

		QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> targetDecoyPairScores;
		QVector<QPair<CandidateScoresDDA*, CandidateScoresDDA*>> targetDecoyPairCandidateScorePntrs;
		e = buildLdaInputData(
			tallyResults,
			features,
			m_weights,
			&m_msCalibratomatic,
			&targetDecoyPairScores,
			&targetDecoyPairCandidateScorePntrs
			); ree;

		QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> targetDecoyPairScoresPntrs;
		targetDecoyPairScoresPntrs.reserve(targetDecoyPairScores.size());
		for (QPair<FeaturesArrayTargets, FeaturesArrayDecoys> &pr : targetDecoyPairScores) {
			targetDecoyPairScoresPntrs.push_back({&pr.first, &pr.second});
		}

		// QVector<float> trainedWeights;
		e = DiscriminantScoretron::trainLDAClassifier(
			targetDecoyPairScoresPntrs,
			false,
			&m_weights
			); ree;

		QVector<FeaturesArray> featuresArray;
		QVector<CandidateScoresDDA*> candidateScoresesPntrs;
		for (int i = 0; i < targetDecoyPairCandidateScorePntrs.size(); i++) {

			const QPair<CandidateScoresDDA*, CandidateScoresDDA*> &csPntrs = targetDecoyPairCandidateScorePntrs[i];

			featuresArray.append({
				csPntrs.first->selectFeaturesArrayFeatures(features),
				csPntrs.second->selectFeaturesArrayFeatures(features)
				});

			e = ErrorUtils::isNotEqual(csPntrs.first->isDecoy, csPntrs.second->isDecoy); ree;
			candidateScoresesPntrs.append({csPntrs.first, csPntrs.second});
		}

		QVector<FeaturesArray*> featuresArrayPntrs;
		std::transform(
			featuresArray.begin(),
			featuresArray.end(),
			std::back_inserter(featuresArrayPntrs),
			[](FeaturesArray &fa){return &fa;}
			);

		QVector<float> discriminantScores;
		e = DiscriminantScoretron::applyWeights(
			m_weights,
			m_parameters.threadCount,
			featuresArrayPntrs,
			&discriminantScores
			); ree;

		e = setCandidateScoresDiscriminantScore(
			discriminantScores,
			&candidateScoresesPntrs
			); ree;

		std::sort(
			candidateScoresesPntrs.rbegin(),
			candidateScoresesPntrs.rend(),
			[](const CandidateScoresDDA *l, const CandidateScoresDDA *r) {
				return l->discriminantScore < r->discriminantScore;
			});

		//Need to speed this up
		e = QValueSettertron::setQValueForCandidates(
			QValueSettertron::QValueScoreType::DiscriminantScore,
			&targetDecoyPairCandidateScorePntrs
			); ree;

		QMap<int, int> fdrVsCounts;
		e = FDRCLassifierNeuralNet::outputFDRResults(
			candidateScoresesPntrs,
			1,
			&fdrVsCounts
			); ree;

		constexpr int fdrKey = 1;
		constexpr int fdrKeyMassCalMS2 = 1;
		constexpr int fdrKeyMassCalMS1 = 1;

		if (bestCount >= fdrVsCounts[fdrKey]) {
			break;
		}
		
		e = honeIRTAndMassCalibration(
			&candidateScoresesPntrs,
			fdrVsCounts[fdrKey],
			fdrVsCounts[fdrKey]
			); ree;

		bestCount = fdrVsCounts[fdrKey];
	}

	ERR_RETURN
}

Err PythiaDDAWorkflow::honeIRTAndMassCalibration(
	QVector<CandidateScoresDDA *> *candidateScoresesPntrs,
	int topNCandidates,
	int topCandidatesMass
	) {

	ERR_INIT

    e = ErrorUtils::isFalse(candidateScoresesPntrs->isEmpty()); ree;

	std::sort(
		candidateScoresesPntrs->rbegin(),
		candidateScoresesPntrs->rend(),
		[](const CandidateScoresDDA *l, const CandidateScoresDDA *r) {
			return l->discriminantScore < r->discriminantScore;
		});

    QVector<CandidateScoresDDA*> candidateScoresVecBatchPntrsResized = *candidateScoresesPntrs;

    constexpr int minTrainingCountTranche = 50;
    candidateScoresVecBatchPntrsResized.resize(std::max(minTrainingCountTranche, topNCandidates));

    if (m_parameters.verbosity > 0) {
        qDebug() << "Using" << candidateScoresVecBatchPntrsResized.size() << "for iRT Estimation";
        qDebug() << candidateScoresVecBatchPntrsResized.size() << "after filtering";
    }

    if (candidateScoresVecBatchPntrsResized.isEmpty()) {
        ERR_RETURN
    }

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
            MSLevelEnum::MS2,
            candidateScoresVecBatchPntrsResized,
            m_parameters.verbosity,
            &msCalibrationReaderRows
            ); ree;

    // for (const CandidateScores *cs : candidateScoresVecBatchPntrsResized) {
    //
    //     if (m_entered.value(cs->targetDecoyCandidatePair) ) {
    //         continue;
    //     }
    //
    //     if (m_excludeDecoys) {
    //         if (cs->isDecoy) {
    //             continue;
    //         }
    //     }
    //
    //     m_targetDecoyCandidatePairsTopScores.push_back(cs->targetDecoyCandidatePair);
    //     m_entered.insert(cs->targetDecoyCandidatePair, true);
    // }

    e = m_msCalibratomatic.buildRTMapper(msCalibrationReaderRows); ree;

	if (m_parameters.verbosity > 0) {
		qDebug()
		<< qPrintable(S_GLOBAL_TIMER.elapsed())
		<< " scanTimeWindowStDev x"
		<< m_parameters.scanTimeWindowStDevs
		<<":"
		<< m_msCalibratomatic.scanTimeStDev(m_parameters.scanTimeWindowStDevs);
	}

    // constexpr int ms2MassRecalCountMin = 200;
    // if (topCandidatesMass > ms2MassRecalCountMin) {
    //
    //     msCalibrationReaderRows.resize(topCandidatesMass);
    //
    //     e = m_msCalibratomatic.setMassCalibrationCoeffs(
    //         msCalibrationReaderRows,
    //         MSLevelEnum::MS2
    //         ); ree;
    // }

    ERR_RETURN
}

namespace {

	Err getTopFrequencyParameters(
			const QVector<QPair<float, int>> &results,
			int verbosity,
			double *ppmSetting
			) {

		ERR_INIT
		e = ErrorUtils::isNotEmpty(results); ree;

		Eigen::MatrixX<double> xyMat(results.size() + 1, 2);
		xyMat.setZero();
		for (int row = 0; row < results.size(); row++) {
			const QPair<float, int> &doeResult = results.at(row);
			xyMat.coeffRef(row + 1, 0) = doeResult.first;
			xyMat.coeffRef(row + 1, 1) = static_cast<double>(doeResult.second);
			if (verbosity > 0) {
				qDebug() << doeResult.first << doeResult.second << xyMat.coeff(row, 1);

			}
		}

		constexpr int polynomialOrder = 2;

		QVector<double> coeffs;
		EigenUtils::fitPolynomialQRDecomposition(xyMat, polynomialOrder, &coeffs);

		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PPM Coeffs" << coeffs;

		QVector<double> xPoints = {results.first().first};
		while (xPoints.back() < results.back().first) {
			constexpr double incrementVal = 0.25;
			xPoints.push_back(xPoints.back() + incrementVal);
		}

		QVector<double> yPoints;
		for (double x : xPoints) {
			double y = 0.0;
			for (int i = 0; i < coeffs.size(); i++) {
				y += coeffs.at(i) * std::pow(x, i);
			}
			if (verbosity > 0) {
				qDebug() << x << y;

			}
			yPoints.push_back(y);
		}

		*ppmSetting = xPoints.at(std::max_element(yPoints.begin(), yPoints.end()) - yPoints.begin());

		ERR_RETURN
	}
}//namespace
Err PythiaDDAWorkflow::optimizeExtractionPPM() {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairsPntrs); ree;

	QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsCalibration;
	for (int i = 0; i < m_targetDecoyCandidatePairsPntrs.size(); i++) {
		if (i % skipCount != 0 ) {
			continue;
		}
		targetDecoyCandidatePairsCalibration.push_back(m_targetDecoyCandidatePairsPntrs[i]);
	}

	int bestCount = 0;
	int lesserCountsCount = 0;
	QVector<QPair<float,int>> allCounts;

	QVector<float> bestWeights = m_weights;

	const QVector<float> ppms = {3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 35, 40, 45, 50};
	for (float ppm : ppms) {

		PythiaParameters paramsPPM = m_parameters;
		paramsPPM.ms2ExtractionWidthPPM = ppm;
		m_msFraggertron.setPythiaParameters(paramsPPM);

		QPair<Err, QVector<CandidateScoresDDATuple>> result
			= m_msFraggertron.performFragging(targetDecoyCandidatePairsCalibration);
		e = result.first; ree;

		QVector<CandidateScoresDDATuple> &tallyResults = result.second;

		const QVector<FeaturesDDA> features = DiscriminantScoretron::featuresCalibrationDDA();

		QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> targetDecoyPairScores;
		QVector<QPair<CandidateScoresDDA*, CandidateScoresDDA*>> targetDecoyPairCandidateScorePntrs;
		e = buildLdaInputData(
			tallyResults,
			features,
			bestWeights,
			&m_msCalibratomatic,
			&targetDecoyPairScores,
			&targetDecoyPairCandidateScorePntrs
			); ree;

		QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> targetDecoyPairScoresPntrs;
		targetDecoyPairScoresPntrs.reserve(targetDecoyPairScores.size());
		for (QPair<FeaturesArrayTargets, FeaturesArrayDecoys> &pr : targetDecoyPairScores) {
			targetDecoyPairScoresPntrs.push_back({&pr.first, &pr.second});
		}

		QVector<float> trainedWeights;
		e = DiscriminantScoretron::trainLDAClassifier(
			targetDecoyPairScoresPntrs,
			false,
			&trainedWeights
			); ree;

		QVector<FeaturesArray> featuresArray;
		QVector<CandidateScoresDDA*> candidateScoresesPntrs;
		for (int i = 0; i < targetDecoyPairCandidateScorePntrs.size(); i++) {

			const QPair<CandidateScoresDDA*, CandidateScoresDDA*> &csPntrs = targetDecoyPairCandidateScorePntrs[i];

			featuresArray.append({
				csPntrs.first->selectFeaturesArrayFeatures(features),
				csPntrs.second->selectFeaturesArrayFeatures(features)
				});

			e = ErrorUtils::isNotEqual(csPntrs.first->isDecoy, csPntrs.second->isDecoy); ree;
			candidateScoresesPntrs.append({csPntrs.first, csPntrs.second});
		}

		QVector<FeaturesArray*> featuresArrayPntrs;
		std::transform(
			featuresArray.begin(),
			featuresArray.end(),
			std::back_inserter(featuresArrayPntrs),
			[](FeaturesArray &fa){return &fa;}
			);

		QVector<float> discriminantScores;
		e = DiscriminantScoretron::applyWeights(
			m_weights,
			m_parameters.threadCount,
			featuresArrayPntrs,
			&discriminantScores
			); ree;

		e = setCandidateScoresDiscriminantScore(
			discriminantScores,
			&candidateScoresesPntrs
			); ree;

		std::sort(
			candidateScoresesPntrs.rbegin(),
			candidateScoresesPntrs.rend(),
			[](const CandidateScoresDDA *l, const CandidateScoresDDA *r) {
				return l->discriminantScore < r->discriminantScore;
			});

		//Need to speed this up
		e = QValueSettertron::setQValueForCandidates(
			QValueSettertron::QValueScoreType::DiscriminantScore,
			&targetDecoyPairCandidateScorePntrs
			); ree;

		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())<< "Testing optimal PPM" << ppm << "ppm";
		QMap<int, int> fdrVsCounts;
		e = FDRCLassifierNeuralNet::outputFDRResults(
			candidateScoresesPntrs,
			1,
			&fdrVsCounts
			); ree;

		constexpr int fdrKey = 1;
		constexpr int fdrKeyMassCalMS2 = 1;
		constexpr int fdrKeyMassCalMS1 = 1;

		const int ppmCount = fdrVsCounts[fdrKey];
		allCounts.push_back({ppm, ppmCount});

		if (ppmCount > bestCount) {
			bestCount = ppmCount;
			lesserCountsCount = 0;
			bestWeights = trainedWeights;

		} else {
			lesserCountsCount++;
		}

		constexpr int lesserCountsCountLimit = 3;
		if (lesserCountsCount >= lesserCountsCountLimit) {
			break;
		}
	}

	m_weights = bestWeights;
	e = getTopFrequencyParameters(
		allCounts, -1,
		&m_parameters.ms2ExtractionWidthPPM
		); ree;
	m_msFraggertron.setPythiaParameters(m_parameters);

	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Optimal PPM" << m_parameters.ms2ExtractionWidthPPM << "ppm";

	ERR_RETURN
}

Err PythiaDDAWorkflow::processAll(QVector<CandidateScoresDDA*> *candidateScoresesPntrs) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairsPntrs); ree;

	QPair<Err, QVector<CandidateScoresDDATuple>> result
			= m_msFraggertron.performFragging(m_targetDecoyCandidatePairsPntrs);
	e = result.first; ree;

	m_tallyResultsProcessAll = result.second;

	const QVector<FeaturesDDA> features = DiscriminantScoretron::featuresCalibrationDDA();

	QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> targetDecoyPairScores;
	m_targetDecoyPairCandidateScorePntrs.clear();
	e = buildLdaInputData(
		m_tallyResultsProcessAll,
		features,
		m_weights,
		&m_msCalibratomatic,
		&targetDecoyPairScores,
		&m_targetDecoyPairCandidateScorePntrs
		); ree;

	QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> targetDecoyPairScoresPntrs;
	targetDecoyPairScoresPntrs.reserve(targetDecoyPairScores.size());
	for (QPair<FeaturesArrayTargets, FeaturesArrayDecoys> &pr : targetDecoyPairScores) {
		targetDecoyPairScoresPntrs.push_back({&pr.first, &pr.second});
	}

	// QVector<float> trainedWeights;
	e = DiscriminantScoretron::trainLDAClassifier(
		targetDecoyPairScoresPntrs,
		false,
		&m_weights
		); ree;

	QVector<FeaturesArray> featuresArray;
	for (int i = 0; i < m_targetDecoyPairCandidateScorePntrs.size(); i++) {

		const QPair<CandidateScoresDDA*, CandidateScoresDDA*> &csPntrs = m_targetDecoyPairCandidateScorePntrs[i];

		featuresArray.append({
			csPntrs.first->selectFeaturesArrayFeatures(features),
			csPntrs.second->selectFeaturesArrayFeatures(features)
			});

		e = ErrorUtils::isNotEqual(csPntrs.first->isDecoy, csPntrs.second->isDecoy); ree;
		candidateScoresesPntrs->append({csPntrs.first, csPntrs.second});
	}

	QVector<FeaturesArray*> featuresArrayPntrs;
	std::transform(
		featuresArray.begin(),
		featuresArray.end(),
		std::back_inserter(featuresArrayPntrs),
		[](FeaturesArray &fa){return &fa;}
		);

	QVector<float> discriminantScores;
	e = DiscriminantScoretron::applyWeights(
		m_weights,
		m_parameters.threadCount,
		featuresArrayPntrs,
		&discriminantScores
		); ree;

	e = setCandidateScoresDiscriminantScore(
		discriminantScores,
		candidateScoresesPntrs
		); ree;

	std::sort(
		candidateScoresesPntrs->rbegin(),
		candidateScoresesPntrs->rend(),
		[](const CandidateScoresDDA *l, const CandidateScoresDDA *r) {
			return l->discriminantScore < r->discriminantScore;
		});

	//Need to speed this up
	e = QValueSettertron::setQValueForCandidates(
		QValueSettertron::QValueScoreType::DiscriminantScore,
		&m_targetDecoyPairCandidateScorePntrs
		); ree;

	const bool sortedScoreDiscScoreDscnd = std::is_sorted(
		candidateScoresesPntrs->rbegin(),
		candidateScoresesPntrs->rend(),
		[](CandidateScoresDDA *l, CandidateScoresDDA *r) {return l->discriminantScore < r->discriminantScore;}
		);
	e = ErrorUtils::isTrue(sortedScoreDiscScoreDscnd); ree;

	ERR_RETURN
}

namespace {

	Err filterScoredCandidatesForNeuralNet(
			int minMs2FragCount,
			QVector<CandidateScoresDDA*> *candidateScoresTargetsAndDecoys
			) {

		ERR_INIT

		e = ErrorUtils::isFalse(candidateScoresTargetsAndDecoys->isEmpty()); ree;

		const auto terminatorLogic = [minMs2FragCount](CandidateScoresDDA *cs) {
			return cs->featuresArray[Occurrences] < static_cast<float>(minMs2FragCount);
		};
		const auto terminator = std::remove_if(
			candidateScoresTargetsAndDecoys->begin(),
			candidateScoresTargetsAndDecoys->end(),
			terminatorLogic
			);
		candidateScoresTargetsAndDecoys->erase(terminator, candidateScoresTargetsAndDecoys->end());

		std::sort(
			candidateScoresTargetsAndDecoys->rbegin(),
			candidateScoresTargetsAndDecoys->rend(),
			[](const CandidateScoresDDA *l, const CandidateScoresDDA *r) {
				return l->discriminantScore < r->discriminantScore;
			});


		int counter = 0;
		for (const CandidateScoresDDA *csp : *candidateScoresTargetsAndDecoys) {
			counter++;
			if (constexpr double fdrTrainingThreshold = 0.85; csp->qValue >= fdrTrainingThreshold && !csp->isDecoy) {
				break;
			}
		}

		candidateScoresTargetsAndDecoys->resize(counter);

		ERR_RETURN
	}

    Err minMaxScaleScores(
            const QVector<KarnnNNTargetDDA> &karnnNNTargets,
            QVector<KarnnNNTargetDDA> *karnnNNTargetsNorm
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(karnnNNTargets); ree;

        QVector<QVector<float>> vecs;
        std::transform(
                karnnNNTargets.begin(),
                karnnNNTargets.end(),
                std::back_inserter(vecs),
                [](const KarnnNNTargetDDA &kt){return kt.scoreVecNormalized;}
                );

        Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(vecs);
        EigenUtils::minMaxScaleMatrix(&mat);
        const QVector<QVector<float>> vecsNorm = EigenUtils::convertEigenMatrixToQVectors(mat);

        e = ErrorUtils::isEqual(vecsNorm.size(), karnnNNTargets.size()); ree;

        for (int i = 0; i < vecsNorm.size(); i++) {
            KarnnNNTargetDDA ktNew = karnnNNTargets.at(i);
            ktNew.scoreVecNormalized = vecsNorm.at(i);
            karnnNNTargetsNorm->push_back(ktNew);
        }

        ERR_RETURN
    }

	Err buildKarnnNNTargetsNormalized(
		const QVector<FeaturesDDA> &neuralNetFeatures,
		const QVector<CandidateScoresDDA*> &candidateScoresTargetsAndDecoysFDRFiltered,
		QVector<KarnnNNTargetDDA> *karnnNNTargetsNorm
		){

		ERR_INIT

		e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoysFDRFiltered); ree;

		QVector<KarnnNNTargetDDA> karnnNNTargets;
		karnnNNTargets.reserve(candidateScoresTargetsAndDecoysFDRFiltered.size());
		for (int i = 0; i < candidateScoresTargetsAndDecoysFDRFiltered.size(); i++) {
			CandidateScoresDDA *cs = candidateScoresTargetsAndDecoysFDRFiltered.at(i);
			KarnnNNTargetDDA karnnNnTarget;
			karnnNnTarget.candidateScoresDDA = cs;


			karnnNnTarget.scoreVecNormalized = CandidateScoresDDA::selectFeaturesArrayFeatures(
				cs->featuresArray,
				neuralNetFeatures
				);

			karnnNNTargets.push_back(karnnNnTarget);
		}

		e = minMaxScaleScores(karnnNNTargets, karnnNNTargetsNorm); ree;

		ERR_RETURN
	}
    Err trainNeuralNetwork(
            const PythiaParameters &pythiaParameters,
            const QVector<KarnnNNTargetDDA> &karnnNNTargetsNorm,
            int seed,
            FDRCLassifierNeuralNet *fdrcLassifierNeuralNet
            ) {

        ERR_INIT

        constexpr int batchSizeMin = 500;
        const int batchSize = std::min(batchSizeMin, std::max(1, static_cast<int>(karnnNNTargetsNorm.size() / 100.0)));

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Batch Size:" << batchSize;

        QVector<QVector<float>> xData;
        QVector<float> yData;
        for (const KarnnNNTargetDDA &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVecNormalized);
            yData.push_back(kt.candidateScoresDDA->isDecoy ? 1.0 : 0.0);
        }

        e = fdrcLassifierNeuralNet->init(
                pythiaParameters.epochs,
                pythiaParameters.baggingSize,
                batchSize,
                pythiaParameters.learningRate,
                pythiaParameters.nodesFraction,
                pythiaParameters.threadCount
        ); ree;

        e = fdrcLassifierNeuralNet->trainClassifier(
                xData,
                yData,
                seed,
                pythiaParameters.verbosity
                ); ree;

        ERR_RETURN
    }

    QPair<Err, QVector<float>> predictClassiferScoresLogic(
        const QVector<KarnnNNTargetDDA> &karnnNNTargetsNorm,
        FDRCLassifierNeuralNet *fdrcLassifierNeuralNet
        ) {

        ERR_INIT

        QVector<float> predictions;

        QVector<QVector<float>> xData;
        xData.reserve(karnnNNTargetsNorm.size());
        QVector<float> yData;
        yData.reserve(karnnNNTargetsNorm.size());
        for (const KarnnNNTargetDDA &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVecNormalized);
            yData.push_back(kt.candidateScoresDDA->isDecoy ? 1.0 : 0.0);
        }

        e = fdrcLassifierNeuralNet->predictBaggedClassifiers(
            xData,
            &predictions
            ); rree;

        return {e, predictions};
    }

    Err predictClassifierScores(
        const QVector<KarnnNNTargetDDA> &karnnNNTargetsNorm,
        FDRCLassifierNeuralNet *fdrcLassifierNeuralNet,
        QVector<float> *predictions
        ) {

        ERR_INIT

        const auto applyWeightsBinder = std::bind(
            predictClassiferScoresLogic,
            std::placeholders::_1,
            fdrcLassifierNeuralNet
            );

        QVector<QVector<KarnnNNTargetDDA>> karnnNNTargetsNormVecsTranched;
        e = ParallelUtils::trancheVectorForParallelizationInOrder(
            karnnNNTargetsNorm,
            ParallelUtils::numberOfAvailableSystemProcessors(),
            0,
            &karnnNNTargetsNormVecsTranched
            ); ree;

        QFuture<QPair<Err, QVector<float>>> future = QtConcurrent::mapped(
            karnnNNTargetsNormVecsTranched,
            applyWeightsBinder
            );
        future.waitForFinished();

        for (const QPair<Err, QVector<float>> &pr : future) {
            e = pr.first; ree;
            predictions->append(pr.second);
        }

        ERR_RETURN
    }

    Err processPredictions(
            const QVector<float> &predictions,
            QVector<KarnnNNTargetDDA> *karnnNNTargetsNorm
            ) {

        ERR_INIT

        e = ErrorUtils::isFalse(karnnNNTargetsNorm->isEmpty()); ree;

        for (int i = 0; i < predictions.size(); i++) {
            (*karnnNNTargetsNorm)[i].candidateScoresDDA->classifierScore = predictions.at(i);
        }

        ERR_RETURN
    }

    Err subsetKarnnNNTargetsForTraining(
        const QVector<KarnnNNTargetDDA> &karnnNNTargetsNorm,
        QVector<KarnnNNTargetDDA> *karnnNNTargetsNormTrain
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(karnnNNTargetsNorm); ree;

        *karnnNNTargetsNormTrain = karnnNNTargetsNorm;
        std::sort(
            karnnNNTargetsNormTrain->begin(),
            karnnNNTargetsNormTrain->end(),
            [](const KarnnNNTargetDDA &l, const KarnnNNTargetDDA &r) {
                if (MathUtils::tSame(l.candidateScoresDDA->discriminantScore, r.candidateScoresDDA->discriminantScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                    return l.candidateScoresDDA->featuresArray[Occurrences] > r.candidateScoresDDA->featuresArray[Occurrences];
                }

                return l.candidateScoresDDA->discriminantScore > r.candidateScoresDDA->discriminantScore;
            });

        int counter = 0;
        for (const KarnnNNTargetDDA &knt : *karnnNNTargetsNormTrain) {
            counter++;
            if (constexpr double fdrTrainingThreshold = 0.65; knt.candidateScoresDDA->qValue >= fdrTrainingThreshold && !knt.candidateScoresDDA->isDecoy) {
                break;
            }
        }
        karnnNNTargetsNormTrain->resize(counter);

        std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
        constexpr int shuffleCount = 3;
        for (int i = 0; i < shuffleCount; i++) {
            std::shuffle(
                    karnnNNTargetsNormTrain->begin(),
                    karnnNNTargetsNormTrain->end(),
                    rng
            );
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDDAWorkflow::applyNeuralNetClassifier(
        const QVector<CandidateScoresDDA*> &candidateScoresTargetsAndDecoys,
        int seed,
        QVector<CandidateScoresDDA*> *candidateScoreClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_targetDecoyPairCandidateScorePntrs); ree;
    e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoys); ree;

     QVector<CandidateScoresDDA*> candidateScoresTargetsAndDecoysNeuralNet = candidateScoresTargetsAndDecoys;
     e = filterScoredCandidatesForNeuralNet(
         m_parameters.minMs2FragCount,
         &candidateScoresTargetsAndDecoysNeuralNet
         ); ree;
     qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Analyzing" << candidateScoresTargetsAndDecoysNeuralNet.size() << "for filtering";

     candidateScoreClassifier->clear();

	QVector<FeaturesDDA> neuralNetFeatures = DiscriminantScoretron::featuresNeuralNetworkDDA();

	QVector<KarnnNNTargetDDA> karnnNNTargetsNorm;
	e = buildKarnnNNTargetsNormalized(
		neuralNetFeatures,
		candidateScoresTargetsAndDecoysNeuralNet,
		&karnnNNTargetsNorm
		); ree;

     QVector<KarnnNNTargetDDA> karnnNNTargetsNormTrain;
     e = subsetKarnnNNTargetsForTraining(
         karnnNNTargetsNorm,
         &karnnNNTargetsNormTrain
         ); ree;

     const int totalCount = karnnNNTargetsNormTrain.size();
     const int decoyCount = static_cast<int>(std::count_if(
             karnnNNTargetsNormTrain.begin(),
             karnnNNTargetsNormTrain.end(),
             [](const KarnnNNTargetDDA &kt){return kt.candidateScoresDDA->isDecoy;}
             ));

     qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
              << "target vs decoy count for neural net training"
              << totalCount - decoyCount << ":" << decoyCount
              << "total" << totalCount;

     FDRCLassifierNeuralNet fdrClassifierNeuralNet;
     e = trainNeuralNetwork(
             m_parameters,
             karnnNNTargetsNormTrain,
             seed,
             &fdrClassifierNeuralNet
             ); ree;

     qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Inference start";
     QVector<float> predictions;
     e = predictClassifierScores(
         karnnNNTargetsNorm,
         &fdrClassifierNeuralNet,
         &predictions
         ); ree;
     qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Inference end";

     e = processPredictions(
             predictions,
             &karnnNNTargetsNorm
             ); ree;

     *candidateScoreClassifier = candidateScoresTargetsAndDecoysNeuralNet;

    e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::NNClassifierScore,
            &m_targetDecoyPairCandidateScorePntrs
            ); ree

    ERR_RETURN
}


