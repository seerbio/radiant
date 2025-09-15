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

	e = processAll(); ree;

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

Err PythiaDDAWorkflow::processAll() {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(m_targetDecoyCandidatePairsPntrs); ree;

	QPair<Err, QVector<CandidateScoresDDATuple>> result
			= m_msFraggertron.performFragging(m_targetDecoyCandidatePairsPntrs);
		e = result.first; ree;

		QVector<CandidateScoresDDATuple> &tallyResults = result.second;

		const QVector<FeaturesDDA> features = DiscriminantScoretron::featuresCalibrationDDA();

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

		const int ppmCount = fdrVsCounts[fdrKey];


	ERR_RETURN
}


