//
// Created by andrewnichols on 8/13/25.
//

#include "CandidateScoresFeatureManager.h"

#include "CandidateScoresV2.h"

	QVector<FTR> CandidateScoresFeatureManager::featuresCalibration() {
		const QVector<FTR> baseFeatures = {
			CosineSimSumMeanCorrelation,
			CosineSimSumStDevCorrelation,
			CosineSimSpectrumOverTimeCubed,
			CosineSimSpectrumOverTimeStDev,
			ScanTimeRelativeDeltaAbs,
			// CosineSimMS1,
			// CosineSimSpectrumCubed,
			// TopBottomRatio,
			// TopBottomRatioNorm,
		};

		return baseFeatures;
	}

Err CandidateScoresFeatureManager::cleanScores(const QVector<CandidateScoresV2 *> &candidateScoresV2s)  {

		ERR_INIT

		e = ErrorUtils::isNotEmpty(candidateScoresV2s); ree;

		QVector<int> nanCounts(FTR::FeaturesSize, 0);
		QVector<float> featuresArrayMaxes(FTR::FeaturesSize, 0);

		for (CandidateScoresV2 *cs : candidateScoresV2s) {

			for (int i = 0; i < FTR::FeaturesSize; i++) {

				if (std::isnan(cs->featuresArray[i])) {
					nanCounts[i]++;
					continue;
				}

				featuresArrayMaxes[i] = std::max(featuresArrayMaxes[i], cs->featuresArray[i]);
			}
		}

		for (CandidateScoresV2 *cs : candidateScoresV2s) {

			for (int i = 0; i < FTR::FeaturesSize; i++) {

				if (!std::isnan(cs->featuresArray[i])) {
					continue;
				}

				cs->featuresArray[i] = featuresArrayMaxes[i];
				nanCounts[i]--;
			}
		}

		qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished cleaning scores";

		ERR_RETURN
	}


