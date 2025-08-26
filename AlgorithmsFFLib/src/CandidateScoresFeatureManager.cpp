//
// Created by andrewnichols on 8/13/25.
//

#include "CandidateScoresFeatureManager.h"

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


