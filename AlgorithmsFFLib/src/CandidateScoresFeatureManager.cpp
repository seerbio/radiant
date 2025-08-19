//
// Created by andrewnichols on 8/13/25.
//

#include "CandidateScoresFeatureManager.h"

	QVector<CandidateScoresFeatureManager::Features> CandidateScoresFeatureManager::featuresCalibration() {
		const QVector<CandidateScoresFeatureManager::Features> baseFeatures = {
			CosineSimSumGreaterThan80,
			CosineSimSpectrumOverTimeCubed,
			CosineSimSpectrumStDev,
			CosineSimMS1,
			CosineSimSpectrumCubed,
			TopBottomRatio,
			TopBottomRatioNorm,
		};

		return baseFeatures;
	}


