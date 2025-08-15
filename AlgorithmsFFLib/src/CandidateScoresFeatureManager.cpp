//
// Created by andrewnichols on 8/13/25.
//

#include "CandidateScoresFeatureManager.h"

	QVector<CandidateScoresFeatureManager::Features> CandidateScoresFeatureManager::featuresCalibration() {
		const QVector<CandidateScoresFeatureManager::Features> baseFeatures = {
			CosineSimSum100GreaterThan80,
			CosineSimSpectrumOverTimeCubed,
			CosineSimSpectrumStDev,
			CosineSim100MS1,
			CosineSimSpectrumCubed,
			TopBottomRatio,
			TopBottomRatioNorm,

		};

		return baseFeatures;
	}


