//
// Created by ubuntu on 9/9/25.
//

#ifndef PYTHIADIACPP_CANDIDATESCORESDDA_H
#define PYTHIADIACPP_CANDIDATESCORESDDA_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "TargetDecoyCandidatePair.h"

using namespace Error;

enum FeaturesDDA {

	Occurences = 0,
	CosineSimilaritySpectrum,
	RelativeIntensityDifferenceAverage,
	IntensitiesSum,
	IndexesFoundY,
	IndexesFoundB,
	SeqTagLongestY,
	SeqTagLongestB,
	ScanTimeDelta,
	ScanTimeDeltaFromMean,

	RelativeIntensityRank0,
	RelativeIntensityRank1,
	RelativeIntensityRank2,
	RelativeIntensityRank3,
	RelativeIntensityRank4,
	RelativeIntensityRank5,
	RelativeIntensityRank6,
	RelativeIntensityRank7,
	RelativeIntensityRank8,
	RelativeIntensityRank9,
	RelativeIntensityRank10,
	RelativeIntensityRank11,

	FeaturesDDASize

};


class ALGORITHMSFFLIB_EXPORTS CandidateScoresDDA {

public:

	CandidateScoresDDA() = default;
	~CandidateScoresDDA() = default;

	TargetDecoyCandidatePair *targetDecoyCandidatePair = nullptr;

	bool isDecoy = false;

	ScanNumber scanNumber = -1;
	ScanTime scanTime = -1.0;
	ScanTime scanTimePredicted = -1.0;

	int foundInWindowsCount = -1;
	int foundInWindowsCountPossible = -1;

	double classifierScore = -1.0;
	double discriminantScore = -1.0;
	double qValue = 1.0;
	double decoyRatio = -1.0;

	QVector<float> featuresArray;

public:

	/**
	* @brief Initializes the features array with default values.
	*/
	void initFeaturesArray();

	QVector<float> selectFeaturesArrayFeatures(const QVector<FeaturesDDA> &enumFeatures);

	static QVector<float> selectFeaturesArrayFeatures(
			const QVector<float> &featureArray,
			const QVector<FeaturesDDA> &enumFeatures
			);

	void printFeatures(const QVector<FeaturesDDA> &featuresToPrint);

};


#endif //PYTHIADIACPP_CANDIDATESCORESDDA_H