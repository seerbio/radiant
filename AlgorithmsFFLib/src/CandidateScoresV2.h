//
// Created by andrewnichols on 8/14/25.
//

#ifndef CANDIDATESCORESV2_H
#define CANDIDATESCORESV2_H

#include "AlgorithmsFFLib_Exports.h"

#include "CandidateScoresFeatureManager.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "TargetDecoyCandidatePair.h"

#include <cstring>

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS CandidateScoresV2 {

public:

	CandidateScoresV2() = default;
	~CandidateScoresV2() = default;

	TargetDecoyCandidatePair *targetDecoyCandidatePair = nullptr;

	MzTargetKey targetKey;
	QString proteinGroup;

	bool isDecoy = false;
	FrameIndex frameIndex = -1;
	FrameIndex frameIndexStart = -1;
	FrameIndex frameIndexEnd = -1;
	ScanNumber scanNumber = -1;
	ScanNumber scanNumberStart = -1;
	ScanNumber scanNumberEnd = -1;
	ScanTime scanTime = -1.0;
	ScanTime scanTimeStart = -1.0;
	ScanTime scanTimeEnd = -1.0;
	ScanTime scanTimePredicted = -1.0;

	double classifierScore = -1.0;
	double discriminantScore = -1.0;
	double qValue = 1.0;

	QVector<float> featuresArray;

	void initFeaturesArray() {
		featuresArray.clear();
		featuresArray.resize(FTR::FeaturesSize);
		std::memset(featuresArray.data(), 0, FTR::FeaturesSize * sizeof(float));
	}

};

using CandidateScoresV2Target = CandidateScoresV2;
using CandidateScoresV2Decoy = CandidateScoresV2;



#endif //CANDIDATESCORESV2_H
