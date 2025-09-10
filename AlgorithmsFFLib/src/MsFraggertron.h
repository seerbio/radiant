//
// Created by andrewnichols on 9/3/25.
//

#ifndef MSFRAGGERTRON_H
#define MSFRAGGERTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "CandidateScoresDDA.h"
#include "MS2Ion.h"
#include "Ms2IonFraggertronManager.h"
#include "MsReaderBase.h"
#include "PrecursorMzArrangetron.h"
#include "TargetDecoyCandidatePair.h"

#include "Error.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"

// struct TallyResult {
// 	ScanNumber scanNumber = -1;
// 	bool isDecoy = false;
// 	Occurrence occurrence = 0;
// 	float cosineSimilarity = -1.0;
// 	float relativeIntensityDifferenceAverage = -1.0;
// 	float intensitiesSum = -1.0;
// 	QVector<int> ranks;
// 	int totalFound = -1;
// 	// int indexesFoundY = 0;
// 	// int indexesFoundB = 0;
// 	// int seqTagLongestY = 0;
// 	// int seqTagLongestB = 0;
// };

struct Tally {
	ScanNumber scanNumber = -1;
	Occurrence occurrence = 0;
	float cosineSimilarity = -1.0;
	float relativeIntensityDifferenceAverage = -1.0;
	float intensitiesSum = -1.0;
	QVector<int> ranks;
	QVector<float> mzEmperical;
	QVector<float> mzTheoretical;
	QVector<float> intensitiesEmperical;
	QVector<float> intensitiesTheoretical;
	// int indexesFoundY = 0;
	// int indexesFoundB = 0;
	// int seqTagLongestY = 0;
	// int seqTagLongestB = 0;
};

struct IonSearchResult2 {
	ScanNumberMzIntensity* msScanPointPntr = nullptr;
	MS2IonLibrary *ms2IonLibraryPntr = nullptr;
};

using namespace Error;
using ScoresTarget = CandidateScoresDDA;
using ScoresDecoy = CandidateScoresDDA;
using CandidateScoresDDATuple = std::tuple<TargetDecoyCandidatePair*, QVector<CandidateScoresDDA>, QVector<CandidateScoresDDA>>;

class ALGORITHMSFFLIB_EXPORTS MsFraggertron {

public:

	MsFraggertron();
	~MsFraggertron() = default;

	Err init(
		const PythiaParameters & params,
		MsReaderPointerAcc *msReaderPtr
		);

	[[nodiscard]] bool isInit() const;

	Err performFragging(const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairsPntrs);

private:

	QPair<Err, QVector<CandidateScoresDDATuple>> processTargetDecoyCandidatePairsPntrsTranch(
		const QVector<TargetDecoyCandidatePair*> &tdcps

		);

private:

	PythiaParameters m_parameters;
	MsReaderPointerAcc *m_msReaderPtr;
	MsCalibratomatic *m_msCalibratomatic;

	PrecursorMzArrangetron m_precursorMzArrangetron;
	QVector<QVector<ScanNumberMzIntensity*>> m_scanNumberMzIntensityTranched;

	QVector<QPair<float, float>> m_mzPrecursorStartVsStopResult;


};


#endif //MSFRAGGERTRON_H
