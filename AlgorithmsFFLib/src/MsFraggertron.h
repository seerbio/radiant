//
// Created by andrewnichols on 9/3/25.
//

#ifndef MSFRAGGERTRON_H
#define MSFRAGGERTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "MS2Ion.h"
#include "Ms2IonFraggertronManager.h"
#include "MsReaderBase.h"
#include "TargetDecoyCandidatePair.h"

#include "Error.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"

struct TallyResult {
	ScanNumber scanNumber = -1;
	bool isDecoy = false;
	Occurrence occurrence = 0;
	float cosineSimilarity = -1.0;
	float intensitiesSum = -1.0;
	QVector<int> ranks;
	// int indexesFoundY = 0;
	// int indexesFoundB = 0;
	// int seqTagLongestY = 0;
	// int seqTagLongestB = 0;
};

struct Tally {
	ScanNumber scanNumber = -1;
	Occurrence occurrence = 0;
	float cosineSimilarity = -1.0;
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

using namespace Error;
using TallyResultTarget = TallyResult;
using TallyResultDecoy = TallyResult;
using TallyResultTuple = std::tuple<TargetDecoyCandidatePair*, QVector<TallyResultTarget>, QVector<TallyResultDecoy>>;

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

	Err buildMsScanPointPntrs();

	Err extractScansParallel(const QVector<MsScanInfo*> &scanInfosPntrs);

	QPair<Err, QVector<TallyResultTuple>> processTargetDecoyCandidatePairsPntrsTranch(
		const QVector<TargetDecoyCandidatePair*> &tdcps,
		const QVector<ProcessingGroup> &processingGroups
		);

private:

	PythiaParameters m_parameters;
	MsReaderPointerAcc *m_msReaderPtr;
	MsCalibratomatic *m_msCalibratomatic;
	QVector<QVector<MsScanPoint>> m_msScanPointsTranched;
	QVector<QVector<MsScanPoint*>> m_msScanPointsPntrsTranched;
	QVector<QPair<float, float>> m_mzPrecursorStartVsStopResult;


};


#endif //MSFRAGGERTRON_H
