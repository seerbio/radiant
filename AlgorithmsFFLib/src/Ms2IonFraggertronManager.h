//
// Created by anichols on 6/7/24.
//

#ifndef MS2IONFRAGGERTRONMANAGER_H
#define MS2IONFRAGGERTRONMANAGER_H

#include "AlgorithmsFFLib_Exports.h"
#include "MS2Ion.h"
#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePair.h"

#include "Error.h"

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
	QVector<float> intensitiesEmperical;
	QVector<float> intensitiesTheoretical;
	// int indexesFoundY = 0;
	// int indexesFoundB = 0;
	// int seqTagLongestY = 0;
	// int seqTagLongestB = 0;
};

using TallyResultTarget = TallyResult;
using TallyResultDecoy = TallyResult;
using TallyResultTuple
	= std::tuple<TargetDecoyCandidatePair*, QVector<TallyResultTarget>, QVector<TallyResultDecoy>>;

using namespace Error;

struct MsScanPoint {
	MsScanInfo *scanInfoPntr = nullptr;
	float mzVal = -1.0;
	float intensityVal = -1.0;
};

struct MS2IonLibrary {
	TargetDecoyCandidatePair *targeDecoyCandidatePairPntr = nullptr;
	MS2Ion *ms2IonPntr = nullptr;
	bool isDecoy = false;
};

struct ProcessingGroup {
	QPair<float, float> mzPrecursorRangeMinMax;
	QVector<MsScanPoint*> *msScanPointsPntr = nullptr;
	QVector<MS2IonLibrary*> ms2IonsLibrary;
};

struct IonSearchResult {
	MsScanPoint* msScanPointPntr = nullptr;
	MS2IonLibrary *ms2IonLibraryPntr = nullptr;
};

class CandidateScores;
class MsCalibratomatic;

class ALGORITHMSFFLIB_EXPORTS Ms2IonFraggertronManager {

    friend class PlayGroundTests;

public:

    Ms2IonFraggertronManager();
    ~Ms2IonFraggertronManager();

    [[nodiscard]] Err init(
    	const PythiaParameters &parameters,
		const QVector<MS2IonLibrary*> &ms2IonsLibrary,
		MsReaderPointerAcc *msReaderPtr
    	);

	Err performFragging();

private:

	Err buildMsScanPointPntrs(MsReaderPointerAcc *msReaderPtr);

	Err extractScansParallel(
		const QVector<MsScanInfo*> &scanInfosPntrs,
		MsReaderPointerAcc *msReaderPtr
		);

	Err extractMs2Points(
		float mzMin,
		float mzMax,
		float iRTMin,
		float iRTMax,
		QVector<MS2IonLibrary*> *tdPeptideFrags
		) const;

    [[nodiscard]] Err initTesting(const QVector<CandidateScores*> &candidateScoresPntrs) const;

private:

	PythiaParameters m_parameters;
	QVector<QVector<MsScanPoint>> m_msScanPointsTranched;
	QVector<QVector<MsScanPoint*>> m_msScanPointsPntrsTranched;
	QVector<QPair<float, float>> m_mzPrecursorStartVsStopResult;

    Q_DISABLE_COPY(Ms2IonFraggertronManager) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //MS2IONFRAGGERTRONMANAGER_H
