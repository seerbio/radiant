//
// Created by anichols on 6/7/24.
//

#ifndef MS2IONFRAGGERTRONMANAGER_H
#define MS2IONFRAGGERTRONMANAGER_H

#include "AlgorithmsFFLib_Exports.h"
#include "MS2Ion.h"
#include "MsReaderBase.h"
#include "TargetDecoyCandidatePair.h"

#include "Error.h"

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
	QVector<MsScanPoint*> msScanPoints;
	QVector<MS2IonLibrary*> ms2IonsLibrary;
};

class CandidateScores;
class MsCalibratomatic;

class ALGORITHMSFFLIB_EXPORTS Ms2IonFraggertronManager {

    friend class PlayGroundTests;

public:

    Ms2IonFraggertronManager();
    ~Ms2IonFraggertronManager();

    [[nodiscard]] Err init(
		const QVector<MS2IonLibrary*> &ms2IonsLibrary
    	) const;

    Err extractMs2Points(
        float mzMin,
        float mzMax,
        float iRTMin,
        float iRTMax,
        QVector<MS2IonLibrary*> *tdPeptideFrags
        ) const;

private:

    [[nodiscard]] Err initTesting(const QVector<CandidateScores*> &candidateScoresPntrs) const;

private:

    Q_DISABLE_COPY(Ms2IonFraggertronManager) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //MS2IONFRAGGERTRONMANAGER_H
