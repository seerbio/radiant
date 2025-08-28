//
// Created by anichols on 6/7/24.
//

#ifndef MS2IONFRAGGERTRONMANAGER_H
#define MS2IONFRAGGERTRONMANAGER_H

#include "AlgorithmsFFLib_Exports.h"
#include "MS2Ion.h"

#include "Error.h"

using namespace Error;

class CandidateScores;
class MsCalibratomatic;
class MsScanInfo;
class TargetDecoyCandidatePair;

struct MS2Frag {
	Id tdcpId = -1;
	float mzVal = -1.f;
	float intensityVal = -1.f;
};

class ALGORITHMSFFLIB_EXPORTS Ms2IonFraggertronManager {

    friend class PlayGroundTests;

public:

    Ms2IonFraggertronManager();
    ~Ms2IonFraggertronManager();

    [[nodiscard]] Err init(
		const QVector<std::tuple<TargetDecoyCandidatePair*, MS2IonsTarget, MS2IonsDecoy>> &ms2Ions
    	) const;

    Err extractMs2Points(
        float mzMin,
        float mzMax,
        float iRTMin,
        float iRTMax,
        QVector<MS2Frag*> *tdPeptideFrags
        );

private:

    [[nodiscard]] Err initTesting(const QVector<CandidateScores*> &candidateScoresPntrs) const;

private:

    Q_DISABLE_COPY(Ms2IonFraggertronManager) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //MS2IONFRAGGERTRONMANAGER_H
