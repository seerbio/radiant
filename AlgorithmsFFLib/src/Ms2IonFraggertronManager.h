//
// Created by anichols on 6/7/24.
//

#ifndef MS2IONFRAGGERTRONMANAGER_H
#define MS2IONFRAGGERTRONMANAGER_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"

using namespace Error;

class CandidateScores;
class MS2Ion;
class MsCalibratomatic;
class MsScanInfo;
class TargetDecoyCandidatePair;

class ALGORITHMSFFLIB_EXPORTS Ms2IonFraggertronManager {

    friend class PlayGroundTests;

public:

    Ms2IonFraggertronManager();
    ~Ms2IonFraggertronManager();

    [[nodiscard]] Err init(const QVector<TargetDecoyCandidatePair*> &tdcpPntrs) const;

    Err extractMs2Points(
        float mzMin,
        float mzMax,
        float scanTimeMin,
        float scanTimeMax,
        QVector<QPair<MS2Ion, CandidateScores*>> *ms2IonsVsCandidateScoresPntrses
        );

private:

    [[nodiscard]] Err initTesting(const QVector<CandidateScores*> &candidateScoresPntrs) const;

private:

    Q_DISABLE_COPY(Ms2IonFraggertronManager) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //MS2IONFRAGGERTRONMANAGER_H
