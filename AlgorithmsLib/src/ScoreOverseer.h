//
// Created by anichols on 10/20/23.
//

#ifndef PYTHIADIACPP_SCOREOVERSEER_H
#define PYTHIADIACPP_SCOREOVERSEER_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "TurboXIC.h"

using namespace Error;

class CandidateScores;
class MS2Ion;
class MsFrame;
class TargetDecoyCandidatePair;


class ALGORITHMSLIB_EXPORTS ScoreOverseer {

public:

    ScoreOverseer(
            int topNMS2Ions,
            double cosineSimToAnchorThreshold,
            double ppmTol,
            const QVector<double> &summedMatVecToVec,
            TurboXIC *turboXICMS1
            );

    ~ScoreOverseer();

    Err buildScores(
            const TargetDecoyCandidatePair* targetDecoyCandidatePair,
            const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints100,
            const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPointsIsotopeShadows,
            MsFrame *msFrame,
            CandidateScores *candidateScores
    );

private:

    Q_DISABLE_COPY(ScoreOverseer) class Private;
    const QScopedPointer<Private> d_ptr;

    TurboXIC *m_turboXICMS1;

};


#endif //PYTHIADIACPP_SCOREOVERSEER_H
