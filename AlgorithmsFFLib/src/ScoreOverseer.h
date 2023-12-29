//
// Created by anichols on 10/20/23.
//

#ifndef PYTHIADIACPP_SCOREOVERSEER_H
#define PYTHIADIACPP_SCOREOVERSEER_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "FeatureFinderHillBuilder.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "PythiaParameterReader.h"
#include "TurboXIC.h"

using namespace Error;


class CandidateScores;
class MS2Ion;
class TargetDecoyCandidatePair;


class ALGORITHMSFFLIB_EXPORTS ScoreOverseer {

public:

    ScoreOverseer();

    Err init(
            const PythiaParameters &pythiaParameters,
            const MzTargetKey &targetKey,
            MsFrame *msFrameTandem,
            TurboXIC *turboXICMS1
            );

    ~ScoreOverseer();

    Err buildScores(
            TargetDecoyCandidatePair* targetDecoyCandidatePair,
            const QPair<float, float> &scanTimeMinMax,
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed , XICPoints> &mzHashedVsXICPoints,
            const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
            const QHash<MzHashed , XICPoints> &mzHashedVsXICPointsShadows,
            CandidateScores *candidateScores
    );

private:

    MzTargetKey m_mzTargetKey;
    MsFrame *m_msFrameTandem;
    TurboXIC *m_turboXICMS1;

    Q_DISABLE_COPY(ScoreOverseer) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_SCOREOVERSEER_H
