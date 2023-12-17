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
            MsFrame *ms1Frame
            );

    ~ScoreOverseer();

    Err buildScores(
            TargetDecoyCandidatePair* targetDecoyCandidatePair,
            const QPair<double, double> &scanTimeMinMax,
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsShadows,
            const MS2Ion &ms2IonUnfragPrecursor,
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &unfragPrecursorVsfeatureFinderHills,
            CandidateScores *candidateScores
    );

    static QPair<FrameIndex, FrameIndex> getMinMaxFrameIndexes(
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills
    );

private:

    TurboXIC m_turboXICMS1;
    MzTargetKey m_mzTargetKey;
    MsFrame *m_ms1Frame;

    Q_DISABLE_COPY(ScoreOverseer) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_SCOREOVERSEER_H
