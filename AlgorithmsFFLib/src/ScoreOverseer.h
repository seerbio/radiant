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
    ~ScoreOverseer();

    /**
    * @brief Initializes the ScoreOverseer with necessary parameters.
    *
    * This function initializes the ScoreOverseer with the provided Pythia parameters, target key, MS frame, and TurboXIC for MS1.
    *
    * @param pythiaParameters The Pythia parameters to use for initialization.
    * @param targetKey The MzTargetKey representing the target.
    * @param msFrameTandem The MSFrame object for tandem MS.
    * @param turboXICMS1 The TurboXIC object for MS1.
    * @return An error code indicating the success or failure of the initialization process.
    */
    Err init(
            const PythiaParameters &pythiaParameters,
            const MzTargetKey &targetKey,
            MsFrame *msFrameTandem,
            TurboXIC *turboXICMS1
            );

    /**
    * @brief Builds scores for a target-decoy candidate pair.
    *
    * This function builds scores for a given target-decoy candidate pair using various input data.
    * It initializes features and calculates alignment matrices, then computes various metrics and populates the features array.
    *
    * @param targetDecoyCandidatePair The target-decoy candidate pair for which to build scores.
    * @param scanTimeMinMax The minimum and maximum scan time range.
    * @param peakIntegrationIndexes The peak integration indexes.
    * @param ms2IonsTheoretical The vector of theoretical MS2 ions.
    * @param mzHashedVsXICPoints The hash map of hashed m/z values to XIC points.
    * @param ms2IonsTheoreticalIsotopeShadows The vector of theoretical MS2 ions isotope shadows.
    * @param mzHashedVsXICPointsShadows The hash map of hashed m/z values to XIC points for isotope shadows.
    * @param candidateScores The candidate scores object to store the computed scores.
    * @return An error code indicating the success or failure of the score building process.
    */
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
