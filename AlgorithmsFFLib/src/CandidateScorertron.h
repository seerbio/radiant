//
// Created by anichols on 10/19/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORERTRON_H
#define PYTHIADIACPP_CANDIDATESCORERTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsFrame.h"
#include "PeakIntegratomatic.h"
#include "PythiaParameterReader.h"
#include "TurboXIC.h"

class CandidateScores;
class FeatureFinderHill;
class FeatureFinderHillBuilder;
class MS2Ion;
class MsFrame;
class TargetDecoyCandidatePair;

using namespace Error;

/**
 * Responsible for scoring an individual TargetDecoyCandidatePair
 */
class ALGORITHMSFFLIB_EXPORTS CandidateScorertron {

public:

    friend class CandidateScorertronTests;

    CandidateScorertron();
    ~CandidateScorertron();

    /**
    * @brief Initializes the CandidateScorertron with the provided parameters.
    *
    * This function initializes the CandidateScorertron with the given parameters, performing necessary checks
    * and initializing internal components. It sets up the PythiaParameters, top N MS2 ions, scan time min-max,
    * mzHashedVsCount, PeakIntegratomatic, MsCalibratomatic, MsFrameTandemScans, TurboXIC, and other internal state.
    *
    * @param diaTargetFrame A map of ScanNumber to ScanPoints representing the DIA target frame.
    * @param scanNumberVsScanTime A map of ScanNumber to ScanTime for associating scan numbers with scan times.
    * @param pythiaParameters The PythiaParameters used for scoring.
    * @param targetKey The MzTargetKey representing the target.
    * @param scanTimeMinMax The min-max range of scan times.
    * @param topNMS2Ions The top N MS2 ions to consider.
    * @param mzHashedVsCount A hash mapping MzHashed to ion counts.
    * @param msCalibratomatic The MsCalibratomatic instance.
    * @param turboXICMS1 The TurboXIC instance for MS1 scans.
    * @return Err The error status indicating the success of the initialization.
    *
    */
    Err init(
            const QMap<ScanNumber, ScanPoints*> &diaTargetFrame,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const PythiaParameters &pythiaParameters,
            const MzTargetKey &targetKey,
            const QPair<double, double> &scanTimeMinMax,
            int topNMS2Ions,
            const QHash<MzHashed, int> &mzHashedVsCount,
            MsCalibratomatic *msCalibratomatic,
            TurboXIC *turboXICMS1
            );

    /**
    * @brief Calculates scores for a given set of MS2 ions and candidate information.
    *
    * This function calculates scores for a given set of MS2 ions and candidate information.
    * It initializes the CandidateScores, predicts the scan time, extracts XICs, builds an integration vector,
    * finds candidate integrations, and processes the results to calculate scores.
    *
    * @param _ms2IonsTheoretical The theoretical MS2 ions to consider.
    * @param targetDecoyCandidatePair The TargetDecoyCandidatePair containing candidate information.
    * @param candidateScores The CandidateScores instance to store the calculated scores.
    * @return Err The error status indicating the success of the score calculation.
    *
    */
    Err calculateScores(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            TargetDecoyCandidatePair* targetDecoyCandidatePair,
            CandidateScores *candidateScores
            );

private:

    Err extractXICs(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            QHash<MzHashed , XICPoints> *mzHashedVsXICPoints
    );

    Err findCandidateIntegrations(
            const QVector<float> &summedMatToVec,
            double filterDeltaScoreValue,
            FrameIndex frameIndexPredictedMin,
            FrameIndex frameIndexPredictedMax,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
    );

    Err simpleIntegrator(
            const QVector<float> &vec,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
    );

    Err processPeakIntegrationIndexes(
            const QVector<QPair<PeakIntegrationIndexes, float>> &peakIntegrationIndexesVsIntensity,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed , XICPoints> &mzHashedVsXICPoints,
            const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
            TargetDecoyCandidatePair* targetDecoyCandidatePair,
            CandidateScores *candidateScores
            );

private:

    PythiaParameters m_pythiaParameters;
    int m_topNMS2Ions;

    PeakIntegratomatic m_peakIntegratomatic;
    MsCalibratomatic *m_msCalibratomatic;

    TurboXIC m_turboXIC;
    MsFrame m_msFrameTandemScans;

    QHash<MzHashed , XICPoints> m_mzHashedVsXICPointsCached;
    QHash<MzHashed, int> m_mzHashedVsCount;

    MzTargetKey m_targetKey;
    QPair<float, float> m_scanTimeMinMax;
    float m_mzMin;
    float m_mzMax;

    Q_DISABLE_COPY(CandidateScorertron) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
