//
// Created by anichols on 10/19/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORERTRON_H
#define PYTHIADIACPP_CANDIDATESCORERTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "CandidateScores.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsFrame.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

#include <QString>

class BestCorrelationResult;
class CandidateScores;
class MatriciesAndVecs;
class MS2Ion;
class MsFrame;
class TargetDecoyCandidatePair;
class TimsMs2IonMobilityIndex;
class XICPeakManager;

using namespace Error;

/**
 * Responsible for scoring an individual TargetDecoyCandidatePair
 */
class ALGORITHMSFFLIB_EXPORTS CandidateScorertron {

public:

    friend class CandidateScorertronTests;

    CandidateScorertron();
    ~CandidateScorertron();

    Err init(
        const PythiaParameters &pythiaParameters,
        const MsCalibratomatic &msCalibratomatic,
        const MzTargetKey &mzTargetKey,
        int topNMS2Ions,
        float minPeakCount,
        float scanTimeRange,
        const QMap<NominalMzMass, QVector<float>> &averagineTable,
        const QVector<Features> &features,
        bool useTopNIntegrationsParameter,
        XICPeakManager *xicPeakManager,
        MsFrame *msFrameMzTarget,
        TurboXIC *turboXicMS1,
        MsFrame *msFrameMS1,
        MsReaderPointerAcc *msReaderPointerAcc,
        TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex
        );

    Err calculateScores(
        const QVector<MS2Ion> &ms2Ions,
        const QVector<float> &weights,
        TargetDecoyCandidatePair* targetDecoyCandidatePair,
        CandidateScores *candidateScores
        ) const;

    [[nodiscard]] QString scoringDiagnosticsSummary() const;
    void printScoringDiagnosticsIfEnabled() const;
    void setUseAdaptiveTimsMobilityCentering(bool useAdaptiveTimsMobilityCentering);


private:

    Err initMatricesdAndVecs(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const QVector<MS2Ion> &ms2Ions,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        MatriciesAndVecs *matriciesAndVecs
        ) const;

    Err setPredictedFrameIndexes(
        float iRT,
        CandidateScores *candidateScores,
        FrameIndex *frameIndexPredictedMin,
        FrameIndex *frameIndexPredictedMax
        ) const;

    [[nodiscard]] float ionMobilityCenter(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair
        ) const;

    Err processIntegrationVectorPeakIntegrations(
        const MatriciesAndVecs &matriciesAndVecs,
        const QVector<QPair<PeakIntegrationIndexes, Intensity>> &peakIntegrationsVsIntensity,
        QVector<BestCorrelationResult> *bestCorrelationResults
        ) const;

    Err setCandidateScores(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const QVector<MS2Ion> &ms2Ions,
        const QVector<BestCorrelationResult> &bestCorrelationResults,
        const QVector<float> &ms1Averagine,
        CandidateScores *candidateScores
        ) const;

    Err setMs1RelatedScores(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const BestCorrelationResult &bestCorrelationResult,
        float ppmTol,
        CandidateScores *candidateScores
        ) const;

    Err setLibraryIonMobilityRelatedScores(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        CandidateScores *candidateScores
        ) const;

    Err setMs2IonMobilityRelatedScores(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const QVector<MS2Ion> &ms2Ions,
        CandidateScores *candidateScores
        ) const;

    Err setFullTheoMs2IonsScores(CandidateScores *candidateScores) const;

private:

    PythiaParameters m_pythiaParameters;
    int m_topNMS2Ions;
    XICPeakManager *m_xicPeakManager;
    MsFrame *m_msFrameMzTarget;
    MsCalibratomatic m_msCalibratomatic;
    TurboXIC *m_turboXicMS1;
    MsFrame *m_msFrameMS1;
    MsReaderPointerAcc *m_msReaderPointerAcc;
    TimsMs2IonMobilityIndex *m_timsMs2IonMobilityIndex;
    MzTargetKey m_mzTargetKey;
    QVector<FrameNumberTIMS> m_ms1FrameNumbersTIMS;
    QVector<QPair<float, IonMobilityIndex>> m_ms1DriftTimeVsIonMobilityIndexTIMS;

    float m_minPeakCount;
    float m_scanTimeRange;

    QMap<NominalMzMass, QVector<float>> m_averagineTable;;

    QVector<Features> m_features;
    bool m_useTopNIntegrationsParam;
    bool m_useAdaptiveTimsMobilityCentering;

    Q_DISABLE_COPY(CandidateScorertron) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
