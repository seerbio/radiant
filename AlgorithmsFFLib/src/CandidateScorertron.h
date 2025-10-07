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
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

class BestCorrelationResult;
class CandidateScores;
class MatriciesAndVecs;
class MS2Ion;
class MsFrame;
class TargetDecoyCandidatePair;
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
        MsFrame *msFrameMS1
        );

    Err calculateScores(
        const QVector<MS2Ion> &ms2Ions,
        const QVector<float> &weights,
        TargetDecoyCandidatePair* targetDecoyCandidatePair,
        CandidateScores *candidateScores
        ) const;


private:

    Err initMatricesdAndVecs(
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

    static Err processIntegrationVectorPeakIntegrations(
        const MatriciesAndVecs &matriciesAndVecs,
        const QVector<QPair<PeakIntegrationIndexes, Intensity>> &peakIntegrationsVsIntensity,
        QVector<BestCorrelationResult> *bestCorrelationResults
        ) ;

    Err setCandidateScores(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
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

    Err setFullTheoMs2IonsScores(CandidateScores *candidateScores) const;

private:

    PythiaParameters m_pythiaParameters;
    int m_topNMS2Ions;
    XICPeakManager *m_xicPeakManager;
    MsFrame *m_msFrameMzTarget;
    MsCalibratomatic m_msCalibratomatic;
    TurboXIC *m_turboXicMS1;
    MsFrame *m_msFrameMS1;
    MzTargetKey m_mzTargetKey;

    float m_minPeakCount;
    float m_scanTimeRange;

    QMap<NominalMzMass, QVector<float>> m_averagineTable;;

    QVector<Features> m_features;
    bool m_useTopNIntegrationsParam;

    Q_DISABLE_COPY(CandidateScorertron) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
