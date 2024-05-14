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
        int topNMS2Ions,
        XICPeakManager *xicPeakManager,
        MsFrame *msFrameMzTarget
        );

    Err calculateScores(
        const QVector<MS2Ion> &ms2Ions,
        TargetDecoyCandidatePair* tdcp,
        CandidateScores *candidateScores
        ) const;


private:

    Err processIntegrationVectorPeakIntegrations(
        const MatriciesAndVecs &matriciesAndVecs,
        const QVector<QPair<PeakIntegrationIndexes, Intensity>> &peakIntegrationsVsIntensity,
        BestCorrelationResult *bestCorrelationResult
        ) const;



private:

    PythiaParameters m_pythiaParameters;
    int m_topNMS2Ions;
    XICPeakManager *m_xicPeakManager;
    MsFrame *m_msFrameMzTarget;
    MsCalibratomatic m_msCalibratomatic;

    Q_DISABLE_COPY(CandidateScorertron) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
