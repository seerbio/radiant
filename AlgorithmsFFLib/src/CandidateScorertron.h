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


class ALGORITHMSFFLIB_EXPORTS CandidateScorertron {

public:

    CandidateScorertron();
    ~CandidateScorertron();

    Err init(
            const QMap<ScanNumber, ScanPoints*> &diaTargetFrame,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const QMap<ScanNumber, ScanPoints> &scanPointsMS1,
            const PythiaParameters &pythiaParameters,
            const MzTargetKey &targetKey,
            const QPair<double, double> &scanTimeMinMax,
            int topNMS2Ions,
            const QHash<MzHashed, int> &mzHashedVsCount,
            MsCalibratomatic *msCalibratomatic
            );

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

    MsFrame m_ms1Frame;
    QMap<ScanNumber, ScanPoints> m_scanNumberVsScanPointsMS1;
    QMap<ScanNumber, ScanPoints*> m_scanNumberVsScanPointsMS1Pntrs;

    QHash<MzHashed , XICPoints> m_mzHashedVsXICPointsCached;
    QHash<MzHashed, int> m_mzHashedVsCount;

    MzTargetKey m_targetKey;
    QPair<double, double> m_scanTimeMinMax;
    float m_frameIndexMin;
    float m_frameIndexMax;
    float m_mzMin;
    float m_mzMax;

    Q_DISABLE_COPY(CandidateScorertron) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
