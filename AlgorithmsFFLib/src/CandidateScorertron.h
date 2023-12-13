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
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const QMap<ScanNumber, ScanPoints> &scanPointsMS1,
            const PythiaParameters &pythiaParameters,
            const MzTargetKey &targetKey,
            int topNMS2Ions,
            MsCalibratomatic *msCalibratomatic,
            FeatureFinderHillBuilder *featureFinderHillsBuilderMS2
            );

    Err calculateScores(
            const TargetDecoyCandidatePair* targetDecoyCandidatePair,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            CandidateScores *candidateScores
            );


private:

    Err extractHills(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            QHash<MzHashed , QVector<FeatureFinderHill*>> *mzHashedVsfeatureFinderHills
    );

    Err findCandidateIntegrations(
            const QVector<float> &summedMatToVec,
            double filterDeltaScoreValue,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
    );

    Err simpleIntegrator(
            const QVector<float> &vec,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
    );

    Err processPeakIntegrationIndexes(
            const TargetDecoyCandidatePair* targetDecoyCandidatePair,
            QVector<QPair<PeakIntegrationIndexes, float>> &peakIntegrationIndexesVsIntensity,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
            CandidateScores *candidateScores
            );

private:

    PythiaParameters m_pythiaParameters;
    int m_topNMS2Ions;

    PeakIntegratomatic m_peakIntegratomatic;
    MsCalibratomatic *m_msCalibratomatic;

    FeatureFinderHillBuilder *m_featureFinderHillsBuilderMS2;
    QHash<MzHashed , QVector<FeatureFinderHill*>> m_mzHashedVsFeatureFinderHillsCached;

    MsFrame m_ms1Frame;
    QMap<ScanNumber, ScanPoints> m_scanNumberVsScanPointsMS1;
    QMap<ScanNumber, ScanPoints*> m_scanNumberVsScanPointsMS1Pntrs;

    MzTargetKey m_targetKey;

    Q_DISABLE_COPY(CandidateScorertron) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
