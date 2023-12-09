//
// Created by anichols on 10/19/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORERTRON_H
#define PYTHIADIACPP_CANDIDATESCORERTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "PeakIntegratomatic.h"
#include "PythiaParameterReader.h"
#include "TurboXIC.h"

class CandidateScores;
class MS2Ion;
class MsFrame;
class TargetDecoyCandidatePair;
class XICPoints;


using namespace Error;


class ALGORITHMSFFLIB_EXPORTS CandidateScorertron {

public:

    CandidateScorertron();
    ~CandidateScorertron() = default;

    Err init(
            const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPointsMS1,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            const PythiaParameters &pythiaParameters,
            int topNMS2Ions
            );

    Err calculateScores(
            const TargetDecoyCandidatePair* targetDecoyCandidatePair,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPointsIsotopeShadows,
            double scanTimePredicted,
            MsFrame *msFrame,
            CandidateScores *candidateScores
            );

private:

    Err findCandidateIntegrations(
            const QVector<double> &summedMatToVec,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
    );

private:

    PythiaParameters m_pythiaParameters;
    int m_topNMS2Ions;

    PeakIntegratomatic m_peakIntegratomatic;

    TurboXIC m_turboXICMS1;

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
