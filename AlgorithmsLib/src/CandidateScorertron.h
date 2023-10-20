//
// Created by anichols on 10/19/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORERTRON_H
#define PYTHIADIACPP_CANDIDATESCORERTRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "PeakIntegratomatic.h"
#include "PythiaParameterReader.h"

class MS2Ion;
class XICPoints;


using namespace Error;


class ALGORITHMSLIB_EXPORTS CandidateScorertron {

public:

    explicit CandidateScorertron();
    ~CandidateScorertron() = default;

    Err init(const PythiaParameters &pythiaParameters, int topNMS2Ions);

    Err calculateScores(
            const QMap<MzHashed, XICPoints> &xicPointMap,
            const QVector<MS2Ion> &ms2Ions,
            ScanTime scanTimePredicted
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

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
