//
// Created by anichols on 10/19/23.
//

#ifndef PYTHIADIACPP_CANDIDATESCORERTRON_H
#define PYTHIADIACPP_CANDIDATESCORERTRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"

class MS2Ion;
class XICPoints;


using namespace Error;


class ALGORITHMSLIB_EXPORTS CandidateScorertron {

public:

    CandidateScorertron(int topNMS2Ions);
    ~CandidateScorertron() = default;

    Err calculateScores(
            const QMap<MzHashed, XICPoints> &xicPointMap,
            const QVector<MS2Ion> &ms2Ions
            );

private:

    int m_topNMS2Ions;

};


#endif //PYTHIADIACPP_CANDIDATESCORERTRON_H
