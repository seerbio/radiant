//
// Created by andrewnichols on 8/12/24.
//

#ifndef QUANTRANSISTIONREFINERTRON_H
#define QUANTRANSISTIONREFINERTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"

using namespace Error;

class CandidateScores;


class ALGORITHMSFFLIB_EXPORTS QuanTransitionRefinertron {

    friend class QuanTransitionRefinertronTests;


public:

    QuanTransitionRefinertron(
        double ppm,
        int frameIndexBuffer
        );
    ~QuanTransitionRefinertron() = default;

    Err refineTransitions(const QVector<CandidateScores*> &candidateScoresPntrs);

    static Err groupTransitionsByMz(
        const QVector<float> &mzVals,
        double ppmThreshold,
        QHash<MzHashed, int> *mzHashedVsCount
        );


private:

    double m_ppmThreshold;
    int m_frameIndexBuffer;

};



#endif //QUANTRANSISTIONREFINERTRON_H
