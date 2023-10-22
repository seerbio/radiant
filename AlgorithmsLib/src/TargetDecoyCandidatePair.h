//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MS2Ion.h"
#include "CandidateScores.h"


using namespace Error;

using TargetDecoyCandidatePairIndex = int;


class ALGORITHMSLIB_EXPORTS TargetDecoyCandidatePair {

public:

    TargetDecoyCandidatePair();
    TargetDecoyCandidatePair(
            const PeptideStringWithMods &peptideStringWithMods,
            const QVector<MS2Ion> &ms2IonsTarget,
            const QVector<MS2Ion> &ms2IonsDecoy,
            int charge,
            double mass,
            double iRt,
            int totalFramentCount,
            TargetDecoyCandidatePairIndex targetDecoyCandidatePairIndex
            );

    ~TargetDecoyCandidatePair() = default;

    CandidateScores* scoresTarget();
    CandidateScores* scoresDecoy();

    [[nodiscard]] PeptideStringWithMods peptideStringWithMods() const;
    [[nodiscard]] QVector<MS2Ion> ms2IonsTarget() const;
    [[nodiscard]] QVector<MS2Ion> ms2IonsDecoy() const;

    [[nodiscard]] double mz() const;
    [[nodiscard]] double charge() const;
    [[nodiscard]] double mass() const;
    [[nodiscard]] double iRt() const;
    [[nodiscard]] int totalFragmentCount() const;
    [[nodiscard]] TargetDecoyCandidatePairIndex targetDecoyCandidatePairIndex() const;

    void setTargetDecoyCandidatePairIndex(TargetDecoyCandidatePairIndex index);

    void clearScores();


private:

    PeptideStringWithMods m_peptideStringWithMods;
    QVector<MS2Ion> m_ms2IonsTarget;
    QVector<MS2Ion> m_ms2IonsDecoy;
    int m_charge;
    double m_mass;
    double m_iRt;
    int m_totalFragmentCount;
    TargetDecoyCandidatePairIndex m_targetDecoyCandidatePairIndex;

    CandidateScores m_scoresTarget;
    CandidateScores m_scoresDecoy;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
