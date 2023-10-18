//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"


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
            double totalFramentCount,
            TargetDecoyCandidatePairIndex targetDecoyCandidatePairIndex
            );

    ~TargetDecoyCandidatePair() = default;

    PeptideStringWithMods peptideStringWithMods();
    QVector<MS2Ion> ms2IonsTarget();
    QVector<MS2Ion> ms2IonsDecoy();

    [[nodiscard]] double mz() const;
    [[nodiscard]] double charge() const;
    [[nodiscard]] double mass() const;
    [[nodiscard]] double iRt() const;
    [[nodiscard]] int totalFragmentCount() const;
    [[nodiscard]] TargetDecoyCandidatePairIndex targetDecoyCandidatePairIndex() const;


private:

    const PeptideStringWithMods m_peptideStringWithMods;

    const QVector<MS2Ion> m_ms2IonsTarget;
    const QVector<MS2Ion> m_ms2IonsDecoy;

    const int m_charge;
    const double m_mass;
    const double m_iRt;
    const int m_totalFragmentCount;

    const TargetDecoyCandidatePairIndex m_targetDecoyCandidatePairIndex;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
