//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MS2Ion.h"


using namespace Error;


class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePair {

public:

    friend class DiscriminantScoretronTest;

    TargetDecoyCandidatePair();

    TargetDecoyCandidatePair(
            const PeptideStringWithMods &peptideStringWithMods,
            const QVector<MS2Ion> &ms2IonsTarget,
            const QVector<MS2Ion> &ms2IonsDecoy,
            int charge,
            float mass,
            float iRt,
            int totalFramentCount
            );

    ~TargetDecoyCandidatePair() = default;

    [[nodiscard]] PeptideStringWithMods peptideStringWithMods() const;
    [[nodiscard]] QVector<MS2Ion> ms2IonsTarget() const;
    [[nodiscard]] QVector<MS2Ion> ms2IonsDecoy() const;
    [[nodiscard]] float mz() const;
    [[nodiscard]] int charge() const;
    [[nodiscard]] float mass() const;
    [[nodiscard]] float iRt() const;
    [[nodiscard]] int totalFragmentCount() const;

private:

    void setPeptideStringWithMods(const PeptideStringWithMods &peptideStringWithMods);
    void setCharge(int charge);

private:

    PeptideStringWithMods m_peptideStringWithMods;
    QVector<MS2Ion> m_ms2IonsTarget;
    QVector<MS2Ion> m_ms2IonsDecoy;
    int m_charge;
    float m_mass;
    float m_iRt;
    int m_totalFragmentCount;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
