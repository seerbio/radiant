//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePair.h"

#include "BiophysicalCalcs.h"

TargetDecoyCandidatePair::TargetDecoyCandidatePair()
: m_peptideStringWithMods("")
, m_ms2IonsTarget({})
, m_ms2IonsDecoy({})
, m_charge(-1)
, m_mass(-1.0)
, m_iRt(-1.0)
, m_totalFragmentCount(-1)
{}

TargetDecoyCandidatePair::TargetDecoyCandidatePair(
        const PeptideStringWithMods &peptideStringWithMods,
        const QVector<MS2Ion> &ms2IonsTarget,
        const QVector<MS2Ion> &ms2IonsDecoy,
        int charge,
        float mass,
        float iRt,
        int totalFramentCount
)
: m_peptideString(peptideStringWithMods.removeUniModChars())
, m_peptideStringWithMods(peptideStringWithMods)
, m_ms2IonsTarget(ms2IonsTarget)
, m_ms2IonsDecoy(ms2IonsDecoy)
, m_charge(charge)
, m_mass(mass)
, m_iRt(iRt)
, m_totalFragmentCount(totalFramentCount)
{}


PeptideString TargetDecoyCandidatePair::peptideString() const {
    return m_peptideString;
}

PeptideStringWithMods TargetDecoyCandidatePair::peptideStringWithMods() const {
    return m_peptideStringWithMods;
}

QVector<MS2Ion> TargetDecoyCandidatePair::ms2IonsTarget() const {
    return m_ms2IonsTarget;
}

QVector<MS2Ion> TargetDecoyCandidatePair::ms2IonsDecoy() const {
    return m_ms2IonsDecoy;
}

float TargetDecoyCandidatePair::mz() const {
    return static_cast<float>(BiophysicalCalcs::calculateThomsonFromMass(m_mass, m_charge));
}

int TargetDecoyCandidatePair::charge() const {
    return m_charge;
}

float TargetDecoyCandidatePair::mass() const {
    return m_mass;
}

float TargetDecoyCandidatePair::iRt() const {
    return m_iRt;
}

int TargetDecoyCandidatePair::totalFragmentCount() const {
    return m_totalFragmentCount;
}

void TargetDecoyCandidatePair::setPeptideStringWithMods(const PeptideStringWithMods &peptideStringWithMods) {
    m_peptideStringWithMods = peptideStringWithMods;
}

void TargetDecoyCandidatePair::setCharge(int charge) {
    m_charge = charge;
}
