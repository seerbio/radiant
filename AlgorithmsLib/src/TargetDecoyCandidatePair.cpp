//
// Created by anichols on 10/17/23.
//

#include "TargetDecoyCandidatePair.h"

#include "ChemConstants.h"

TargetDecoyCandidatePair::TargetDecoyCandidatePair()
: m_peptideStringWithMods("")
, m_ms2IonsTarget({})
, m_ms2IonsDecoy({})
, m_charge(-1)
, m_mass(-1.0)
, m_iRt(-1.0)
, m_totalFragmentCount(-1)
, m_targetDecoyCandidatePairIndex(-1)
{}

TargetDecoyCandidatePair::TargetDecoyCandidatePair(
        const PeptideStringWithMods &peptideStringWithMods,
        const QVector<MS2Ion> &ms2IonsTarget,
        const QVector<MS2Ion> &ms2IonsDecoy,
        int charge,
        double mass,
        double iRt,
        int totalFramentCount,
        TargetDecoyCandidatePairIndex targetDecoyCandidatePairIndex
        )
        : m_peptideStringWithMods(peptideStringWithMods)
        , m_ms2IonsTarget(ms2IonsTarget)
        , m_ms2IonsDecoy(ms2IonsDecoy)
        , m_charge(charge)
        , m_mass(mass)
        , m_iRt(iRt)
        , m_totalFragmentCount(totalFramentCount)
        , m_targetDecoyCandidatePairIndex(targetDecoyCandidatePairIndex)
{}

PeptideStringWithMods TargetDecoyCandidatePair::peptideStringWithMods() const {
    return m_peptideStringWithMods;
}

QVector<MS2Ion> TargetDecoyCandidatePair::ms2IonsTarget() const {
    return m_ms2IonsTarget;
}

QVector<MS2Ion> TargetDecoyCandidatePair::ms2IonsDecoy() const {
    return m_ms2IonsDecoy;
}

double TargetDecoyCandidatePair::mz() const {
    return (m_mass + (ChemConstants::PROTON * m_charge)) / m_charge;
}

double TargetDecoyCandidatePair::charge() const {
    return m_charge;
}

double TargetDecoyCandidatePair::mass() const {
    return m_mass;
}

double TargetDecoyCandidatePair::iRt() const {
    return m_iRt;
}

int TargetDecoyCandidatePair::totalFragmentCount() const {
    return m_totalFragmentCount;
}

TargetDecoyCandidatePairIndex TargetDecoyCandidatePair::targetDecoyCandidatePairIndex() const {
    return m_targetDecoyCandidatePairIndex;
}

void TargetDecoyCandidatePair::setTargetDecoyCandidatePairIndex(TargetDecoyCandidatePairIndex index) {
    m_targetDecoyCandidatePairIndex = index;
}

CandidateScores* TargetDecoyCandidatePair::scoresTarget() {
    return &m_scoresTarget;
}

CandidateScores* TargetDecoyCandidatePair::scoresDecoy() {
    return &m_scoresDecoy;
}

void TargetDecoyCandidatePair::clearScores() {
    m_scoresDecoy.clear();
    m_scoresTarget.clear();
}
