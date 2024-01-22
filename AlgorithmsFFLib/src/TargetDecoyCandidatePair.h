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
    friend class QValueSettertronTests;

    TargetDecoyCandidatePair();

    TargetDecoyCandidatePair(
            const PeptideString &peptideString,
            const PeptideStringWithMods &peptideStringWithMods,
            const QVector<MS2Ion> &ms2IonsTarget,
            const QVector<MS2Ion> &ms2IonsDecoy,
            int charge,
            float mass,
            float iRt,
            int totalFramentCount
            );

    ~TargetDecoyCandidatePair() = default;

    /**
    * @brief Gets the peptide string without modifications for the target-decoy candidate pair.
    *
    * @return PeptideString representing the peptide string without modifications.
    */
    [[nodiscard]] PeptideString peptideString() const;

    /**
    * @brief Gets the peptide string with modifications for the target-decoy candidate pair.
    *
    * @return PeptideStringWithMods representing the peptide string with modifications.
    */
    [[nodiscard]] PeptideStringWithMods peptideStringWithMods() const;

    /**
    * @brief Gets the MS2 ions for the target-decoy candidate pair.
    *
    * @return QVector<MS2Ion> representing the MS2 ions for the target.
    */
    [[nodiscard]] QVector<MS2Ion> ms2IonsTarget() const;

    /**
    * @brief Gets the MS2 ions for the decoy of the target-decoy candidate pair.
    *
    * @return QVector<MS2Ion> representing the MS2 ions for the decoy.
    */
    [[nodiscard]] QVector<MS2Ion> ms2IonsDecoy() const;

    /**
    * @brief Calculates and gets the m/z (mass-to-charge ratio) for the target-decoy candidate pair.
    *
    * @return float representing the calculated m/z.
    */
    [[nodiscard]] float mz() const;

    /**
    * @brief Gets the charge of the target-decoy candidate pair.
    *
    * @return int representing the charge of the candidate pair.
    */
    [[nodiscard]] int charge() const;

    /**
    * @brief Gets the mass of the target-decoy candidate pair.
    *
    * @return float representing the mass of the candidate pair.
    */
    [[nodiscard]] float mass() const;

    /**
    * @brief Gets the iRT (indexed retention time) value of the target-decoy candidate pair.
    *
    * @return float representing the iRT value.
    */
    [[nodiscard]] float iRt() const;


    /**
    * @brief Gets the total fragment count of the target-decoy candidate pair.
    *
    * @return int representing the total fragment count.
    */
    [[nodiscard]] int totalFragmentCount() const;

private:

    void setPeptideStringWithMods(const PeptideStringWithMods &peptideStringWithMods);
    void setCharge(int charge);

private:

    PeptideString m_peptideString;
    PeptideStringWithMods m_peptideStringWithMods;
    QVector<MS2Ion> m_ms2IonsTarget;
    QVector<MS2Ion> m_ms2IonsDecoy;
    int m_charge;
    float m_mass;
    float m_iRt;
    int m_totalFragmentCount;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
