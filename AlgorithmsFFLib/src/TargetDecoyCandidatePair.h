//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "FragLibReaderRow.h"
#include "GlobalSettings.h"
#include "MS2Ion.h"
#include "PeptideStringWithMods.h"


using namespace Error;


/*
 * Decoy fragment m/z shifts can be interpreted in two modes.
 * Both use the same decoy-sequence mutation pattern (mutating the penultimate residues at both termini),
 * but they differ in how far the corresponding residue-level mass deltas are propagated.
 */
enum class DecoyFragmentShiftMode {
	/*
	 * Chemistry-aligned mode (current default): shifts are propagated only to fragments that
	 * physically contain the mutated residues, preserving sequence-derived fragment chemistry.
	 */
	ShiftPenultimate,

	/*
	 * Terminus-shift mode: in addition to chemistry-aligned propagation, terminal fragments inherit
	 * the penultimate-residue mutation mass shift.
	 * Example (VATVSLPR): based on the A->L mutation (from our mutation table), the mass delta (+42.0469 Da)
	 * is applied to short V-derived fragments (e.g., 100.0757 -> 142.1226 m/z). The shifted mass remains
	 * chemically valid by atom-count arithmetic (delta C3H6; C5H10NO+ + C3H6 = C8H16NO+), but it is not a
	 * canonical sequence-derived peptide fragment. Entrapment analyses indicated no material degradation
	 * of decoy-based FDR calibration under this mode.
	 */
	ShiftTerminalByPenultimate,
};

class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePair {

public:

    friend class DiscriminantScoretronTests;
    friend class QValueSettertronTests;
    friend class TargetDecoyCandidatePairTests;

    TargetDecoyCandidatePair();

    TargetDecoyCandidatePair(
            const PeptideStringWithMods &peptideStringWithMods,
            float decoyMassDelta
            );

    ~TargetDecoyCandidatePair() = default;

    void setFragLibReaderRowPntr(FragLibReaderRow *fragLibReaderRowPntr);


    [[nodiscard]] QString proteinGroups() const;

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
    [[nodiscard]] float mz(bool isDecoy) const;

	[[nodiscard]] bool isDecoy() const;

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
    [[nodiscard]] float iRt(bool isDecoy) const;
    [[nodiscard]] float iIM() const;


    /**
    * @brief Gets the total fragment count of the target-decoy candidate pair.
    *
    * @return int representing the total fragment count.
    */
    [[nodiscard]] int totalFragmentCount() const;

    static void mangleMs2IonsDecoy(QVector<MS2Ion> *ms2Ions); //TODO document
    void decoySharesSequenceWithOtherTarget(bool val);
    void setDecoyFragmentShiftMode(DecoyFragmentShiftMode mode);


private:

    void setPeptideStringWithMods(const PeptideStringWithMods &peptideStringWithMods);

	static QVector<MS2Ion> mutateCandidatePeptideTargetTestAccess(
		const PeptideStringWithMods &peptideStringWithMods,
		const QVector<MS2Ion> &ms2IonTarget,
		const DecoyFragmentShiftMode shiftMode
		);

private:

    PeptideStringWithMods m_peptideStringWithMods;
    float m_decoyMassDelta;

    FragLibReaderRow *m_fragLibReaderRowPntr;
    bool m_decoySharesSequenceWithOtherTarget;
    DecoyFragmentShiftMode m_decoyFragmentShiftMode = DecoyFragmentShiftMode::ShiftPenultimate;
};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
