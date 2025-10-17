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

struct ResidueMutation {
	int index = -1;
	QChar residueMutatedTo;
	float residueDeltaMass = 0;

	friend QDebug operator<<(QDebug dbg, const ResidueMutation& obj) {
		dbg.nospace() << "ResidueMutation(" << obj.index << ", " << obj.residueMutatedTo << ", " << obj.residueDeltaMass << ") ";
		return dbg;
	}
};


class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePair {

public:

    friend class DiscriminantScoretronTests;
    friend class QValueSettertronTests;
    friend class TargetDecoyCandidatePairTests;

    TargetDecoyCandidatePair();

    ~TargetDecoyCandidatePair() = default;

	Err init();

    void setFragLibReaderRowPntr(FragLibReaderRow *fragLibReaderRowPntr);


    [[nodiscard]] QString proteinGroups() const;

    /**
    * @brief Gets the peptide string without modifications for the target-decoy candidate pair.
    *
    * @return PeptideString representing the peptide string without modifications.
    */
    [[nodiscard]] PeptideString peptideString() const;
	[[nodiscard]] PeptideString peptideStringDecoy() const;

    /**
    * @brief Gets the peptide string with modifications for the target-decoy candidate pair.
    *
    * @return PeptideStringWithMods representing the peptide string with modifications.
    */
    [[nodiscard]] PeptideStringWithMods peptideStringWithMods() const;
    [[nodiscard]] PeptideStringWithMods peptideStringWithModsDecoy() const;

    /**
    * @brief Gets the MS2 ions for the target-decoy candidate pair.
    *
    * @return QVector<MS2Ion> representing the MS2 ions for the target.
    */
    [[nodiscard]] QVector<MS2Ion> ms2IonsTarget();
    /**
    * @brief Gets the MS2 ions for the decoy of the target-decoy candidate pair.
    *
    * @return QVector<MS2Ion> representing the MS2 ions for the decoy.
    */
    [[nodiscard]] QVector<MS2Ion> ms2IonsDecoy();

	QVector<MS2Ion> mutateCandidatePeptideTarget();

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
    [[nodiscard]] float mass(bool isDecoy) const;

    /**
    * @brief Gets the iRT (indexed retention time) value of the target-decoy candidate pair.
    *
    * @return float representing the iRT value.
    */
    [[nodiscard]] float iRt(bool isDecoy) const;
    [[nodiscard]] float iIM() const;
    [[nodiscard]] bool isInit() const;

    /**
    * @brief Gets the total fragment count of the target-decoy candidate pair.
    *
    * @return int representing the total fragment count.
    */
    [[nodiscard]] int totalFragmentCount() const;

    static void mangleMs2IonsDecoy(QVector<MS2Ion> *ms2Ions); //TODO document
    void decoySharesSequenceWithOtherTarget(bool val);


private:

	Err buildResidueMutations();

	[[nodiscard]] float calculateDecoyMass() const;

	static void mutateCandidatePeptideTargetTestAccess(
		const PeptideStringWithMods &peptideStringWithMods,
		const QVector<MS2Ion> &ms2IonTarget
		);

private:

    PeptideString m_peptideString;
    PeptideStringWithMods m_peptideStringWithMods;
    float m_decoyMassDelta;

    FragLibReaderRow *m_fragLibReaderRowPntr;
    bool m_decoySharesSequenceWithOtherTarget;

	QVector<MS2Ion> m_ms2Targets;
	QVector<ResidueMutation> m_residueMutations;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIR_H
