#ifndef AMINOACIDS_H
#define AMINOACIDS_H

#include "ChemLib_Exports.h"
#include "ChemConstants.h"
#include "Error.h"
#include "Molecule.h"
#include "MolecularFormula.h"

#include <QMap>

using namespace Error;

class CHEMLIB_EXPORTS AminoAcids
{

public:

    AminoAcids();
    ~AminoAcids() = default;

    /**
    * @brief Retrieves the Molecule representation of an amino acid.
    *
    * This method returns the Molecule representation of the specified amino acid.
    *
    * @param aminoAcid The one-letter code of the amino acid.
    * @return The Molecule representation of the amino acid.
    */
    [[nodiscard]] Molecule aminoAcid(QChar aminoAcid) const;

    /**
    * @brief Adds a fixed modification to the specified amino acid.
    *
    * This method adds a fixed modification to the specified amino acid, updating its molecular formula
    * and maintaining the modifications information.
    *
    * @param aminoAcid The one-letter code of the amino acid.
    * @param modificationFormula The molecular formula of the fixed modification.
    * @return An error code indicating the success or failure of the operation.
    */
    Err addFixedModification(QChar aminoAcid, MolecularFormula modificationFormula);

    /**
    * @brief Returns a list of allowed amino acids.
    *
    * This method retrieves a list of allowed amino acids based on the available amino acid data.
    *
    * @return A list of allowed amino acids.
    */
    static QList<QChar> allowedAminoAcids();

    /**
    * @brief Returns a map of amino acids and their corresponding molecular formulas.
    *
    * This method provides a map of amino acids and their associated molecular formulas.
    *
    * @return A QMap<QChar, Molecule> containing amino acids and their molecular formulas.
    */
    static QMap<QChar, Molecule> aminoAcids();

    /**
    * @brief Validates a peptide sequence by checking if all characters are allowed amino acid residues.
    *
    * This method checks whether all characters in the given peptide sequence are valid amino acid residues.
    *
    * @param sequence The peptide sequence to be validated.
    * @return True if the peptide sequence is valid; otherwise, false.
    */
    static bool validPeptideSequence(const QString &sequence);

    /**
    * @brief Returns the map of fixed modifications for amino acids.
    *
    * This method provides access to the map containing fixed modifications for each amino acid.
    *
    * @return The map of fixed modifications (amino acid -> molecular formula).
    */
    [[nodiscard]] const QMap<QChar, MolecularFormula> &fixedModifications() const;

    /**
    * @brief Returns the map of mass differences for DIANN amino acid mutations.
    *
    * This method provides a map of mass differences for DIANN mutations from one amino acid to another.
    * The map is represented as amino acid -> mass difference.
    *
    * @return The map of DIANN mutations (amino acid -> mass difference).
    */
    static QMap<QChar, double> diannMutateAminoAcidTo();

private:

    QMap<QChar, Molecule> m_aminoAcids;
    QMap<QChar, MolecularFormula> m_fixedModifications;


};

#endif // AMINOACIDS_H
