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

    //TODO Write unit test
    Molecule aminoAcid(QChar aminoAcid) const;

    //TODO Write unit test
    Err addFixedModification(QChar aminoAcid, MolecularFormula modificationFormula);

    static QList<QChar> allowedAminoAcids();

    static QMap<QChar, Molecule> aminoAcids();

    static bool validPeptideSequence(const QString &sequence);

    [[nodiscard]] const QMap<QChar, MolecularFormula> &fixedModifications() const;

private:

    QMap<QChar, Molecule> m_aminoAcids;
    QMap<QChar, MolecularFormula> m_fixedModifications;


};

#endif // AMINOACIDS_H
