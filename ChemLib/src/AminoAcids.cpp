#include "AminoAcids.h"

#include "ErrorUtils.h"


using namespace MolecularFormulas;


QMap<QChar, Molecule> AminoAcids::aminoAcids()
{
    QMap<QChar, Molecule> aa = {
       {'A', Molecule(alanineFormula)},
       {'C', Molecule(cysteineFormula)},
       {'D', Molecule(asparticAcidFormula)},
       {'E', Molecule(glutamicAcidFormula)},
       {'F', Molecule(phenylalanineFormula)},
       {'G', Molecule(glycineFormula)},
       {'H', Molecule(histidineFormula)},
       {'I', Molecule(isoleucineFormula)},
       {'K', Molecule(lysineFormula)},
       {'L', Molecule(leucineFormula)},
       {'M', Molecule(methionineFormula)},
       {'N', Molecule(asparagineFormula)},
       {'P', Molecule(prolineFormula)},
       {'Q', Molecule(glutamineFormula)},
       {'R', Molecule(arginineFormula)},
       {'S', Molecule(serineFormula)},
       {'T', Molecule(threonineFormula)},
       {'V', Molecule(valineFormula)},
       {'W', Molecule(tryptophanFormula)},
       {'Y', Molecule(tyrosineFormula)},
       {'X', Molecule(leucineFormula)}
    };

    return aa;
}


AminoAcids::AminoAcids() {
    m_aminoAcids = aminoAcids();
}


Molecule AminoAcids::aminoAcid(QChar aminoAcid) const {

    return m_aminoAcids.value(aminoAcid);
}


Err AminoAcids::addFixedModification(
        QChar aminoAcid,
        MolecularFormula modificationFormula
        ) {

    ERR_INIT

    const MolecularFormula newMolecularFormula
            = m_aminoAcids.value(aminoAcid).molecularFormula() + modificationFormula;

    const Molecule newAminoAcidWithFixedModification(newMolecularFormula);
    m_aminoAcids[aminoAcid] = newAminoAcidWithFixedModification;
    m_fixedModifications.insert(aminoAcid, modificationFormula);

    ERR_RETURN
}

QList<QChar> AminoAcids::allowedAminoAcids() {
    return aminoAcids().keys();
}


bool AminoAcids::validPeptideSequence(const QString &sequence) {

    const QList<QChar> allowedResidues =aminoAcids().keys();
    const auto anyLogic = [allowedResidues](const QChar &c){return allowedResidues.contains(c);};
    return std::all_of(sequence.begin(),  sequence.end(), anyLogic);
}


const QMap<QChar, MolecularFormula> &AminoAcids::fixedModifications() const {
    return m_fixedModifications;
}

