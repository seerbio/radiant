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

        // J is ambiguous but has a single mass / formula:
        {'J', Molecule(leucineFormula)},

        // O is non-standard but has a single mass / formula:
        {'O', Molecule(pyrrolysineFormula)},

        // U is non-standard, has a single mass / formula, but contains non-CHONPS elements, so it is unsupported:
        // {'U', Molecule(selenocysteineFormula)},
        // {'U', Molecule(cysteicAcidFormula + cysteineFormula)}, //FROM DIANN

        // These are ambiguous and have two or more masses / formulae; they are not supported:
        // B (D/N)
        // Z (E/Q)
        // X
        // {'X', Molecule(methionineFormula + cysteineFormula)} //FROM DIANN
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

QMap<QChar, double> AminoAcids::diannMutateAminoAcidToMass() {

    QMap<QChar, double> diannMutateAminoAcidToMass = {
            {'G', Molecule(leucineFormula).monoisotopicMass() - Molecule(glycineFormula).monoisotopicMass()},
            {'A', Molecule(leucineFormula).monoisotopicMass() - Molecule(alanineFormula).monoisotopicMass()},
            {'V', Molecule(leucineFormula).monoisotopicMass() - Molecule(valineFormula).monoisotopicMass()},
            {'L', Molecule(valineFormula).monoisotopicMass() - Molecule(leucineFormula).monoisotopicMass()},
            {'X', Molecule(valineFormula).monoisotopicMass() - Molecule(leucineFormula + prolineFormula).monoisotopicMass()},
            {'I', Molecule(valineFormula).monoisotopicMass() - Molecule(isoleucineFormula).monoisotopicMass()},
            {'F', Molecule(leucineFormula).monoisotopicMass() - Molecule(phenylalanineFormula).monoisotopicMass()},
            {'M', Molecule(leucineFormula).monoisotopicMass() - Molecule(methionineFormula).monoisotopicMass()},
            {'P', Molecule(leucineFormula).monoisotopicMass() - Molecule(prolineFormula).monoisotopicMass()},
            {'W', Molecule(leucineFormula).monoisotopicMass() - Molecule(tryptophanFormula).monoisotopicMass()},
            {'S', Molecule(threonineFormula).monoisotopicMass() - Molecule(serineFormula).monoisotopicMass()},
            {'C', Molecule(serineFormula).monoisotopicMass() - Molecule(cysteineFormula).monoisotopicMass()},
            {'U', Molecule(serineFormula).monoisotopicMass() - Molecule(cysteineFormula + cysteicAcidFormula).monoisotopicMass()},
            {'T', Molecule(serineFormula).monoisotopicMass() - Molecule(threonineFormula).monoisotopicMass()},
            {'Y', Molecule(serineFormula).monoisotopicMass() - Molecule(tyrosineFormula).monoisotopicMass()},
            {'H', Molecule(serineFormula).monoisotopicMass() - Molecule(histidineFormula).monoisotopicMass()},
            {'K', Molecule(leucineFormula).monoisotopicMass() - Molecule(lysineFormula).monoisotopicMass()},
            {'R', Molecule(leucineFormula).monoisotopicMass() - Molecule(arginineFormula).monoisotopicMass()},
            {'Q', Molecule(asparagineFormula).monoisotopicMass() - Molecule(glutamineFormula).monoisotopicMass()},
            {'E', Molecule(asparticAcidFormula).monoisotopicMass() - Molecule(glutamicAcidFormula).monoisotopicMass()},
            {'N', Molecule(glutamineFormula).monoisotopicMass() - Molecule(asparagineFormula).monoisotopicMass()},
            {'D', Molecule(glutamicAcidFormula).monoisotopicMass() - Molecule(asparticAcidFormula).monoisotopicMass()}
    };

    return diannMutateAminoAcidToMass;
}

QMap<QChar, QChar> AminoAcids::diannMutateAminoAcidToResidue() {

    QMap<QChar, QChar> diannMutateAminoAcidToMass = {
            {'G', 'L'},
            {'A', 'L'},
            {'V', 'L'},
            {'L', 'V'},
            {'X', 'V'},
            {'I', 'V'},
            {'F', 'L'},
            {'M', 'L'},
            {'P', 'L'},
            {'W', 'L'},
            {'S', 'T'},
            {'C', 'S'},
            {'U', 'S'},
            {'T', 'S'},
            {'Y', 'S'},
            {'H', 'S'},
            {'K', 'L'},
            {'R', 'L'},
            {'Q', 'N'},
            {'E', 'D'},
            {'N', 'Q'},
            {'D', 'E'}
    };

    return diannMutateAminoAcidToMass;
}

PeptideStringWithMods AminoAcids::mutatePenultimatePeptideResidues(const PeptideStringWithMods &peptideStringWithMods) {

    const QMap<QChar, QChar> diannMutateAminoAcidToResidues = AminoAcids::diannMutateAminoAcidToResidue();

    const int peptideLength = peptideStringWithMods.sizeNoMods();

    QString moddedPeptide;

    bool unimodOn = false;

    int residueCounter = 0;
    for (const QChar &c : peptideStringWithMods) {

        if (c == '[' || c == '(') {
            unimodOn = true;
        }

        else if (c == ']' || c == ')') {
            moddedPeptide += c;
            unimodOn = false;
            continue;
        }

        if (unimodOn) {
            moddedPeptide += c;
            continue;
        }

        if (residueCounter == 1 || residueCounter == peptideLength - 2) {
            moddedPeptide += diannMutateAminoAcidToResidues.value(c);
            residueCounter++;
            continue;
        }

        moddedPeptide += c;
        residueCounter++;
    }

    return PeptideStringWithMods(moddedPeptide);
}
