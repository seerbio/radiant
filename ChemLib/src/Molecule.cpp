#include "Molecule.h"
#include "MolecularFormula.h"

Molecule::Molecule(const MolecularFormula &molecularFormula)
    : m_molecularFormula(molecularFormula) {
    computeMoleculeMass();
}


double Molecule::averageMass() const {
    return m_averageMass;
}


double Molecule::monoisotopicMass() const {
    return m_monoisotopicMass;
}


MolecularFormula Molecule::molecularFormula() const {
    return m_molecularFormula;
}


void Molecule::computeMoleculeMass() {
    m_averageMass += m_molecularFormula.carbonCount * AverageMass::CARBON;
    m_averageMass += m_molecularFormula.hydrogenCount * AverageMass::HYDROGEN;
    m_averageMass += m_molecularFormula.nitrogenCount * AverageMass::NITROGEN;
    m_averageMass += m_molecularFormula.oxygenCount * AverageMass::OXYGEN;
    m_averageMass += m_molecularFormula.phosphorusCount * AverageMass::PHOSPHORUS;
    m_averageMass += m_molecularFormula.sulfurCount * AverageMass::SULFUR;

    m_monoisotopicMass += m_molecularFormula.carbonCount * MonoIsotopicMass::CARBON;
    m_monoisotopicMass += m_molecularFormula.hydrogenCount * MonoIsotopicMass::HYDROGEN;
    m_monoisotopicMass += m_molecularFormula.nitrogenCount * MonoIsotopicMass::NITROGEN;
    m_monoisotopicMass += m_molecularFormula.oxygenCount * MonoIsotopicMass::OXYGEN;
    m_monoisotopicMass += m_molecularFormula.phosphorusCount * MonoIsotopicMass::PHOSPHORUS;
    m_monoisotopicMass += m_molecularFormula.sulfurCount * MonoIsotopicMass::SULFUR;
}


namespace CommonMolecules {
    const Molecule H2O(MolecularFormulas::waterFormula);
}