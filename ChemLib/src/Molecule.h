#ifndef MOLECULE_H
#define MOLECULE_H

#include "ChemLib_Exports.h"
#include "ChemConstants.h"
#include "MolecularFormula.h"


using namespace ChemConstants;
using namespace MolecularFormulas;

class CHEMLIB_EXPORTS Molecule {

public:

    explicit Molecule(const MolecularFormula &molecularFormula);

    Molecule() = default;
    ~Molecule() = default;

    [[nodiscard]] double averageMass() const;
    [[nodiscard]] double monoisotopicMass() const;
    [[nodiscard]] MolecularFormula molecularFormula() const;


private:
    void computeMoleculeMass();


private:
    double m_monoisotopicMass = 0.0;
    double m_averageMass = 0.0;
    MolecularFormula m_molecularFormula;

};


namespace CommonMolecules {

    const extern Molecule CHEMLIB_EXPORTS H2O;

}


#endif // MOLECULE_H
