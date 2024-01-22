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

    /**
    * @brief Returns the average mass of the molecule.
    *
    * This function returns the average mass of the molecule, which is precalculated and stored
    * during the initialization of the Molecule object.
    *
    * @return The average mass of the molecule.
    */
    [[nodiscard]] double averageMass() const;

    /**
    * @brief Returns the monoisotopic mass of the molecule.
    *
    * This function returns the monoisotopic mass of the molecule, which is precalculated and stored
    * during the initialization of the Molecule object.
    *
    * @return The monoisotopic mass of the molecule.
    */
    [[nodiscard]] double monoisotopicMass() const;

    /**
    * @brief Returns the molecular formula of the molecule.
    *
    * This function returns the molecular formula of the molecule, which is precalculated and stored
    * during the initialization of the Molecule object.
    *
    * @return The molecular formula of the molecule.
    */
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
