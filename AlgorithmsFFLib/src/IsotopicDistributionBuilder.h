#ifndef ISOTOPICDISTRIBUTIONBUILDER_H
#define ISOTOPICDISTRIBUTIONBUILDER_H

#include "AlgorithmsFFLib_Exports.h"
#include "Molecule.h"

#include <QMap>
#include <QVector>
#include <QPoint>



class ALGORITHMSFFLIB_EXPORTS IsotopicDistributionBuilder {


public:

    /**
    * @brief Gets the isotopic distribution for a given bio-polymer sequence.
    *
    * This function calculates and returns the isotopic distribution for a given bio-polymer sequence.
    *
    * @param sequence A QStringList representing the bio-polymer sequence.
    * @return A QVector<double> containing the isotopic distribution.
    */
    static QVector<double> getIsotopicDistribution(const QStringList &sequence);

    /**
    * @brief Gets the isotopic distribution for a given peptide mass.
    *
    * This function calculates and returns the isotopic distribution for a given peptide mass.
    *
    * @param mass A double representing the peptide mass.
    * @return A QVector<double> containing the isotopic distribution.
    */
    static QVector<double> getIsotopicDistribution(double mass);

    /**
    * @brief Gets the molecular formula of a bio-polymer sequence.
    *
    * This function calculates and returns the molecular formula of a bio-polymer sequence.
    *
    * @param sequence A QStringList representing the bio-polymer sequence.
    * @return The MolecularFormula of the bio-polymer.
    */
    static MolecularFormula getMolecularFormulaOfBioPolymer(const QStringList &sequence);

};

#endif // ISOTOPICDISTRIBUTIONBUILDER_H
