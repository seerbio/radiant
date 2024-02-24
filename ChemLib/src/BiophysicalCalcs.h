//
// Created by anichols on 2/17/23.
//

#ifndef PYTHIADIACPP_BIOPHYSICALCALCS_H
#define PYTHIADIACPP_BIOPHYSICALCALCS_H

#include "ChemLib_Exports.h"
#include "AminoAcids.h"
#include "GlobalSettings.h"

#include <QHash>

class AminoAcids;
class PeptideStringWithMods;

class CHEMLIB_EXPORTS BiophysicalCalcs {

public:

    /**
    * @brief Calculates the mass of a peptide sequence with modifications.
    *
    * This function calculates the mass of a peptide sequence specified by the `sequence`
    * parameter, considering the monoisotopic masses of amino acids provided in the `aminoAcids`
    * parameter. Additionally, it incorporates modification masses specified in the `mods` parameter.
    * The mass calculation includes the monoisotopic mass of water for each amino acid in the sequence.
    * The function utilizes lambda function `accumulateLogic` for accumulating amino acid masses and
    * then calculates the total mass by summing up the masses of amino acids and modification masses.
    *
    * @param sequence The peptide sequence for which the mass is to be calculated.
    * @param aminoAcids A collection of amino acids and their monoisotopic masses.
    * @param mods A hash map containing modification masses for specific residue indices.
    * @return double The calculated mass of the peptide with modifications.
    *
    */
    static double calculatePeptideMass(
            const QString &sequence,
            const AminoAcids &aminoAcids,
            const QHash<ResidueIndex, ModificationMass> &mods = {}
    );

    /**
    * @brief Calculates masses for a collection of peptide sequences with modifications.
    *
    * This function calculates the masses for a collection of peptide sequences along with their
    * modifications. The input is a vector of pairs, where each pair consists of a peptide sequence
    * and a hash map of modification masses for specific residue indices. The function iterates through
    * each pair, calculates the mass using the `calculatePeptideMass` function, and stores the result
    * in a vector of pairs containing the peptide sequence and its corresponding mass. The final
    * collection of peptide masses is returned.
    *
    * @param sequenceAndMods A vector of pairs containing peptide sequences and modification masses.
    * @param aminoAcids A collection of amino acids and their monoisotopic masses.
    * @return QVector<QPair<PeptideString, double>> The calculated masses for each peptide sequence.
    *
    */
    static QVector<QPair<PeptideString, double>> calculatePeptideMasses(
            const QVector<QPair<PeptideString, QHash<ResidueIndex, ModificationMass>>> &sequenceAndMods,
            const AminoAcids &aminoAcids
    );

    /**
    * @brief Calculates the mass from Thomson units for a given m/z value and charge.
    *
    * This function calculates the mass in Thomson units based on the provided m/z value (`mz`),
    * charge state (`charge`), and a monoisotopic offset (`monoOffset`). The formula used considers
    * the product of m/z and charge, subtracts the product of charge and the proton mass, and further
    * subtracts the product of monoOffset and the neutron mass.
    *
    * @param mz The m/z value for which the mass is to be calculated.
    * @param charge The charge state of the molecule.
    * @param monoOffset The monoisotopic offset.
    * @return double The calculated mass in Thomson units.
    *
    */
    static double calculateMassFromThomson(
            double mz,
            int charge,
            int monoOffset
            );

    /**
    * @brief Calculates the Thomson units from mass for a given mass, charge, and monoisotopic offset.
    *
    * This templated static function calculates the mass in Thomson units based on the provided mass (`mass`),
    * charge state (`charge`), and an optional monoisotopic offset (`monoOffset`). The formula used considers
    * the sum of mass, the product of charge and the proton mass, and the product of monoOffset and the neutron mass.
    * The result is divided by the charge.
    *
    * @tparam T The type of the input mass (e.g., double, float).
    * @param mass The mass for which the Thomson units are to be calculated.
    * @param charge The charge state of the molecule.
    * @param monoOffset The optional monoisotopic offset (default is 0).
    * @return T The calculated Thomson units.
    *
    */
    template<typename T>
    static T calculateThomsonFromMass(
            T mass,
            int charge,
            int monoOffset = 0
    ) {
        return static_cast<T>((mass + (charge * ChemConstants::PROTON) + (monoOffset * ChemConstants::NEUTRON) ) / charge) ;
    }

    enum class FragmentSeriesType {
            bSeries,
            ySeries,
            aSeries
        };

    /**
    * @brief Builds tandem fragment masses for a given peptide sequence, charge, and start mass.
    *
    * This function calculates tandem fragment masses for a peptide sequence (`seq`) based on the provided
    * charge state (`charge`), start mass (`startMass`), and maximum length (`maxLength`). The amino acid masses
    * are retrieved from the `aa` parameter, and the cumulative sum of amino acid masses is computed. The resulting
    * values are adjusted for charge and NULL_VALUE is used for padding beyond the sequence length.
    * The final result is returned as a QVector<double>.
    *
    * @param seq The peptide sequence for which tandem fragment masses are to be calculated.
    * @param charge The charge state of the molecule.
    * @param startMass The starting mass for the calculation.
    * @param maxLength The maximum length of tandem fragment masses to be calculated.
    * @param aa A collection of amino acids and their monoisotopic masses.
    * @return QVector<double> The calculated tandem fragment masses.
    *
    */
    static QVector<double> buildTandemFragmentMasses(
            const PeptideStringWithMods &peptideStringWithMods,
            const FragmentSeriesType &fragmentSeriesType,
            int charge,
            const AminoAcids &aa
            );

    /**
    * @brief Calculates Thomson units for a given peptide sequence, amino acids, and charge.
    *
    * This function calculates Thomson units for a peptide sequence (`sequence`) based on the provided
    * amino acids (`aminoAcids`) and charge state (`charge`). The peptide mass is first calculated
    * using the `calculatePeptideMass` function, and the Thomson units are then determined by adding
    * the product of charge and proton mass to the peptide mass and dividing by the charge.
    *
    * @param sequence The peptide sequence for which Thomson units are to be calculated.
    * @param aminoAcids A collection of amino acids and their monoisotopic masses.
    * @param charge The charge state of the molecule.
    * @return double The calculated Thomson units.
    *
    */
    static double calculateThomson(
            const QString &sequence,
            const AminoAcids &aminoAcids,
            int charge
            );

    /**
    * @brief Calculates the charge state for a given peptide sequence, amino acids, and precursor mass-to-charge ratio.
    *
    * This function calculates the charge state for a peptide sequence (`peptideSequence`) based on the provided
    * amino acids (`aminoAcids`) and precursor mass-to-charge ratio (`mz`). The peptide mass is determined using
    * the `calculatePeptideMass` function with an empty modification set, and the charge state is calculated by
    * rounding the ratio of the peptide mass to the precursor mass-to-charge ratio.
    *
    * @param peptideSequence The peptide sequence for which the charge state is to be calculated.
    * @param aminoAcids A collection of amino acids and their monoisotopic masses.
    * @param mz The precursor mass-to-charge ratio.
    * @return int The calculated charge state.
    *
    */
    static int calculateChargeFromSequence(
            const QString &peptideSequence,
            const AminoAcids &aminoAcids,
            double mz
            );

    /**
    * @brief Calculates isotopic peaks for a given precursor mass-to-charge ratio, charge state, and range of ion counts.
    *
    * This function calculates isotopic peaks for a precursor mass-to-charge ratio (`mz`), charge state (`charge`),
    * and specified range of ion counts. The isotopic peaks are determined by adding multiples of the isotopic
    * difference (`S_GLOBAL_SETTINGS.ISO_DIFF`) to the precursor mass-to-charge ratio based on the given ion count range.
    *
    * @param mz The precursor mass-to-charge ratio.
    * @param charge The charge state of the molecule.
    * @param ionCountLeft The number of isotopic peaks to the left of the main peak.
    * @param ionCountRight The number of isotopic peaks to the right of the main peak.
    * @return QVector<double> The calculated isotopic peaks.
    *
    */
    static QVector<double> calculateIsotopesFromMz(
            double mz,
            double charge,
            int ionCountLeft,
            int ionCountRight
            );

};


#endif //PYTHIADIACPP_BIOPHYSICALCALCS_H
