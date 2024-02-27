//
// Created by anichols on 2/17/23.
//

#include "BiophysicalCalcs.h"

#include "PeptideStringWithMods.h"

#include <Eigen/Dense>

#include <iostream>
#include <numeric>


double BiophysicalCalcs::calculatePeptideMass(
        const QString &sequence,
        const AminoAcids &aminoAcids,
        const QHash<ResidueIndex, ModificationMass> &mods
) {

    const auto accumulateLogic
            = [&](double sum, const QChar &aa){return sum + aminoAcids.aminoAcid(aa).monoisotopicMass(); };

    const double mass = std::accumulate(sequence.begin(),  sequence.end(), CommonMolecules::H2O.monoisotopicMass(), accumulateLogic);

    const QList<double> &modMasses = mods.values();
    const double modMassesSum = std::accumulate(modMasses.begin(), modMasses.end(), 0.0);

    return mass + modMassesSum;
}

QVector<QPair<PeptideString, double>> BiophysicalCalcs::calculatePeptideMasses(
        const QVector<QPair<PeptideString, QHash<ResidueIndex, ModificationMass>>> &sequenceAndMods,
        const AminoAcids &aminoAcids
        ) {

    QVector<QPair<PeptideString, double>> output;
    for (const QPair<PeptideString, QHash<ResidueIndex, ModificationMass>> &pr : sequenceAndMods) {

        const double mass = calculatePeptideMass(
                pr.first,
                aminoAcids,
                pr.second
                );

        output.push_back({pr.first, mass});
    }

    return output;
}

double BiophysicalCalcs::calculateMassFromThomson(
        double mz,
        int charge,
        int monoOffset
) {
    return (mz * charge) - (charge * ChemConstants::PROTON) - (monoOffset * ChemConstants::NEUTRON);
}

QVector<double> BiophysicalCalcs::buildTandemFragmentMasses(
        const PeptideStringWithMods &peptideStringWithMods,
        const FragmentSeriesType &fragmentSeriesType,
        int charge,
        const AminoAcids &aa
        ) {

    const double NULL_VALUE = -1.0;

    const PeptideString peptideString = peptideStringWithMods.removeUniModChars();
    const QMap<Index, Mass> modsMap = peptideStringWithMods.modificationsMap();

    const int pepLength = peptideString.size();

    Eigen::VectorXd aminoAcidsMassValues(pepLength);
    aminoAcidsMassValues.setZero();
    aminoAcidsMassValues.coeffRef(0) = fragmentSeriesType == FragmentSeriesType::bSeries
            ? (charge * ChemConstants::PROTON)
            : 0.0;

    aminoAcidsMassValues.coeffRef(pepLength - 1) = fragmentSeriesType == FragmentSeriesType::bSeries
            ? CommonMolecules::H2O.monoisotopicMass()
            : CommonMolecules::H2O.monoisotopicMass() + (charge * ChemConstants::PROTON);

    for(int i = 0; i < pepLength; i++){

        if (i >= aminoAcidsMassValues.size()){
            break;
        }

        const QChar aminoAcid = peptideString.at(i);

        aminoAcidsMassValues.coeffRef(i) += modsMap.contains(i)
                ? aa.aminoAcid(aminoAcid).monoisotopicMass() + modsMap.value(i)
                : aa.aminoAcid(aminoAcid).monoisotopicMass();
    }

    Eigen::VectorXd cumSumAminoAcidsMassValues(pepLength);

    if (fragmentSeriesType == FragmentSeriesType::ySeries) {
        std::reverse(aminoAcidsMassValues.begin(), aminoAcidsMassValues.end());
    }

    std::partial_sum(
            aminoAcidsMassValues.begin(),
            aminoAcidsMassValues.begin() + pepLength,
            cumSumAminoAcidsMassValues.begin(),
            std::plus<double>()
                    );

    cumSumAminoAcidsMassValues /= charge;

    for(int i = pepLength; i < cumSumAminoAcidsMassValues.size(); i++){
        cumSumAminoAcidsMassValues.coeffRef(i) = NULL_VALUE;
    }

    std::vector<double> vec(cumSumAminoAcidsMassValues.data(), cumSumAminoAcidsMassValues.data() + pepLength);
    return {vec.begin(), vec.end()};
}

double BiophysicalCalcs::calculateThomson(
        const QString &sequence,
        const AminoAcids &aminoAcids,
        int charge
        ) {

    const double peptideMass = calculatePeptideMass(sequence, aminoAcids);
    return ((charge * ChemConstants::PROTON) + peptideMass) / charge;
}

int BiophysicalCalcs::calculateChargeFromSequence(
        const QString &peptideSequence,
        const AminoAcids &aminoAcids,
        double mz
        ) {

    const double mass = calculatePeptideMass(
            peptideSequence,
            aminoAcids,
            {}
            );

    return static_cast<int>(std::round(mass / mz));

}

QVector<double> BiophysicalCalcs::calculateIsotopesFromMz(
        double mz,
        double charge,
        int ionCountLeft,
        int ionCountRight
        ) {

    ionCountLeft = std::max(ionCountLeft, 0);
    ionCountRight = std::max(ionCountRight, 0);

    QVector<double> isotopes;
    const double chargeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / charge;

    for (int i = -ionCountLeft; i <= ionCountRight; i++) {
        const double isotopeIsoDiff = i * chargeDistance;
        isotopes.push_back(mz + isotopeIsoDiff);
    }

    return isotopes;
}
