//
// Created by anichols on 2/17/23.
//

#include "BiophysicalCalcs.h"

#include "AminoAcids.h"

#include <Eigen/Dense>

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

double BiophysicalCalcs::calculateMassFromThomson(
        double mz,
        int charge,
        int monoOffset
        ) {
    return (mz * charge) - (charge * ChemConstants::PROTON) - (monoOffset * ChemConstants::NEUTRON) ;
}

double BiophysicalCalcs::calculateThomsonFromMass(
        double mass,
        int charge
        ) {
    return (mass + (charge * ChemConstants::PROTON)) / charge ;
}

QVector<double>
BiophysicalCalcs::buildTandemFragmentMasses(
        const QString &seq,
        int charge,
        double startMass,
        int maxLength,
        const AminoAcids &aa
        ) {

    const double NULL_VALUE = -1.0;

    Eigen::VectorXd aminoAcidsMassValues(maxLength);
    aminoAcidsMassValues.setZero();

    for(int i = 0; i < seq.size(); i++){

        if (i >= aminoAcidsMassValues.size()){
            break;
        }

        const QChar aminoAcid = seq.at(i);
        aminoAcidsMassValues.coeffRef(i) = aa.aminoAcid(aminoAcid).monoisotopicMass();
    }

    Eigen::VectorXd cumSumAminoAcidsMassValues(maxLength);

    std::partial_sum(aminoAcidsMassValues.begin(),
                     aminoAcidsMassValues.begin() + maxLength,
                     cumSumAminoAcidsMassValues.begin(),
                     std::plus<double>());

    cumSumAminoAcidsMassValues = cumSumAminoAcidsMassValues.array() + startMass;
    cumSumAminoAcidsMassValues /= charge;

    for(int i = seq.length(); i < cumSumAminoAcidsMassValues.size(); i++){
        cumSumAminoAcidsMassValues.coeffRef(i) = NULL_VALUE;
    }

    std::vector<double> vec(cumSumAminoAcidsMassValues.data(), cumSumAminoAcidsMassValues.data() + maxLength);
    return QVector<double>::fromStdVector(vec);
}

double BiophysicalCalcs::calculateThomson(
        const QString &sequence,
        const AminoAcids &aminoAcids,
        int charge
        ) {

    const double peptideMass = calculatePeptideMass(sequence, aminoAcids);
    return ((charge * ChemConstants::PROTON) + peptideMass) / charge;
}
