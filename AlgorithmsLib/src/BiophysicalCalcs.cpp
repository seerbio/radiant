//
// Created by anichols on 2/17/23.
//

#include "BiophysicalCalcs.h"

#include "AminoAcids.h"


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