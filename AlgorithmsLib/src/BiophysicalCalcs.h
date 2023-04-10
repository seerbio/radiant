//
// Created by anichols on 2/17/23.
//

#ifndef PYTHIADIACPP_BIOPHYSICALCALCS_H
#define PYTHIADIACPP_BIOPHYSICALCALCS_H

#include "AlgorithmsLib_Exports.h"
#include "GlobalSettings.h"

#include <QHash>

class AminoAcids;

class ALGORITHMSLIB_EXPORTS BiophysicalCalcs {

public:

    static double calculatePeptideMass(
            const QString &sequence,
            const AminoAcids &aminoAcids,
            const QHash<ResidueIndex, ModificationMass> &mods = {}
    );

    static double calculateMassFromThomson(
            double mz,
            int charge,
            int monoOffset
            );

    static double calculateThomsonFromMass(
            double mass,
            int charge
            );

    static QVector<double> buildTandemFragmentMasses(
            const QString &seq,
            int charge,
            double startMass,
            int maxLength,
            const AminoAcids &aa
            );

    static double calculateThomson(
            const QString &sequence,
            const AminoAcids &aminoAcids,
            int charge
            );
};


#endif //PYTHIADIACPP_BIOPHYSICALCALCS_H
