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
};


#endif //PYTHIADIACPP_BIOPHYSICALCALCS_H
