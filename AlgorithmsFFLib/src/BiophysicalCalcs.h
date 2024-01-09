//
// Created by anichols on 2/17/23.
//

#ifndef PYTHIADIACPP_BIOPHYSICALCALCS_H
#define PYTHIADIACPP_BIOPHYSICALCALCS_H

#include "AlgorithmsFFLib_Exports.h"
#include "AminoAcids.h"
#include "GlobalSettings.h"

#include <QHash>

class AminoAcids;

class ALGORITHMSFFLIB_EXPORTS BiophysicalCalcs {

public:

    static double calculatePeptideMass(
            const QString &sequence,
            const AminoAcids &aminoAcids,
            const QHash<ResidueIndex, ModificationMass> &mods = {}
    );

    static QVector<QPair<PeptideString, double>> calculatePeptideMasses(
            const QVector<QPair<PeptideString,QHash<ResidueIndex, ModificationMass>>> &sequenceAndMods,
            const AminoAcids &aminoAcids
    );

    static double calculateMassFromThomson(
            double mz,
            int charge,
            int monoOffset
            );

    template<typename T>
    static T calculateThomsonFromMass(
            T mass,
            int charge,
            int monoOffset = 0
            ) {
            return static_cast<T>((mass + (charge * ChemConstants::PROTON) + (monoOffset * ChemConstants::NEUTRON) ) / charge) ;
    }

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

    static int calculateChargeFromSequence(
            const QString &peptideSequence,
            const AminoAcids &aminoAcids,
            double mz
            );

    static QVector<double> calculateIsotopesFromMz(
            double mz,
            double charge,
            int ionCountLeft,
            int ionCountRight
            );

};


#endif //PYTHIADIACPP_BIOPHYSICALCALCS_H
