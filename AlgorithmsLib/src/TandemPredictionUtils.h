//
// Created by Drucifer on 12/7/2021.
//

#ifndef TANDEM_PREDICTION_UTILS_H
#define TANDEM_PREDICTION_UTILS_H

#include "AlgorithmsLib_Exports.h"
#include "AminoAcids.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"

#include <QVector>
#include <utility>


struct ALGORITHMSLIB_EXPORTS FragmentIon{

    FragmentIon() = default;

    FragmentIon(double mz, double intensity, QString ionLabel)
    : mz(mz)
    , intensity(intensity)
    , ionLabel(std::move(ionLabel)){}

    QString ionLabel;
    double mz = -1.0;
    double intensity = -1.0;
};


class ALGORITHMSLIB_EXPORTS TandemPredictionUtils  : public QObject {

private:

    Q_OBJECT
    Q_CLASSINFO("RegisterEnumClassesUnscoped", "false")

public:

    enum class IonType {
        a
        , b
        , c
        , x
        , y
        , z
    };
    Q_ENUM(IonType)

    enum class IonModifier {
        None
        , NH3
        , H2O
        , Chrg2
        , Chrg3
    };
    Q_ENUM(IonModifier)


    static bool isValidPeptideForPrediction(const QString &peptideSeq, int charge);

    static int getMaxLengthForPredictionInputByCharge(int charge);

    static QMap<QChar, int> aminoAcidIndex();

    static QVector<int> peptideSequenceToArray(const QString &sequence, int maxSequenceLength);

    static QStringList buildIonLabels(int charge);

    static void filterFragmentIonVectorByIntensity(double intensityThreshold, QVector<FragmentIon> *fragIons);

    static int calculateNormalizedCollisionEnergy(double mz, int charge, int collisionEnergy);

    static void sortPredictionByIonLabel(QVector<FragmentIon> *frags);

    static QVector<double> buildTheoreticalMzListBYOnly(
            const QString &seq,
            double mzMin,
            double mzMax,
            const AminoAcids &aa
            );

    static QPair<QVector<double>, QVector<double>> buildTheoreticalMzListBYSeparated(
            const QString &seq,
            double mzMin,
            double mzMax,
            const AminoAcids &aa,
            int charge,
            const QHash<ResidueIndex, ModificationMass> &modifications
    );

    static Err calculateMzValuesForPrediction(
            const QString &peptideSequence,
            const QString &modificatonString,
            const AminoAcids &aminoAcids,
            int charge,
            QVector<double> *vec
            );

    static Err extractSequenceAndChargeFromPeptideSequenceChargeKey(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            PeptideString *peptideString,
            int *charge
            );

};


#endif //TANDEM_PREDICTION_UTILS_H
