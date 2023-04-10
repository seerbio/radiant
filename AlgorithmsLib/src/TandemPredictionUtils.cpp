//
// Created by Drucifer on 12/7/2021.
//

#include "TandemPredictionUtils.h"

#include "BiophysicalCalcs.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "StringUtils.h"

#include <Eigen/Dense>
#include <QVector>

#include <iostream>


int TandemPredictionUtils::getMaxLengthForPredictionInputByCharge(int charge) {

    switch (charge) {
            case 1:
                return 22;
            case 2:
                return 40;
            case 3:
            case 4:
                return 50;
            default:
                return -1;
        }
}

QMap<QChar, int> TandemPredictionUtils::aminoAcidIndex() {

    return {
            {'D',1},
            {'E',2},
            {'P',3},
            {'M',4},
            {'C',5},
            {'N',6},
            {'Q',7},
            {'S',8},
            {'T',9},
            {'G',10},
            {'A',11},
            {'I',12},
            {'L',13},
            {'V',14},
            {'F',15},
            {'W',16},
            {'Y',17},
            {'H',18},
            {'K',19},
            {'R',20},
            {'X',12},
    };
}

QVector<int> TandemPredictionUtils::peptideSequenceToArray(
        const QString &sequence,
        int maxSequenceLength
        ) {

    // TODO replace w/ error check
    if (sequence.size() > maxSequenceLength) {
        return {};
    }

    const QMap<QChar, int> aaIndexes = aminoAcidIndex();

    QVector<int> arr(maxSequenceLength, 0);

    for (int i = 0; i < sequence.size(); i++) {
        const QChar aa = sequence.at(i);
        arr[i] =  aaIndexes.value(aa);
    }

    return arr;
}

namespace {

    struct IonLengths {
        int lengthCharge1Ions = -1;
        int lengthCharge2Ions = -1;
        int lengthAIons = -1;
    };

    IonLengths ionLengthsByCharge(int charge) {

        IonLengths ionLen;

        switch (charge) {

            case 1:
                ionLen.lengthCharge1Ions = 22;
                ionLen.lengthAIons = 22;
                return ionLen;

            case 2:
                ionLen.lengthCharge1Ions = 40;
                ionLen.lengthCharge2Ions = 40;
                ionLen.lengthAIons = 40;
                return ionLen;

            case 3:
            case 4:
                ionLen.lengthCharge1Ions = 50;
                ionLen.lengthCharge2Ions = 50;
                ionLen.lengthAIons = 50;
                return ionLen;

            default:
                return ionLen;
        }
    }

    QStringList buildIonList(
            int length,
            const TandemPredictionUtils::IonType &ionType,
            const TandemPredictionUtils::IonModifier &ionModifier
            ) {

        const QString stringTemplate = QStringLiteral("%1%2%3");

        QStringList ions;
        for (int i = 0; i < length; i++){

            const QString ionTypeName =  StringUtils::enumClassToString(ionType);
            const QString ionModifierName = StringUtils::enumClassToString(ionModifier);

            QString ionModifierNameTemplate = ionModifier == TandemPredictionUtils::IonModifier::None
                                              ? QString()
                                              : "-" + ionModifierName;

            ionModifierNameTemplate = ionModifier == TandemPredictionUtils::IonModifier::Chrg2 || ionModifier == TandemPredictionUtils::IonModifier::Chrg3
                                      ? ionModifierNameTemplate.replace(QStringLiteral("-Chrg"), "^")
                                      : ionModifierNameTemplate;

            const QString ionName = stringTemplate.arg(ionTypeName).arg(i + 1).arg(ionModifierNameTemplate);
            ions.push_back(ionName);
        }

        return ions;
    }

}//NAMESPACE
QStringList TandemPredictionUtils::buildIonLabels(int charge) {

    const IonLengths il = ionLengthsByCharge(charge);

    QStringList ionList;

    if (charge < 3) {
        ionList.append({"p",  "p-H2O", "p-NH3"});
    }

    ionList.append(buildIonList(il.lengthCharge1Ions,
                                TandemPredictionUtils::IonType::y,
                                TandemPredictionUtils::IonModifier::None
                                ));
    if (charge > 1) {
        ionList.append(buildIonList(il.lengthCharge2Ions,
                                    TandemPredictionUtils::IonType::y,
                                    TandemPredictionUtils::IonModifier::Chrg2
                                    ));
    }
    ionList.append(buildIonList(il.lengthCharge1Ions,
                                TandemPredictionUtils::IonType::y,
                                TandemPredictionUtils::IonModifier::NH3));

    ionList.append(buildIonList(il.lengthCharge1Ions,
                                TandemPredictionUtils::IonType::y,
                                TandemPredictionUtils::IonModifier::H2O
                                ));

    ionList.append(buildIonList(il.lengthCharge1Ions,
                                TandemPredictionUtils::IonType::b,
                                TandemPredictionUtils::IonModifier::None
                                ));
    if (charge > 1) {
        ionList.append(buildIonList(il.lengthCharge2Ions,
                                    TandemPredictionUtils::IonType::b,
                                    TandemPredictionUtils::IonModifier::Chrg2
                                    ));
    }
    ionList.append(buildIonList(il.lengthCharge1Ions,
                                TandemPredictionUtils::IonType::b,
                                TandemPredictionUtils::IonModifier::NH3
                                ));

    ionList.append(buildIonList(il.lengthCharge1Ions,
                                TandemPredictionUtils::IonType::b,
                                TandemPredictionUtils::IonModifier::H2O
                                ));


    ionList.append(buildIonList(il.lengthAIons,
                                TandemPredictionUtils::IonType::a,
                                TandemPredictionUtils::IonModifier::None
                                ));

    return ionList;
}


void TandemPredictionUtils::filterFragmentIonVectorByIntensity(
        double intensityThreshold,
        QVector<FragmentIon> *fragIons
        ) {

    const auto filterLogic = [intensityThreshold](const FragmentIon &fi){
        return fi.intensity < intensityThreshold;
    };

    const auto terminator = std::remove_if(fragIons->begin(),  fragIons->end(), filterLogic);
    fragIons->erase(terminator, fragIons->end());
}


bool TandemPredictionUtils::isValidPeptideForPrediction(
        const QString &peptideSeq,
        int charge
        ) {
    const int maxPredictionLengthForCharge = TandemPredictionUtils::getMaxLengthForPredictionInputByCharge(charge);
    return peptideSeq.size() <= maxPredictionLengthForCharge;
}


int TandemPredictionUtils::calculateNormalizedCollisionEnergy(
        double mz,
        int charge,
        int collisionEnergy
        ) {

//    http://proteomicsnews.blogspot.com/2014/06/normalized-collision-energy-calculation.html
//    Absolute energy (eV) = (settling NCE) x (Isolation center) / (500 m/z) x (charge factor)
//    Charge factors are as so:
//        1/1.0
//        2/0.9
//        3/0.85
//        4/0.8
//        5/0.75
//        >5/0.75

    const double thermoMzNormalizer = 500.0;

    double chargeFactor = 1.0;

    switch (charge) {

        case 1:
            break;
        case 2:
            chargeFactor = 0.9;
            break;
        case 3:
            chargeFactor = 0.85;
            break;
        case 4:
            chargeFactor = 0.8;
            break;
        default:
            chargeFactor = 0.75;

    }

    return static_cast<int>(collisionEnergy / ( (mz / thermoMzNormalizer) * chargeFactor ));
}

void TandemPredictionUtils::sortPredictionByIonLabel(QVector<FragmentIon> *frags) {

    const auto sortLogic = [](const FragmentIon &l, const FragmentIon &r){
        return l.ionLabel < r.ionLabel;
    };

    std::sort(frags->begin(), frags->end(), sortLogic);
}


namespace {

    Err modificationsMapFromString(
            const QString &modString,
            QHash<ResidueIndex, ModificationMass> *modifications
            ) {

        ERR_INIT

        const int expectedModSplitSize = 2;

        const QStringList modStringSplit
                = modString.split(S_GLOBAL_SETTINGS.SEPARATOR, QString::SkipEmptyParts);

        for (const QString &mod : modStringSplit) {

            const QStringList modSplit
                    = mod.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP, QString::SkipEmptyParts);

            e = ErrorUtils::isEqual(modSplit.size(), expectedModSplitSize); ree;

            ResidueIndex residueIndex;
            e = ErrorUtils::toInt(modSplit.front(), &residueIndex); ree;

            ModificationMass modificationMass;
            e = ErrorUtils::toDouble(modSplit.back(), &modificationMass); ree;

            modifications->insert(residueIndex, modificationMass);
        }

        ERR_RETURN
    }


    Eigen::VectorXd buildModificationVector(
            int pepLength,
            const QHash<ResidueIndex, ModificationMass> &modifications,
            bool isYSeries
            ) {

        Eigen::VectorXd modVec(pepLength);
        modVec.setZero();

        if (isYSeries) {

            for (auto it = modifications.begin(); it != modifications.end(); it++) {

                const int modIndex = it.key();
                const double modMass = it.value();

                for (int i = modIndex; i >= 0; --i) {
                    modVec.coeffRef(i) = modMass;
                }
            }

            return modVec.reverse();
        }

        for (auto it = modifications.begin(); it != modifications.end(); it++) {

            const int modIndex = it.key();
            const double modMass = it.value();

            for (int i = modIndex; i < pepLength; ++i) {
                modVec.coeffRef(i) = modMass;
            }
        }
        return modVec;
    }


    void addModificatonSeriesToVector(const QHash<ResidueIndex, ModificationMass> &modifications,
                                      int peptideSequenceSize,
                                      bool isYSeries,
                                      Charge charge,
                                      QVector<double> *vec) {

        Eigen::VectorXd seriesModVec = buildModificationVector(
                peptideSequenceSize,
                modifications,
                isYSeries
                );

        seriesModVec /= charge;

        for (int i = 0; i < seriesModVec.size(); ++i) {
            (*vec)[i] += seriesModVec.coeff(i);
        }
    }

}//namespace
Err TandemPredictionUtils::calculateMzValuesForPrediction(
        const QString &peptideSequence,
        const QString &modificatonString,
        const AminoAcids &aminoAcids,
        int charge,
        QVector<double> *vec
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(peptideSequence); ree;

    QHash<ResidueIndex, ModificationMass> modifications;
    e = modificationsMapFromString(modificatonString, &modifications); ree;

    const double waterMonoisotopicMass = Molecule(MolecularFormulas::waterFormula).monoisotopicMass();
    const double nh3MonoisotopicMass = Molecule(MolecularFormulas::ammoniaFormula).monoisotopicMass();
    const IonLengths il = ionLengthsByCharge(charge);
    const int maxPeptideLength = getMaxLengthForPredictionInputByCharge(charge);
    const int peptideSequenceSize = peptideSequence.size();
    QString peptideSequenceReversed = peptideSequence;
    std::reverse(peptideSequenceReversed.begin(),  peptideSequenceReversed.end());

    if (peptideSequence.size() > maxPeptideLength) {
        return Error::eError;
    }

    QVector<double> bySeries;

    if (charge < 3) {
        const double precursorMass = BiophysicalCalcs::calculateThomson(
                peptideSequence,
                aminoAcids,
                charge
                );

        const double precursorMassDesH2O = precursorMass - waterMonoisotopicMass;
        const double precursorMassDesNH3 = precursorMass - nh3MonoisotopicMass;

        bySeries.append(precursorMass);
        bySeries.append(precursorMassDesH2O);
        bySeries.append(precursorMassDesNH3);
    }

    const Charge charge1 = 1;
    const Charge charge2 = 2;

    QVector<double> ySeries = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideSequenceReversed,
            charge1,
            PROTON + waterMonoisotopicMass,
            il.lengthCharge1Ions,
            aminoAcids
            );

    addModificatonSeriesToVector(modifications, peptideSequenceSize, true, charge1, &ySeries);
    bySeries.append(ySeries);

    if (charge > 1) {

        QVector<double> ySeriesCharge2 = BiophysicalCalcs::buildTandemFragmentMasses(
                peptideSequenceReversed,
                charge2,
                (PROTON * charge2) +
                waterMonoisotopicMass,
                il.lengthCharge2Ions,
                aminoAcids
                );
        addModificatonSeriesToVector(
                modifications,
                peptideSequenceSize,
                true,
                charge2,
                &ySeriesCharge2
                );
        bySeries.append(ySeriesCharge2);
    }

    QVector<double> ySeriesNH3 = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideSequenceReversed,
            charge1,
            PROTON + waterMonoisotopicMass - nh3MonoisotopicMass,
            il.lengthCharge1Ions,
            aminoAcids
            );

    addModificatonSeriesToVector(
            modifications,
            peptideSequenceSize,
            true,
            charge1,
            &ySeriesNH3
            );

    bySeries.append(ySeriesNH3);


    QVector<double> ySeriesH2O = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideSequenceReversed,
            charge1,
            PROTON,
            il.lengthCharge1Ions,
            aminoAcids
            );

    addModificatonSeriesToVector(modifications, peptideSequenceSize, true, charge1, &ySeriesH2O);
    bySeries.append(ySeriesH2O);


    ///////////////
    QVector<double> bSeries = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideSequence,
            charge1,
            PROTON,
            il.lengthCharge1Ions,
            aminoAcids
            );

    const int lastResidueIndex = peptideSequence.size() - 1;

    if (lastResidueIndex < bSeries.size()) {
        bSeries[lastResidueIndex] += waterMonoisotopicMass;
    }

    addModificatonSeriesToVector(modifications, peptideSequenceSize, false, charge1, &bSeries);
    bySeries.append(bSeries);

    if (charge > 1) {

        QVector<double> bSeriesCharge2 = BiophysicalCalcs::buildTandemFragmentMasses(
                peptideSequence,
                charge2,
                (PROTON * charge2),
                il.lengthCharge2Ions,
                aminoAcids
                );
        if (lastResidueIndex < bSeries.size()) {
            bSeriesCharge2[peptideSequence.size() - 1] += (waterMonoisotopicMass / charge2);
        }
        addModificatonSeriesToVector(
                modifications,
                peptideSequenceSize,
                false,
                charge2,
                &bSeriesCharge2
                );
        bySeries.append(bSeriesCharge2);
    }

    QVector<double> bSeriesNH3 = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideSequence,
            charge1,
            PROTON - nh3MonoisotopicMass,
            il.lengthCharge1Ions,
            aminoAcids
            );
    if (lastResidueIndex < bSeries.size()) {
        bSeriesNH3[peptideSequence.size() - 1] += waterMonoisotopicMass;
    }
    addModificatonSeriesToVector(
            modifications,
            peptideSequenceSize,
            false,
            charge1,
            &bSeriesNH3
            );
    bySeries.append(bSeriesNH3);

    QVector<double> bSeriesH2O = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideSequence,
            charge1,
            PROTON - waterMonoisotopicMass,
            il.lengthCharge1Ions,
            aminoAcids
            );
    if (lastResidueIndex < bSeries.size()) {
        bSeriesH2O[peptideSequence.size() - 1] += waterMonoisotopicMass;
    }
    addModificatonSeriesToVector(
            modifications,
            peptideSequenceSize,
            false,
            charge1,
            &bSeriesH2O
            );
    bySeries.append(bSeriesH2O);

    //////////////
    QVector<double> aSeries = BiophysicalCalcs::buildTandemFragmentMasses(
            peptideSequence,
            charge1,
            PROTON - (MonoIsotopicMass::CARBON + MonoIsotopicMass::OXYGEN),
            il.lengthAIons,
            aminoAcids
            );

    addModificatonSeriesToVector(
            modifications,
            peptideSequenceSize,
            false,
            charge1,
            &aSeries
            );
    bySeries.append(aSeries);

    *vec = bySeries;
    ERR_RETURN
}


namespace {

    template <typename T>
    void removeMzMassesOutOfRange(
            double mzMin,
            double mzMax,
            T *bySeries
            ) {

        const auto terminatorLogic = [mzMin, mzMax](double mz) {
            return mz < mzMin || mz > mzMax;
        };

        const auto terminator = std::remove_if(bySeries->begin(),  bySeries->end(), terminatorLogic);
        bySeries->erase(terminator, bySeries->end());
    }

    Eigen::VectorXd buildAminoAcidMassesVec(
            const QString &seq,
            const AminoAcids &aa
            ) {

        const int sequenceLength = seq.size();

        Eigen::VectorXd aminoAcidsMassValues(sequenceLength);
        aminoAcidsMassValues.setZero();

        for(int i = 0; i < seq.size(); i++){
            const QChar aminoAcid = seq.at(i);
            aminoAcidsMassValues.coeffRef(i) = aa.aminoAcid(aminoAcid).monoisotopicMass();
        }

        return aminoAcidsMassValues;
    }


    Eigen::VectorXd buildThomsonValues(
            const QString &peptideSequence,
            const AminoAcids &aminoAcids,
            bool reverseSequence
            ) {

        Eigen::VectorXd aminoAcidsMassValues = buildAminoAcidMassesVec(peptideSequence, aminoAcids);
        if (reverseSequence) {
            std::reverse(aminoAcidsMassValues.begin(), aminoAcidsMassValues.end());
        }

        Eigen::VectorXd cumSumSeries(peptideSequence.size());
        std::partial_sum(aminoAcidsMassValues.begin(), aminoAcidsMassValues.end(), cumSumSeries.begin(), std::plus<double>());

        return cumSumSeries;
    }

} //END NAMESPACE
QVector<double> TandemPredictionUtils::buildTheoreticalMzListBYOnly(
        const QString &seq,
        double mzMin,
        double mzMax,
        const AminoAcids &aa
        ) {

    QVector<double> bySeries;
    //TODO see if you can delete this method.  It is only used in ScanCalibratron
    const double waterMonoisotopicMass = Molecule(MolecularFormulas::waterFormula).monoisotopicMass();

    Eigen::VectorXd cumSumBSeries = buildThomsonValues(seq, aa, false);
    cumSumBSeries = cumSumBSeries.array() + PROTON;
    cumSumBSeries.coeffRef(cumSumBSeries.size() - 1) += waterMonoisotopicMass; //TODO comment this line out and uncomment other lines when TandemScanLibrary is implemented
    std::vector<double> bSeriesVec(cumSumBSeries.data(), cumSumBSeries.data() + cumSumBSeries.size());
//    bSeriesVec.pop_back();
    bySeries.append(QVector<double>::fromStdVector(bSeriesVec));

    Eigen::VectorXd cumSumYSeries = buildThomsonValues(seq, aa, true);
    cumSumYSeries = cumSumYSeries.array() + PROTON + waterMonoisotopicMass;
    std::vector<double> ySeriesVec(cumSumYSeries.data(), cumSumYSeries.data() + cumSumYSeries.size());
//    ySeriesVec.pop_front();
    bySeries.append(QVector<double>::fromStdVector(ySeriesVec));

//#define USE_CHARGE2_IONS_FRAGGING
#ifdef USE_CHARGE2_IONS_FRAGGING
    Eigen::VectorXd cumSumBSeriesCharge2 = (cumSumBSeries.array() + PROTON) / 2.0;
    QVector<double> bSeriesVecCharge2(cumSumBSeriesCharge2.data(), cumSumBSeriesCharge2.data() + cumSumBSeriesCharge2.size());
    bSeriesVecCharge2.pop_back();
    bySeries.append(bSeriesVecCharge2);

    Eigen::VectorXd cumSumYSeriesCharge2 = (cumSumYSeries.array() + PROTON) / 2.0;
    QVector<double> ySeriesVecCharge2(cumSumYSeriesCharge2.data(), cumSumYSeriesCharge2.data() + cumSumYSeriesCharge2.size());
    bySeries.append(ySeriesVecCharge2);
#endif

    removeMzMassesOutOfRange(mzMin, mzMax, &bySeries);
    return bySeries;
}

QPair<QVector<double>, QVector<double>> TandemPredictionUtils::buildTheoreticalMzListBYSeparated(
        const QString &seq,
        double mzMin,
        double mzMax,
        const AminoAcids &aa,
        int charge,
        const QHash<ResidueIndex, ModificationMass> &modifications
) {

    const double waterMonoisotopicMass = Molecule(MolecularFormulas::waterFormula).monoisotopicMass();
    charge = std::max(1, charge);

    Eigen::VectorXd bSeriesModVec = buildModificationVector(seq.size(), modifications, false);
    Eigen::VectorXd ySeriesModVec = buildModificationVector(seq.size(), modifications, true);

    Eigen::VectorXd cumSumBSeries = buildThomsonValues(seq, aa, false);
    cumSumBSeries = (cumSumBSeries.array() + (PROTON * charge) + bSeriesModVec.array()) / charge;

    std::vector<double> bSeriesStdVec(cumSumBSeries.data(), cumSumBSeries.data() + cumSumBSeries.size());
    bSeriesStdVec.pop_back();

    Eigen::VectorXd cumSumYSeries = buildThomsonValues(seq, aa, true);
    cumSumYSeries = (cumSumYSeries.array() + (PROTON * charge) + waterMonoisotopicMass + ySeriesModVec.array()) / charge;

    std::vector<double> ySeriesStdVec(cumSumYSeries.data(), cumSumYSeries.data() + cumSumYSeries.size());
    ySeriesStdVec.pop_back();

    removeMzMassesOutOfRange(mzMin, mzMax, &bSeriesStdVec);
    removeMzMassesOutOfRange(mzMin, mzMax, &ySeriesStdVec);

    const QVector<double> bSeriesVec = QVector<double>::fromStdVector(bSeriesStdVec);
    const QVector<double> ySeriesVec = QVector<double>::fromStdVector(ySeriesStdVec);

    return {bSeriesVec, ySeriesVec};
}

Err TandemPredictionUtils::extractSequenceAndChargeFromPeptideSequenceChargeKey(
        const PeptideSequenceChargeKey &peptideSequenceChargeKey,
        PeptideString *peptideString,
        int *charge
        ) {

    ERR_INIT

    const QStringList peptideSequenceChargeKeySplit
            = peptideSequenceChargeKey.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP);

    const int expectedSplitSize = 2;
    e = ErrorUtils::isEqual(
            peptideSequenceChargeKeySplit.size(),
            expectedSplitSize
            ); ree;

    *peptideString = peptideSequenceChargeKeySplit.first();

    e = ErrorUtils::toInt(
            peptideSequenceChargeKeySplit.back(),
            charge
            ); ree;

    ERR_RETURN
}
