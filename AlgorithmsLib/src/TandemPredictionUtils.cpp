//
// Created by Drucifer on 12/7/2021.
//

#include "TandemPredictionUtils.h"

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

QVector<int> TandemPredictionUtils::peptideSequenceToArray(const QString &sequence, int maxSequenceLength) {

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

