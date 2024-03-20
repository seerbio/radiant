//
// Created by anichols on 2/23/24.
//

#include "PeptideStringWithMods.h"

#include "AminoAcids.h"
#include "BiophysicalCalcs.h"
#include "MolecularFormula.h"

PeptideStringWithMods::PeptideStringWithMods(const QString &seq)
: QString(seq) {}

PeptideStringWithMods::PeptideStringWithMods() {}

int PeptideStringWithMods::sizeNoMods() const {
    const PeptideString peptideString = this->removeUniModChars();
    return peptideString.size();
}

PeptideString PeptideStringWithMods::removeUniModChars() const {

    PeptideString peptideString;

    bool unimodOn = false;
    for (const QChar &c : *this) {
        if (c == "(" || c == "[") {
            unimodOn = true;
            continue;
        }
        else if (c == ")" || c == "]" ) {
            unimodOn = false;
            continue;
        }
        else if (unimodOn) {
            continue;
        }

        peptideString += c;
    }

    return peptideString;
}


QMap<Index, Mass> PeptideStringWithMods::modificationsMap() const {

    QMap<Index, Mass> modMap;

    bool unimodOn = false;
    QString modBuilder;

    Index index = 0;

    for (int i = 0; i < this->size(); i++) {

        const QChar &c = this->at(i);

        if (c == "(" || c == "[") {
            unimodOn = true;
            index = index - 1;
            continue;
        }
        else if (c == ")" || c == "]" ) {
            unimodOn = false;

            bool isNumber;
            const double num = modBuilder.toDouble(&isNumber);

            if (isNumber) {
                modMap.insert(index, num);
            }

            else {
                 if (UniModNamespace::uniModNameVsModificationMass.contains(modBuilder)) {
                     modMap.insert(index, UniModNamespace::uniModNameVsModificationMass.value(modBuilder));
                 }
            }

            modBuilder.clear();
            index++;

            continue;
        }
        else if (unimodOn) {
            modBuilder += c;
        }
        else {
            index++;
        }
    }

    return modMap;

}

QVector<double> PeptideStringWithMods::bSeries(
        int charge,
        const AminoAcids &aminoAcids
        ) const {
    return BiophysicalCalcs::buildTandemFragmentMasses(
            *this,
            BiophysicalCalcs::FragmentSeriesType::bSeries,
            charge,
            aminoAcids
            );
}

QVector<double> PeptideStringWithMods::ySeries(
        int charge,
        const AminoAcids &aminoAcids
        ) const {

    return BiophysicalCalcs::buildTandemFragmentMasses(
            *this,
            BiophysicalCalcs::FragmentSeriesType::ySeries,
            charge,
            aminoAcids
            );
}

QStringList PeptideStringWithMods::bSeriesIonLabels(const QString &modifier) const {
    QStringList bLabels;
    for (int i = 0; i < this->sizeNoMods(); ++i) {
        const QString ionLabel = 'b' + QString::number(i + 1) + modifier;
        bLabels.push_back(ionLabel);
    }
    return bLabels;
}

QStringList PeptideStringWithMods::ySeriesIonLabels(const QString &modifier) const {
    QStringList yLabels;
    for (int i = 0; i < this->sizeNoMods(); ++i) {
        const QString ionLabel = 'y' + QString::number(i + 1) + modifier;
        yLabels.push_back(ionLabel);
    }
    return yLabels;
}
