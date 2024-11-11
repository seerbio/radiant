//
// Created by anichols on 2/23/24.
//

#ifndef PYTHIADIACPP_PEPTIDESTRINGWITHMODS_H
#define PYTHIADIACPP_PEPTIDESTRINGWITHMODS_H

#include "ChemLib_Exports.h"
#include "GlobalSettings.h"

#include <QString>

class AminoAcids;

class  CHEMLIB_EXPORTS PeptideStringWithMods : public QString {

public:

    PeptideStringWithMods();

    [[nodiscard]] int sizeNoMods() const;

    explicit PeptideStringWithMods(const QString &seq);

    ~PeptideStringWithMods() = default;

    [[nodiscard]] PeptideString removeUniModChars() const;

    [[nodiscard]] QMap<Index, double> modificationsMap() const;

    [[nodiscard]] QVector<double> bSeries(
            int charge,
            const AminoAcids &aminoAcids
            ) const;

    [[nodiscard]] QVector<double> ySeries(
            int charge,
            const AminoAcids &aminoAcids
            ) const;

    [[nodiscard]] QStringList bSeriesIonLabels(const QString &modifier) const;

    [[nodiscard]] QStringList ySeriesIonLabels(const QString &modifier) const;


private:


};


#endif //PYTHIADIACPP_PEPTIDESTRINGWITHMODS_H
