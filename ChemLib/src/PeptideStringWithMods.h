//
// Created by anichols on 2/23/24.
//

#ifndef PYTHIADIACPP_PEPTIDESTRINGWITHMODS_H
#define PYTHIADIACPP_PEPTIDESTRINGWITHMODS_H

#include "ChemLib_Exports.h"
#include "GlobalSettings.h"

#include <QString>


class  CHEMLIB_EXPORTS PeptideStringWithMods : public QString {

public:

    PeptideStringWithMods();

    [[nodiscard]] int sizeNoMods() const;

    explicit PeptideStringWithMods(const QString &seq);

    ~PeptideStringWithMods() = default;

    [[nodiscard]] PeptideString removeUniModChars() const;

    QMap<Index, Mass> modificationsMap() const;





private:


};


#endif //PYTHIADIACPP_PEPTIDESTRINGWITHMODS_H
