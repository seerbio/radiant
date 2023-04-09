//
// Created by anichols on 4/3/23.
//

#ifndef PYTHIADIACPP_DIAMZTARGETSREADER_H
#define PYTHIADIACPP_DIAMZTARGETSREADER_H

#include "CSVReader.h"
#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"

using namespace Error;

namespace DIAMzTargetsReaderNamespace {

    const QString TARGET_MZ = QStringLiteral("targetMz");
    const QString ISO_WIN_LO = QStringLiteral("isoWindowLo");
    const QString ISO_WIN_HI = QStringLiteral("isoWindowHi");
    const QString COLLISION_ENERGY = QStringLiteral("collisionEnergy");

    const QStringList keysToCheck = {
            TARGET_MZ,
            ISO_WIN_LO,
            ISO_WIN_HI,
            COLLISION_ENERGY
    };
}

struct DIAMzTargetsReaderRow : public CSVReaderInputBase {

    double targetMz = -1.0;
    double isoWindowLo = -1.0;
    double isoWindowHi = -1.0;
    double collisionEnergy = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace DIAMzTargetsReaderNamespace;

        return {
            {TARGET_MZ, targetMz},
            {ISO_WIN_LO, isoWindowLo},
            {ISO_WIN_HI, isoWindowHi},
            {COLLISION_ENERGY, collisionEnergy}
        };
    }

    Err initFromRead(const CSVReaderInputBase &row) override {

        using namespace DIAMzTargetsReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = CSVReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        targetMz = dataMap.value(TARGET_MZ).toDouble();
        collisionEnergy = dataMap.value(COLLISION_ENERGY).toDouble();
        isoWindowLo = dataMap.value(ISO_WIN_LO).toDouble();
        isoWindowHi = dataMap.value(ISO_WIN_HI).toDouble();

        ERR_RETURN
    }
    
};


class FILEREADERSLIB_EXPORTS DIAMzTargetsReader {

public:

    DIAMzTargetsReader() = default;
    ~DIAMzTargetsReader() = default;

    Err write(
            const QVector<DIAMzTargetsReaderRow> &diaMzTargetReaderRows,
            const QString &outputFilePath
            );

    Err read(
            const QString &fileURI,
            QVector<DIAMzTargetsReaderRow> *diaMzTargetReaderRows
            );

};

#endif //PYTHIADIACPP_DIAMZTARGETSREADER_H
