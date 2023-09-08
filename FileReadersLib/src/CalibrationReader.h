//
// Created by anichols on 7/24/23.
//

#ifndef PYTHIADIACPP_CALIBRATIONREADER_H
#define PYTHIADIACPP_CALIBRATIONREADER_H


#include "CSVReader.h"
#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"

using namespace Error;

namespace CalibarationReaderNamespace {

    const QString PEP_STR_MODS = QStringLiteral("peptideStringWithMods");
    const QString CHARGE = QStringLiteral("charge");
    const QString SCAN_TIME = QStringLiteral("scanTime");

    const QStringList keysToCheck = {
            PEP_STR_MODS,
            CHARGE,
            SCAN_TIME
    };
}

struct CalibarationReaderRow: public CSVReaderInputBase {

    PeptideStringWithMods peptideStringWithMods;
    Charge charge = -1;
    ScanTime scanTime = -1.0;

    QMap<QString, QVariant> map() override {

        using namespace CalibarationReaderNamespace;

        return {
                {PEP_STR_MODS, peptideStringWithMods},
                {CHARGE, charge},
                {SCAN_TIME, scanTime}
        };
    }

    Err initFromRead(const CSVReaderInputBase &row) override {

        using namespace CalibarationReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = CSVReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        if (!allKeysPresent) {
            qDebug() << dataMap;
            qDebug() << keysToCheck;
        }

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        peptideStringWithMods = dataMap.value(PEP_STR_MODS).toString();
        charge = dataMap.value(CHARGE).toInt();
        scanTime = dataMap.value(SCAN_TIME).toDouble();

        ERR_RETURN
    }

};





#endif //PYTHIADIACPP_CALIBRATIONREADER_H
