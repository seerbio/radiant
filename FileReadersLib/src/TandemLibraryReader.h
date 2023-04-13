//
// Created by anichols on 4/3/23.
//

#ifndef PYTHIADIACPP_TANDEMLIBRARYREADER_H
#define PYTHIADIACPP_TANDEMLIBRARYREADER_H

#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

using namespace Error;

namespace TandemLibraryReaderNamespace {

    const QString PEP_SEQ_CHRG_KEY = QStringLiteral("peptideSequenceChargeKey");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");
    const QString ION_LABELS = QStringLiteral("ionLabels");

    const QStringList keysToCheck = {
            PEP_SEQ_CHRG_KEY,
            INTENSITY_VALS,
            ION_LABELS
    };
}

struct FILEREADERSLIB_EXPORTS TandemLibraryReaderRow : public ParquetReaderInputBase {

    PeptideSequenceChargeKey peptideSequenceChargeKey;
    QVector<double> intensityVals;
    QStringList ionLabels;

    QMap<QString, QVariant> map() override {

        using namespace TandemLibraryReaderNamespace;

        return {
            {PEP_SEQ_CHRG_KEY, QVariant(peptideSequenceChargeKey)},
            {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
            {ION_LABELS, QVariant(joinQStringList(ionLabels))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace TandemLibraryReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
                );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        peptideSequenceChargeKey = dataMap.value(PEP_SEQ_CHRG_KEY).toString();
        intensityVals = bytesArrayToQVector<double>(dataMap.value(INTENSITY_VALS).toByteArray());
        ionLabels = dataMap.value(ION_LABELS).toString().split(S_GLOBAL_SETTINGS.SEPARATOR);

        ERR_RETURN
    }

};


#endif //PYTHIADIACPP_TANDEMLIBRARYREADER_H
