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
    const QString MZ_VALS = QStringLiteral("mzVals");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");
    const QString ION_LABELS = QStringLiteral("ionLabels");
    const QString IS_DECOY = QStringLiteral("isDecoy");

    const QStringList keysToCheck = {
            PEP_SEQ_CHRG_KEY,
            MZ_VALS,
            INTENSITY_VALS,
            ION_LABELS,
            IS_DECOY
    };
}

struct FILEREADERSLIB_EXPORTS TandemLibraryReaderRow : public ParquetReaderInputBase {

    PeptideSequenceChargeKey peptideSequenceChargeKey;
    QVector<double> mzVals;
    QVector<double> intensityVals;
    QStringList ionLabels;
    bool isDecoy = false;

    QMap<QString, QVariant> map() override {

        using namespace TandemLibraryReaderNamespace;

        return {
            {PEP_SEQ_CHRG_KEY, QVariant(peptideSequenceChargeKey)},
            {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
            {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
            {ION_LABELS, QVariant(joinQStringList(ionLabels))},
            {IS_DECOY, QVariant(isDecoy)},
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
        mzVals = bytesArrayToQVector<double>(dataMap.value(MZ_VALS).toByteArray());
        intensityVals = bytesArrayToQVector<double>(dataMap.value(INTENSITY_VALS).toByteArray());
        ionLabels = dataMap.value(ION_LABELS).toString().split(S_GLOBAL_SETTINGS.SEPARATOR);
        isDecoy = dataMap.value(IS_DECOY).toBool();

        ERR_RETURN
    }

};


#endif //PYTHIADIACPP_TANDEMLIBRARYREADER_H
