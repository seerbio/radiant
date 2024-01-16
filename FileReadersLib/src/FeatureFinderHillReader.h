//
// Created by anichols on 5/3/23.
//

#ifndef PYTHIADIACPP_FEATUREFINDERHILLREADER_H
#define PYTHIADIACPP_FEATUREFINDERHILLREADER_H

#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

using namespace Error;

namespace FeatureFinderHillReaderNamespace {

    const QString MZ_VALS = QStringLiteral("mzVals");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");
    const QString SCAN_NUBMERS = QStringLiteral("scanNumbers");
    const QString SCAN_NUMBER_INDEXES = QStringLiteral("scanNumberIndexes");

    const QStringList keysToCheck = {
            MZ_VALS,
            INTENSITY_VALS,
            SCAN_NUBMERS,
            SCAN_NUMBER_INDEXES
    };
}

/**
 * See ParquetReaderInputBase for documentation
 */
struct FILEREADERSLIB_EXPORTS FeatureFinderHillReaderRow  : public ParquetReaderInputBase {

    QVector<double> mzVals;
    QVector<double> intensityVals;
    QVector<int> scanNumbers;
    QVector<int> scanNumberIndexes;

    QMap<QString, QVariant> map() override {

        using namespace FeatureFinderHillReaderNamespace;

        return {
                {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
                {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
                {SCAN_NUBMERS, QVariant(qVectorToQByteArray(scanNumbers))},
                {SCAN_NUMBER_INDEXES, QVariant(qVectorToQByteArray(scanNumberIndexes))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace FeatureFinderHillReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        mzVals = bytesArrayToQVector<double>(dataMap.value(MZ_VALS).toByteArray());
        intensityVals = bytesArrayToQVector<double>(dataMap.value(INTENSITY_VALS).toByteArray());
        scanNumbers = bytesArrayToQVector<int>(dataMap.value(SCAN_NUBMERS).toByteArray());
        scanNumberIndexes = bytesArrayToQVector<int>(dataMap.value(SCAN_NUMBER_INDEXES).toByteArray());

        ERR_RETURN
    }

};



#endif //PYTHIADIACPP_FEATUREFINDERHILLREADER_H
