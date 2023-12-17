//
// Created by anichols on 7/24/23.
//

#ifndef PYTHIADIACPP_MSCALIBRATIONREADER_H
#define PYTHIADIACPP_MSCALIBRATIONREADER_H


#include "ParquetReader.h"
#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"

using namespace Error;

namespace MsCalibarationReaderNamespace {

    const QString PEP_STR_MODS = QStringLiteral("peptideStringWithMods");
    const QString IRT_PRED = QStringLiteral("iRTPredicted");
    const QString SCAN_TIME = QStringLiteral("scanTime");
    const QString MZ_SRCH_V = QStringLiteral("mzSearchedVec");
    const QString MZ_FND_MEAN_V = QStringLiteral("mzFoundMeanVec");
    const QString MZ_FND_STDEV_V = QStringLiteral("mzFoundStDevVec");

    const QStringList keysToCheck = {
            PEP_STR_MODS,
            IRT_PRED,
            SCAN_TIME,
            MZ_SRCH_V,
            MZ_FND_MEAN_V,
            MZ_FND_STDEV_V
    };
}

struct MsCalibarationReaderRow: public ParquetReaderInputBase {

    PeptideStringWithMods peptideStringWithMods;
    IRT iRTPredicted = -1.0;
    ScanTime scanTime = -1.0;
    QVector<float> mzSearchedVec;
    QVector<float> mzFoundMeanVec;
    QVector<float> mzFoundStDevVec;

    QMap<QString, QVariant> map() override {

        using namespace MsCalibarationReaderNamespace;

        return {
                {PEP_STR_MODS, QVariant(peptideStringWithMods)},
                {IRT_PRED, QVariant(iRTPredicted)},
                {SCAN_TIME, QVariant(scanTime)},
                {MZ_SRCH_V, QVariant(qVectorToQByteArray(mzSearchedVec))},
                {MZ_FND_MEAN_V, QVariant(qVectorToQByteArray(mzFoundMeanVec))},
                {MZ_FND_STDEV_V, QVariant(qVectorToQByteArray(mzFoundStDevVec))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace MsCalibarationReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        if (!allKeysPresent) {
            qDebug() << dataMap;
            qDebug() << keysToCheck;
        }

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        peptideStringWithMods = dataMap.value(PEP_STR_MODS).toString();
        iRTPredicted = dataMap.value(IRT_PRED).toDouble();
        scanTime = dataMap.value(SCAN_TIME).toDouble();
        mzSearchedVec = bytesArrayToQVector<float>(dataMap.value(MZ_SRCH_V).toByteArray());
        mzFoundMeanVec = bytesArrayToQVector<float>(dataMap.value(MZ_FND_MEAN_V).toByteArray());
        mzFoundStDevVec = bytesArrayToQVector<float>(dataMap.value(MZ_FND_STDEV_V).toByteArray());

        ERR_RETURN
    }

};

#endif //PYTHIADIACPP_MSCALIBRATIONREADER_H
