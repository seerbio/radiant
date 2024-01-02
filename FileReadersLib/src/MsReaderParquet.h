//
// Created by anichols on 2/6/23.
//

#ifndef SPARKDIA_PARQUETREADER_H
#define SPARKDIA_PARQUETREADER_H

#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "ParquetReader.h"

#include <QMap>
#include <vector>

using namespace Error;

namespace MsParquetReaderNamespace {

    const QString PRQ_DF_SUFFIX = QStringLiteral("prqDF");

    const QString MS_LEVEL = QStringLiteral("msLevel");
    const QString SCAN_NUMBER = QStringLiteral("scanNumber");
    const QString SCAN_TIME = QStringLiteral("scanTime");
    const QString COLLISION_ENERGY = QStringLiteral("collisionEnergy");
    const QString PERCURSOR_TARGET_MZ = QStringLiteral("precursorTargetMz");
    const QString ISO_WINDOW_LOWER = QStringLiteral("isoWindowLower");
    const QString ISO_WINDOW_UPPER = QStringLiteral("isoWindowUpper");
    const QString IM_DRIFT_TIME = QStringLiteral("ionMobilityDriftTime");
    const QString IM_IND = QStringLiteral("ionMobilityIndex");
    const QString MZ_VALS = QStringLiteral("mzVals");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");
    const QString TARGET_KEY = QStringLiteral("targetKey");

    const QStringList keysToCheck = {
        MS_LEVEL,
        SCAN_NUMBER ,
        SCAN_TIME,
        COLLISION_ENERGY,
        PERCURSOR_TARGET_MZ,
        ISO_WINDOW_LOWER ,
        ISO_WINDOW_UPPER ,
        IM_DRIFT_TIME ,
        IM_IND ,
        MZ_VALS ,
        INTENSITY_VALS,
        TARGET_KEY
    };
}

struct FILEREADERSLIB_EXPORTS MsParquetReaderRow : public ParquetReaderInputBase {

    int msLevel = -1;
    ScanNumber scanNumber = 1;
    float scanTime = -1.0;
    float collisionEnergy = -1.0;
    float precursorTargetMz = -1.0;
    float isoWindowLower = -1.0;
    float isoWindowUpper = -1.0;
    float ionMobilityDriftTime = -1.0;
    IonMobilityIndex ionMobilityIndex = -1;
    QVector<float> mzVals;
    QVector<float> intensityVals;
    QString targetKey;

    QMap<QString, QVariant> map() override {

        using namespace MsParquetReaderNamespace;

        return {
            {MS_LEVEL, QVariant(msLevel)},
            {SCAN_NUMBER, QVariant(scanNumber)},
            {SCAN_TIME, QVariant(scanTime)},
            {COLLISION_ENERGY, QVariant(collisionEnergy)},
            {PERCURSOR_TARGET_MZ, QVariant(precursorTargetMz)},
            {ISO_WINDOW_LOWER, QVariant(isoWindowLower)},
            {ISO_WINDOW_UPPER, QVariant(isoWindowUpper)},
            {IM_DRIFT_TIME, QVariant(ionMobilityDriftTime)},
            {IM_IND, QVariant(ionMobilityIndex)},
            {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
            {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))},
            {TARGET_KEY, QVariant(QVariant(targetKey))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace MsParquetReaderNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        msLevel = dataMap.value(MS_LEVEL).toInt();
        scanNumber = dataMap.value(SCAN_NUMBER).toInt();
        scanTime = dataMap.value(SCAN_TIME).toFloat();
        collisionEnergy = dataMap.value(COLLISION_ENERGY).toFloat();
        precursorTargetMz = dataMap.value(PERCURSOR_TARGET_MZ).toFloat();
        isoWindowLower = dataMap.value(ISO_WINDOW_LOWER).toFloat();
        isoWindowUpper = dataMap.value(ISO_WINDOW_UPPER).toFloat();
        ionMobilityDriftTime = dataMap.value(IM_DRIFT_TIME).toFloat();
        ionMobilityIndex = dataMap.value(IM_IND).toInt();
        mzVals = bytesArrayToQVector<float>(dataMap.value(MZ_VALS).toByteArray());
        intensityVals = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS).toByteArray());
        targetKey = dataMap.value(TARGET_KEY).toString();

        ERR_RETURN
    }

};


class MsReaderPointerAcc;

class FILEREADERSLIB_EXPORTS MsReaderParquet : public MsReaderBase {

public:

    MsReaderParquet() = default;
    ~MsReaderParquet() = default;

    Err openFile(const QString &filePath) override;

    Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy,
            const QPair<double, double> &filterRange
            ) override;

    Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy
    ) override;

    Err closeFile() override;

    static Err writeMsReaderToParquet(
            const QString &outputFilePath,
            const QSharedPointer<MsReaderBase> &sharedMsReaderBase
            );

};

#endif //SPARKDIA_PARQUETREADER_H
