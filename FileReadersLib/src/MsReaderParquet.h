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

    const QString PRQ_FF_SUFFIX = QStringLiteral("prqFF");

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

    /**
    * @brief Opens a Parquet file using MsReaderParquet.
    *
    * This function opens a Parquet file specified by the `filePath` parameter with the
    * MsReaderParquet instance. It checks if the file exists, validates the file extension,
    * reads the data from the Parquet file into MsParquetReaderRow objects, converts the
    * data for member variables, and prints file information. Any encountered errors during
    * the process are indicated by the returned Err code.
    *
    * @param filePath The file path of the Parquet file to be opened.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err openFile(const QString &filePath) override;

    /**
    * @brief Opens a Parquet file using MsReaderParquet with optional column filtering.
    *
    * This function opens a Parquet file specified by the `filePath` parameter with the
    * MsReaderParquet instance. It checks if the file exists, validates the file extension,
    * reads the data from the Parquet file into MsParquetReaderRow objects, filters the data
    * based on the specified `columnToFilterBy` and `filterRange`, converts the data for
    * member variables, and prints file information. Any encountered errors during the process
    * are indicated by the returned Err code.
    *
    * @param filePath The file path of the Parquet file to be opened.
    * @param columnToFilterBy The name of the column to filter data by (optional).
    * @param filterRange The numeric range within which data is filtered (optional).
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy,
            const QPair<double, double> &filterRange
            ) override;

    /**
    * @brief Opens a Parquet file using MsReaderParquet with optional column filtering.
    *
    * This function opens a Parquet file specified by the `filePath` parameter with the
    * MsReaderParquet instance. It checks if the file exists, validates the file extension,
    * reads the data from the Parquet file into MsParquetReaderRow objects, filters the data
    * based on the specified `columnToFilterBy`, converts the data for member variables, and
    * prints file information. Any encountered errors during the process are indicated by
    * the returned Err code.
    *
    * @param filePath The file path of the Parquet file to be opened.
    * @param columnToFilterBy The name of the column to filter data by (optional).
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy
    ) override;

    /**
    * @brief Closes the Parquet file associated with MsReaderParquet.
    *
    * This function clears the member variables `m_msScanInfo` and `m_scanPoints`,
    * ensuring that they are empty after the file is closed. It then checks if both
    * member variables are empty to confirm a successful closure. Any encountered errors
    * during the process are indicated by the returned Err code.
    *
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err closeFile() override;

    /**
    * @brief Writes MsReader data to a Parquet file.
    *
    * This function writes data from a specified MsReader instance (`sharedMsReaderBase`)
    * to a Parquet file specified by the `outputFilePath`. It builds the rows to be written
    * using the MsReader data, checks if the rows are not empty, and then writes the data
    * to the Parquet file using ParquetReader. Any encountered errors during the process
    * are indicated by the returned Err code.
    *
    * @param outputFilePath The file path where the Parquet file will be written.
    * @param sharedMsReaderBase A shared pointer to the MsReader instance containing the data to be written.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    static Err writeMsReaderToParquet(
            const QString &outputFilePath,
            const QSharedPointer<MsReaderBase> &sharedMsReaderBase
            );

};

#endif //SPARKDIA_PARQUETREADER_H
