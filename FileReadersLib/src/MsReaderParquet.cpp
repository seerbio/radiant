//
// Created by anichols on 2/6/23.
//

#include "MsReaderParquet.h"

#include "ErrorUtils.h"


Err MsReaderParquet::openFile(const QString &filePath) {

    ERR_INIT

    ParquetReader reader;

    e = ErrorUtils::isNotEmpty(filePath); ree;
    m_filePath = filePath;

    QVector<ParquetReaderInputBase> ptrsRead;
    e = reader.readDataFromParquet(
            filePath,
            &ptrsRead
    ); ree;

    QVector<MsParquetReaderRow> msParquetReaderRows;
    e = ParquetReaderInputBase::convertSharedPointersToInputStruct(
            ptrsRead,
            &msParquetReaderRows
    ); ree;

    for (const MsParquetReaderRow &row : msParquetReaderRows) {

        MsScanInfo msScanInfo;

        msScanInfo.msLevel = row.msLevel;
        msScanInfo.scanNumber = row.scanNumber;
        msScanInfo.scanTime = row.scanTime;
        msScanInfo.collisionEnergy = row.collisionEnergy;
        msScanInfo.precursorTargetMz = row.precursorTargetMz;
        msScanInfo.isoWindowLower = row.isoWindowLower;
        msScanInfo.isoWindowUpper = row.isoWindowUpper;
        msScanInfo.ionMobilityDriftTime = row.ionMobilityDriftTime;
        msScanInfo.ionMobilityIndex = row.ionMobilityIndex;

        ScanPoints scanPoints;
        e = MsReaderBase::zipScanPoints(
                row.mzVals,
                row.intensityVals,
                &scanPoints
                ); ree;

        e = ErrorUtils::isFalse(m_msScanInfo.contains(msScanInfo.scanNumber)); ree;
        e = ErrorUtils::isFalse(m_scanPoints.contains(msScanInfo.scanNumber)); ree;

        m_msScanInfo.insert(msScanInfo.scanNumber, msScanInfo);
        m_scanPoints.insert(msScanInfo.scanNumber, scanPoints);
    }

    e = printFileInfo(); ree;

    ERR_RETURN
}

Err MsReaderParquet::closeFile() {

    ERR_INIT

    m_msScanInfo.clear();
    m_scanPoints.clear();

    e = ErrorUtils::isTrue(m_msScanInfo.isEmpty()); ree;
    e = ErrorUtils::isTrue(m_scanPoints.isEmpty()); ree;

    ERR_RETURN
}

namespace {

    Err buildRowsToWrite(
            const QSharedPointer<MsReaderBase> &sharedMsReaderBase,
            QVector<QSharedPointer<ParquetReaderInputBase>> *ptrs
            ) {

        ERR_INIT

        const QMap<ScanNumber, MsScanInfo> &msScanInfos = sharedMsReaderBase->getMsScanInfos();
        const QMap<ScanNumber, ScanPoints> &scanPoints = sharedMsReaderBase->getScanPoints();

        QVector<MsParquetReaderRow> rowsToWrite;
        for (auto it = msScanInfos.begin(); it != msScanInfos.end(); it++) {

            const ScanNumber &scanNumber = it.key();
            const MsScanInfo &msScanInfo = it.value();

            e = ErrorUtils::isTrue(scanPoints.contains(scanNumber)); ree;

            const ScanPoints &scanPointsVec = scanPoints.value(scanNumber);

            QVector<double> mzVals;
            QVector<double> intensityVals;
            e = MsReaderBase::splitScanPoints(
                    scanPointsVec,
                    &mzVals,
                    &intensityVals
            ); ree;

            MsParquetReaderRow row;
            row.msLevel = msScanInfo.msLevel;
            row.scanNumber = scanNumber;
            row.scanTime = msScanInfo.scanTime;
            row.collisionEnergy = msScanInfo.collisionEnergy;
            row.precursorTargetMz = msScanInfo.precursorTargetMz;
            row.isoWindowLower = msScanInfo.isoWindowLower;
            row.isoWindowUpper = msScanInfo.isoWindowUpper;
            row.ionMobilityDriftTime = msScanInfo.ionMobilityDriftTime;
            row.ionMobilityIndex = msScanInfo.ionMobilityIndex;
            row.mzVals = mzVals;
            row.intensityVals = intensityVals;

            rowsToWrite.push_back(row);
        }

        *ptrs = ParquetReaderInputBase::convertInputStructToSharedPointers(rowsToWrite);

        ERR_RETURN
    }

}//namespace
Err MsReaderParquet::writeMsReaderToParquet(
        const QString &outputFilePath,
        const QSharedPointer<MsReaderBase> &sharedMsReaderBase
        ) {

    ERR_INIT

    QVector<QSharedPointer<ParquetReaderInputBase>> ptrs;
    e = buildRowsToWrite(sharedMsReaderBase, &ptrs); ree;

    qDebug() << "SDFS" << ptrs.size();
    e = ErrorUtils::isNotEmpty(ptrs); ree;

    ParquetReader reader;
    e = reader.writeDataToParquet(
            outputFilePath,
            ptrs
    ); ree;

    ERR_RETURN
}
