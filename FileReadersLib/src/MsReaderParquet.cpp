//
// Created by anichols on 2/6/23.
//

#include "MsReaderParquet.h"

#include "ErrorUtils.h"
#include "MsReaderPointerAcc.h"


namespace {

    Err convertForMemberVars(
            const QVector<MsParquetReaderRow> &msParquetReaderRows,
            QMap<ScanNumber, MsScanInfo> *memberMsScanInfo,
            QMap<ScanNumber, ScanPoints>  *memberScanPoints
            ) {

        ERR_INIT

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

            e = ErrorUtils::isFalse(memberMsScanInfo->contains(msScanInfo.scanNumber)); ree;
            e = ErrorUtils::isFalse(memberScanPoints->contains(msScanInfo.scanNumber)); ree;

            memberMsScanInfo->insert(msScanInfo.scanNumber, msScanInfo);
            memberScanPoints->insert(msScanInfo.scanNumber, scanPoints);
        }

        ERR_RETURN
    }

}
Err MsReaderParquet::openFile(const QString &filePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(filePath); ree;
    m_filePath = filePath;

    QVector<MsParquetReaderRow> msParquetReaderRows;
    e = ParquetReader::read(
            filePath,
            &msParquetReaderRows
    ); ree;

    e = convertForMemberVars(
            msParquetReaderRows,
            &m_msScanInfo,
            &m_scanPoints
            ); ree;

    e = printFileInfo(); ree;

    ERR_RETURN
}

Err MsReaderParquet::openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
) {

    ERR_INIT

    e = ErrorUtils::fileExists(filePath); ree;
    m_filePath = filePath;

    QVector<MsParquetReaderRow> msParquetReaderRows;
    e = ParquetReader::read(
            filePath,
            columnToFilterBy,
            filterRange,
            &msParquetReaderRows
    ); ree;

    e = convertForMemberVars(
            msParquetReaderRows,
            &m_msScanInfo,
            &m_scanPoints
    ); ree;

    ERR_RETURN
}

Err MsReaderParquet::openFile(
        const QString &filePath,
        const QString &columnToFilterBy
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(filePath); ree;
    m_filePath = filePath;

    QVector<MsParquetReaderRow> msParquetReaderRows;
    e = ParquetReader::read(
            filePath,
            columnToFilterBy,
            &msParquetReaderRows
    ); ree;

    e = convertForMemberVars(
            msParquetReaderRows,
            &m_msScanInfo,
            &m_scanPoints
    ); ree;

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
            const QMap<ScanNumber, MsScanInfo> &msScanInfos,
            const QMap<ScanNumber, ScanPoints> &scanPoints,
            QVector<QSharedPointer<ParquetReaderInputBase>> *ptrs
            ) {

        ERR_INIT

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
            row.targetKey = QString::number(static_cast<int>(std::round(row.precursorTargetMz *1000)));

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
    e = buildRowsToWrite(
            sharedMsReaderBase->getMsScanInfos(),
            sharedMsReaderBase->getScanPoints(),
            &ptrs
            ); ree;

    e = ErrorUtils::isNotEmpty(ptrs); ree;

    ParquetReader reader;
    e = reader.writeDataToParquet(
            outputFilePath,
            ptrs
    ); ree;

    ERR_RETURN
}
