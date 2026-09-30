//
// Created by anichols on 2/6/23.
//

#include "MsReaderParquet.h"

#include "ErrorUtils.h"
#include "MsReaderPointerAcc.h"

namespace {

    Err prepareForOpen(
            const QString &filePath,
            QString *normalizedFilePath
            ) {

        ERR_INIT

        e = ErrorUtils::fileExists(filePath); ree;

        QFileInfo fi(filePath);
        const QString fileSuffix = fi.suffix();

        e = ErrorUtils::isTrue(
                MsParquetReaderNamespace::PRQ_FF_SUFFIX == fileSuffix,
                eFileIncorrectTypeError
        ); ree;

        *normalizedFilePath = filePath;

        ERR_RETURN
    }

    Err convertForMemberVars(
            const QVector<MsParquetReaderRow> &msParquetReaderRows,
            QMap<ScanNumber, MsScanInfo> *memberMsScanInfo,
            QMap<ScanNumber, ScanPoints>  *memberScanPoints,
            QMap<FrameIndex, double> *frameIndexVsDriftTime,
            bool *hasIonMobility
            ) {

        ERR_INIT

        for (const MsParquetReaderRow &row : msParquetReaderRows) {

            if (row.rowType != MsParquetReaderNamespace::ROW_TYPE_SPECTRUM) {
                continue;
            }

            MsScanInfo msScanInfo;

            msScanInfo.msLevel = row.msLevel;
            msScanInfo.scanNumber = row.scanNumber;
            msScanInfo.scanTime = row.scanTime;
            msScanInfo.collisionEnergy = row.collisionEnergy;
            msScanInfo.precursorTargetMz = row.precursorTargetMz;
            msScanInfo.isoWindowLower = row.isoWindowLower;
            msScanInfo.isoWindowUpper = row.isoWindowUpper;
            msScanInfo.ionMobilityDriftTime = row.ionMobilityDriftTime;
            msScanInfo.ionMobilityWindowLower = row.ionMobilityWindowLower;
            msScanInfo.ionMobilityWindowUpper = row.ionMobilityWindowUpper;
            msScanInfo.ionMobilityIndex = row.ionMobilityIndex;
            msScanInfo.nativeFrameNumber = row.nativeFrameNumber;
            msScanInfo.nativeScanNumber = row.nativeScanNumber;

            ScanPoints scanPoints;
            e = MsReaderBase::zipScanPoints(
                    row.mzVals,
                    row.intensityVals,
                    &scanPoints
            ); ree;

            e = ErrorUtils::doesNotContain(msScanInfo.scanNumber, *memberMsScanInfo); ree;
            e = ErrorUtils::doesNotContain(msScanInfo.scanNumber, *memberScanPoints); ree;

            memberMsScanInfo->insert(msScanInfo.scanNumber, msScanInfo);
            memberScanPoints->insert(msScanInfo.scanNumber, scanPoints);
            if (row.ionMobilityIndex >= 0 && row.ionMobilityDriftTime >= 0) {
                frameIndexVsDriftTime->insert(row.ionMobilityIndex, row.ionMobilityDriftTime);
            }
            if (row.ionMobilityIndex >= 0
                || row.ionMobilityDriftTime >= 0.0f
                || row.ionMobilityWindowLower >= 0.0f
                || row.ionMobilityWindowUpper >= 0.0f) {
                *hasIonMobility = true;
            }
        }

        ERR_RETURN
    }

}//namespace

Err MsReaderParquet::openFile(const QString &filePath) {

    ERR_INIT

    e = prepareForOpen(filePath, &m_filePath); ree;
    e = closeFile(); ree;

    QVector<MsParquetReaderRow> msParquetReaderRows;
    e = ParquetReader::read(
            filePath,
            &msParquetReaderRows
    ); ree;

    e = convertForMemberVars(
            msParquetReaderRows,
            &m_msScanInfo,
            &m_scanPoints,
            &m_frameIndexVsDriftTime,
            &m_hasIonMobility
            ); ree;

    e = printFileInfo(); ree;

    QVector<MsParquetReaderRow>().swap(msParquetReaderRows);

    ERR_RETURN
}

Err MsReaderParquet::openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
) {

    ERR_INIT

    e = prepareForOpen(filePath, &m_filePath); ree;
    e = closeFile(); ree;

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
            &m_scanPoints,
            &m_frameIndexVsDriftTime,
            &m_hasIonMobility
    ); ree;

    QVector<MsParquetReaderRow>().swap(msParquetReaderRows);

    ERR_RETURN
}

Err MsReaderParquet::openFile(
        const QString &filePath,
        const QString &columnToFilterBy
        ) {

    ERR_INIT

    e = prepareForOpen(filePath, &m_filePath); ree;
    e = closeFile(); ree;

    QVector<MsParquetReaderRow> msParquetReaderRows;
    e = ParquetReader::read(
            filePath,
            columnToFilterBy,
            &msParquetReaderRows
    ); ree;

    e = convertForMemberVars(
            msParquetReaderRows,
            &m_msScanInfo,
            &m_scanPoints,
            &m_frameIndexVsDriftTime,
            &m_hasIonMobility
    ); ree;

    QVector<MsParquetReaderRow>().swap(msParquetReaderRows);

    ERR_RETURN
}

Err MsReaderParquet::closeFile() {

    ERR_INIT

    m_msScanInfo.clear();
    m_scanPoints.clear();
    m_frameIndexVsDriftTime.clear();
    m_hasIonMobility = false;

    e = ErrorUtils::isTrue(m_msScanInfo.isEmpty()); ree;
    e = ErrorUtils::isTrue(m_scanPoints.isEmpty()); ree;

    ERR_RETURN
}

namespace {

    Err buildRowsToWrite(
            const QMap<ScanNumber, MsScanInfo> &msScanInfos,
            const QMap<ScanNumber, ScanPoints*> &scanPoints,
            QVector<QSharedPointer<ParquetReaderInputBase>> *ptrs
            ) {

        ERR_INIT

        QVector<MsParquetReaderRow> rowsToWrite;
        for (auto it = msScanInfos.begin(); it != msScanInfos.end(); it++) {

            const ScanNumber &scanNumber = it.key();
            const MsScanInfo &msScanInfo = it.value();

            e = ErrorUtils::contains(scanNumber, scanPoints); ree;

            const ScanPoints *scanPointsVec = scanPoints.value(scanNumber);

            QVector<float> mzVals;
            QVector<float> intensityVals;
            e = MsReaderBase::splitScanPoints(
                    *scanPointsVec,
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
            row.ionMobilityWindowLower = msScanInfo.ionMobilityWindowLower;
            row.ionMobilityWindowUpper = msScanInfo.ionMobilityWindowUpper;
            row.ionMobilityIndex = msScanInfo.ionMobilityIndex;
            row.nativeFrameNumber = msScanInfo.nativeFrameNumber;
            row.nativeScanNumber = msScanInfo.nativeScanNumber;
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
            sharedMsReaderBase->getScanPointsPntrs(),
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
