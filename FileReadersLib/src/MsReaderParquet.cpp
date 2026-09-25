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
            QMap<ScanNumber, ScanPoints>  *memberScanPoints,
            QMap<FrameNumberTIMS, Ms1FrameTIMS> *frameNumberVsMS1FrameTIMS,
            MzTargetKeyVsMs2FrameTIMS *mzTargetKeyVsFrameNumberVsMS2FrameTIMS,
            QMap<FrameIndex, double> *frameIndexVsDriftTime,
            bool *isTIMS
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
        }

        for (const MsParquetReaderRow &row : msParquetReaderRows) {

            if (row.rowType != MsParquetReaderNamespace::ROW_TYPE_TIMS_MS1_SCAN) {
                continue;
            }

            ScanPoints scanPoints;
            e = MsReaderBase::zipScanPoints(
                    row.mzVals,
                    row.intensityVals,
                    &scanPoints
                    ); ree;

            (*frameNumberVsMS1FrameTIMS)[row.scanNumber].insert(row.ionMobilityIndex, scanPoints);

            if (row.ionMobilityIndex >= 0 && row.ionMobilityDriftTime >= 0) {
                frameIndexVsDriftTime->insert(row.ionMobilityIndex, row.ionMobilityDriftTime);
            }
        }

        for (const MsParquetReaderRow &row : msParquetReaderRows) {

            if (row.rowType != MsParquetReaderNamespace::ROW_TYPE_TIMS_MS2_SCAN) {
                continue;
            }

            ScanPoints scanPoints;
            e = MsReaderBase::zipScanPoints(
                    row.mzVals,
                    row.intensityVals,
                    &scanPoints
                    ); ree;

            (*mzTargetKeyVsFrameNumberVsMS2FrameTIMS)[row.targetKey][row.scanNumber].insert(row.ionMobilityIndex, scanPoints);

            if (row.ionMobilityIndex >= 0 && row.ionMobilityDriftTime >= 0) {
                frameIndexVsDriftTime->insert(row.ionMobilityIndex, row.ionMobilityDriftTime);
            }
        }

        *isTIMS = !frameNumberVsMS1FrameTIMS->isEmpty() || !mzTargetKeyVsFrameNumberVsMS2FrameTIMS->isEmpty();

        ERR_RETURN
    }

}//namespace

Err MsReaderParquet::openFile(const QString &filePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(filePath); ree;
    m_filePath = filePath;

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    e = ErrorUtils::isTrue(
            MsParquetReaderNamespace::PRQ_FF_SUFFIX == fileSuffix,
            eFileIncorrectTypeError
    ); ree;

    QVector<MsParquetReaderRow> msParquetReaderRows;
    e = ParquetReader::read(
            filePath,
            &msParquetReaderRows
    ); ree;

    e = convertForMemberVars(
            msParquetReaderRows,
            &m_msScanInfo,
            &m_scanPoints,
            &m_frameNumberVsMS1FrameTIMS,
            &m_mzTargetKeyVsFrameNumberVsMS2FrameTIMS,
            &m_frameIndexVsDriftTime,
            &m_isTIMS
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

    e = ErrorUtils::fileExists(filePath); ree;
    m_filePath = filePath;

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    e = ErrorUtils::isTrue(
            MsParquetReaderNamespace::PRQ_FF_SUFFIX == fileSuffix,
            eFileIncorrectTypeError
    ); ree;

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
            &m_frameNumberVsMS1FrameTIMS,
            &m_mzTargetKeyVsFrameNumberVsMS2FrameTIMS,
            &m_frameIndexVsDriftTime,
            &m_isTIMS
    ); ree;

    QVector<MsParquetReaderRow>().swap(msParquetReaderRows);

    ERR_RETURN
}

Err MsReaderParquet::openFile(
        const QString &filePath,
        const QString &columnToFilterBy
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(filePath); ree;
    m_filePath = filePath;

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    e = ErrorUtils::isTrue(
            MsParquetReaderNamespace::PRQ_FF_SUFFIX == fileSuffix,
            eFileIncorrectTypeError
    ); ree;

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
            &m_frameNumberVsMS1FrameTIMS,
            &m_mzTargetKeyVsFrameNumberVsMS2FrameTIMS,
            &m_frameIndexVsDriftTime,
            &m_isTIMS
    ); ree;

    QVector<MsParquetReaderRow>().swap(msParquetReaderRows);

    ERR_RETURN
}

Err MsReaderParquet::closeFile() {

    ERR_INIT

    m_msScanInfo.clear();
    m_scanPoints.clear();
    m_frameNumberVsMS1FrameTIMS.clear();
    m_mzTargetKeyVsFrameNumberVsMS2FrameTIMS.clear();
    m_frameIndexVsDriftTime.clear();
    m_isTIMS = false;

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

    if (sharedMsReaderBase->isTIMS()) {
        QMap<FrameNumberTIMS, Ms1FrameTIMS> *ms1Frames = sharedMsReaderBase->frameNumberVsMS1FrameTIMSPntr();

        for (auto frameIt = ms1Frames->begin(); frameIt != ms1Frames->end(); ++frameIt) {

            MsScanInfo msScanInfo;
            e = sharedMsReaderBase->getMsScanInfo(frameIt.key(), &msScanInfo); ree;

            const Ms1FrameTIMS &ms1FrameTims = frameIt.value();
            for (auto mobilityIt = ms1FrameTims.begin(); mobilityIt != ms1FrameTims.end(); ++mobilityIt) {

                QVector<float> mzVals;
                QVector<float> intensityVals;
                e = MsReaderBase::splitScanPoints(
                        mobilityIt.value(),
                        &mzVals,
                        &intensityVals
                        ); ree;

                double driftTime = -1.0;
                e = sharedMsReaderBase->driftTimeFromIonMobilityIndex(
                        mobilityIt.key(),
                        &driftTime
                        ); eee_absorb;

                MsParquetReaderRow row;
                row.rowType = MsParquetReaderNamespace::ROW_TYPE_TIMS_MS1_SCAN;
                row.msLevel = 1;
                row.scanNumber = frameIt.key();
                row.scanTime = msScanInfo.scanTime;
                row.ionMobilityIndex = mobilityIt.key();
                row.ionMobilityDriftTime = static_cast<float>(driftTime);
                row.nativeFrameNumber = msScanInfo.nativeFrameNumber;
                row.nativeScanNumber = mobilityIt.key();
                row.mzVals = mzVals;
                row.intensityVals = intensityVals;
                row.targetKey = S_GLOBAL_SETTINGS.MS1Key;

                ptrs.push_back(QSharedPointer<ParquetReaderInputBase>(new MsParquetReaderRow(row)));
            }
        }

        MzTargetKeyVsMs2FrameTIMS *ms2Frames = sharedMsReaderBase->mzTargetKeyVsFrameNumberVsMS2FrameTIMSPntr();

        for (auto targetIt = ms2Frames->begin(); targetIt != ms2Frames->end(); ++targetIt) {

            const MzTargetKey &targetKey = targetIt.key();
            const QMap<FrameNumberTIMS, Ms2FrameTIMS> &frameNumberVsMs2FrameTims = targetIt.value();

            for (auto frameIt = frameNumberVsMs2FrameTims.begin(); frameIt != frameNumberVsMs2FrameTims.end(); ++frameIt) {

                MsScanInfo msScanInfo;
                e = sharedMsReaderBase->getMsScanInfo(frameIt.key(), &msScanInfo); ree;

                const Ms2FrameTIMS &ms2FrameTims = frameIt.value();
                for (auto mobilityIt = ms2FrameTims.begin(); mobilityIt != ms2FrameTims.end(); ++mobilityIt) {

                    QVector<float> mzVals;
                    QVector<float> intensityVals;
                    e = MsReaderBase::splitScanPoints(
                            mobilityIt.value(),
                            &mzVals,
                            &intensityVals
                            ); ree;

                    double driftTime = -1.0;
                    e = sharedMsReaderBase->driftTimeFromIonMobilityIndex(
                            mobilityIt.key(),
                            &driftTime
                            ); eee_absorb;

                    MsParquetReaderRow row;
                    row.rowType = MsParquetReaderNamespace::ROW_TYPE_TIMS_MS2_SCAN;
                    row.msLevel = 2;
                    row.scanNumber = frameIt.key();
                    row.scanTime = msScanInfo.scanTime;
                    row.collisionEnergy = msScanInfo.collisionEnergy;
                    row.precursorTargetMz = msScanInfo.precursorTargetMz;
                    row.isoWindowLower = msScanInfo.isoWindowLower;
                    row.isoWindowUpper = msScanInfo.isoWindowUpper;
                    row.ionMobilityIndex = mobilityIt.key();
                    row.ionMobilityDriftTime = static_cast<float>(driftTime);
                    row.ionMobilityWindowLower = msScanInfo.ionMobilityWindowLower;
                    row.ionMobilityWindowUpper = msScanInfo.ionMobilityWindowUpper;
                    row.nativeFrameNumber = msScanInfo.nativeFrameNumber;
                    row.nativeScanNumber = mobilityIt.key();
                    row.mzVals = mzVals;
                    row.intensityVals = intensityVals;
                    row.targetKey = targetKey;

                    ptrs.push_back(QSharedPointer<ParquetReaderInputBase>(new MsParquetReaderRow(row)));
                }
            }
        }
    }

    e = ErrorUtils::isNotEmpty(ptrs); ree;

    ParquetReader reader;
    e = reader.writeDataToParquet(
            outputFilePath,
            ptrs
    ); ree;

    ERR_RETURN
}
