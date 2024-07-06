//
// Created by anichols on 11/07/2021.
//

#include "MsReaderParquet.h"

#include "MsReaderBase.h"

#include <QtTest/QtTest>

class MsReaderParquetTests : public QObject
{
    Q_OBJECT

public:
    MsReaderParquetTests() = default;
    ~MsReaderParquetTests() override = default;

private Q_SLOTS:

    static void saveMsReaderToParquetOpenFileCombinedTests();
    static void readFilteredTest();
    static void readUniqueTest();

};

void MsReaderParquetTests::saveMsReaderToParquetOpenFileCombinedTests() {

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("SoLetItBeWritten.prqFF");

    MsScanInfo msScanInfo;
    msScanInfo.msLevel = 2;
    msScanInfo.scanNumber = 666;
    msScanInfo.scanTime = 66.0;
    msScanInfo.collisionEnergy = 6.0;
    msScanInfo.precursorTargetMz = 666.0;
    msScanInfo.isoWindowLower = 6.0;
    msScanInfo.isoWindowUpper = 6.0;
    msScanInfo.ionMobilityDriftTime = 61.0;
    msScanInfo.ionMobilityIndex = 666;

    const ScanPoints scanPoints = {{6.6, 666.6}, {666.6, 666666.6}, {66.6, 6666666666.6}};

    MsReaderBase msReaderBase;
    msReaderBase.m_msScanInfo.insert(msScanInfo.scanNumber, msScanInfo);
    msReaderBase.m_scanPoints.insert(msScanInfo.scanNumber, scanPoints);

    ERR_INIT

    e = MsReaderParquet::writeMsReaderToParquet(
            msParquetFilePath,
            QSharedPointer<MsReaderBase>(new MsReaderBase(msReaderBase))
            );
    QCOMPARE(e, eNoError);

    QFileInfo checkFile(msParquetFilePath);
    QCOMPARE(checkFile.exists(), true);

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(msParquetFilePath);
    QCOMPARE(e, eNoError);

    MsScanInfo readMsScanInfo;
    e = msReaderParquet.getMsScanInfo(
            msScanInfo.scanNumber,
            &readMsScanInfo
            );
    QCOMPARE(e, eNoError);

    QCOMPARE(readMsScanInfo.scanNumber, msScanInfo.scanNumber);
    QCOMPARE(readMsScanInfo.msLevel, msScanInfo.msLevel);
    QCOMPARE(readMsScanInfo.collisionEnergy, msScanInfo.collisionEnergy);
    QCOMPARE(readMsScanInfo.scanTime, msScanInfo.scanTime);
    QCOMPARE(readMsScanInfo.precursorTargetMz, msScanInfo.precursorTargetMz);
    QCOMPARE(readMsScanInfo.isoWindowLower, msScanInfo.isoWindowLower);
    QCOMPARE(readMsScanInfo.isoWindowUpper, msScanInfo.isoWindowUpper);
    QCOMPARE(readMsScanInfo.ionMobilityDriftTime, msScanInfo.ionMobilityDriftTime);
    QCOMPARE(readMsScanInfo.ionMobilityIndex, msScanInfo.ionMobilityIndex);


    QPair<Err, ScanPoints*> scanPointsResult = msReaderParquet.getScanPoints(readMsScanInfo.scanNumber);
    e = scanPointsResult.first;
    QCOMPARE(e, eNoError);

    ScanPoints* scanPointsRead = scanPointsResult.second;
    QCOMPARE(*scanPointsRead, scanPoints);

    QFile::remove(msParquetFilePath);
    QFileInfo checkFile2(msParquetFilePath);
    QCOMPARE(checkFile2.exists(), false);
}

void MsReaderParquetTests::readFilteredTest() {

    ERR_INIT

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

    QVERIFY2(QFileInfo::exists(msParquetFilePath), qPrintable(QString("File not found: %1").arg(msParquetFilePath)));

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(
            msParquetFilePath,
            "targetKey",
            {404934, 404934}
            );

    QMap<ScanNumber, MsScanInfo> scanInfos = msReaderParquet.getMsScanInfos();
    qDebug() << "ScanInfos size filtered" << scanInfos.size();
    QCOMPARE(scanInfos.size(), 80);
}

void MsReaderParquetTests::readUniqueTest() {

    ERR_INIT

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

    QVERIFY2(QFileInfo::exists(msParquetFilePath), qPrintable(QString("File not found: %1").arg(msParquetFilePath)));

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(
            msParquetFilePath,
            "targetKey"
    );
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, MsScanInfo> scanInfos = msReaderParquet.getMsScanInfos();
    qDebug() << "ScanInfos unique size filtered" << scanInfos.size();
    QCOMPARE(scanInfos.size(), 63);
}


QTEST_MAIN(MsReaderParquetTests)
#include "MsReaderParquetTests.moc"
