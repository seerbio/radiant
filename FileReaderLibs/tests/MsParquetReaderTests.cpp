//
// Created by anichols on 11/07/2021.
//

#include "MsParquetReader.h"

#include <QtTest/QtTest>

class MsParquetReaderTests : public QObject
{
    Q_OBJECT

public:
    MsParquetReaderTests() = default;
    ~MsParquetReaderTests() override = default;

private Q_SLOTS:

    void checkParquetStatusTest();
    void readFileTest();
    void getScansTest();
    void getScansInfoTest();
    void writeFileTest();

};

void MsParquetReaderTests::checkParquetStatusTest() {

    const bool isOk = MsParquetReader::checkParquetStatus();
    QCOMPARE(isOk, true);
}

void MsParquetReaderTests::readFileTest() {

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.parquet");

    MsParquetReader reader;

    bool e;
    e = reader.readFile(msParquetFilePath.toStdString());
    QCOMPARE(e, true);

}

void MsParquetReaderTests::getScansTest() {

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.parquet");

    MsParquetReader reader;

    bool e;
    e = reader.readFile(msParquetFilePath.toStdString());
    QCOMPARE(e, true);

    QMap<ScanNumber, ScanPoints> scanPoints;
    e = reader.getScans(2, &scanPoints);
    QCOMPARE(e, true);

    QCOMPARE(scanPoints.firstKey(), 2);
    QCOMPARE(scanPoints.first().size(), 51);

    const QPointF &qp = scanPoints.value(2).at(0);
    QCOMPARE(std::round(qp.x()), std::round(147.064));
    QCOMPARE(std::round(qp.y()), std::round(2569));
}

void MsParquetReaderTests::getScansInfoTest() {

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.parquet");

    MsParquetReader reader;

    bool e;
    e = reader.readFile(msParquetFilePath.toStdString());
    QCOMPARE(e, true);

    QMap<ScanNumber, ScanInfo> scansInfoMap;
    e = reader.getScansInfo(&scansInfoMap);
    QCOMPARE(e, true);
    QCOMPARE(scansInfoMap.size(), 2);

    const ScanInfo &scanInfo = scansInfoMap.value(2);
    QCOMPARE(scanInfo.msLevel, 2);
    QCOMPARE(scanInfo.scanNumber, 2);
    QCOMPARE(scanInfo.retentionTime, 0.006662306341);
    QCOMPARE(scanInfo.collisionEnergy, 28);
    QCOMPARE(scanInfo.mzTarget, 385.213684082);
    QCOMPARE(scanInfo.isoWindowLower, 5.5);
    QCOMPARE(scanInfo.isoWindowUpper, 5.5);

}

void MsParquetReaderTests::writeFileTest() {

    const QString &msParquetFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.parquet");

    MsParquetReader reader;

    bool e;
    e = reader.readFile(msParquetFilePath.toStdString());
    QCOMPARE(e, true);

    const QString writeFileName = QStringLiteral("SoLetItBeWritten.parquet");
    e = reader.writeFile(writeFileName.toStdString());
    QCOMPARE(e, true);

    QFileInfo checkFile(writeFileName);
    QCOMPARE(checkFile.exists(), true);

    QFile::remove(writeFileName);
    QFileInfo checkFile2(writeFileName);
    QCOMPARE(checkFile2.exists(), false);
}


QTEST_MAIN(MsParquetReaderTests)
#include "MsParquetReaderTests.moc"
