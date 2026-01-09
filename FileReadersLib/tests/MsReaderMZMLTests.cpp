//
// Created by anichols on 11/07/2021.
//

#include "MsReaderMzML.h"
#include "FastaReader.h"

#include <QString>
#include <QtTest/QtTest>
#include <QXmlStreamReader>


class MsReaderMZMLTests : public QObject
{
    Q_OBJECT

public:
    MsReaderMZMLTests() = default;
    ~MsReaderMZMLTests() override = default;


//TODO Add more test coverage.
private Q_SLOTS:

    static void openFileTest();
    static void troubleShoot();

};


void MsReaderMZMLTests::openFileTest() {

    const QString &msParquetFilePath = QDir(qApp->applicationDirPath()).filePath("1min.mzML");

    ERR_INIT

    QVERIFY2(QFileInfo::exists(msParquetFilePath), qPrintable(QString("File not found: %1").arg(msParquetFilePath)));

    MsReaderMzML reader;
    e = reader.openFile(msParquetFilePath);
    QCOMPARE(e, Error::eNoError);
    QCOMPARE(reader.m_msScanInfo.size(), 372);
    QCOMPARE(reader.m_scanPoints.size(), 372);

    MsScanInfo msScanInfo;
    e = reader.getMsScanInfo(372, &msScanInfo);
    QCOMPARE(e, eNoError);
    QCOMPARE(msScanInfo.scanNumber, 372);
    QCOMPARE(msScanInfo.msLevel, 1);
    QCOMPARE(msScanInfo.collisionEnergy, -1);
    QVERIFY(std::abs(msScanInfo.scanTime - 60.0 * 10.041118 < 0.0001));  // parsed value should be converted to seconds; compare with tolerance
    QCOMPARE(QString::number(msScanInfo.precursorTargetMz), "-1");
    QCOMPARE(msScanInfo.isoWindowLower, -1);
    QCOMPARE(msScanInfo.isoWindowUpper, -1);

    QPair<Err, ScanPoints*> scanPointsResult = reader.getScanPoints(msScanInfo.scanNumber);
    e = scanPointsResult.first;
    QCOMPARE(e, eNoError);

    const ScanPoints scanPoints = *scanPointsResult.second;

    QCOMPARE(int(scanPoints.first().x()), 150);
    QCOMPARE(int(scanPoints.first().y()), 814);
    QCOMPARE(int(scanPoints.last().x()), 1250);
    QCOMPARE(int(scanPoints.last().y()), 2646);
    QCOMPARE(scanPoints.size(), 2092);

}

void MsReaderMZMLTests::troubleShoot() {

    QSKIP("This should be skipped in production");

    ERR_INIT

    const QString brukerMzMlFilePath = QStringLiteral("/home/anichols/Desktop/EXP23140_2023ms1194X42_A_BB6_1_884.d.mzML");
    MsReaderMzML msReaderMzMl;
    e = msReaderMzMl.openFile(brukerMzMlFilePath);
    QCOMPARE(e, eNoError);


}


QTEST_MAIN(MsReaderMZMLTests)
#include "MsReaderMZMLTests.moc"
