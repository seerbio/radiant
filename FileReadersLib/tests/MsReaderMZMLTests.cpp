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

    void openFileTest();

private:

    //TODO use proper path procedures after finding small file.
    const QString m_filepath
            = QStringLiteral("/home/anichols/Downloads/EXP22092_2022ms0742X32_A.raw.mzML");

};


void MsReaderMZMLTests::openFileTest() {
    QSKIP("TODO: enable with internal test data");

    QSKIP("activate when proper pathing is used");

    QSKIP("activate when proper pathing is used");

//    const QString &msParquetFilePath
//            = QDir(qApp->applicationDirPath()).filePath("SoLetItBeWritten.prq");

    ERR_INIT

    QVERIFY2(QFileInfo::exists(m_filepath), qPrintable(QString("File not found: %1").arg(m_filepath)));

    MsReaderMzML reader;
    e = reader.openFile(m_filepath);
    QCOMPARE(e, Error::eNoError);
    QCOMPARE(reader.m_msScanInfo.size(), 26010);
    QCOMPARE(reader.m_scanPoints.size(), 26010);

    MsScanInfo msScanInfo;
    e = reader.getMsScanInfo(666, &msScanInfo);
    QCOMPARE(e, eNoError);
    QCOMPARE(msScanInfo.scanNumber, 666);
    QCOMPARE(msScanInfo.msLevel, 2);
    QCOMPARE(msScanInfo.collisionEnergy, 28);
    QCOMPARE(QString::number(msScanInfo.scanTime), "0.826693");
    QCOMPARE(QString::number(msScanInfo.precursorTargetMz), "725.079");
    QCOMPARE(msScanInfo.isoWindowLower, 5.5);
    QCOMPARE(msScanInfo.isoWindowUpper, 5.5);

    QPair<Err, ScanPoints*> scanPointsResult = reader.getScanPoints(msScanInfo.scanNumber);
    e = scanPointsResult.first;
    QCOMPARE(e, eNoError);

    const ScanPoints scanPoints = *scanPointsResult.second;

    QCOMPARE(int(149.045), int(scanPoints.first().x()));
    QCOMPARE(int(3574.26), int(scanPoints.first().y()));
    QCOMPARE(int(1166.42), int(scanPoints.last().x()));
    QCOMPARE(int(1360.38), int(scanPoints.last().y()));
    qDebug() << scanPoints.last();
    QCOMPARE(scanPoints.size(), 35);

}

QTEST_MAIN(MsReaderMZMLTests)
#include "MsReaderMZMLTests.moc"
