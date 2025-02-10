//
// Created by anichols on 11/07/2021.
//

#include "MsReaderMzMLMapped.h"
#include "FastaReader.h"

#include <QString>
#include <QtTest/QtTest>
#include <QXmlStreamReader>


class MsReaderMZMLMappedTests : public QObject
{
    Q_OBJECT

public:
    MsReaderMZMLMappedTests() = default;
    ~MsReaderMZMLMappedTests() override = default;


//TODO Add more test coverage.
private Q_SLOTS:

    static void openFileTest();
    static void troubleShoot();

};


void MsReaderMZMLMappedTests::openFileTest() {

    QSKIP("Turn on after complete and making tests.");

    // const QString &msParquetFilePath = QDir(qApp->applicationDirPath()).filePath("1min.mzML");
    //
    // ERR_INIT
    //
    // QVERIFY2(QFileInfo::exists(msParquetFilePath), qPrintable(QString("File not found: %1").arg(msParquetFilePath)));
    //
    // MsReaderMzML reader;
    // e = reader.openFile(msParquetFilePath);
    // QCOMPARE(e, Error::eNoError);
    // QCOMPARE(reader.m_msScanInfo.size(), 372);
    // QCOMPARE(reader.m_scanPoints.size(), 372);
    //
    // MsScanInfo msScanInfo;
    // e = reader.getMsScanInfo(372, &msScanInfo);
    // QCOMPARE(e, eNoError);
    // QCOMPARE(msScanInfo.scanNumber, 372);
    // QCOMPARE(msScanInfo.msLevel, 1);
    // QCOMPARE(msScanInfo.collisionEnergy, -1);
    // QCOMPARE(QString::number(msScanInfo.scanTime), "10.0411");
    // QCOMPARE(QString::number(msScanInfo.precursorTargetMz), "-1");
    // QCOMPARE(msScanInfo.isoWindowLower, -1);
    // QCOMPARE(msScanInfo.isoWindowUpper, -1);
    //
    // QPair<Err, ScanPoints*> scanPointsResult = reader.getScanPoints(msScanInfo.scanNumber);
    // e = scanPointsResult.first;
    // QCOMPARE(e, eNoError);
    //
    // const ScanPoints scanPoints = *scanPointsResult.second;
    //
    // QCOMPARE(int(scanPoints.first().x()), 150);
    // QCOMPARE(int(scanPoints.first().y()), 814);
    // QCOMPARE(int(scanPoints.last().x()), 1250);
    // QCOMPARE(int(scanPoints.last().y()), 2646);
    // QCOMPARE(scanPoints.size(), 2092);

}

//TODO make more tests.
void MsReaderMZMLMappedTests::troubleShoot() {

    ERR_INIT

    QSKIP("Turn on after complete and making tests.");

    // const QString filename  = "/home/anichols/Desktop/Data/MsData/EXP23111_2023ms0979bX45_A.raw.mzML";
    // const QString filename = "/home/anichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML";
    // const QString filename = "/home/anichols/Desktop/Data/MsData/EXP23140_2023ms1194X42_A_BB6_1_884.d.mzML";
    const QString filename = "/home/andrewnichols/Desktop/Data/ExplorisDiscreps/TESTING_20240213RC8_30minE2A22.raw.mzML";

    // const QString &filename = QDir(qApp->applicationDirPath()).filePath("1min.mzML");

    MsReaderMzMLMapped msReaderMzMlMapped;
    e = msReaderMzMlMapped.openFile(filename);
    QCOMPARE(e, eNoError);

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    e = msReaderMzMlMapped.collateMS2MzTargetFrames(&diaTargetFrames);
    QCOMPARE(e, eNoError);

    for (auto it = diaTargetFrames.begin(); it != diaTargetFrames.end(); ++it) {

        const QMap<ScanNumber, ScanPoints*> &scanPoints = it.value();
        int scanPointsInFrame = std::accumulate(
            scanPoints.begin(),
            scanPoints.end(),
            0, [](int sum, const ScanPoints* scanPoints) {return scanPoints->size() + sum;}
            );

        qDebug() << it.key() << scanPointsInFrame;
    }

}


QTEST_MAIN(MsReaderMZMLMappedTests)
#include "MsReaderMZMLMappedTests.moc"
