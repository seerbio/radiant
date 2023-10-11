//
// Created by anichols on 11/07/2021.
//

#include "MsReaderPointerAcc.h"

#include <QString>
#include <QtTest/QtTest>


class MsReaderPointerAccTests : public QObject
{
    Q_OBJECT

public:

    MsReaderPointerAccTests() = default;
    ~MsReaderPointerAccTests() override = default;


//TODO Add more test coverage.
private Q_SLOTS:

    void openFileTest();

    void openFileTest2();

private:

    const QString m_filepath
            = QStringLiteral("/home/anichols/Desktop/RawData/EXP21155_2021ms0609X7_A.raw.mzML");

};

void MsReaderPointerAccTests::openFileTest() {

    ERR_INIT

    MsReaderPointerAcc reader;
    e = reader.openFile(m_filepath);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, MsScanInfo> msScanInfos = reader.ptr->getMsScanInfos();
    QCOMPARE(e, eNoError);
    QCOMPARE(msScanInfos.size(), 33429);

    const MsScanInfo &testInfo = msScanInfos.value(30815);

    QCOMPARE(testInfo.scanNumber, 30815);
    QCOMPARE(testInfo.collisionEnergy, 30);
}


void MsReaderPointerAccTests::openFileTest2() {

    ERR_INIT

    const QString prqFile
        = QStringLiteral("/home/anichols/Desktop/PythiaDIAData/EXP22092_2022ms0742X32_A.raw.mzML.prq");


    MsReaderPointerAcc reader;
    e = reader.openFile(prqFile);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, MsScanInfo> msScanInfos = reader.ptr->getMsScanInfos();
    QCOMPARE(e, eNoError);
    QCOMPARE(msScanInfos.size(), 26010);

    const MsScanInfo &testInfo = msScanInfos.value(10101);
    qDebug() << testInfo.scanNumber << testInfo.collisionEnergy;

    QCOMPARE(testInfo.scanNumber, 10101);
    QCOMPARE(testInfo.collisionEnergy, 28);

}


QTEST_MAIN(MsReaderPointerAccTests)
#include "MsReaderPointerAccTests.moc"
