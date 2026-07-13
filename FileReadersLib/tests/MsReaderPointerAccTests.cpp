//
// Created by anichols on 11/07/2021.
//

#include "MsReaderPointerAcc.h"
#include "MsReaderTimsbukIndex.h"

#include <QDir>
#include <QFile>
#include <QTemporaryDir>
#include <QString>
#include <QtTest/QtTest>


class MsReaderPointerAccTests : public QObject
{
    Q_OBJECT

public:

    MsReaderPointerAccTests() = default;
    ~MsReaderPointerAccTests() override = default;


private Q_SLOTS:

    static void openFileTest();
    static void openFileTest2();
    static void openFileTest3();
    static void openFileTest4();
    static void openFileTest5();
    static void openFileTest6();



};

void MsReaderPointerAccTests::openFileTest() {

    ERR_INIT

    const QString mzMLFilePath = QDir(qApp->applicationDirPath()).filePath("1min.mzML");
    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(mzMLFilePath);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderPointerAcc.ptr->getMsScanInfos();
    QCOMPARE(e, eNoError);
    QCOMPARE(msScanInfos.size(), 372);

    const MsScanInfo &testInfo = msScanInfos.value(372);
    qDebug() << testInfo.msLevel;

    QCOMPARE(testInfo.scanNumber, 372);
    QCOMPARE(testInfo.msLevel, 1);
}


void MsReaderPointerAccTests::openFileTest2() {

    ERR_INIT

    const QString prqFFFilePath
        = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(prqFFFilePath);
    QCOMPARE(e, eNoError);

    QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderPointerAcc.ptr->getMsScanInfos();
    QCOMPARE(e, eNoError);
    QCOMPARE(msScanInfos.size(), 5000);

    const MsScanInfo &testInfo = msScanInfos.value(5000);
    QCOMPARE(testInfo.scanNumber, 5000);
    QCOMPARE(testInfo.collisionEnergy, 28);
}

void MsReaderPointerAccTests::openFileTest3() {

    ERR_INIT

    const QString prqFFFilePath
        = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

    const QString column = QStringLiteral("scanNumber");

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(
            prqFFFilePath,
            column,
            {1, 20}
            );
    QCOMPARE(e, eNoError);

    const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderPointerAcc.ptr->getMsScanInfos();
    QCOMPARE(msScanInfos.size(), 20);
    QCOMPARE(msScanInfos.first().scanNumber, 1);
    QCOMPARE(msScanInfos.last().scanNumber, 20);

}

void MsReaderPointerAccTests::openFileTest4() {

    ERR_INIT

    const QString prqFFFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

    const QString column = QStringLiteral("scanNumber");

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(
            prqFFFilePath,
            column
    );
    QCOMPARE(e, eFunctionNotImplemented);

    // const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderPointerAcc.ptr->getMsScanInfos();
    // QCOMPARE(msScanInfos.first().scanNumber, 1);
    // QCOMPARE(msScanInfos.first().msLevel, -1);

}

void MsReaderPointerAccTests::openFileTest5() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString sidecarRootPath = QDir(temporaryDir.path()).filePath("run.d.idx");
    QVERIFY(QDir().mkpath(sidecarRootPath));

    QFile metadataFile(QDir(sidecarRootPath).filePath("metadata.json"));
    QVERIFY(metadataFile.open(QIODevice::WriteOnly | QIODevice::Text));
    metadataFile.write("{}");
    metadataFile.close();

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(sidecarRootPath);
    QCOMPARE(e, eFunctionNotImplemented);
    QVERIFY(dynamic_cast<MsReaderTimsbukIndex*>(msReaderPointerAcc.ptr.data()) != nullptr);
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), QDir::cleanPath(sidecarRootPath));
}

void MsReaderPointerAccTests::openFileTest6() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString brukerPath = QDir(temporaryDir.path()).filePath("run.d");
    QVERIFY(QDir().mkpath(brukerPath));

    const QString sidecarRootPath = brukerPath + QStringLiteral(".idx");
    QVERIFY(QDir().mkpath(sidecarRootPath));

    QFile metadataFile(QDir(sidecarRootPath).filePath("metadata.json"));
    QVERIFY(metadataFile.open(QIODevice::WriteOnly | QIODevice::Text));
    metadataFile.write("{}");
    metadataFile.close();

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(brukerPath);
    QCOMPARE(e, eFunctionNotImplemented);
    QVERIFY(dynamic_cast<MsReaderTimsbukIndex*>(msReaderPointerAcc.ptr.data()) != nullptr);
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), QDir::cleanPath(sidecarRootPath));
}


QTEST_MAIN(MsReaderPointerAccTests)
#include "MsReaderPointerAccTests.moc"
