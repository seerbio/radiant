//
// Created by anichols on 11/07/2021.
//

#include "MsReaderPointerAcc.h"
#include "MsReaderTimsreader.h"

#include <QDir>
#include <QTemporaryDir>
#include <QString>
#include <QtTest/QtTest>

namespace {

struct TestTimsbukInputPaths {
    QString brukerPath;
    QString sidecarRootPath;
};

TestTimsbukInputPaths createMinimalTimsbukInputPaths(const QString &temporaryPath) {
    TestTimsbukInputPaths paths;
    paths.brukerPath = QDir(temporaryPath).filePath("run.d");
    paths.sidecarRootPath = paths.brukerPath + QStringLiteral(".idx");

    QDir().mkpath(paths.brukerPath);
    return paths;
}

} // namespace

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
    static void openFileTest7();
    static void openFileTest8();
    static void openFileTest9();
    static void openFileTest10();
    static void openFileTest11();

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

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(sidecarRootPath);
    QCOMPARE(e, eFileIncorrectTypeError);
}

void MsReaderPointerAccTests::openFileTest6() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString brukerPath = QDir(temporaryDir.path()).filePath("run.d");
    QVERIFY(QDir().mkpath(brukerPath));

    const QString sidecarRootPath = brukerPath + QStringLiteral(".idx");
    QVERIFY(QDir().mkpath(sidecarRootPath));

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(brukerPath);
    QVERIFY(e != eNoError);
    QVERIFY(dynamic_cast<MsReaderTimsreader*>(msReaderPointerAcc.ptr.data()) != nullptr);
}

void MsReaderPointerAccTests::openFileTest7() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString sidecarRootPath = QDir(temporaryDir.path()).filePath("run.d.idx");

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(sidecarRootPath + QStringLiteral("/"));
    QCOMPARE(e, eFileIncorrectTypeError);
}

void MsReaderPointerAccTests::openFileTest8() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString brukerPath = QDir(temporaryDir.path()).filePath("run.d");
    QVERIFY(QDir().mkpath(brukerPath));

    const QString sidecarRootPath = brukerPath + QStringLiteral(".idx");
    QVERIFY(QDir().mkpath(sidecarRootPath));

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(brukerPath + QStringLiteral("/"));
    QVERIFY(e != eNoError);
    QVERIFY(dynamic_cast<MsReaderTimsreader*>(msReaderPointerAcc.ptr.data()) != nullptr);
}

void MsReaderPointerAccTests::openFileTest9() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString brukerPath = QDir(temporaryDir.path()).filePath("run.d");
    QVERIFY(QDir().mkpath(brukerPath));

    MsReaderPointerAcc msReaderPointerAcc;
    msReaderPointerAcc.setImHandlingMode(ImHandlingMode::Raw4D);
    e = msReaderPointerAcc.openFile(brukerPath);
    QCOMPARE(e, eFunctionNotImplemented);
    QVERIFY(dynamic_cast<MsReaderTimsreader*>(msReaderPointerAcc.ptr.data()) != nullptr);
}

void MsReaderPointerAccTests::openFileTest10() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const TestTimsbukInputPaths paths = createMinimalTimsbukInputPaths(temporaryDir.path());
    QVERIFY(QDir().mkpath(paths.sidecarRootPath));

    const QStringList inputPaths = {
        paths.brukerPath,
        paths.brukerPath + QStringLiteral("/")
    };

    for (const QString &inputPath : inputPaths) {
        MsReaderPointerAcc msReaderPointerAcc;
        e = msReaderPointerAcc.openFile(
            inputPath,
            QStringLiteral("scanTime"),
            {0.0, 0.02}
            );
        QVERIFY(e != eNoError);
        QVERIFY(dynamic_cast<MsReaderTimsreader*>(msReaderPointerAcc.ptr.data()) != nullptr);
    }
}

void MsReaderPointerAccTests::openFileTest11() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const TestTimsbukInputPaths paths = createMinimalTimsbukInputPaths(temporaryDir.path());
    QVERIFY(QDir().mkpath(paths.sidecarRootPath));

    MsReaderPointerAcc msReaderPointerAcc;
    msReaderPointerAcc.setUseLazyLoading(true);
    msReaderPointerAcc.setImHandlingMode(ImHandlingMode::Summed);
    e = msReaderPointerAcc.openFile(paths.brukerPath);
    QVERIFY(e != eNoError);
    QVERIFY(dynamic_cast<MsReaderTimsreader*>(msReaderPointerAcc.ptr.data()) != nullptr);
    QCOMPARE(msReaderPointerAcc.useLazyLoading(), false);
}


QTEST_MAIN(MsReaderPointerAccTests)
#include "MsReaderPointerAccTests.moc"
