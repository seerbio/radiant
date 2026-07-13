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

namespace {

bool writeFile(const QString &filePath, const QByteArray &content = {}) {
    QFile file(filePath);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        return false;
    }

    if (!content.isEmpty() && file.write(content) != content.size()) {
        return false;
    }

    file.close();
    return true;
}

bool createMinimalTimsbukSidecar(const QString &sidecarRootPath) {
    if (!QDir().mkpath(sidecarRootPath)) {
        return false;
    }

    const QString ms2DirectoryPath = QDir(sidecarRootPath).filePath("ms2");
    if (!QDir().mkpath(ms2DirectoryPath)) {
        return false;
    }

    if (!writeFile(QDir(sidecarRootPath).filePath("ms1.parquet"))) {
        return false;
    }

    if (!writeFile(QDir(ms2DirectoryPath).filePath("group_0.parquet"))) {
        return false;
    }

    static const QByteArray metadataJson = R"json(
{
  "version": "2.0",
  "created_at": "2026-07-13T22:44:50.135679776+00:00",
  "ms1_peaks": {
    "relative_path": "ms1.parquet",
    "cycle_to_rt_ms": [0, 811, 1770],
    "bucket_size": 4096
  },
  "ms2_window_groups": [
    {
      "id": 0,
      "quadrupole_isolation": [
        {
          "Aabb": {
            "mz": [400.0, 425.0],
            "im": [0.6408845279309161, 0.830315627111824]
          }
        }
      ],
      "group_info": {
        "relative_path": "ms2/group_0.parquet",
        "cycle_to_rt_ms": [0, 918, 1877],
        "bucket_size": 4096
      }
    }
  ]
}
)json";

    return writeFile(QDir(sidecarRootPath).filePath("metadata.json"), metadataJson);
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
    QVERIFY(createMinimalTimsbukSidecar(sidecarRootPath));

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(sidecarRootPath);
    QCOMPARE(e, eNoError);
    QVERIFY(dynamic_cast<MsReaderTimsbukIndex*>(msReaderPointerAcc.ptr.data()) != nullptr);
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), QDir::cleanPath(sidecarRootPath));
    QCOMPARE(msReaderPointerAcc.ptr->isTIMS(), false);
}

void MsReaderPointerAccTests::openFileTest6() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString brukerPath = QDir(temporaryDir.path()).filePath("run.d");
    QVERIFY(QDir().mkpath(brukerPath));

    const QString sidecarRootPath = brukerPath + QStringLiteral(".idx");
    QVERIFY(createMinimalTimsbukSidecar(sidecarRootPath));

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(brukerPath);
    QCOMPARE(e, eNoError);
    QVERIFY(dynamic_cast<MsReaderTimsbukIndex*>(msReaderPointerAcc.ptr.data()) != nullptr);
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), QDir::cleanPath(sidecarRootPath));
    QCOMPARE(msReaderPointerAcc.ptr->isTIMS(), false);
}

void MsReaderPointerAccTests::openFileTest7() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString sidecarRootPath = QDir(temporaryDir.path()).filePath("run.d.idx");
    QVERIFY(createMinimalTimsbukSidecar(sidecarRootPath));

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(sidecarRootPath + QStringLiteral("/"));
    QCOMPARE(e, eNoError);
    QVERIFY(dynamic_cast<MsReaderTimsbukIndex*>(msReaderPointerAcc.ptr.data()) != nullptr);
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), QDir::cleanPath(sidecarRootPath));
    QCOMPARE(msReaderPointerAcc.ptr->isTIMS(), false);
}

void MsReaderPointerAccTests::openFileTest8() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString brukerPath = QDir(temporaryDir.path()).filePath("run.d");
    QVERIFY(QDir().mkpath(brukerPath));

    const QString sidecarRootPath = brukerPath + QStringLiteral(".idx");
    QVERIFY(createMinimalTimsbukSidecar(sidecarRootPath));

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(brukerPath + QStringLiteral("/"));
    QCOMPARE(e, eNoError);
    QVERIFY(dynamic_cast<MsReaderTimsbukIndex*>(msReaderPointerAcc.ptr.data()) != nullptr);
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), QDir::cleanPath(sidecarRootPath));
    QCOMPARE(msReaderPointerAcc.ptr->isTIMS(), false);
}

void MsReaderPointerAccTests::openFileTest9() {

    ERR_INIT

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString sidecarRootPath = QDir(temporaryDir.path()).filePath("run.d.idx");
    QVERIFY(QDir().mkpath(sidecarRootPath));
    QVERIFY(writeFile(QDir(sidecarRootPath).filePath("ms1.parquet")));
    QVERIFY(QDir().mkpath(QDir(sidecarRootPath).filePath("ms2")));
    QVERIFY(writeFile(QDir(sidecarRootPath).filePath("ms2/group_0.parquet")));
    QVERIFY(writeFile(QDir(sidecarRootPath).filePath("metadata.json"), QByteArrayLiteral("{}")));

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(sidecarRootPath);
    QCOMPARE(e, eFileError);
}


QTEST_MAIN(MsReaderPointerAccTests)
#include "MsReaderPointerAccTests.moc"
