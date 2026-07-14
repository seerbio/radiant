//
// Created by anichols on 11/07/2021.
//

#include "MsReaderPointerAcc.h"
#include "MsReaderTimsbukIndex.h"
#include "ParquetReader.h"

#include <QDir>
#include <QFile>
#include <QTemporaryDir>
#include <QString>
#include <QtTest/QtTest>

namespace {

struct TestTimsbukPeakRow : public ParquetReaderInputBase {
    float mz = -1.0f;
    float intensity = -1.0f;
    float mobilityOok0 = -1.0f;
    int cycleIndex = -1;

    QMap<QString, QVariant> map() override {
        return {
            {QStringLiteral("mz"), QVariant(mz)},
            {QStringLiteral("intensity"), QVariant(intensity)},
            {QStringLiteral("mobility_ook0"), QVariant(mobilityOok0)},
            {QStringLiteral("cycle_index"), QVariant(cycleIndex)}
        };
    }
};

TestTimsbukPeakRow makeTestTimsbukPeakRow(
    float mz,
    float intensity,
    float mobilityOok0,
    int cycleIndex
    ) {

    TestTimsbukPeakRow row;
    row.mz = mz;
    row.intensity = intensity;
    row.mobilityOok0 = mobilityOok0;
    row.cycleIndex = cycleIndex;
    return row;
}

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

    const QVector<TestTimsbukPeakRow> ms2Rows = {
        makeTestTimsbukPeakRow(150.0f, 10.0f, 0.90f, 1),
        makeTestTimsbukPeakRow(151.0f, 11.0f, 0.70f, 1),
        makeTestTimsbukPeakRow(152.0f, 12.0f, 0.90f, 2),
        makeTestTimsbukPeakRow(153.0f, 13.0f, 0.70f, 2)
    };
    if (ParquetReader::write(ms2Rows, QDir(ms2DirectoryPath).filePath("group_0.parquet")) != eNoError) {
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
            "mz": [600.0, 625.0],
            "im": [0.830315627111824, 1.0103231398002852]
          }
        },
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
    QCOMPARE(
        msReaderPointerAcc.ptr->filePath(),
        QDir(temporaryDir.path()).filePath("run.d")
        );
    QCOMPARE(msReaderPointerAcc.ptr->isTIMS(), false);

    const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderPointerAcc.ptr->getMsScanInfos();
    QCOMPARE(msScanInfos.size(), 4);

    const QMap<ScanNumber, ScanPoints> scanPoints = msReaderPointerAcc.ptr->getScanPoints();
    QCOMPARE(scanPoints.size(), 4);

    const MzTargetKey lowTargetKey = MsScanInfo::targetKey(400.0f, 425.0f);
    const MzTargetKey highTargetKey = MsScanInfo::targetKey(600.0f, 625.0f);

    QMap<ScanNumber, ScanPoints> lowTargetScanPoints;
    e = msReaderPointerAcc.ptr->getMzTargetScanPoints(lowTargetKey, &lowTargetScanPoints);
    QCOMPARE(e, eNoError);
    QCOMPARE(lowTargetScanPoints.size(), 2);
    QCOMPARE(lowTargetScanPoints.first().size(), 1);
    QCOMPARE(lowTargetScanPoints.first().first().x(), 151.0f);
    QCOMPARE(lowTargetScanPoints.last().first().x(), 153.0f);

    QMap<ScanNumber, ScanPoints> highTargetScanPoints;
    e = msReaderPointerAcc.ptr->getMzTargetScanPoints(highTargetKey, &highTargetScanPoints);
    QCOMPARE(e, eNoError);
    QCOMPARE(highTargetScanPoints.size(), 2);
    QCOMPARE(highTargetScanPoints.first().first().x(), 150.0f);
    QCOMPARE(highTargetScanPoints.last().first().x(), 152.0f);
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
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), QDir::cleanPath(brukerPath));
    QCOMPARE(msReaderPointerAcc.ptr->isTIMS(), false);
    QCOMPARE(msReaderPointerAcc.ptr->getMsScanInfos().size(), 4);
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
    QCOMPARE(
        msReaderPointerAcc.ptr->filePath(),
        QDir(temporaryDir.path()).filePath("run.d")
        );
    QCOMPARE(msReaderPointerAcc.ptr->isTIMS(), false);
    QCOMPARE(msReaderPointerAcc.ptr->getMsScanInfos().size(), 4);
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
    QCOMPARE(msReaderPointerAcc.ptr->filePath(), QDir::cleanPath(brukerPath));
    QCOMPARE(msReaderPointerAcc.ptr->isTIMS(), false);
    QCOMPARE(msReaderPointerAcc.ptr->getMsScanInfos().size(), 4);
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
