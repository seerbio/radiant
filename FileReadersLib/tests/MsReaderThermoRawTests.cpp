#include "MsReaderMzMLMapped.h"
#include "MsReaderPointerAcc.h"
#include "MsReaderThermoRaw.h"

#include <QCryptographicHash>
#include <QDataStream>
#include <QDir>
#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QTemporaryDir>
#include <QtTest/QtTest>

#include <cstring>

namespace {

    QByteArray makeStream() {
        QByteArray bytes;
        QDataStream stream(&bytes, QIODevice::WriteOnly);
        stream.setByteOrder(QDataStream::LittleEndian);
        stream.writeRawData("RDRRAW01", 8);
        stream << quint32(1) << quint32(2);
        for (qint32 scan = 1; scan <= 2; ++scan) {
            stream << quint32(36 + 3 * 16);
            stream << scan << qint32(scan + 10) << scan;
            stream.setFloatingPointPrecision(QDataStream::SinglePrecision);
            stream << float(scan * 0.5) << float(scan == 1 ? -1 : 30)
                   << float(scan == 1 ? -1 : 500)
                   << float(scan == 1 ? -1 : 10)
                   << float(scan == 1 ? -1 : 10);
            stream << quint32(3);
            stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
            stream << double(200.125) << double(100.5)
                   << double(300.25) << double(0)
                   << double(400.5) << double(250.75);
        }
        stream.writeRawData("RAWDONE1", 8);
        stream << quint32(2);
        return bytes;
    }

    bool writeFile(const QString &path, const QByteArray &bytes) {
        QFile file(path);
        return file.open(QIODevice::WriteOnly) && file.write(bytes) == bytes.size();
    }

    QString makeHelper(const QTemporaryDir &directory, int exitCode = 0, bool delay = false) {
        const QString path = QDir(directory.path()).filePath("native reader helper");
        const QByteArray script =
            "#!/bin/sh\n"
            + (delay ? QByteArray("sleep 1\n") : QByteArray())
            + "cat \"$(dirname \"$0\")/payload.bin\"\nexit "
            + QByteArray::number(exitCode) + "\n";
        if (!writeFile(path, script)
            || !QFile::setPermissions(path, QFile::ReadOwner | QFile::WriteOwner | QFile::ExeOwner)) {
            return {};
        }
        return path;
    }

    quint32 floatBits(float value) {
        quint32 bits;
        std::memcpy(&bits, &value, sizeof(bits));
        return bits;
    }

    void hashScan(QCryptographicHash *hash, const MsScanInfo &info, const ScanPoints &points) {
        QByteArray bytes;
        QDataStream stream(&bytes, QIODevice::WriteOnly);
        stream.setByteOrder(QDataStream::LittleEndian);
        stream.setFloatingPointPrecision(QDataStream::SinglePrecision);
        stream << qint32(info.scanNumber) << qint32(info.msLevel)
               << info.scanTime << info.collisionEnergy << info.precursorTargetMz
               << info.isoWindowLower << info.isoWindowUpper << quint32(points.size());
        for (const ScanPoint &point : points) {
            stream << point.x() << point.y();
        }
        hash->addData(bytes);
    }

} // namespace

class MsReaderThermoRawTests : public QObject {
    Q_OBJECT

private Q_SLOTS:
    void nativeStreamAndRouting();
    void filteredReadPath();
    void timesOutWaitingForReader();
    void rejectsIncompleteOrFailedReader_data();
    void rejectsIncompleteOrFailedReader();
    void failedReopenClearsState();
    void realRawMzMLExactEquality();
};

void MsReaderThermoRawTests::nativeStreamAndRouting() {
#ifdef Q_OS_WIN
    QSKIP("The fake helper uses a POSIX shell.");
#endif
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString rawPath = QDir(directory.path()).filePath("test input.RAW");
    QVERIFY(writeFile(rawPath, "raw fixture"));
    QVERIFY(writeFile(QDir(directory.path()).filePath("payload.bin"), makeStream()));
    const QString helper = makeHelper(directory);
    QVERIFY(!helper.isEmpty());

    MsReaderThermoRaw reader(helper);
    QCOMPARE(reader.openFile(rawPath), eNoError);
    QCOMPARE(reader.filePath(), rawPath);
    QVERIFY(reader.isInit());
    const auto infos = reader.getMsScanInfos();
    const auto points = reader.getScanPoints();
    QCOMPARE(infos.size(), 2);
    QCOMPARE(infos[1].msLevel, 1);
    QCOMPARE(infos[1].scanTime, 30.0f);
    QCOMPARE(infos[2].msLevel, 2);
    QCOMPARE(infos[2].scanTime, 60.0f);
    QCOMPARE(infos[2].collisionEnergy, 30.0f);
    QCOMPARE(infos[2].precursorTargetMz, 500.0f);
    QCOMPARE(infos[2].isoWindowLower, 10.0f);
    QCOMPARE(infos[2].isoWindowUpper, 10.0f);
    QCOMPARE(points[1].size(), 2);
    QCOMPARE(points[1][0].x(), 200.125f);
    QCOMPARE(points[1][0].y(), 100.5f);
    QCOMPARE(points[1][1].x(), 400.5f);
    QCOMPARE(points[1][1].y(), 250.75f);

    const QByteArray previous = qgetenv("RADIANT_THERMO_RAW_READER");
    qputenv("RADIANT_THERMO_RAW_READER", helper.toUtf8());
    MsReaderPointerAcc accessor;
    accessor.setUseLazyLoading(true);
    const Err result = accessor.openFile(rawPath);
    if (previous.isNull()) {
        qunsetenv("RADIANT_THERMO_RAW_READER");
    } else {
        qputenv("RADIANT_THERMO_RAW_READER", previous);
    }
    QCOMPARE(result, eNoError);
    QVERIFY(dynamic_cast<MsReaderThermoRaw *>(accessor.ptr.data()) != nullptr);
    QVERIFY(!accessor.useLazyLoading());
    QCOMPARE(accessor.ptr->filePath(), rawPath);
    QCOMPARE(reader.closeFile(), eNoError);
    QVERIFY(!reader.isInit());
    QVERIFY(reader.getScanPoints().isEmpty());
}

void MsReaderThermoRawTests::filteredReadPath() {
#ifdef Q_OS_WIN
    QSKIP("The fake helper uses a POSIX shell.");
#endif
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString rawPath = QDir(directory.path()).filePath("test.raw");
    QVERIFY(writeFile(rawPath, "raw fixture"));
    QVERIFY(writeFile(QDir(directory.path()).filePath("payload.bin"), makeStream()));
    const QString helper = makeHelper(directory);
    QVERIFY(!helper.isEmpty());

    MsReaderThermoRaw reader(helper);
    QCOMPARE(reader.openFile(rawPath, QStringLiteral("scanTime"), qMakePair(29.0, 31.0)), eNoError);
    QCOMPARE(reader.getMsScanInfos().size(), 1);
    QCOMPARE(reader.getMsScanInfos().first().scanTime, 30.0f);
    QCOMPARE(reader.getScanPoints().size(), 1);
    QVERIFY(reader.getScanPoints().contains(1));
    QVERIFY(!reader.getScanPoints().contains(2));
}

void MsReaderThermoRawTests::timesOutWaitingForReader() {
#ifdef Q_OS_WIN
    QSKIP("The fake helper uses a POSIX shell.");
#endif
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString rawPath = QDir(directory.path()).filePath("test.raw");
    QVERIFY(writeFile(rawPath, "raw fixture"));
    QVERIFY(writeFile(QDir(directory.path()).filePath("payload.bin"), makeStream()));
    const QString helper = makeHelper(directory, 0, true);
    QVERIFY(!helper.isEmpty());

    MsReaderThermoRaw reader(helper);
    QVERIFY(reader.openFile(rawPath) != eNoError);
    QVERIFY(!reader.isInit());
}

void MsReaderThermoRawTests::rejectsIncompleteOrFailedReader_data() {
    QTest::addColumn<QByteArray>("payload");
    QTest::addColumn<int>("exitCode");
    QTest::newRow("empty") << QByteArray() << 0;
    QTest::newRow("wrong-magic") << QByteArray("not a native RAW stream") << 0;
    QTest::newRow("truncated-header") << makeStream().left(12) << 0;
    QTest::newRow("truncated-record") << makeStream().left(75) << 0;
    QTest::newRow("missing-footer") << makeStream().chopped(12) << 0;
    QTest::newRow("extra-output") << (makeStream() + "unexpected") << 0;
    QTest::newRow("helper-failure") << makeStream() << 7;
    QByteArray wrongCount = makeStream();
    wrongCount[12] = 3;
    QTest::newRow("count-mismatch") << wrongCount << 0;
    QByteArray wrongFrame = makeStream();
    wrongFrame[16] = 1;
    QTest::newRow("invalid-frame-size") << wrongFrame << 0;
    QByteArray wrongVersion = makeStream();
    wrongVersion[8] = 2;
    QTest::newRow("unsupported-version") << wrongVersion << 0;
}

void MsReaderThermoRawTests::rejectsIncompleteOrFailedReader() {
#ifdef Q_OS_WIN
    QSKIP("The fake helper uses a POSIX shell.");
#endif
    QFETCH(QByteArray, payload);
    QFETCH(int, exitCode);
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString rawPath = QDir(directory.path()).filePath("test.raw");
    QVERIFY(writeFile(rawPath, "raw fixture"));
    QVERIFY(writeFile(QDir(directory.path()).filePath("payload.bin"), payload));
    const QString helper = makeHelper(directory, exitCode);
    QVERIFY(!helper.isEmpty());
    MsReaderThermoRaw reader(helper);
    QVERIFY(reader.openFile(rawPath) != eNoError);
    QVERIFY(!reader.isInit());
    QVERIFY(reader.getMsScanInfos().isEmpty());
    QVERIFY(reader.getScanPoints().isEmpty());
}

void MsReaderThermoRawTests::failedReopenClearsState() {
#ifdef Q_OS_WIN
    QSKIP("The fake helper uses a POSIX shell.");
#endif
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString rawPath = QDir(directory.path()).filePath("test.raw");
    const QString payloadPath = QDir(directory.path()).filePath("payload.bin");
    QVERIFY(writeFile(rawPath, "raw fixture"));
    QVERIFY(writeFile(payloadPath, makeStream()));
    const QString helper = makeHelper(directory);
    QVERIFY(!helper.isEmpty());
    MsReaderThermoRaw reader(helper);
    QCOMPARE(reader.openFile(rawPath), eNoError);
    QVERIFY(reader.isInit());
    QCOMPARE(reader.restrictScanTimeRange(0.0f, 100.0f), eNoError);
    QVERIFY(writeFile(payloadPath, makeStream().left(75)));
    QVERIFY(reader.openFile(rawPath) != eNoError);
    QVERIFY(!reader.isInit());
    QVERIFY(reader.getScanPoints().isEmpty());
    QCOMPARE(reader.openFile(QDir(directory.path()).filePath("missing.raw")) != eNoError, true);
}

void MsReaderThermoRawTests::realRawMzMLExactEquality() {
    const QString rawPath = qEnvironmentVariable("RADIANT_TEST_THERMO_RAW");
    const QString mzmlPath = qEnvironmentVariable("RADIANT_TEST_THERMO_MZML");
    if (rawPath.isEmpty() || mzmlPath.isEmpty()) {
        QSKIP("Set RADIANT_TEST_THERMO_RAW and RADIANT_TEST_THERMO_MZML for a real paired run.");
    }

    MsReaderThermoRaw raw;
    MsReaderMzMLMapped mzml;
    QCOMPARE(raw.openFile(rawPath), eNoError);
    QCOMPARE(mzml.openFile(mzmlPath), eNoError);
    const auto rawInfos = raw.getMsScanInfos();
    const auto mzmlInfos = mzml.getMsScanInfos();
    const auto rawPoints = raw.getScanPoints();
    const auto mzmlPoints = mzml.getScanPoints();
    QCOMPARE(rawInfos.keys(), mzmlInfos.keys());
    QCOMPARE(rawPoints.keys(), mzmlPoints.keys());
    QCryptographicHash rawHash(QCryptographicHash::Sha256);
    QCryptographicHash mzmlHash(QCryptographicHash::Sha256);
    qint64 totalPoints = 0;
    int ms1 = 0;
    int ms2 = 0;
    for (auto it = rawInfos.cbegin(); it != rawInfos.cend(); ++it) {
        const int scan = it.key();
        const MsScanInfo &a = it.value();
        const MsScanInfo &b = mzmlInfos[scan];
        QCOMPARE(a.msLevel, b.msLevel);
        QCOMPARE(a.scanNumber, b.scanNumber);
        QCOMPARE(a.ionMobilityIndex, b.ionMobilityIndex);
        QCOMPARE(a.nativeFrameNumber, b.nativeFrameNumber);
        QCOMPARE(a.nativeScanNumber, b.nativeScanNumber);
        QCOMPARE(a.scanOffsetStart, b.scanOffsetStart);
        QCOMPARE(a.scanOffsetEnd, b.scanOffsetEnd);
        const float rawValues[] = {a.scanTime, a.collisionEnergy, a.precursorTargetMz,
                                   a.isoWindowLower, a.isoWindowUpper, a.ionMobilityDriftTime,
                                   a.ionMobilityWindowLower, a.ionMobilityWindowUpper};
        const float mzmlValues[] = {b.scanTime, b.collisionEnergy, b.precursorTargetMz,
                                    b.isoWindowLower, b.isoWindowUpper, b.ionMobilityDriftTime,
                                    b.ionMobilityWindowLower, b.ionMobilityWindowUpper};
        for (int field = 0; field < 8; ++field) {
            QVERIFY2(floatBits(rawValues[field]) == floatBits(mzmlValues[field]),
                     qPrintable(QString("Scan %1 metadata field %2 differs: RAW=%3 mzML=%4")
                                    .arg(scan).arg(field)
                                    .arg(rawValues[field], 0, 'g', 9)
                                    .arg(mzmlValues[field], 0, 'g', 9)));
        }
        const ScanPoints &aPoints = rawPoints[scan];
        const ScanPoints &bPoints = mzmlPoints[scan];
        QCOMPARE(aPoints.size(), bPoints.size());
        for (int point = 0; point < aPoints.size(); ++point) {
            QVERIFY2(floatBits(aPoints[point].x()) == floatBits(bPoints[point].x())
                         && floatBits(aPoints[point].y()) == floatBits(bPoints[point].y()),
                     qPrintable(QString("Scan %1 peak %2 differs").arg(scan).arg(point)));
        }
        totalPoints += aPoints.size();
        ms1 += a.msLevel == 1;
        ms2 += a.msLevel == 2;
        hashScan(&rawHash, a, aPoints);
        hashScan(&mzmlHash, b, bPoints);
    }
    const QString rawDigest = rawHash.result().toHex();
    const QString mzmlDigest = mzmlHash.result().toHex();
    QCOMPARE(rawDigest, mzmlDigest);
    const QJsonObject report{
        {"raw", rawPath}, {"mzml", mzmlPath}, {"exact_equal", true},
        {"spectra", rawInfos.size()}, {"ms1_spectra", ms1}, {"ms2_spectra", ms2},
        {"retained_peaks", double(totalPoints)}, {"raw_spectra_sha256", rawDigest},
        {"mzml_spectra_sha256", mzmlDigest}
    };
    qInfo().noquote() << QJsonDocument(report).toJson(QJsonDocument::Compact);
    const QString reportPath = qEnvironmentVariable("RADIANT_TEST_THERMO_REPORT");
    if (!reportPath.isEmpty()) {
        QVERIFY(writeFile(reportPath, QJsonDocument(report).toJson(QJsonDocument::Indented)));
    }
}

QTEST_GUILESS_MAIN(MsReaderThermoRawTests)
#include "MsReaderThermoRawTests.moc"
