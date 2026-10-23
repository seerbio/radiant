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
#include <limits>

namespace {

    struct NativeScanFixture {
        qint32 nativeScanNumber;
        qint32 msLevel;
        float retentionTimeMinutes;
        float collisionEnergy;
        float precursorTargetMz;
        float isoWindowLower;
        float isoWindowUpper;
        QVector<QPair<float, float>> peaks;
    };

    void appendLittleEndianUInt32(QByteArray *bytes, quint32 value) {
        for (int byte = 0; byte < 4; ++byte) {
            bytes->append(static_cast<char>((value >> (byte * 8)) & 0xff));
        }
    }

    quint32 littleEndianUInt32(const QByteArray &bytes, int offset) {
        quint32 value = 0;
        for (int byte = 0; byte < 4; ++byte) {
            value |= static_cast<quint32>(static_cast<unsigned char>(bytes[offset + byte]))
                << (byte * 8);
        }
        return value;
    }

    void setLittleEndianUInt32(QByteArray *bytes, int offset, quint32 value) {
        for (int byte = 0; byte < 4; ++byte) {
            (*bytes)[offset + byte] = static_cast<char>((value >> (byte * 8)) & 0xff);
        }
    }

    NativeScanFixture makeNativeScan(
            qint32 nativeScanNumber,
            qint32 msLevel,
            float retentionTimeMinutes,
            float mz,
            float intensity) {
        const bool isMs2 = msLevel == 2;
        NativeScanFixture scan{
            nativeScanNumber,
            msLevel,
            retentionTimeMinutes,
            isMs2 ? 30.0f : -1.0f,
            isMs2 ? 500.0f : -1.0f,
            isMs2 ? 10.0f : -1.0f,
            isMs2 ? 10.0f : -1.0f,
            {}
        };
        scan.peaks.append(qMakePair(mz, intensity));
        return scan;
    }

    QByteArray makeNativeFrame(const NativeScanFixture &scan) {
        QByteArray frame;
        QDataStream stream(&frame, QIODevice::WriteOnly);
        stream.setByteOrder(QDataStream::LittleEndian);
        stream.setFloatingPointPrecision(QDataStream::SinglePrecision);
        stream << quint32(36 + scan.peaks.size() * 8)
               << scan.nativeScanNumber << scan.nativeScanNumber << scan.msLevel
               << scan.retentionTimeMinutes << scan.collisionEnergy
               << scan.precursorTargetMz << scan.isoWindowLower << scan.isoWindowUpper
               << quint32(scan.peaks.size());
        for (const auto &peak : scan.peaks) {
            stream << peak.first << peak.second;
        }
        return frame;
    }

    QByteArray makeNativeBatch(const QList<NativeScanFixture> &scans) {
        QByteArray batchFrames;
        for (const NativeScanFixture &scan : scans) {
            batchFrames.append(makeNativeFrame(scan));
        }

        QByteArray batch;
        appendLittleEndianUInt32(&batch, static_cast<quint32>(batchFrames.size()));
        appendLittleEndianUInt32(&batch, static_cast<quint32>(scans.size()));
        batch.append(batchFrames);
        return batch;
    }

    QByteArray makeNativeKeyedStream(const QList<QList<NativeScanFixture>> &batches) {
        QByteArray bytes("RDRRAW01", 8);
        appendLittleEndianUInt32(&bytes, 4);
        appendLittleEndianUInt32(&bytes, std::numeric_limits<quint32>::max());

        quint32 scanCount = 0;
        for (const QList<NativeScanFixture> &batch : batches) {
            bytes.append(makeNativeBatch(batch));
            scanCount += static_cast<quint32>(batch.size());
        }
        bytes.append("RAWDONE1", 8);
        appendLittleEndianUInt32(&bytes, scanCount);
        return bytes;
    }

    QByteArray makeUnsupportedVersionHeader(quint32 version, quint32 spectrumCount) {
        QByteArray bytes("RDRRAW01", 8);
        appendLittleEndianUInt32(&bytes, version);
        appendLittleEndianUInt32(&bytes, spectrumCount);
        return bytes;
    }

    QByteArray makeStream() {
        NativeScanFixture ms1 = makeNativeScan(11, 1, 0.5f, 400.5f, 250.75f);
        ms1.peaks.append(qMakePair(300.25f, 0.0f));
        ms1.peaks.append(qMakePair(200.125f, 100.5f));
        NativeScanFixture ms2 = makeNativeScan(12, 2, 1.0f, 400.5f, 250.75f);
        ms2.peaks.append(qMakePair(300.25f, 0.0f));
        ms2.peaks.append(qMakePair(200.125f, 100.5f));
        return makeNativeKeyedStream({{ms1, ms2}});
    }

    QByteArray makeOutOfOrderMultiBatchStream() {
        QList<QList<NativeScanFixture>> batches;
        batches.append({
            makeNativeScan(90, 2, 9.0f, 900.25f, 9000.5f),
            makeNativeScan(10, 1, 1.0f, 100.25f, 1000.5f)
        });
        batches.append({makeNativeScan(70, 2, 7.0f, 700.25f, 7000.5f)});
        batches.append({
            makeNativeScan(30, 1, 3.0f, 300.25f, 3000.5f),
            makeNativeScan(50, 2, 5.0f, 500.25f, 5000.5f)
        });
        return makeNativeKeyedStream(batches);
    }

    QByteArray makeDuplicateNativeScanStream() {
        QList<QList<NativeScanFixture>> batches;
        batches.append({makeNativeScan(42, 1, 1.0f, 420.25f, 4200.5f)});
        batches.append({makeNativeScan(42, 2, 2.0f, 421.25f, 4210.5f)});
        return makeNativeKeyedStream(batches);
    }

    bool writeFile(const QString &path, const QByteArray &bytes) {
        QFile file(path);
        return file.open(QIODevice::WriteOnly) && file.write(bytes) == bytes.size();
    }

    QString makeHelper(const QTemporaryDir &directory, int exitCode = 0) {
        const QString path = QDir(directory.path()).filePath("native reader helper");
        const QByteArray script =
            QByteArray("#!/bin/sh\n")
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
    void protocol4BatchesPreserveNativeIdsAndPeakAssociations();
    void filteredReadPath();
    void rejectsIncompleteOrFailedReader_data();
    void rejectsIncompleteOrFailedReader();
    void rejectsNonEmittedProtocolVersions_data();
    void rejectsNonEmittedProtocolVersions();
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
    QCOMPARE(infos[1].nativeScanNumber, 11);
    QCOMPARE(infos[1].scanTime, 30.0f);
    QCOMPARE(infos[2].msLevel, 2);
    QCOMPARE(infos[2].nativeScanNumber, 12);
    QCOMPARE(infos[2].scanTime, 60.0f);
    QCOMPARE(infos[2].collisionEnergy, 30.0f);
    QCOMPARE(infos[2].precursorTargetMz, 500.0f);
    QCOMPARE(infos[2].isoWindowLower, 10.0f);
    QCOMPARE(infos[2].isoWindowUpper, 10.0f);
    QCOMPARE(points[1].size(), 2);
    QCOMPARE(points[1][0].x(), 400.5f);
    QCOMPARE(points[1][0].y(), 250.75f);
    QCOMPARE(points[1][1].x(), 200.125f);
    QCOMPARE(points[1][1].y(), 100.5f);

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

void MsReaderThermoRawTests::protocol4BatchesPreserveNativeIdsAndPeakAssociations() {
#ifdef Q_OS_WIN
    QSKIP("The fake helper uses a POSIX shell.");
#endif
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString rawPath = QDir(directory.path()).filePath("test.raw");
    QVERIFY(writeFile(rawPath, "raw fixture"));
    QVERIFY(writeFile(
        QDir(directory.path()).filePath("payload.bin"),
        makeOutOfOrderMultiBatchStream()));
    const QString helper = makeHelper(directory);
    QVERIFY(!helper.isEmpty());

    MsReaderThermoRaw reader(helper);
    QCOMPARE(reader.openFile(rawPath), eNoError);
    const auto infos = reader.getMsScanInfos();
    const auto points = reader.getScanPoints();
    const qint32 expectedNativeScanNumbers[] = {10, 30, 50, 70, 90};
    const float expectedMz[] = {100.25f, 300.25f, 500.25f, 700.25f, 900.25f};
    const float expectedIntensity[] = {1000.5f, 3000.5f, 5000.5f, 7000.5f, 9000.5f};
    const float expectedScanTimes[] = {60.0f, 180.0f, 300.0f, 420.0f, 540.0f};
    QCOMPARE(infos.size(), 5);
    QCOMPARE(points.size(), 5);
    for (int index = 0; index < 5; ++index) {
        const ScanNumber sequentialScanNumber = index + 1;
        const MsScanInfo &scanInfo = infos[sequentialScanNumber];
        QCOMPARE(scanInfo.scanNumber, sequentialScanNumber);
        QCOMPARE(scanInfo.nativeScanNumber, expectedNativeScanNumbers[index]);
        QCOMPARE(scanInfo.scanTime, expectedScanTimes[index]);
        QCOMPARE(points[sequentialScanNumber].size(), 1);
        QCOMPARE(points[sequentialScanNumber][0].x(), expectedMz[index]);
        QCOMPARE(points[sequentialScanNumber][0].y(), expectedIntensity[index]);
    }
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
    QByteArray zeroBatchCount = makeStream();
    setLittleEndianUInt32(&zeroBatchCount, 20, 0);
    QTest::newRow("zero-batch-count") << zeroBatchCount << 0;
    QByteArray excessiveBatchCount = makeStream();
    setLittleEndianUInt32(&excessiveBatchCount, 20, 3);
    QTest::newRow("batch-count-exceeds-frames") << excessiveBatchCount << 0;
    QByteArray shortBatchLength = makeStream();
    setLittleEndianUInt32(
        &shortBatchLength, 16,
        littleEndianUInt32(shortBatchLength, 16) - 1);
    QTest::newRow("truncated-batch-length") << shortBatchLength << 0;
    QByteArray longBatchLength = makeStream();
    setLittleEndianUInt32(
        &longBatchLength, 16,
        littleEndianUInt32(longBatchLength, 16) + 1);
    QTest::newRow("batch-length-includes-extra-byte") << longBatchLength << 0;
    QByteArray understatedBatchCount = makeStream();
    setLittleEndianUInt32(&understatedBatchCount, 20, 1);
    QTest::newRow("batch-count-understates-payload") << understatedBatchCount << 0;
    QByteArray invalidBatchLength = makeStream();
    setLittleEndianUInt32(&invalidBatchLength, 16, 1);
    QTest::newRow("invalid-batch-length") << invalidBatchLength << 0;
    QTest::newRow("protocol4-fixed-count-header")
        << makeUnsupportedVersionHeader(4, 2) << 0;
    QTest::newRow("duplicate-native-scan") << makeDuplicateNativeScanStream() << 0;
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

void MsReaderThermoRawTests::rejectsNonEmittedProtocolVersions_data() {
    QTest::addColumn<QByteArray>("payload");
    QTest::newRow("protocol-2") << makeUnsupportedVersionHeader(2, 1);
    QTest::newRow("protocol-3")
        << makeUnsupportedVersionHeader(3, std::numeric_limits<quint32>::max());
    QTest::newRow("unknown-protocol")
        << makeUnsupportedVersionHeader(5, std::numeric_limits<quint32>::max());
}

void MsReaderThermoRawTests::rejectsNonEmittedProtocolVersions() {
#ifdef Q_OS_WIN
    QSKIP("The fake helper uses a POSIX shell.");
#endif
    QFETCH(QByteArray, payload);
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString rawPath = QDir(directory.path()).filePath("test.raw");
    QVERIFY(writeFile(rawPath, "raw fixture"));
    QVERIFY(writeFile(QDir(directory.path()).filePath("payload.bin"), payload));
    const QString helper = makeHelper(directory);
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
        //QCOMPARE(a.nativeScanNumber, b.nativeScanNumber);
        QCOMPARE(a.scanOffsetStart, b.scanOffsetStart);
        QCOMPARE(a.scanOffsetEnd, b.scanOffsetEnd);
        const float rawValues[] = {a.scanTime, a.collisionEnergy, a.precursorTargetMz,
                                   a.isoWindowLower, a.isoWindowUpper, a.ionMobilityDriftTime,
                                   a.ionMobilityWindowLower, a.ionMobilityWindowUpper
                                   };
        const float mzmlValues[] = {b.scanTime, b.collisionEnergy, b.precursorTargetMz,
                                    b.isoWindowLower, b.isoWindowUpper, b.ionMobilityDriftTime,
                                    b.ionMobilityWindowLower, b.ionMobilityWindowUpper
                                    };
        for (int field = 0; field < std::size(rawValues); ++field) {
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
