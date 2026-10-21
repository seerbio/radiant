#include "MsReaderThermoRaw.h"

#include <QCoreApplication>
#include <QDir>
#include <QElapsedTimer>
#include <QFileInfo>
#include <QProcess>
#include <QStandardPaths>
#include <QtEndian>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <utility>

namespace {

constexpr int streamBufferSize = 64 * 1024;
constexpr int diagnosticLimit = 64 * 1024;
constexpr quint32 streamVersion = 2;
constexpr quint32 metadataByteCount = 36;
constexpr quint32 pointByteCount = 8;
constexpr quint32 maxPointsPerScan = 10'000'000;
constexpr int thermoRawReaderTimeoutMs = 30'000;

template<typename T>
T readLittleEndian(const char *data) {
    T value;
    std::memcpy(&value, data, sizeof(T));
#if Q_BYTE_ORDER == Q_BIG_ENDIAN
    auto *bytes = reinterpret_cast<char *>(&value);
    std::reverse(bytes, bytes + sizeof(T));
#endif
    return value;
}

struct ReaderInvocation {
    QString program;
    QStringList arguments;
    QString error;
};

QString resolveReaderPath(const QString &configuredPath) {
    const QFileInfo configuredInfo(configuredPath);
    if (configuredInfo.isDir()) {
        for (const QString &name : {
                QStringLiteral("RadiantThermoRawReader"),
                QStringLiteral("RadiantThermoRawReader.exe"),
                QStringLiteral("RadiantThermoRawReader.dll")
            }) {
            const QString candidate = QDir(configuredPath).filePath(name);
            if (QFileInfo(candidate).isFile()) {
                return candidate;
            }
        }
        return QString();
    }
    if (configuredInfo.isFile()) {
        return configuredInfo.absoluteFilePath();
    }
    return QStandardPaths::findExecutable(configuredPath);
}

ReaderInvocation readerInvocation(const QString &explicitPath) {
    ReaderInvocation invocation;
    QString configuredPath = explicitPath;
    if (configuredPath.isEmpty()) {
        configuredPath = qEnvironmentVariable("RADIANT_THERMO_RAW_READER");
    }
    QString readerPath;
    if (!configuredPath.isEmpty()) {
        readerPath = resolveReaderPath(configuredPath);
        if (readerPath.isEmpty()) {
            invocation.error = QStringLiteral("Thermo RAW reader not found: %1").arg(configuredPath);
            return invocation;
        }
    }
    else {
        const QString applicationDirectory = QCoreApplication::applicationDirPath();
        for (const QString &candidate : {
                QDir(applicationDirectory).filePath(QStringLiteral("RadiantThermoRawReader")),
                QDir(applicationDirectory).filePath(QStringLiteral("RadiantThermoRawReader.exe")),
                QDir(applicationDirectory).filePath(QStringLiteral("RadiantThermoRawReader.dll")),
                QDir(applicationDirectory).filePath(QStringLiteral("thermoraw")),
                QStringLiteral("RadiantThermoRawReader")
            }) {
            readerPath = resolveReaderPath(candidate);
            if (!readerPath.isEmpty()) {
                break;
            }
        }
    }

    if (readerPath.isEmpty()) {
        invocation.error = QStringLiteral(
            "Thermo RAW input requires RadiantThermoRawReader. Set "
            "RADIANT_THERMO_RAW_READER or use --thermo-raw-reader to select it."
            );
        return invocation;
    }

    const QString suffix = QFileInfo(readerPath).suffix();
    if (suffix.compare(QStringLiteral("dll"), Qt::CaseInsensitive) == 0) {
        invocation.program = QStandardPaths::findExecutable(QStringLiteral("dotnet"));
        if (invocation.program.isEmpty()) {
            invocation.error = QStringLiteral("The .NET runtime (dotnet) is required for %1").arg(readerPath);
            return invocation;
        }
        invocation.arguments.append(readerPath);
    }
#ifndef Q_OS_WIN
    else if (suffix.compare(QStringLiteral("exe"), Qt::CaseInsensitive) == 0) {
        invocation.program = QStandardPaths::findExecutable(QStringLiteral("mono"));
        if (invocation.program.isEmpty()) {
            invocation.error = QStringLiteral("Mono is required for the managed executable %1").arg(readerPath);
            return invocation;
        }
        invocation.arguments.append(readerPath);
    }
#endif
    else {
        invocation.program = readerPath;
    }
    return invocation;
}

class BinaryProcessReader {
public:
    BinaryProcessReader(QProcess &process, int timeoutMs)
        : m_process(process), m_timeoutMs(timeoutMs) {
        QObject::connect(
            &m_process,
            &QProcess::readyReadStandardError,
            &m_process,
            [this]() { collectDiagnostics(); }
            );
    }

    bool readBytes(int byteCount, QByteArray *bytes) {
        QElapsedTimer timeout;
        timeout.start();
        while (m_buffer.size() - m_offset < byteCount) {
            compactBuffer();
            collectDiagnostics();
            const QByteArray incoming = m_process.read(streamBufferSize);
            if (!incoming.isEmpty()) {
                m_buffer.append(incoming);
                continue;
            }
            if (m_process.state() == QProcess::NotRunning) {
                collectDiagnostics();
                return false;
            }
            const qint64 remainingMs = m_timeoutMs - timeout.elapsed();
            if (remainingMs <= 0) {
                m_timedOut = true;
                return false;
            }
            m_process.waitForReadyRead(static_cast<int>(std::min<qint64>(remainingMs, 1000)));
        }
        *bytes = m_buffer.mid(m_offset, byteCount);
        m_offset += byteCount;
        return true;
    }

    bool finish() {
        QElapsedTimer timeout;
        timeout.start();
        while (m_process.state() != QProcess::NotRunning) {
            const qint64 remainingMs = m_timeoutMs - timeout.elapsed();
            if (remainingMs <= 0) {
                m_timedOut = true;
                return false;
            }
            m_process.waitForFinished(static_cast<int>(std::min<qint64>(remainingMs, 1000)));
            collectDiagnostics();
            if (!m_process.read(streamBufferSize).isEmpty()) {
                return false;
            }
        }
        collectDiagnostics();
        return m_buffer.size() == m_offset
            && m_process.readAllStandardOutput().isEmpty()
            && m_process.exitStatus() == QProcess::NormalExit
            && m_process.exitCode() == 0;
    }

    bool timedOut() const {
        return m_timedOut;
    }

    QString diagnostics() {
        collectDiagnostics();
        return QString::fromUtf8(m_diagnostics).trimmed();
    }

private:
    void compactBuffer() {
        if (m_offset != 0) {
            m_buffer.remove(0, m_offset);
            m_offset = 0;
        }
    }

    void collectDiagnostics() {
        m_diagnostics.append(m_process.readAllStandardError());
        if (m_diagnostics.size() > diagnosticLimit) {
            m_diagnostics = m_diagnostics.right(diagnosticLimit);
        }
    }

    QProcess &m_process;
    QByteArray m_buffer;
    QByteArray m_diagnostics;
    int m_timeoutMs;
    int m_offset = 0;
    bool m_timedOut = false;
};

void stopProcess(QProcess *process) {
    if (process->state() != QProcess::NotRunning) {
        process->terminate();
        if (!process->waitForFinished(1000)) {
            process->kill();
            process->waitForFinished(1000);
        }
    }
}

} // namespace

MsReaderThermoRaw::MsReaderThermoRaw(const QString &readerPath)
    : m_readerPath(readerPath) {}

Err MsReaderThermoRaw::closeFile() {
    MsReaderBase::closeFile();
    m_filePath.clear();
    m_fileIsCalibrated = false;
    m_mzMs1Min = std::numeric_limits<float>::max();
    m_mzMs1Max = -1.0f;
    m_mzMs2Min = std::numeric_limits<float>::max();
    m_mzMs2Max = -1.0f;
    return eNoError;
}

Err MsReaderThermoRaw::openFile(const QString &filePath) {
    closeFile();
    const QFileInfo inputInfo(filePath);
    if (inputInfo.suffix().compare(S_GLOBAL_SETTINGS.THERMO_RAW_FILE_EXTENSION, Qt::CaseInsensitive) != 0) {
        return eFileIncorrectTypeError;
    }
    if (!inputInfo.isFile() || !inputInfo.isReadable()) {
        qCritical() << "Cannot read Thermo RAW file:" << filePath;
        return eFileError;
    }

    const ReaderInvocation invocation = readerInvocation(m_readerPath);
    if (!invocation.error.isEmpty()) {
        qCritical().noquote() << invocation.error;
        return eFileError;
    }

    QProcess process;
    process.setProcessChannelMode(QProcess::SeparateChannels);
    QStringList arguments = invocation.arguments;
    arguments.append(inputInfo.absoluteFilePath());
    BinaryProcessReader reader(process, thermoRawReaderTimeoutMs);
    process.start(invocation.program, arguments, QIODevice::ReadOnly);
    if (!process.waitForStarted(10000)) {
        qCritical() << "Cannot start Thermo RAW reader:" << process.errorString();
        return eFileError;
    }

    const auto fail = [&](const QString &message) {
        stopProcess(&process);
        const QString failure = reader.timedOut()
            ? QStringLiteral("Timed out waiting for the reader process.")
            : message;
        qCritical().noquote() << "Thermo RAW reader:" << failure;
        const QString diagnostics = reader.diagnostics();
        if (!diagnostics.isEmpty()) {
            qCritical().noquote() << diagnostics;
        }
        closeFile();
        return eFileError;
    };

    QByteArray bytes;
    if (!reader.readBytes(16, &bytes)
            || bytes.left(8) != QByteArrayLiteral("RDRRAW01")
            || readLittleEndian<quint32>(bytes.constData() + 8) != streamVersion) {
        return fail(QStringLiteral("Missing or unsupported binary stream header."));
    }
    const quint32 spectrumCount = readLittleEndian<quint32>(bytes.constData() + 12);
    if (spectrumCount == 0 || spectrumCount > static_cast<quint32>(std::numeric_limits<int>::max())) {
        return fail(QStringLiteral("Invalid spectrum count in binary stream."));
    }

    for (quint32 spectrumIndex = 0; spectrumIndex < spectrumCount; ++spectrumIndex) {
        if (!reader.readBytes(4, &bytes)) {
            return fail(QStringLiteral("Truncated spectrum frame."));
        }
        const quint32 frameByteCount = readLittleEndian<quint32>(bytes.constData());
        if (frameByteCount < metadataByteCount || !reader.readBytes(metadataByteCount, &bytes)) {
            return fail(QStringLiteral("Invalid or truncated spectrum metadata."));
        }

        MsScanInfo scanInfo;
        scanInfo.scanNumber = readLittleEndian<qint32>(bytes.constData());
        const qint32 nativeScanNumber = readLittleEndian<qint32>(bytes.constData() + 4);
        scanInfo.nativeScanNumber = nativeScanNumber;
        scanInfo.msLevel = readLittleEndian<qint32>(bytes.constData() + 8);
        const float retentionTimeMinutes = readLittleEndian<float>(bytes.constData() + 12);
        scanInfo.scanTime = retentionTimeMinutes * 60.0f;
        scanInfo.collisionEnergy = readLittleEndian<float>(bytes.constData() + 16);
        scanInfo.precursorTargetMz = readLittleEndian<float>(bytes.constData() + 20);
        scanInfo.isoWindowLower = readLittleEndian<float>(bytes.constData() + 24);
        scanInfo.isoWindowUpper = readLittleEndian<float>(bytes.constData() + 28);
        const quint32 pointCount = readLittleEndian<quint32>(bytes.constData() + 32);
        const quint64 expectedFrameBytes = metadataByteCount + static_cast<quint64>(pointCount) * pointByteCount;
        if (frameByteCount != expectedFrameBytes
                || pointCount > maxPointsPerScan
                || scanInfo.scanNumber != static_cast<int>(spectrumIndex + 1)
                || nativeScanNumber <= 0
                || (scanInfo.msLevel != 1 && scanInfo.msLevel != 2)
                || !std::isfinite(scanInfo.scanTime)
                || retentionTimeMinutes < 0.0f
                || !std::isfinite(scanInfo.collisionEnergy)
                || !std::isfinite(scanInfo.precursorTargetMz)
                || !std::isfinite(scanInfo.isoWindowLower)
                || !std::isfinite(scanInfo.isoWindowUpper)) {
            return fail(QStringLiteral("Invalid spectrum metadata or frame size."));
        }

        // Match the mzML reader: cast each value to float before filtering zeros.
        // Read bounded groups of points instead of retaining a second whole scan.
        ScanPoints scanPoints;
        scanPoints.reserve(static_cast<int>(pointCount));
        constexpr quint32 pointsPerBuffer = streamBufferSize / pointByteCount;
        for (quint32 pointIndex = 0; pointIndex < pointCount;) {
            const quint32 pointsToRead = std::min(pointsPerBuffer, pointCount - pointIndex);
            if (!reader.readBytes(static_cast<int>(pointsToRead * pointByteCount), &bytes)) {
                return fail(QStringLiteral("Truncated spectrum point data."));
            }
            for (quint32 i = 0; i < pointsToRead; ++i) {
                const char *point = bytes.constData() + i * pointByteCount;
                const float mz = readLittleEndian<float>(point);
                const float intensity = readLittleEndian<float>(point + 4);
                if (!std::isfinite(mz) || !std::isfinite(intensity)) {
                    return fail(QStringLiteral("Non-finite value in spectrum point data."));
                }
                if (!MathUtils::tZero(intensity)) {
                    scanPoints.append(ScanPoint(mz, intensity));
                }
            }
            pointIndex += pointsToRead;
        }
        m_msScanInfo.insert(scanInfo.scanNumber, scanInfo);
        m_scanPoints.insert(scanInfo.scanNumber, std::move(scanPoints));
    }

    if (!reader.readBytes(12, &bytes)
            || bytes.left(8) != QByteArrayLiteral("RAWDONE1")
            || readLittleEndian<quint32>(bytes.constData() + 8) != spectrumCount) {
        return fail(QStringLiteral("Missing or incorrect binary stream footer."));
    }
    if (!reader.finish()) {
        return fail(QStringLiteral("Reader failed or emitted data after its footer (exit code %1).").arg(process.exitCode()));
    }

    m_filePath = filePath;
    const bool hasEmptyScan = std::any_of(
        m_scanPoints.cbegin(),
        m_scanPoints.cend(),
        [](const ScanPoints &points) { return points.isEmpty(); }
        );
    if (!hasEmptyScan && !getUniqueTandemMsScanInfos().isEmpty()) {
        return printFileInfo();
    }
    // The shared detailed reporter requires populated MS2 scans and nonempty
    // point arrays. Valid MS1-only or empty spectra still remain readable.
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Thermo RAW scan count" << m_msScanInfo.size();
    return eNoError;
}

Err MsReaderThermoRaw::openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
        ) {

    ERR_INIT

    if (!columnToFilterBy.isEmpty()
        && columnToFilterBy.compare(QStringLiteral("scanTime"), Qt::CaseInsensitive) != 0) {
        rrr(eFunctionNotImplemented);
    }

    e = openFile(filePath); ree;
    e = restrictScanTimeRange(
        static_cast<ScanTime>(filterRange.first),
        static_cast<ScanTime>(filterRange.second)
        ); ree;

    ERR_RETURN
}
