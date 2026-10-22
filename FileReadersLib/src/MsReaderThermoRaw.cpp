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
#include <condition_variable>
#include <deque>
#include <functional>
#include <future>
#include <limits>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace {

constexpr int streamBufferSize = 64 * 1024;
constexpr int diagnosticLimit = 64 * 1024;
constexpr quint32 fixedCountStreamVersion = 2;
constexpr quint32 streamingStreamVersion = 3;
constexpr quint32 unknownSpectrumCount = std::numeric_limits<quint32>::max();
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
        bytes->resize(byteCount);
        return readBytes(byteCount, bytes->data());
    }

    bool readBytes(int byteCount, char *destination) {
        QElapsedTimer timeout;
        timeout.start();
        int bytesRead = 0;
        while (bytesRead < byteCount) {
            collectDiagnostics();
            const qint64 incoming = m_process.read(
                destination + bytesRead,
                byteCount - bytesRead
                );
            if (incoming > 0) {
                bytesRead += static_cast<int>(incoming);
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
            if (!m_process.readAllStandardOutput().isEmpty()) {
                return false;
            }
        }
        collectDiagnostics();
        return m_process.readAllStandardOutput().isEmpty()
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
    void collectDiagnostics() {
        m_diagnostics.append(m_process.readAllStandardError());
        if (m_diagnostics.size() > diagnosticLimit) {
            m_diagnostics = m_diagnostics.right(diagnosticLimit);
        }
    }

    QProcess &m_process;
    QByteArray m_diagnostics;
    int m_timeoutMs;
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

class ThermoRawFrameQueue {
public:
    explicit ThermoRawFrameQueue(size_t capacity)
        : m_capacity(capacity) {}

    bool push(QByteArray frame) {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_notFull.wait(lock, [this]() {
            return m_frames.size() < m_capacity || m_aborted;
        });
        if (m_aborted) {
            return false;
        }
        m_frames.push_back(std::move(frame));
        m_notEmpty.notify_one();
        return true;
    }

    bool pop(QByteArray *frame) {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_notEmpty.wait(lock, [this]() {
            return !m_frames.empty() || m_closed;
        });
        if (m_frames.empty()) {
            return false;
        }
        *frame = std::move(m_frames.front());
        m_frames.pop_front();
        m_notFull.notify_one();
        return true;
    }

    void close(bool success, QString error, QString diagnostics, quint32 emittedCount) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_success = success;
        m_error = std::move(error);
        m_diagnostics = std::move(diagnostics);
        m_emittedCount = emittedCount;
        m_closed = true;
        m_notEmpty.notify_all();
        m_notFull.notify_all();
    }

    void abort() {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_aborted = true;
        m_notEmpty.notify_all();
        m_notFull.notify_all();
    }

    bool success() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_success;
    }

    QString error() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_error;
    }

    QString diagnostics() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_diagnostics;
    }

    quint32 emittedCount() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_emittedCount;
    }

private:
    const size_t m_capacity;
    mutable std::mutex m_mutex;
    std::condition_variable m_notEmpty;
    std::condition_variable m_notFull;
    std::deque<QByteArray> m_frames;
    QString m_error;
    QString m_diagnostics;
    quint32 m_emittedCount = 0;
    bool m_closed = false;
    bool m_aborted = false;
    bool m_success = false;
};

struct ParsedThermoRawFrame {
    MsScanInfo scanInfo;
    ScanPoints scanPoints;
};

struct ThermoRawFrameParseResult {
    QString error;
    ParsedThermoRawFrame frame;

    bool success() const {
        return error.isEmpty();
    }
};

ThermoRawFrameParseResult parseThermoRawFrame(
        const QByteArray &frame,
        quint32 emittedSpectrumCount) {
    if (frame.size() < static_cast<int>(4 + metadataByteCount)) {
        return {QStringLiteral("Invalid or truncated spectrum metadata."), {}};
    }
    const quint32 frameByteCount = readLittleEndian<quint32>(frame.constData());
    if (frameByteCount < metadataByteCount
        || frame.size() != static_cast<int>(4 + frameByteCount)) {
        return {QStringLiteral("Invalid or truncated spectrum frame."), {}};
    }

    const char *metadata = frame.constData() + 4;
    MsScanInfo scanInfo;
    scanInfo.scanNumber = readLittleEndian<qint32>(metadata);
    const qint32 nativeScanNumber = readLittleEndian<qint32>(metadata + 4);
    scanInfo.nativeScanNumber = nativeScanNumber;
    scanInfo.msLevel = readLittleEndian<qint32>(metadata + 8);
    const float retentionTimeMinutes = readLittleEndian<float>(metadata + 12);
    scanInfo.scanTime = retentionTimeMinutes * 60.0f;
    scanInfo.collisionEnergy = readLittleEndian<float>(metadata + 16);
    scanInfo.precursorTargetMz = readLittleEndian<float>(metadata + 20);
    scanInfo.isoWindowLower = readLittleEndian<float>(metadata + 24);
    scanInfo.isoWindowUpper = readLittleEndian<float>(metadata + 28);
    const quint32 pointCount = readLittleEndian<quint32>(metadata + 32);
    const quint64 expectedFrameBytes = metadataByteCount
        + static_cast<quint64>(pointCount) * pointByteCount;
    if (frameByteCount != expectedFrameBytes
        || pointCount > maxPointsPerScan
        || scanInfo.scanNumber != static_cast<int>(emittedSpectrumCount + 1)
        || nativeScanNumber <= 0
        || (scanInfo.msLevel != 1 && scanInfo.msLevel != 2)
        || !std::isfinite(scanInfo.scanTime)
        || retentionTimeMinutes < 0.0f
        || !std::isfinite(scanInfo.collisionEnergy)
        || !std::isfinite(scanInfo.precursorTargetMz)
        || !std::isfinite(scanInfo.isoWindowLower)
        || !std::isfinite(scanInfo.isoWindowUpper)) {
        return {QStringLiteral("Invalid spectrum metadata or frame size."), {}};
    }

    ScanPoints scanPoints;
    scanPoints.reserve(static_cast<int>(pointCount));
    const char *points = metadata + metadataByteCount;
    for (quint32 pointIndex = 0; pointIndex < pointCount; ++pointIndex) {
        const char *point = points + pointIndex * pointByteCount;
        const float mz = readLittleEndian<float>(point);
        const float intensity = readLittleEndian<float>(point + 4);
        if (!std::isfinite(mz) || !std::isfinite(intensity)) {
            return {QStringLiteral("Non-finite value in spectrum point data."), {}};
        }
        if (!MathUtils::tZero(intensity)) {
            scanPoints.append(ScanPoint(mz, intensity));
        }
    }
    return {{}, {scanInfo, std::move(scanPoints)}};
}

class ThermoRawFrameParserPool {
public:
    explicit ThermoRawFrameParserPool(size_t workerCount) {
        m_workers.reserve(workerCount);
        try {
            for (size_t i = 0; i < workerCount; ++i) {
                m_workers.emplace_back([this]() { processTasks(); });
            }
        }
        catch (...) {
            {
                std::lock_guard<std::mutex> lock(m_mutex);
                m_stopping = true;
            }
            m_taskAvailable.notify_all();
            for (std::thread &worker : m_workers) {
                if (worker.joinable()) {
                    worker.join();
                }
            }
            throw;
        }
    }

    ~ThermoRawFrameParserPool() {
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_stopping = true;
        }
        m_taskAvailable.notify_all();
        for (std::thread &worker : m_workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }

    std::future<ThermoRawFrameParseResult> submit(QByteArray frame, quint32 frameNumber) {
        auto task = std::make_shared<std::packaged_task<ThermoRawFrameParseResult()>>(
            [frame = std::move(frame), frameNumber]() {
                return parseThermoRawFrame(frame, frameNumber);
            }
            );
        std::future<ThermoRawFrameParseResult> result = task->get_future();
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            m_tasks.emplace_back([task]() { (*task)(); });
        }
        m_taskAvailable.notify_one();
        return result;
    }

private:
    void processTasks() {
        while (true) {
            std::function<void()> task;
            {
                std::unique_lock<std::mutex> lock(m_mutex);
                m_taskAvailable.wait(lock, [this]() {
                    return m_stopping || !m_tasks.empty();
                });
                if (m_stopping && m_tasks.empty()) {
                    return;
                }
                task = std::move(m_tasks.front());
                m_tasks.pop_front();
            }
            task();
        }
    }

    std::mutex m_mutex;
    std::condition_variable m_taskAvailable;
    std::deque<std::function<void()>> m_tasks;
    std::vector<std::thread> m_workers;
    bool m_stopping = false;
};

void produceThermoRawFrames(
        const ReaderInvocation &invocation,
        const QString &filePath,
        ThermoRawFrameQueue *queue) {
    QProcess process;
    process.setProcessChannelMode(QProcess::SeparateChannels);
    QStringList arguments = invocation.arguments;
    arguments.append(filePath);
    BinaryProcessReader reader(process, thermoRawReaderTimeoutMs);
    process.start(invocation.program, arguments, QIODevice::ReadOnly);
    if (!process.waitForStarted(10000)) {
        queue->close(false, QStringLiteral("Cannot start Thermo RAW reader: %1").arg(process.errorString()), {}, 0);
        return;
    }

    const auto fail = [&](const QString &message) {
        stopProcess(&process);
        const QString failure = reader.timedOut()
            ? QStringLiteral("Timed out waiting for the reader process.")
            : message;
        queue->close(false, failure, reader.diagnostics(), 0);
    };

    QByteArray bytes;
    if (!reader.readBytes(16, &bytes)
            || bytes.left(8) != QByteArrayLiteral("RDRRAW01")) {
        fail(QStringLiteral("Missing or unsupported binary stream header."));
        return;
    }
    const quint32 version = readLittleEndian<quint32>(bytes.constData() + 8);
    const quint32 spectrumCount = readLittleEndian<quint32>(bytes.constData() + 12);
    const bool streaming = version == streamingStreamVersion
        && spectrumCount == unknownSpectrumCount;
    if ((version != fixedCountStreamVersion && version != streamingStreamVersion)
        || (!streaming && (spectrumCount == 0
                           || spectrumCount > static_cast<quint32>(std::numeric_limits<int>::max())))
        || (version == fixedCountStreamVersion && spectrumCount == unknownSpectrumCount)) {
        fail(QStringLiteral("Invalid spectrum count in binary stream."));
        return;
    }

    quint32 emittedCount = 0;
    bool footerSeen = false;
    while (streaming || emittedCount < spectrumCount) {
        if (!reader.readBytes(4, &bytes)) {
            fail(QStringLiteral("Truncated spectrum frame."));
            return;
        }
        if (streaming && bytes == QByteArrayLiteral("RAWD")) {
            QByteArray footerBytes;
            if (!reader.readBytes(8, &footerBytes)
                || footerBytes.left(4) != QByteArrayLiteral("ONE1")
                || readLittleEndian<quint32>(footerBytes.constData() + 4) != emittedCount) {
                fail(QStringLiteral("Missing or incorrect binary stream footer."));
                return;
            }
            footerSeen = true;
            break;
        }

        const quint32 frameByteCount = readLittleEndian<quint32>(bytes.constData());
        if (frameByteCount < metadataByteCount
            || frameByteCount > metadataByteCount + maxPointsPerScan * pointByteCount) {
            fail(QStringLiteral("Invalid spectrum frame size."));
            return;
        }
        QByteArray frame;
        frame.resize(static_cast<int>(4 + frameByteCount));
        std::memcpy(frame.data(), bytes.constData(), 4);
        if (!reader.readBytes(static_cast<int>(frameByteCount), frame.data() + 4)) {
            fail(QStringLiteral("Truncated spectrum frame."));
            return;
        }
        if (!queue->push(std::move(frame))) {
            stopProcess(&process);
            return;
        }
        ++emittedCount;
    }

    if ((!streaming && (!reader.readBytes(12, &bytes)
                        || bytes.left(8) != QByteArrayLiteral("RAWDONE1")
                        || readLittleEndian<quint32>(bytes.constData() + 8) != spectrumCount))
        || (streaming && !footerSeen)) {
        fail(QStringLiteral("Missing or incorrect binary stream footer."));
        return;
    }
    if (!reader.finish()) {
        fail(QStringLiteral("Reader failed or emitted data after its footer (exit code %1).").arg(process.exitCode()));
        return;
    }
    queue->close(true, {}, {}, emittedCount);
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

    ThermoRawFrameQueue frameQueue(8);
    std::thread producer([&frameQueue, &invocation, inputPath = inputInfo.absoluteFilePath()]() {
        produceThermoRawFrames(invocation, inputPath, &frameQueue);
    });
    const auto fail = [&](const QString &message) {
        frameQueue.abort();
        if (producer.joinable()) {
            producer.join();
        }
        qCritical().noquote() << "Thermo RAW reader:" << message;
        const QString diagnostics = frameQueue.diagnostics();
        if (!diagnostics.isEmpty()) {
            qCritical().noquote() << diagnostics;
        }
        closeFile();
        return eFileError;
    };

    QByteArray frame;
    quint32 emittedSpectrumCount = 0;
    const unsigned int hardwareThreads = std::thread::hardware_concurrency();
    const size_t workerCount = std::max<size_t>(
        1,
        std::min<size_t>(4, hardwareThreads == 0 ? 1 : hardwareThreads)
        );
    ThermoRawFrameParserPool parserPool(workerCount);
    std::deque<std::future<ThermoRawFrameParseResult>> pendingFrames;
    const auto consumeFrame = [&](std::future<ThermoRawFrameParseResult> future) {
        ThermoRawFrameParseResult result = future.get();
        if (!result.success()) {
            return fail(result.error);
        }
        m_msScanInfo.insert(result.frame.scanInfo.scanNumber, result.frame.scanInfo);
        m_scanPoints.insert(result.frame.scanInfo.scanNumber, std::move(result.frame.scanPoints));
        ++emittedSpectrumCount;
        return eNoError;
    };
    while (frameQueue.pop(&frame)) {
        const quint32 frameNumber = emittedSpectrumCount
            + static_cast<quint32>(pendingFrames.size());
        pendingFrames.push_back(parserPool.submit(std::move(frame), frameNumber));
        if (pendingFrames.size() >= workerCount) {
            if (consumeFrame(std::move(pendingFrames.front())) != eNoError) {
                return eFileError;
            }
            pendingFrames.pop_front();
        }
    }
    while (!pendingFrames.empty()) {
        if (consumeFrame(std::move(pendingFrames.front())) != eNoError) {
            return eFileError;
        }
        pendingFrames.pop_front();
    }

    if (producer.joinable()) {
        producer.join();
    }
    if (!frameQueue.success()) {
        return fail(frameQueue.error().isEmpty()
            ? QStringLiteral("Thermo RAW reader failed.")
            : frameQueue.error());
    }
    if (emittedSpectrumCount != frameQueue.emittedCount()) {
        return fail(QStringLiteral("Thermo RAW reader frame count mismatch."));
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
