#include "PerFileLog.h"

#include <QCoreApplication>
#include <QFile>
#include <QSet>
#include <QTemporaryDir>
#include <QtTest/QtTest>

#include <atomic>
#include <thread>
#include <vector>

namespace {
std::atomic<int> forwardedCount{0};

void countForwardedMessage(QtMsgType, const QMessageLogContext &, const QString &) {
    ++forwardedCount;
}

class ScopedCountingHandler
{
public:
    ScopedCountingHandler()
        : m_previous(qInstallMessageHandler(countForwardedMessage)) {
        forwardedCount = 0;
    }

    ~ScopedCountingHandler() {
        qInstallMessageHandler(m_previous);
    }

private:
    QtMessageHandler m_previous;
};

QByteArray readLog(const QString &path) {
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) {
        return {};
    }
    return file.readAll();
}
}

class RadiantPerFileLogTests : public QObject
{
    Q_OBJECT

private Q_SLOTS:
    void init();
    void workerMessagesRemainIntactAndAreForwarded();
    void sequentialLogsAreIsolatedAndRestorePreviousHandler();
    void destructorClosesScopeAndRestoresPreviousHandler();
    void refusesASecondActiveLogWithoutDisruptingTheFirst();
    void refusesDirectoryAsLogFile();
    void reportsWriteFailure();
    void preservesSeverityAndMessagePattern();
    void suppressesEmptyFormattedMessages();
};

void RadiantPerFileLogTests::init() {
    qSetMessagePattern(QStringLiteral("%{message}"));
}

void RadiantPerFileLogTests::workerMessagesRemainIntactAndAreForwarded() {
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString path = directory.filePath(QStringLiteral("workers.log"));
    ScopedCountingHandler console;
    PerFileLog log;
    QVERIFY2(log.open(path), qPrintable(log.errorString()));

    constexpr int workerCount = 4;
    constexpr int messagesPerWorker = 100;
    std::vector<std::thread> workers;
    for (int worker = 0; worker < workerCount; ++worker) {
        workers.emplace_back([worker] {
            for (int message = 0; message < messagesPerWorker; ++message) {
                qInfo().noquote() << QStringLiteral("worker-%1-line-%2").arg(worker).arg(message);
            }
        });
    }
    for (std::thread &worker : workers) {
        worker.join();
    }

    QVERIFY2(log.close(), qPrintable(log.errorString()));
    QCOMPARE(forwardedCount.load(), workerCount * messagesPerWorker);
    const QByteArray contents = readLog(path);
    QVERIFY(contents.endsWith('\n'));
    const QList<QByteArray> lines = contents.trimmed().split('\n');
    QCOMPARE(lines.size(), workerCount * messagesPerWorker);
    QSet<QByteArray> uniqueLines;
    for (const QByteArray &line : lines) {
        uniqueLines.insert(line);
    }
    QCOMPARE(uniqueLines.size(), workerCount * messagesPerWorker);
    for (int worker = 0; worker < workerCount; ++worker) {
        for (int message = 0; message < messagesPerWorker; ++message) {
            const QByteArray expected = QStringLiteral("worker-%1-line-%2")
                .arg(worker).arg(message).toUtf8();
            QVERIFY2(uniqueLines.contains(expected), expected.constData());
        }
    }
}

void RadiantPerFileLogTests::sequentialLogsAreIsolatedAndRestorePreviousHandler() {
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString firstPath = directory.filePath(QStringLiteral("first.log"));
    const QString secondPath = directory.filePath(QStringLiteral("second.log"));
    ScopedCountingHandler console;
    PerFileLog log;

    QVERIFY2(log.open(firstPath), qPrintable(log.errorString()));
    qInfo().noquote() << "first-input";
    QVERIFY2(log.close(), qPrintable(log.errorString()));
    qInfo().noquote() << "between-inputs";
    QVERIFY2(log.open(secondPath), qPrintable(log.errorString()));
    qWarning().noquote() << "second-input";
    QVERIFY2(log.close(), qPrintable(log.errorString()));
    qInfo().noquote() << "batch-summary";

    QCOMPARE(readLog(firstPath), QByteArray("first-input\n"));
    QCOMPARE(readLog(secondPath), QByteArray("second-input\n"));
    QCOMPARE(forwardedCount.load(), 4);
    const QtMessageHandler restored = qInstallMessageHandler(countForwardedMessage);
    QVERIFY(restored == &countForwardedMessage);
}

void RadiantPerFileLogTests::destructorClosesScopeAndRestoresPreviousHandler() {
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString path = directory.filePath(QStringLiteral("scope.log"));
    ScopedCountingHandler console;
    {
        PerFileLog log;
        QVERIFY2(log.open(path), qPrintable(log.errorString()));
        qInfo().noquote() << "inside-scope";
    }
    qInfo().noquote() << "after-scope";

    QCOMPARE(readLog(path), QByteArray("inside-scope\n"));
    QCOMPARE(forwardedCount.load(), 2);
    const QtMessageHandler restored = qInstallMessageHandler(countForwardedMessage);
    QVERIFY(restored == &countForwardedMessage);
}

void RadiantPerFileLogTests::refusesASecondActiveLogWithoutDisruptingTheFirst() {
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString firstPath = directory.filePath(QStringLiteral("first.log"));
    const QString secondPath = directory.filePath(QStringLiteral("second.log"));
    ScopedCountingHandler console;
    PerFileLog first;
    PerFileLog second;

    QVERIFY2(first.open(firstPath), qPrintable(first.errorString()));
    QVERIFY(!second.open(secondPath));
    QVERIFY(!second.errorString().isEmpty());
    QVERIFY(!second.close());
    qInfo().noquote() << "first-still-active";
    QVERIFY2(first.close(), qPrintable(first.errorString()));

    QCOMPARE(readLog(firstPath), QByteArray("first-still-active\n"));
    QVERIFY(!QFile::exists(secondPath));
    QCOMPARE(forwardedCount.load(), 1);
}

void RadiantPerFileLogTests::refusesDirectoryAsLogFile() {
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    ScopedCountingHandler console;
    PerFileLog log;

    QVERIFY(!log.open(directory.path()));
    QVERIFY(!log.errorString().isEmpty());
    QVERIFY(!log.close());
    qInfo().noquote() << "console-after-open-failure";
    QCOMPARE(forwardedCount.load(), 1);
    const QtMessageHandler restored = qInstallMessageHandler(countForwardedMessage);
    QVERIFY(restored == &countForwardedMessage);
}

void RadiantPerFileLogTests::reportsWriteFailure() {
#if defined(Q_OS_LINUX)
    if (!QFile::exists(QStringLiteral("/dev/full"))) {
        QSKIP("/dev/full is unavailable");
    }
    ScopedCountingHandler console;
    PerFileLog log;
    QVERIFY2(log.open(QStringLiteral("/dev/full")), qPrintable(log.errorString()));
    qInfo().noquote() << "cannot-be-written";
    QVERIFY(!log.errorString().isEmpty());
    QVERIFY(!log.close());
    qInfo().noquote() << "console-after-write-failure";
    QCOMPARE(forwardedCount.load(), 2);
    const QtMessageHandler restored = qInstallMessageHandler(countForwardedMessage);
    QVERIFY(restored == &countForwardedMessage);
#else
    QSKIP("Write failure check requires Linux /dev/full");
#endif
}

void RadiantPerFileLogTests::preservesSeverityAndMessagePattern() {
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString path = directory.filePath(QStringLiteral("formatted.log"));
    qSetMessagePattern(QStringLiteral("severity=%{type}|%{message}"));
    ScopedCountingHandler console;
    PerFileLog log;
    QVERIFY2(log.open(path), qPrintable(log.errorString()));
    qDebug().noquote() << "debug-detail";
    qInfo().noquote() << "progress";
    qWarning().noquote() << "warning-detail";
    qCritical().noquote() << "failure-detail";
    QVERIFY2(log.close(), qPrintable(log.errorString()));

    QCOMPARE(readLog(path), QByteArray(
        "severity=debug|debug-detail\n"
        "severity=info|progress\n"
        "severity=warning|warning-detail\n"
        "severity=critical|failure-detail\n"));
    QCOMPARE(forwardedCount.load(), 4);
}

void RadiantPerFileLogTests::suppressesEmptyFormattedMessages() {
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString path = directory.filePath(QStringLiteral("suppressed.log"));
    qSetMessagePattern(QStringLiteral("%{if-warning}%{message}%{endif}"));
    ScopedCountingHandler console;
    PerFileLog log;
    QVERIFY2(log.open(path), qPrintable(log.errorString()));
    qInfo().noquote() << "suppressed-progress";
    qWarning().noquote() << "visible-warning";
    QVERIFY2(log.close(), qPrintable(log.errorString()));

    QCOMPARE(readLog(path), QByteArray("visible-warning\n"));
    QCOMPARE(forwardedCount.load(), 2);
}

int main(int argc, char *argv[]) {
    // Isolate these formatting checks from the caller's console configuration.
    qunsetenv("QT_MESSAGE_PATTERN");
    QCoreApplication application(argc, argv);
    RadiantPerFileLogTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "PerFileLogTests.moc"
