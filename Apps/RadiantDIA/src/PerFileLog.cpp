#include "PerFileLog.h"

#include <QMutex>
#include <QMutexLocker>

#include <cstdio>

namespace {
QMutex logMutex;
PerFileLog *activeLog = nullptr;
QtMessageHandler previousHandler = nullptr;
}

PerFileLog::~PerFileLog() {
    close();
}

bool PerFileLog::open(const QString &path) {
    QMutexLocker lock(&logMutex);
    if (activeLog != nullptr) {
        m_error = QStringLiteral("Another per-file log is already open");
        return false;
    }
    m_error.clear();
    m_file.setFileName(path);
    if (!m_file.open(QIODevice::WriteOnly | QIODevice::Truncate)) {
        m_error = m_file.errorString();
        return false;
    }
    activeLog = this;
    previousHandler = qInstallMessageHandler(messageHandler);
    return true;
}

bool PerFileLog::close() {
    QMutexLocker lock(&logMutex);
    if (activeLog == this) {
        qInstallMessageHandler(previousHandler);
        activeLog = nullptr;
        if (!m_file.flush() && m_error.isEmpty()) {
            m_error = m_file.errorString();
        }
        m_file.close();
    }
    return m_error.isEmpty();
}

QString PerFileLog::errorString() const {
    QMutexLocker lock(&logMutex);
    return m_error;
}

void PerFileLog::messageHandler(QtMsgType type, const QMessageLogContext &context,
                               const QString &message) {
    const QString formatted = qFormatLogMessage(type, context, message);
    const QByteArray bytes = formatted.isEmpty()
        ? QByteArray{} : formatted.toUtf8() + '\n';
    QtMessageHandler consoleHandler;
    {
        QMutexLocker lock(&logMutex);
        if (activeLog != nullptr && !bytes.isEmpty()) {
            if (activeLog->m_file.write(bytes) != bytes.size()
                || !activeLog->m_file.flush()) {
                // Logging here would recursively invoke this handler.
                if (activeLog->m_error.isEmpty()) {
                    activeLog->m_error = activeLog->m_file.errorString();
                }
            }
        }
        consoleHandler = previousHandler;
    }
    // A previous handler may itself log, so call it without holding our mutex.
    if (consoleHandler != nullptr) {
        consoleHandler(type, context, message);
    } else if (!bytes.isEmpty()) {
        std::fwrite(bytes.constData(), 1, static_cast<size_t>(bytes.size()), stderr);
    }
}
