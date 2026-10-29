#ifndef RADIANT_PER_FILE_LOG_H
#define RADIANT_PER_FILE_LOG_H

#include <QFile>
#include <QtGlobal>

// A scoped copy of Radiant's Qt diagnostics. Console output remains available.
class PerFileLog
{
public:
    PerFileLog() = default;
    ~PerFileLog();
    PerFileLog(const PerFileLog &) = delete;
    PerFileLog &operator=(const PerFileLog &) = delete;

    bool open(const QString &path);
    bool close();
    QString errorString() const;

private:
    static void messageHandler(QtMsgType type, const QMessageLogContext &context,
                               const QString &message);

    QFile m_file;
    QString m_error;
};

#endif // RADIANT_PER_FILE_LOG_H
