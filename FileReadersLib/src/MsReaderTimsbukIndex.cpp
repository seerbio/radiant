//
// Created by Codex on 7/13/26.
//

#include "MsReaderTimsbukIndex.h"

#include "ErrorUtils.h"

#include <QDebug>
#include <QFileInfo>


namespace {

    Err validateInputPath(const QString &filePath) {

        ERR_INIT

        e = ErrorUtils::fileExists(filePath); ree;

        const QFileInfo fileInfo(filePath);
        e = ErrorUtils::isTrue(fileInfo.isDir(), eFileIncorrectTypeError); ree;

        ERR_RETURN
    }

}//namespace

Err MsReaderTimsbukIndex::openFile(const QString &filePath) {

    ERR_INIT

    e = closeFile(); ree;
    e = validateInputPath(filePath); ree;

    m_filePath = filePath;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "MsReaderTimsbukIndex sidecar reader path active"
             << filePath;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "MsReaderTimsbukIndex parsing not implemented yet";

    rrr(eFunctionNotImplemented);
}

Err MsReaderTimsbukIndex::openFile(
    const QString &filePath,
    const QString &columnToFilterBy,
    const QPair<double, double> &filterRange
    ) {

    Q_UNUSED(columnToFilterBy)
    Q_UNUSED(filterRange)

    return openFile(filePath);
}

Err MsReaderTimsbukIndex::closeFile() {
    return MsReaderBase::closeFile();
}
