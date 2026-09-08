//
// Created by anichols on 8/24/22.
//

#include "MsReaderPointerAcc.h"

#include "GlobalSettings.h"
#include "MsReaderParquet.h"
#include "MsReaderTimsreader.h"
#include "MsReaderMzMLLazyLoad.h"
#include "MsReaderMzMLMapped.h"
#include "StringUtils.h"

#include <QByteArray>
#include <QCryptographicHash>
#include <QDataStream>
#include <QDir>
#include <QFile>
#include <QFileInfo>

namespace {

    bool isBrukerDirectoryInputPath(const QString &filePath) {
        const QString normalizedPath = QDir::cleanPath(filePath);
        const QFileInfo fileInfo(normalizedPath);
        if (!fileInfo.isDir()) {
            return false;
        }

        return StringUtils::stringsMatch(
            fileInfo.suffix(),
            S_GLOBAL_SETTINGS.BRUKER_FILE_EXTENSION,
            false
            );
    }

    Err openBrukerDirectoryReader(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> *filterRange,
        MsReaderPointerAcc *msReaderPointerAcc
        ) {

        ERR_INIT

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Cannot use lazy loading w/ timsreader Bruker inputs";
        msReaderPointerAcc->setUseLazyLoading(false);

        QSharedPointer<MsReaderBase> msReader(new MsReaderTimsreader(msReaderPointerAcc->imHandlingMode()));
        msReaderPointerAcc->ptr = msReader;

        if (filterRange == nullptr) {
            e = msReaderPointerAcc->ptr->openFile(filePath); ree;
        }
        else {
            e = msReaderPointerAcc->ptr->openFile(filePath, columnToFilterBy, *filterRange); ree;
        }
        ERR_RETURN
    }

}//namespace

MsReaderPointerAcc::MsReaderPointerAcc()
    : m_useLazyLoading(false)
    , m_imHandlingMode(ImHandlingMode::Centroid) {}

void MsReaderPointerAcc::setUseLazyLoading(bool useLazyLoading) {
    m_useLazyLoading = useLazyLoading;
}

bool MsReaderPointerAcc::useLazyLoading() const {
    return m_useLazyLoading;
}

void MsReaderPointerAcc::setImHandlingMode(ImHandlingMode imHandlingMode) {
    m_imHandlingMode = imHandlingMode;
}

ImHandlingMode MsReaderPointerAcc::imHandlingMode() const {
    return m_imHandlingMode;
}

Err MsReaderPointerAcc::openFile(const QString &filePath) {
    ERR_INIT
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Reading MsFile";
    e = setMsReaderPointer(filePath); ree;
    ERR_RETURN
}

Err MsReaderPointerAcc::setMsReaderPointer(const QString &filePath) {

    ERR_INIT

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION, false) && fi.isFile()) {
        if (m_useLazyLoading) {
            QSharedPointer<MsReaderBase> msReader(new MsReaderMzMLLazyLoad);
            ptr = msReader;
            e = ptr->openFile(filePath); ree;
        }
        else {
            QSharedPointer<MsReaderBase> msReader(new MsReaderMzMLMapped);
            ptr = msReader;
            e = ptr->openFile(filePath); ree;
        }

    }

    else if (
            (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION, false)
                || StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.CACHED_FILE_EXTENSION, false))
            && fi.isFile()) {

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Cannot use lazy loading w/ .prq files";
        m_useLazyLoading = false;

        QSharedPointer<MsReaderBase> msReader(new MsReaderParquet);
        ptr = msReader;
        e = ptr->openFile(filePath); ree;
    }

    else if (isBrukerDirectoryInputPath(filePath)) {
        e = openBrukerDirectoryReader(filePath, QString(), nullptr, this); ree;
    }

    else {
        qDebug() << "Filepath" << filePath;
        qDebug() << "Suffix" << fileSuffix;
        rrr(eFileIncorrectTypeError);
    }

    const QString msReaderType = typeid(*ptr).name();
    const bool isMsReaderBase = msReaderType.contains(QStringLiteral("MsReaderBase"));
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsReader Derived Type" << msReaderType << isMsReaderBase;

    ERR_RETURN
}

Err MsReaderPointerAcc::openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
        ) {

    ERR_INIT

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Reading MsFile";

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION, false) && fi.isFile()) {

        QSharedPointer<MsReaderBase> msReader(new MsReaderMzMLMapped);
        ptr = msReader;
        e = ptr->openFile(filePath, columnToFilterBy, filterRange); ree;
    }

    else if (
            (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION, false)
             || StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.CACHED_FILE_EXTENSION, false))
            && fi.isFile()) {

        QSharedPointer<MsReaderBase> msReader(new MsReaderParquet);
        ptr = msReader;
        e = ptr->openFile(filePath, columnToFilterBy, filterRange); ree;
    }

    else if (isBrukerDirectoryInputPath(filePath)) {
        e = openBrukerDirectoryReader(filePath, columnToFilterBy, &filterRange, this); ree;
    }

    else {
        qDebug() << "Filepath" << filePath;
        qDebug() << "Suffix" << fileSuffix;
        rrr(eFileIncorrectTypeError);
    }

    const QString msReaderType = typeid(*ptr).name();
    const bool isMsReaderBase = msReaderType.contains(QStringLiteral("MsReaderBase"));
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsReader Derived Type" << msReaderType << isMsReaderBase;

    ERR_RETURN
}

Err MsReaderPointerAcc::openFile(const QString &filePath, const QString &columnToFilterBy) {
    return eFunctionNotImplemented;
}

bool MsReaderPointerAcc::isInit() {

    if (ptr) {
        return ptr->isInit();
    }

    return false;
}
