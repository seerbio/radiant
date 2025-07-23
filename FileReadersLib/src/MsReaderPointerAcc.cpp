//
// Created by anichols on 8/24/22.
//

#include "MsReaderPointerAcc.h"

#include "GlobalSettings.h"
// #include "MsReaderBrukerTims.h"
// #include "MsReaderParquet.h"
#include "MsReaderMzMLLazyLoad.h"
// #include "MsReaderMzMLMapped.h"
#include "StringUtils.h"

#include <QByteArray>
#include <QCryptographicHash>
#include <QDataStream>
#include <QFile>
#include <QFileInfo>


MsReaderPointerAcc::MsReaderPointerAcc() : m_useLazyLoading(false) {}

void MsReaderPointerAcc::setUseLazyLoading(bool useLazyLoading) {
    m_useLazyLoading = useLazyLoading;
}

bool MsReaderPointerAcc::useLazyLoading() const {
    return m_useLazyLoading;
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
        // if (m_useLazyLoading) {
            QSharedPointer<MsReaderBase> msReader(new MsReaderMzMLLazyLoad);
            ptr = msReader;
            e = ptr->openFile(filePath); ree;
        }
        // else {
        //     QSharedPointer<MsReaderBase> msReader(new MsReaderMzMLMapped);
        //     ptr = msReader;
        //     e = ptr->openFile(filePath); ree;
        // }

    // }

    // else if (
    //         (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION, false)
    //             || StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.CACHED_FILE_EXTENSION, false))
    //         && fi.isFile()) {
    //
    //     qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Cannot use lazy loading w/ .prq files";
    //     m_useLazyLoading = false;
    //
    //     QSharedPointer<MsReaderBase> msReader(new MsReaderParquet);
    //     ptr = msReader;
    //     e = ptr->openFile(filePath); ree;
    // }
    //
    // else if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.BRUKER_FILE_EXTENSION, false) && fi.isDir()) {
    //     rrr(eFunctionNotImplemented);
    //     // // //TODO uncomment code when ARM is sorted out.
    //     //
    //     // qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Cannot use lazy loading w/ Bruker files";
    //     // m_useLazyLoading = false;
    //     //
    //     // QSharedPointer<MsReaderBase> msReader(new MsReaderBrukerTims);
    //     // ptr = msReader;
    //     // e = ptr->openFile(filePath); ree;
    // }

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

// Err MsReaderPointerAcc::openFile(
//         const QString &filePath,
//         const QString &columnToFilterBy,
//         const QPair<double, double> &filterRange
//         ) {
//
//     ERR_INIT
//
//     QFileInfo fi(filePath);
//     const QString fileSuffix = fi.suffix();
//
//     if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION, false) && fi.isFile()) {
//
//         QSharedPointer<MsReaderBase> msReader(new MsReaderMzMLMapped);
//         ptr = msReader;
//         e = ptr->openFile(filePath, columnToFilterBy, filterRange); ree;
//     }
//
//     else if (
//             (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION, false)
//              || StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.CACHED_FILE_EXTENSION, false))
//             && fi.isFile()) {
//
//         QSharedPointer<MsReaderBase> msReader(new MsReaderParquet);
//         ptr = msReader;
//         e = ptr->openFile(filePath, columnToFilterBy, filterRange); ree;
//     }
//
//     else {
//         qDebug() << "Filepath" << filePath;
//         qDebug() << "Suffix" << fileSuffix;
//         rrr(eFileIncorrectTypeError);
//     }
//
//     const QString msReaderType = typeid(*ptr).name();
//     const bool isMsReaderBase = msReaderType.contains(QStringLiteral("MsReaderBase"));
//     qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsReader Derived Type" << msReaderType << isMsReaderBase;
//
//     ERR_RETURN
// }

Err MsReaderPointerAcc::openFile(const QString &filePath, const QString &columnToFilterBy) {
    return eFunctionNotImplemented;
}

bool MsReaderPointerAcc::isInit() {

    if (ptr) {
        return ptr->isInit();
    }

    return false;
}
