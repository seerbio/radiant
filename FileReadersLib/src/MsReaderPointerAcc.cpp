//
// Created by anichols on 8/24/22.
//

#include "MsReaderPointerAcc.h"

#include "GlobalSettings.h"
#include "MsReaderParquet.h"
#include "MsReaderMzML.h"
#include "StringUtils.h"

#include <QByteArray>
#include <QCryptographicHash>
#include <QDataStream>
#include <QFile>
#include <QFileInfo>


Err MsReaderPointerAcc::openFile(const QString &filePath) {
    ERR_INIT
    e = setMsReaderPointer(filePath); ree;
    ERR_RETURN
}

Err MsReaderPointerAcc::setMsReaderPointer(const QString &filePath) {

    ERR_INIT

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION, false) && fi.isFile()) {
        QSharedPointer<MsReaderBase> msReader(new MsReaderMzML);
        ptr = msReader;
        e = ptr->openFile(filePath); ree;
    }

    else if (
            (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION, false)
                || StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.CACHED_FILE_EXTENSION, false))
            && fi.isFile()) {

        QSharedPointer<MsReaderBase> msReader(new MsReaderParquet);
        ptr = msReader;
        e = ptr->openFile(filePath); ree;
    }

    else {
        qDebug() << "Filepath" << filePath;
        qDebug() << "Suffix" << fileSuffix;
        rrr(eFileIncorrectTypeError);
    }

    const QString msReaderType = typeid(*ptr).name();
    const bool isMsReaderBase = msReaderType.contains(QStringLiteral("MsReaderBase"));
    qDebug() << "MsReader Derived Type" << msReaderType << isMsReaderBase;

    ERR_RETURN
}

namespace {

    QPair<Err, QString> buildFrameCacheFileLogic(
            const FilePath &filePath,
            const MzTargetKey &mzTargetKey,
            const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPointsPntrs,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scanNumberVsScanPointsPntrs); rree;

        const QString filePathTargetKey = filePath + mzTargetKey;
        const QByteArray hash = QCryptographicHash::hash(filePath.toUtf8(), QCryptographicHash::Sha256);
        const QString hashHex = hash.toHex().toUpper() + ".pythiaCache";

        QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
        for (auto it = scanNumberVsScanPointsPntrs.begin(); it != scanNumberVsScanPointsPntrs.end(); it++) {
            scanNumberVsScanPoints.insert(it.key(), *it.value());
        }

        QFile file(hashHex);
        if(file.open(QIODevice::WriteOnly)){
            QDataStream out(&file);
            out.setVersion(QDataStream::Qt_5_12);
            out << scanNumberVsScanTime;
            out << scanNumberVsScanPoints;
        }
        else {
            return {eFileError, {}};
        }

        file.close();

        return {e, hashHex};
    }

}//namespace
Err MsReaderPointerAcc::openFileWithCache(const QString &filePath) {

    ERR_INIT

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    QSharedPointer<MsReaderBase> localPointer;

    if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION, false) && fi.isFile()) {
        QSharedPointer<MsReaderBase> msReader(new MsReaderMzML);
        ptr = msReader;
        e = localPointer->openFile(filePath); ree;
    }

    else if (
            (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION, false)
             || StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.CACHED_FILE_EXTENSION, false))
            && fi.isFile()) {

        QSharedPointer<MsReaderBase> msReader(new MsReaderParquet);
        ptr = msReader;
        e = localPointer->openFile(filePath); ree;
    }

    else {
        qDebug() << "Filepath" << filePath;
        qDebug() << "Suffix" << fileSuffix;
        rrr(eFileIncorrectTypeError);
    }

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    e = localPointer->collateMS2MzTargetFrames(
            &diaTargetFrames
            ); ree;

    for (auto it = diaTargetFrames.begin(); it != diaTargetFrames.end(); it++) {

        const MzTargetKey &mzTargetKey = it.key();
        const QMap<ScanNumber, ScanPoints*> &mzTargetVsScanPointPntr = it.value();

        const QPair<Err, FilePath> result = buildFrameCacheFileLogic(
                filePath,
                mzTargetKey,
                mzTargetVsScanPointPntr,
                localPointer->getScanNumberVsScanTime()
                ); ree;

        e = result.first; ree;
        m_mzTargetKeyVsFilePathCache.insert(mzTargetKey, result.second);
    }

    const QString msReaderType = typeid(*localPointer).name();
    const bool isMsReaderBase = msReaderType.contains(QStringLiteral("MsReaderBase"));
    qDebug() << "MsReader Derived Type" << msReaderType << isMsReaderBase;

    ERR_RETURN
}

Err MsReaderPointerAcc::readFrameCache(
        const QString &cachedFilePath,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints,
        QMap<ScanNumber, ScanTime> *scanNumberVsScanTime
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(cachedFilePath); ree;

    QFile file(cachedFilePath);
    if(file.open(QIODevice::ReadOnly)) {
        QDataStream in(&file);
        in.setVersion(QDataStream::Qt_5_12);
        in >> *scanNumberVsScanTime;
        in >> *scanNumberVsScanPoints;
    }
    else {
        rrr(eFileError);
    }

    file.close();

    ERR_RETURN
}

