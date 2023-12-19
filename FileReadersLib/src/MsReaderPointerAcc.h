//
// Created by anichols on 8/24/22.
//

#ifndef PYTHIACPP_MSREADERPOINTERACC_H
#define PYTHIACPP_MSREADERPOINTERACC_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "MsReaderBase.h"

#include <QSharedPointer>


using namespace Error;

class FILEREADERSLIB_EXPORTS MsReaderPointerAcc {

public:

    MsReaderPointerAcc() = default;

    ~MsReaderPointerAcc() = default;

    Err openFile(const QString &filePath);

    Err openFileWithCache(const QString &filePath);

    bool usingCache();

    static Err readFrameCache(
            const QString &cachedFilePath,
            QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints,
            QMap<ScanNumber, ScanTime> *scanNumberVsScanTime
            );

    QMap<MzTargetKey, FilePath> mzTargetKeyVsFilePathCache();

    QSharedPointer<MsReaderBase> ptr;


private:

    Err setMsReaderPointer(const QString &filePath);

private:

    QString m_cachedFilePath;
    QMap<MzTargetKey, FilePath> m_mzTargetKeyVsFilePathCache;

};


#endif //PYTHIACPP_MSREADERPOINTERACC_H
