//
// Created by anichols on 4/7/23.
//

#include "ConvertMzMLToParquetWorkFlow.h"

#include "GlobalSettings.h"

#include <QFileInfo>

Err ConvertMzMLToParquetWorkFlow::convertMzMLToParquet(
        const QString &msDataFilePath,
        QString *outputFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;

    MsReaderPointerAcc msReaderPointerAcc;
    msReaderPointerAcc.setUseLazyLoading(false);
    e = msReaderPointerAcc.openFile(msDataFilePath); ree;

    const QFileInfo fileInfo(msDataFilePath);
    *outputFilePath = fileInfo.absoluteFilePath() + S_GLOBAL_SETTINGS.DOT_PRQ_FF;

    e = MsReaderParquet::writeMsReaderToParquet(
            *outputFilePath,
            msReaderPointerAcc.ptr
            ); ree;

    ERR_RETURN
}
