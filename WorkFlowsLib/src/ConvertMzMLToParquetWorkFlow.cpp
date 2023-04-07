//
// Created by anichols on 4/7/23.
//

#include "ConvertMzMLToParquetWorkFlow.h"

#include "GlobalSettings.h"

Err ConvertMzMLToParquetWorkFlow::convertMzMLToParquet(const QString &mzmlFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(mzmlFilePath); ree;

    MsReaderMzML msReaderMzMl;
    e = msReaderMzMl.openFile(mzmlFilePath); ree;

    const QString prqFileName = mzmlFilePath + S_GLOBAL_SETTINGS.DOT_PRQ;

    e = MsReaderParquet::writeMsReaderToParquet(
            prqFileName,
            QSharedPointer<MsReaderBase>(new MsReaderBase(msReaderMzMl.msReaderBase()))
            ); ree;

    ERR_RETURN
}
