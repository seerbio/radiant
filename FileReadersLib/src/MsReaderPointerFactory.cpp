//
// Created by anichols on 8/24/22.
//

#include "MsReaderPointerFactory.h"

#include "GlobalSettings.h"
#include "MsReaderMzML.h"
#include "MsReaderParquet.h"
#include "StringUtils.h"

#include <QFileInfo>


QPair<Err, MsReaderPointer> MsReaderPointerFactory::createInstance(const QString &filePath) {

    ERR_INIT

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION, false) &&
        fi.isFile()) {

        QSharedPointer<MsReaderBase> msReader(new MsReaderMzML);
        e = msReader->openFile(filePath);
        if (e != eNoError) {
            return {eFileIncorrectTypeError, nullptr};
        }

        const QString msReaderType = typeid(msReader).name();
        qDebug() << "MsReader Derived Type" << msReaderType;

        return {e, msReader};

    }

    else if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION, false)
            && fi.isFile()) {

        QSharedPointer<MsReaderBase> msReader(new MsReaderParquet);
        e = msReader->openFile(filePath);
        if (e != eNoError) {
            return {eFileIncorrectTypeError, nullptr};
        }

        const QString msReaderType = typeid(msReader).name();
        qDebug() << "MsReader Derived Type" << msReaderType;
        return {e, msReader};
    }

    else {
        return {eFileIncorrectTypeError, nullptr};
    }
}
