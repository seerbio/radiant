//
// Created by anichols on 8/24/22.
//

#include "MsReaderPointerAcc.h"

#include "GlobalSettings.h"
#include "MsReaderParquet.h"
#include "MsReaderMzML.h"
#include "StringUtils.h"

#include <QFileInfo>


Err MsReaderPointerAcc::openFile(const QString &filePath) {

    ERR_INIT

    qDebug() << "Open File 2";

    const double defaultIsolationWindow = 1.0;
    const int defaultCollisionEnergy = 30;

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

    else if (StringUtils::stringsMatch(fileSuffix, S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION, false) && fi.isFile()) {

        QSharedPointer<MsReaderBase> msReader(new MsReaderParquet);
        ptr = msReader;
        e = ptr->openFile(filePath); ree;
    }

    else {
        rrr(eFileIncorrectTypeError);
    }

    const QString msReaderType = typeid(*ptr).name();
    const bool isMsReaderBase = msReaderType.contains(QStringLiteral("MsReaderBase"));
    qDebug() << "MsReader Derived Type" << msReaderType << isMsReaderBase;

    ERR_RETURN
}
