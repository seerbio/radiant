// Created by anichols on 11/07/2021.

#include "ErrorUtils.h"

#include <QFileInfo>

Err ErrorUtils::fileExists(const QString &filePath) {

    QFileInfo fi(filePath);

    ERR_INIT

    e = ErrorUtils::isNotEmpty(filePath); ree;

    if (!fi.exists() || !fi.isFile()) {
        qDebug() << filePath << "not found";
        rrr(eFileError);
    }

    return eNoError;
}
