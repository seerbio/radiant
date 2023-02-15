// Created by anichols on 11/07/2021.

#include "ErrorUtils.h"

#include <QFileInfo>

Err ErrorUtils::fileExists(const QString &filePath, Err e) {

    QFileInfo fi(filePath);

    if (!fi.exists() || !fi.isFile()) {
        qDebug() << filePath << "not found";
        rrr(e);
    }

    return eNoError;
}
