//
// Created by anichols on 5/17/23.
//

#include "CommandLineParserUtils.h"

#include "StringUtils.h"

#include <QFileInfo>



bool CommandLineParserUtils::checkFileNameExtension(
        const QString &filePath,
        const QString &expectedFileExtension
) {

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();

    qDebug() << "drewholio" << fileSuffix << expectedFileExtension;
    qDebug() << StringUtils::stringsMatch(fileSuffix, expectedFileExtension, false);
    qDebug() << fi.isFile();
    qDebug() << fi.isDir();

    if (StringUtils::stringsMatch(fileSuffix, expectedFileExtension, false) && fi.isFile()) {
        return true;
    }

    else if (fi.isDir()) {
        return true;
    }

    return false;
}