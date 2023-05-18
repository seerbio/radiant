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

    qDebug() << "drewholio" << filePath << expectedFileExtension << fi.isFile() << fi.isDir();

    if (fi.isFile()) {
        return StringUtils::stringsMatch(fileSuffix, expectedFileExtension, false);
    }

    else if (fi.isDir()) {
        return true;
    }

    return false;
}