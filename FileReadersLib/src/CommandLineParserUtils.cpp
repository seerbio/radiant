//
// Created by anichols on 5/17/23.
//

#include "CommandLineParserUtils.h"

#include "StringUtils.h"

#include <QFileInfo>
#include <QDirIterator>


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

Err CommandLineParserUtils::getDataFilesFromDirectory(
        const QString &dataFilesDirectory,
        QStringList *dataFiles
) {

    ERR_INIT

    QDirIterator it(dataFilesDirectory);
    while (it.hasNext()) {

        const QString &dataFilePath = it.next();
        const QFileInfo fi(dataFilePath);
        const QString fileExtension = fi.suffix();

        const QString prqSuffix = S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION;

        if (StringUtils::stringsMatch(prqSuffix, fileExtension, false)) {
            dataFiles->push_back(dataFilePath);
        }
    }

    ERR_RETURN
}