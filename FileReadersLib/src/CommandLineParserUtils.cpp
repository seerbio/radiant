//
// Created by anichols on 5/17/23.
//

#include "CommandLineParserUtils.h"

#include "StringUtils.h"

#include <QFileInfo>
#include <QDirIterator>


bool CommandLineParserUtils::checkFileNameExtensions(
        const QString &filePath,
        const QStringList &expectedFileExtensions
) {

    QFileInfo fi(filePath);
    const QString fileSuffix = fi.suffix();
    qDebug() << "Expected extensions:" << expectedFileExtensions << "Got:" << fileSuffix;
    return expectedFileExtensions.contains(fileSuffix);
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