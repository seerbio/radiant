//
// Created by anichols on 5/17/23.
//

#include "CommandLineParserUtils.h"

#include "MsReaderTimsbukIndex.h"
#include "StringUtils.h"

#include <QDir>
#include <QFileInfo>
#include <QDirIterator>


bool CommandLineParserUtils::checkFileNameExtensions(
        const QString &filePath,
        const QStringList &expectedFileExtensions
) {

    const QString normalizedPath = QDir::cleanPath(filePath);
    QFileInfo fi(normalizedPath);
    const QString fileSuffix = fi.suffix().toLower();

    return std::any_of(
            expectedFileExtensions.begin(),
            expectedFileExtensions.end(),
            [fileSuffix](const QString &s){return s.toLower() == fileSuffix;}
            );
}

bool CommandLineParserUtils::isMassSpectrometryDataPath(const QString &filePath) {

    return checkFileNameExtensions(
               filePath,
               {
                   S_GLOBAL_SETTINGS.PRQ_FILE_EXTENSION,
                   S_GLOBAL_SETTINGS.MZML_FILE_EXTENSION,
                   S_GLOBAL_SETTINGS.BRUKER_FILE_EXTENSION,
                   QStringLiteral("idx")
               }
               )
        || MsReaderTimsbukIndex::isDirectIndexRootPath(filePath);
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
