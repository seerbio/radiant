//
// Created by anichols on 5/17/23.
//

#ifndef PYTHIADIACPP_COMMANDLINEPARSERUTILS_H
#define PYTHIADIACPP_COMMANDLINEPARSERUTILS_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"


using namespace Error;


class FILEREADERSLIB_EXPORTS CommandLineParserUtils {

public:

    static bool checkFileNameExtension(
            const QString &filePath,
            const QString &expectedFileExtension
    );

    static Err getDataFilesFromDirectory(
            const QString &dataFilesDirectory,
            QStringList *dataFiles
    );


};


#endif //PYTHIADIACPP_COMMANDLINEPARSERUTILS_H
