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

    /**
    * @brief Checks if the provided file path has the expected file extension.
    *
    * This method is part of the CommandLineParserUtils class and is used to ensure whether
    * a certain file has the desired extension or not.
    *
    * @param filePath The full path to the file to check.
    * @param expectedFileExtension The extension that the file is expected to have, not including the leading dot.
    *
    * @return Returns true if the file has the expected extension, false otherwise.
    */
    static bool checkFileNameExtension(
            const QString &filePath,
            const QString &expectedFileExtension
    );

    /**
    * @brief Retrieve data files with a specific extension from a directory.
    *
    * This is a method of CommandLineParserUtils class that is intended to fetch all files
    * with a particular extension from a specified directory.
    * It uses QDirIterator to recursively retrieve the files.
    * Files with matching extensions are added to the dataFiles QStringList.
    *
    * @param dataFilesDirectory A QString representing the directory from which the data files should be fetched.
    * @param dataFiles A pointer to a QStringList. This is where the paths of the data files,
    * with matching extensions, fetched from the directory will be stored.
    *
    * @return Returns an Err object which indicates the success or failure of the operation. If the operation is successful,
    * an Err object initialized with a success state is returned.
    */
    static Err getDataFilesFromDirectory(
            const QString &dataFilesDirectory,
            QStringList *dataFiles
    );

};


#endif //PYTHIADIACPP_COMMANDLINEPARSERUTILS_H
