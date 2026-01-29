//
// Created by andrewnichols on 11/15/24.
//

#ifndef SPECLIBREADER_H
#define SPECLIBREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"

using namespace Error;

class FragLibReaderRow;

class FILEREADERSLIB_EXPORTS SpecLibReader {

public:

    SpecLibReader() = default;
    ~SpecLibReader() = default;

    static Err getFragLibReaerRows(
        const QString &fragLibFilePath,
        QList<FragLibReaderRow> *fragLibReaderRows
        );

};



#endif //SPECLIBREADER_H
