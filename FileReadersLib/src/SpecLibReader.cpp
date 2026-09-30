//
// Created by andrewnichols on 11/15/24.
//

#include "SpecLibReader.h"

#include "ErrorUtils.h"
#include "FragLibReaderRow.h"

#include "SpecLibSrc/Library.h"

Err SpecLibReader::getFragLibReaerRows(
    const QString &fragLibFilePath,
    QList<FragLibReaderRow> *fragLibReaderRows
    ) {

    ERR_INIT

    e = ErrorUtils::fileExists(fragLibFilePath); ree;

    Library library;
    MappedFileInput mappedInput(fragLibFilePath);
    e = ErrorUtils::isTrue(mappedInput.isOpen(), eFileError); ree;
    e = library.read(mappedInput, fragLibReaderRows); ree;

    ERR_RETURN
}
