//
// Created by andrewnichols on 11/15/24.
//

#include "SpecLibReader.h"

#include "ErrorUtils.h"
#include "FragLibReaderRow.h"

#include "SpecLibSrc/Library.h"

Err SpecLibReader::getFragLibReaerRows(
    const QString &fragLibFilePath,
    QVector<FragLibReaderRow> *fragLibReaderRows
    ) {

    ERR_INIT

    e = ErrorUtils::fileExists(fragLibFilePath); ree;

    std::ifstream speclibStream(fragLibFilePath.toStdString(), std::ifstream::binary);

    Library library;
    library.read(speclibStream, fragLibReaderRows);

    ERR_RETURN
}

