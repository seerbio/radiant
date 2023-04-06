//
// Created by anichols on 4/3/23.
//

#include "TandemLibraryReader.h"

#include "ParquetReader.h"

Err TandemLibraryReader::writeTandemPredictions(
        const QVector<TandemLibraryReaderRow> &tandemLibraryReaderRows,
        const QString &outputFilePath
        ) {

    ERR_INIT

    const QVector<QSharedPointer<ParquetReaderInputBase>> ptrs
            = ParquetReaderInputBase::convertInputStructToSharedPointers(tandemLibraryReaderRows);

    ParquetReader reader;
    e = reader.writeDataToParquet(
            outputFilePath,
            ptrs
            ); ree;

    ERR_RETURN
}

Err TandemLibraryReader::readTandemPredictions(
        const QString &fileURI,
        QVector<TandemLibraryReaderRow> *tandemLibraryReaderRows
        ) {

    ERR_INIT

    ParquetReader reader;

    QVector<ParquetReaderInputBase> ptrsRead;
    e = reader.readDataFromParquet(
            fileURI,
            &ptrsRead
    ); ree;


    e = ParquetReaderInputBase::convertSharedPointersToInputStruct(
            ptrsRead,
            tandemLibraryReaderRows
            ); ree;

    ERR_RETURN
}
