//
// Created by anichols on 4/3/23.
//

#include "DIAMzTargetsReader.h"

Err DIAMzTargetsReader::write(
        const QVector<DIAMzTargetsReaderRow> &diaMzTargetsReaderRow,
        const QString &outputFilePath
        ) {

    ERR_INIT

    const QVector<QSharedPointer<CSVReaderInputBase>> ptrs
            = CSVReaderInputBase::convertInputStructToSharedPointers(diaMzTargetsReaderRow);

    CSVReader reader;
    e = reader.writeDataToCSV(
            outputFilePath,
            ptrs
            ); ree;

    ERR_RETURN
}

Err DIAMzTargetsReader::read(
        const QString &fileURI,
        QVector<DIAMzTargetsReaderRow> *diaMzTargetsReaderRow
        ) {

    ERR_INIT

    CSVReader reader;

    QVector<CSVReaderInputBase> ptrsRead;
    e = reader.readDataFromCSV(
            fileURI,
            &ptrsRead
    ); ree;

    e = CSVReaderInputBase::convertSharedPointersToInputStruct(
            ptrsRead,
            diaMzTargetsReaderRow
            ); ree;

    ERR_RETURN
}
