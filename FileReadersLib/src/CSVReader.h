//
// Created by anichols on 4/4/23.
//

#ifndef PYTHIADIACPP_CSVREADER_H
#define PYTHIADIACPP_CSVREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"


using namespace Error;


struct CSVReaderBase {

public:

    virtual QMap<QString, QVariant> map() = 0;

};


class FILEREADERSLIB_EXPORTS CSVReader {


public:

    static Err writeDataToCSV(
            const QString &outputFilePath,
            const QVector<QSharedPointer<CSVReaderBase>> &rowsToWrite
            );



};


#endif //PYTHIADIACPP_CSVREADER_H
