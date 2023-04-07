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

    virtual QMap<QString, QVariant> map() {return {};};

    template<typename T>
    static QVector<QSharedPointer<CSVReaderBase>> convertInputStructToSharedPointers(
            const QVector<T> &csvReaderInputBaseDerivedStruct
            ) {
        QVector<QSharedPointer<CSVReaderBase>> ptrs;
        for (const T &tr : csvReaderInputBaseDerivedStruct) {
            QSharedPointer<CSVReaderBase> ptr(new T(tr));
            ptrs.push_back(ptr);
        }

        return ptrs;
    }

};


class FILEREADERSLIB_EXPORTS CSVReader {


public:

    static Err writeDataToCSV(
            const QString &outputFilePath,
            const QVector<QSharedPointer<CSVReaderBase>> &rowsToWrite
            );

    static Err readDataFromCSV(
            const QString &outputFilePath,
            QVector<CSVReaderBase> &readRows
            );



};


#endif //PYTHIADIACPP_CSVREADER_H
