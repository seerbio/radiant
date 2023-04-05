//
// Created by anichols on 4/5/23.
//

#ifndef PYTHIADIACPP_PARQUETREADER_H
#define PYTHIADIACPP_PARQUETREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "SqlUtils.h"

#include <QDebug>

using namespace Error;

struct FILEREADERSLIB_EXPORTS ParquetReaderInputBase {

public:

    virtual QMap<QString, QVariant> map() = 0;

};

class ParquetReader {

public:

    ParquetReader();
    ~ParquetReader();

    Err writeDataToParquet(
            const QString &outputFilePath,
            const QVector<QSharedPointer<ParquetReaderInputBase>> &rowsToWrite
            );



private:

    Q_DISABLE_COPY(ParquetReader) class Private;
    const QScopedPointer<Private> d_ptr;




};


#endif //PYTHIADIACPP_PARQUETREADER_H
