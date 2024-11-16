//
// Created by anichols on 4/3/23.
//

#ifndef PYTHIADIACPP_FRAGLIBREADER_H
#define PYTHIADIACPP_FRAGLIBREADER_H

#include <utility>

#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "FragLibReaderRow.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

using namespace Error;

class FragLibReaderRow;


class FILEREADERSLIB_EXPORTS FragLibReader {

public:


    static Err getFragLibReaderRows(
            const QString &fragLibFilePath,
            QVector<FragLibReaderRow> *fragLibReaderRows
    );

};




#endif //PYTHIADIACPP_FRAGLIBREADER_H
