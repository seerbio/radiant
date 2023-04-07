//
// Created by anichols on 2/6/23.
//

#ifndef SPARKDIA_PARQUETREADER_H
#define SPARKDIA_PARQUETREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"

#include <QMap>
#include <vector>

using namespace Error;


class FILEREADERSLIB_EXPORTS MsReaderParquet : public MsReaderBase {

public:

    MsReaderParquet() = default;
    ~MsReaderParquet() = default;

    Err openFile(const QString &filePath) override;

    Err closeFile() override;

};


#endif //SPARKDIA_PARQUETREADER_H
