//
// Created by Codex on 7/13/26.
//

#ifndef MSREADERTIMSBUKINDEX_H
#define MSREADERTIMSBUKINDEX_H

#include "FileReadersLib_Exports.h"

#include "MsReaderBase.h"

using namespace Error;

class FILEREADERSLIB_EXPORTS MsReaderTimsbukIndex : public MsReaderBase {

public:

    MsReaderTimsbukIndex() = default;
    ~MsReaderTimsbukIndex() override = default;

    Err openFile(const QString &filePath) override;

    Err openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
        ) override;

    Err closeFile() override;
};

#endif // MSREADERTIMSBUKINDEX_H
