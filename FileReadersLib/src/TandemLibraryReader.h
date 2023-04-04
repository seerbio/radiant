//
// Created by anichols on 4/3/23.
//

#ifndef PYTHIADIACPP_TANDEMLIBRARYREADER_H
#define PYTHIADIACPP_TANDEMLIBRARYREADER_H

#include "Error.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReaderBase.h"

using namespace Error;


struct FILEREADERSLIB_EXPORTS TandemLibraryReaderRow {
    PeptideString peptideString;
    std::vector<double> intensityVals;
    std::vector<std::string> ionLabels;
};


class FILEREADERSLIB_EXPORTS TandemLibraryReader : public ParquetReaderBase {

public:

    TandemLibraryReader() = default;
    ~TandemLibraryReader() = default;

    Err writeTandemPrediction(
            const QVector<TandemLibraryReaderRow> &tandemLibraryReaderRows,
            const QString &outputFilePath
            );





private:


};


#endif //PYTHIADIACPP_TANDEMLIBRARYREADER_H
