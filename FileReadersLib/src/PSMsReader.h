//
// Created by anichols on 4/26/23.
//

#ifndef PYTHIADIACPP_PSMSREADER_H
#define PYTHIADIACPP_PSMSREADER_H

#include "Error.h"
#include "ErrorUtils.h"
#include "FileReadersLib_Exports.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"

#include <QMap>
#include <vector>

using namespace Error;


struct FILEREADERSLIB_EXPORTS PSMsReaderRow : public ParquetReaderInputBase {

    FrameIndex frameIndex = -1;
    ScanNumber scanNumber = -1;
    Charge charge = -1;
    UniqueMsInfoScanKey uniqueMsInfoScanKey;

    PeptideStringWithMods peptideStringWithMods;

    Score score = -1.0;
    int frameRankScore = -1;

    DiscScore discScore = -1.0;
    int frameRankDiscScore = -1;

    bool isDecoy = false;

    QVector<double> mzFound;
    QVector<double> mzTheo;
    QVector<double> intensityFound;
    QVector<double> intensityTheo;

};



class FILEREADERSLIB_EXPORTS PSMsReader {


public:

    PSMsReader() = default;
    ~PSMsReader() = default;





};


#endif //PYTHIADIACPP_PSMSREADER_H
