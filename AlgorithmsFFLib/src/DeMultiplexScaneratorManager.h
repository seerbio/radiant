//
// Created by anichols on 2/9/24.
//

#ifndef PYTHIADIACPP_DEMULTIPLEXSCANERATORMANAGER_H
#define PYTHIADIACPP_DEMULTIPLEXSCANERATORMANAGER_H


#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"

using namespace Error;

class MsReaderPointerAcc;

class ALGORITHMSFFLIB_EXPORTS DeMultiplexScaneratorManager {

public:

    DeMultiplexScaneratorManager() = default;
    ~DeMultiplexScaneratorManager() = default;

    Err demultiplexMsReader(MsReaderPointerAcc *msReaderPointerAcc);
};


#endif //PYTHIADIACPP_DEMULTIPLEXSCANERATORMANAGER_H
