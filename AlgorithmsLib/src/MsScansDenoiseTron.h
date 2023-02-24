//
// Created by anichols on 2/23/23.
//

#ifndef PYTHIADIACPP_MSSCANSDENOISETRON_H
#define PYTHIADIACPP_MSSCANSDENOISETRON_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"

using namespace Error;


class ALGORITHMSLIB_EXPORTS MsScansDenoiseTron {

public:

    MsScansDenoiseTron() = default;
    ~MsScansDenoiseTron() = default;

    Err init(const PythiaParameters &pythiaParameter);

    Err denoiseScansFrame(const QMap<ScanNumber, ScanPoints> &scansFrame);

private:

    PythiaParameters m_pythiaParameters;


};


#endif //PYTHIADIACPP_MSSCANSDENOISETRON_H
