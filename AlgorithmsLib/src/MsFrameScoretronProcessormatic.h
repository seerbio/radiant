//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
#define PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H


#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"


using namespace Error;


class ALGORITHMSLIB_EXPORTS MsFrameScoretronProcessormatic {

public:

    MsFrameScoretronProcessormatic() = default;
    ~MsFrameScoretronProcessormatic() = default;

    Err init(const PythiaParameters &pythiaParameters);

private:

    PythiaParameters m_params;

};


#endif //PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
