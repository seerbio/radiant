//
// Created by anichols on 4/7/23.
//

#ifndef PYTHIADIACPP_MSFILEPREPROCESSORWORKFLOW_H
#define PYTHIADIACPP_MSFILEPREPROCESSORWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsReaderPointerFactory.h"
#include "PythiaParameterReader.h"

using namespace Error;

class WORKFLOWSLIB_EXPORTS MsFilePreProcessorWorkFlow {

public:

    MsFilePreProcessorWorkFlow() = default;
    ~MsFilePreProcessorWorkFlow() = default;

    Err init(const PythiaParameters &pythiaParameters);

private:

    PythiaParameters m_params;

};


#endif //PYTHIADIACPP_MSFILEPREPROCESSORWORKFLOW_H
