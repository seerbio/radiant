//
// Created by anichols on 2/17/23.
//

#ifndef PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H
#define PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "PythiaParameterReader.h"


using namespace Error;


class WORKFLOWSLIB_EXPORTS LibraryBuilderWorkFlow {

public:

    LibraryBuilderWorkFlow() = default;
    ~LibraryBuilderWorkFlow() = default;

    Err exec(
            const PythiaParameters &pythiaParameters,
            const QString &filePath
            );

};


#endif //PYTHIADIACPP_LIBRARYBUILDERWORKFLOW_H
