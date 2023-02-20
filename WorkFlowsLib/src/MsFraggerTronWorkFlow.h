//
// Created by anichols on 2/19/23.
//

#ifndef PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H
#define PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"


using namespace Error;


class WORKFLOWSLIB_EXPORTS MsFraggerTronWorkFlow {

public:

    explicit MsFraggerTronWorkFlow(const PythiaParameters &pythiaParameters);

    ~MsFraggerTronWorkFlow() = default;

    Err exec(const QString &mzmLFileURI);


private:

    PythiaParameters m_pythiaParameters;

};


#endif //PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H
