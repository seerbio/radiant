//
// Created by anichols on 2/19/23.
//

#ifndef PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H
#define PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibraryTron.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"


using namespace Error;


class WORKFLOWSLIB_EXPORTS MsFraggerTronWorkFlow {

public:

    MsFraggerTronWorkFlow() = default;
    ~MsFraggerTronWorkFlow();

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri
            );

    Err processFile(const QString &mzmLFileURI);


private:


private:

    PythiaParameters m_pythiaParameters;
    FragLibraryTron *m_fragLibraryTron;

};


#endif //PYTHIADIACPP_MSFRAGGERTRONWORKFLOW_H
