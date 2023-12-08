//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "PythiaParameterReader.h"

using namespace Error;




class WORKFLOWSLIB_EXPORTS PythiaDIAFFWorkflow {

public:

    PythiaDIAFFWorkflow();
    ~PythiaDIAFFWorkflow() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &fastaUri
            );

    Err processFile(const QString &msDataFilePath);

private:



private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_fastaUri;


    
};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
