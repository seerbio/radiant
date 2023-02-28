//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"

#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "ReCalibratomatic.h"


using namespace Error;


class WORKFLOWSLIB_EXPORTS PythiaDIAWorkflow {

public:

    PythiaDIAWorkflow() = default;
    ~PythiaDIAWorkflow() = default;


    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &pepLibUri
            );

    Err processFile(const QString &mzmlFilePath);


private:

    Err runFirstPassMsFraggerTronWorkFlow(
            const QString &mzmlFilePath,
            QString *firstPassPSMsFilePath
            );

    Err initReCalibratomatic(const QString &firstPassPSMsFilePath);


private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_pepLibUri;

    ReCalibratomatic m_reCalibratomatic;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
