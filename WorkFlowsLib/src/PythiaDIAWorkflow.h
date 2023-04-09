//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"

#include "GlobalSettings.h"
#include "MsFrame.h"
#include "MsReaderPointerFactory.h"
#include "PythiaParameterReader.h"


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

    Err processFile(const QString &msDatalFilePath);


private:

    Err preprocessDIAFrames(
            const MsReaderPointer &msReaderPointer,
            QVector<MsFrame> *msFrames
            );

    Err scoreCandidatesPerFrame(const QVector<MsFrame> &msFrames);


private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_pepLibUri;
    

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
