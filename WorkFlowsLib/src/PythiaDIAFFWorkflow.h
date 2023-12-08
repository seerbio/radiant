//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

using namespace Error;


class MsScanInfo;
class TargetDecoyCandidatePair;


class WORKFLOWSLIB_EXPORTS PythiaDIAFFWorkflow {

public:

    friend class PythiaDIAFFWorkflowTests;

    PythiaDIAFFWorkflow();
    ~PythiaDIAFFWorkflow() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &fastaUri
            );

    Err processFile(const QString &msDataFilePath);

private:

    Err buildCalibration();

    Err buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            const QVector<MsScanInfo> &msScanInfos,
            double selectionFraction,
            QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *uniqueInfoScanKeyVsTargetDecoyCandidatePointers
            );

private:

    TargetDecoyCandidatePairManager m_targetDecoyCandidatePairManager;

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_fastaUri;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
