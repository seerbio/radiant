//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsFrame.h"
#include "MsReaderPointerFactory.h"
#include "ProteinDigestomatic.h"
#include "PythiaParameterReader.h"


using namespace Error;
using ScoreVecFilePath = QString;
using ExtractsFilePath = QString;

class FrameParallelInput;
class PSMsReaderRow;

class WORKFLOWSLIB_EXPORTS PythiaDIAWorkflow {

public:

    PythiaDIAWorkflow() = default;
    ~PythiaDIAWorkflow() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri
            );

    Err processFile(const QString &msDataFilePath);


private:

    Err buildPSMResultsForCalibrationFile(
            const QVector<FrameParallelInput> &frameParallelInputs,
            QVector<QPair<ScoreVecFilePath, ExtractsFilePath>> *frameScoreVectorsAndExtractFilePaths
            );

    static Err buildFrameScoreVectors(
            const QVector<FrameParallelInput> &frameParallelInputs,
            QVector<QPair<ScoreVecFilePath, ExtractsFilePath>> *frameScoreVectorsAndExtractFilePaths
            );

    static Err processFrameScoreVectors(
            const QStringList &frameScoreVectorsFilePaths
            );

private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;


};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
