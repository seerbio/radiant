//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibraryTronDIA.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "MsReaderPointerFactory.h"
#include "ProteinDigestomatic.h"
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

    Err preprocessDIAFramesParallel(
            const MsReaderPointer &msReaderPointer,
            QVector<MsFrame> *msFrames
            );

    Err scoreCandidatesPerFrameParallelWrite(
            const QVector<MsFrame> &msFrames,
            const QString &msDataFilePath,
            QMap<UniqueMsInfoScanKey, QString> *uniqueMsInfoScanKeyVsScoredFrameFilePaths
    );

    Err buildTargetCandidatesForFrame(
            const QVector<MsFrame> &msFrames,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, QVector<MS2Ion>>> *framePredictions
    );


private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_pepLibUri;

    FragLibraryTronDIA m_fragLibraryTronDia;
    QMap<PeptideStringWithMods, PeptideSequence> m_peptidesWithModsVsPeptideSequences;


};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
