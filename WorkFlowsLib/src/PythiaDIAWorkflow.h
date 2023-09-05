//
// Created by anichols on 2/28/23.
//

#ifndef PYTHIADIACPP_PYTHIADIAWORKFLOW_H
#define PYTHIADIACPP_PYTHIADIAWORKFLOW_H

#include "WorkFlowsLib_Exports.h"

#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsFrame.h"
#include "MsReaderPointerFactory.h"
#include "ProteinDigestomatic.h"
#include "PythiaParameterReader.h"


using namespace Error;

class PSMsReaderRow;
class ScoredCandidate;
class FeatureFinderHill;

class WORKFLOWSLIB_EXPORTS PythiaDIAWorkflow {

public:

    PythiaDIAWorkflow() = default;
    ~PythiaDIAWorkflow() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri
            );

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibUri,
            const QString &iRTReCalFilePath
    );

    Err processFile(const QString &msDataFilePath);


private:

    Err buildCalibration(const QString &msDataFilePath);

    Err extractTargetDecoyData(
            const QString &msDataFilePath,
            int topNMs2Ions,
            double selectionFraction,
            QVector<ScoredCandidate> *combinedResults
    );

    Err buildCandidates(
            const QString &msDataFilePath,
            int topNMs2Ions,
            double selectionListFraction,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide
    );

    static Err buildUniqueMsInfoScanKeyVsScanPoints(
            const QString &msDataFilePath,
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
            QMap<ScanNumber, ScanTime> *scanNumberVsScanTime,
            QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
            );

private:

    PythiaParameters m_pythiaParameters;
    QString m_fragLibUri;
    QString m_iRTReCalFilePath;

    FragLibReader m_fragLibReader;
    QVector<MsScanInfo> m_msScanInfos;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
