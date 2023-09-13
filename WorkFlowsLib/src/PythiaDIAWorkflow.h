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

class MsReaderParquet;
class PSMsReaderRow;
class ScoredCandidate;

class WORKFLOWSLIB_EXPORTS PythiaDIAWorkflow {

public:

    PythiaDIAWorkflow();
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

    Err deisotopeScans(MsReaderParquet *msReaderParquet);

    Err buildCalibration(MsReaderParquet *msReaderParquet);

    Err extractTargetDecoyData(
            const PythiaParameters &pythiaParameters,
            int topNMs2Ions,
            double selectionFraction,
            MsReaderParquet *msReaderParquet,
            QVector<ScoredCandidate> *combinedResults
    );

    Err buildCandidates(
            int topNMs2Ions,
            double selectionListFraction,
            MsReaderParquet *msReaderParquet,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide
    );

    Err optimizeParameters(MsReaderParquet *msReaderParquet);

    Err mainAnalysis(MsReaderParquet *msReaderParquet);

    static Err buildUniqueMsInfoScanKeyVsScanPoints(
            MsReaderParquet *msReaderParquet,
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

    MsCalibratomatic m_msCalibratomatic;

    const int m_minTopNMs2Ions;

};


#endif //PYTHIADIACPP_PYTHIADIAWORKFLOW_H
