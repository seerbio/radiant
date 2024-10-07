//
// Created by andrewnichols on 9/28/24.
//

#ifndef PYTHIADIAFFWORKFLOWSHAREDMETHODS_H
#define PYTHIADIAFFWORKFLOWSHAREDMETHODS_H

#include "WorkFlowsLib_Exports.h"

#include "CandidateScores.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibratomatic.h"
#include "MsReaderBase.h"
#include "MsReaderPointerAcc.h"
#include "PythiaParameterReader.h"

using namespace Error;

class TurboXIC;

class WORKFLOWSLIB_EXPORTS PythiaDIAFFWorkflowSharedMethods {

public:

    static Err buildUniqueMsScanInfosForProcessing(
        const QVector<MsScanInfo> &uniqueMsScanInfos,
        int numberOfUniqueScanInfosForPurpose,
        int offset,
        QVector<MsScanInfo> *uniqueMsScanInfosPurpose
        );

    static Err buildMzTargetKeyVsTurboXicPntrs(
        const QVector<MsScanInfo> &uniqueMsScanInfos,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<MzTargetKey, TurboXIC*> *mzTargetKeyVsTurboXicPntr
        );

    static Err buildCandidateScoresPtrs(
        QVector<CandidateScores> &candidateScores,
        QVector<CandidateScores*> *candidateScoresPntrs
        );

    static void sortCandidatePointersDiscScoreDesc(
        QVector<CandidateScores*> *candidateScoresPntrs
        );

    static Err processBatch(
            QVector<CandidateScores> &candidateScoresVecBatch,
            const PythiaParameters &pythiaParameters,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            QVector<CandidateScores*> *candidateScoresVecBatchPntrs,
            QMap<int, int> *fdrVsCounts
        );

    static Err buildPeptideKeyVsTargetDecoyCandidateScoresPntrs(
        const QVector<CandidateScores*> &candidateScores,
        QMap<PeptideSequenceWithModsChargeAndTargetKey , QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *peptideKeyVsTargetDecoyCandidateScoresPntrs
        );

    static Err buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
        const PythiaParameters &pythiaParameters,
        const QVector<MsScanInfo> &msScanInfos,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
        );

    static Err buildMsCalibrationReaderRows(
        const MSLevelEnum &msLevel,
        const QVector<CandidateScores*> &_candidateScores,
        int verbosity,
        QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
        );

};


#endif //PYTHIADIAFFWORKFLOWSHAREDMETHODS_H
