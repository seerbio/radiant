//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
#define PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H


#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "MsFrameScoreVectorReader.h"
#include "MsUtils.h"
#include "PythiaParameterReader.h"


using namespace Error;


struct PeptideStringWithModsScoreResult {
    PeptideStringWithMods peptideStringWithMods;
    Score score = -1.0;
    FrameIndex frameIndex = -1;
    Charge charge = -1;
};


struct ScoredPSM : public ParquetReaderInputBase {
    FrameIndex frameIndex = -1;
    PeptideStringWithMods peptideStringWithMods;
    ExtractPoints extractPoints;

};


class FrameIndexCandidate;

class ALGORITHMSLIB_EXPORTS MsFrameScoretronProcessormatic {

public:

    friend class MsFrameScoretronProcessormaticTests;

    static Err processLogicForFrameScores(
            const QString &scoreVectorsFilePath,
            const QString &msFrameScansFilePath,
            int topNPSMs,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
    );

private:

    static Err getTopNCandidatesPerFrameIndex(
            const QVector<MsFrameScoreVectorReaderRow> &scoreVectors,
            const QMap<FrameIndex, ScanPoints> &frameIndexVsScanPoints,
            int topNPSMs,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
    );

};


#endif //PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
