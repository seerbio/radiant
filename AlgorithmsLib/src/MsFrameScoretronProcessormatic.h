//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
#define PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H


#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FragLibReader.h"
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


class ALGORITHMSLIB_EXPORTS MsFrameScoretronProcessormatic {

public:

    friend class MsFrameScoretronProcessormaticTests;

    MsFrameScoretronProcessormatic() = default;
    ~MsFrameScoretronProcessormatic() = default;

    Err init(
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &fragPreds,
            const MsFrame &msFrame,
            const PythiaParameters &params,
            const QString &scoreVectorsFilePath
            );

    Err processLogicForFrameScores(
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> *topCansInFrameIndexVsDiscScore
    );


private:

    Err getTopNCandidatesPerFrameIndex(
            const QVector<MsFrameScoreVectorReaderRow> &scoreVectors,
            const QMap<FrameIndex, ScanPoints> &frameIndexVsScanPoints,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
    );

    Err calculateDiscriminateScoreForFrameIndexes(
            const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> &topCansInFrameIndex,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore >>> *topCansInFrameIndexVsDiscScore
    );

    Err calculateDiscriminateScoreForFrame(
            const QVector<QPair<PeptideStringWithMods, Score>> &peptideStringWithModsScore,
            const ScanPoints &scanPoints,
            FrameIndex frameIndex,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> *frameIndexVsPeptideStringWithModsDiscScore
    );

private:

    QMap<PeptideStringWithMods, QVector<MS2Ion>> m_fragPreds;
    MsFrame m_msFrame;
    PythiaParameters m_params;
    QString m_scoreVectorsFilePath;

};


#endif //PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
