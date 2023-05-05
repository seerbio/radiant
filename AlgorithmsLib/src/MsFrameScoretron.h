//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRON_H
#define PYTHIADIACPP_MSFRAMESCORETRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "PSMsReader.h"
#include "PythiaParameterReader.h"

using namespace Error;

class FrameIndexScoreResultOfTarget;

class ALGORITHMSLIB_EXPORTS MsFrameScoretron {

public:

    MsFrameScoretron() = default;
    ~MsFrameScoretron() = default;

    QPair<Err, QVector<PSMsReaderRow>> scoreCandidates(
            const PythiaParameters &params,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop,
            bool applySmooth2D
    );

private:

    Err calculateDiscriminateScoreForFrameIndexes(
            const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> &topCansInFrameIndex,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore >>> *topCansInFrameIndexVsDiscScore
            );

    Err calculateDiscriminateScoreForFrame(
            const QVector<QPair<PeptideStringWithMods, Score>> &peptideStringWithModsScore,
            FrameIndex frameIndex,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> *frameIndexVsPeptideStringWithModsDiscScore
            );

    Err buildPSMsReaderRows(
            const QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> &pepStrWModsVsFrameIndexScoreResultOfTargets,
            const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> &topCansInFrameIndex,
            const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> &topCansInFrameIndexVsDiscScore,
            QVector<PSMsReaderRow> *psmsReaderRows
    );

private:

    QMap<PeptideStringWithMods, QVector<MS2Ion>> m_fragPreds;
    QMap<PeptideStringWithMods, bool> m_fragPredsIsDecoy;
    MsFrame m_msFrame;
    PythiaParameters m_params;

};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
