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

struct FrameIndexScoreResultOfTarget {
    QVector<double> cosineSimPerFrameIndexOfTargetVec;
    QVector<double> intensityPerFrameIndexOfTargetVec;
    QVector<int> foundIonsPerFrameIndexOfTargetVec;
    QVector<double> scorePerFrameIndexOfTargetVec;
    PeakIntegrationIndexes bestScorePeakLimits = {-1, -1};
    int charge = -1;
};

class ALGORITHMSLIB_EXPORTS MsFrameScoretron {

public:

    friend class MissingPeptideManualTroubleshooter;
    friend class MsFrameScoretronProcessormaticTests;

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

    Err buildPSMsReaderRows(
            const QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> &pepStrWModsVsFrameIndexScoreResultOfTargets,
            const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> &topCansInFrameIndex,
            const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> &topCansInFrameIndexVsDiscScore,
            QVector<PSMsReaderRow> *psmsReaderRows
    );

    static Err buildFragIonLibForTargetMz(
            const PythiaParameters &params,
            const QString &fragLibUri,
            const QPair<double, double> &mzTargetStartStop,
            QMap<PeptideStringWithMods, QVector<MS2Ion>> *outputIons,
            QMap<PeptideStringWithMods, bool> *outputDecoy
    );

    static Err buildMsFrame(
            const QString &msDataFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const PythiaParameters &params,
            QPair<double, double> mzTargetStartStop,
            bool applySmooth2D,
            MsFrame *msFrame
    );

    static Err processFrameLogic(
            const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk,
            const PythiaParameters &params,
            QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
    );

    static void filterByFoundMzCount(
            int minFoundMzCount,
            QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
    );

    static Err writeFrameTargetScoreVectors(
            const QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> &pepStrWModsVsFrameIndexScoreResultOfTargets,
            const QString &outputFilePath
    );

private:

    QMap<PeptideStringWithMods, QVector<MS2Ion>> m_fragPreds;
    QMap<PeptideStringWithMods, bool> m_fragPredsIsDecoy;
    MsFrame m_msFrame;
    PythiaParameters m_params;

};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
