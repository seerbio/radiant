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

class TandemDeconvolverResult;

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

    Err init(
            const PythiaParameters &params,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
            );

    Err buildFrameScoreVectors(QString *frameScoreVectorsFilePath);

    Err buildAllExtractedTheoriticalPointsFromTargetKeyFrame(QString *frameExtractedPointsFilePath);

private:

    Err buildFragIonLibForTargetMz(const QString &fragLibUri);

    Err buildMsFrame();

    Err processFrameLogic(
            const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk,
            const PythiaParameters &params,
            QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
    );

    void filterByFoundMzCount(
            int minFoundMzCount,
            QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
    );

    Err writeFrameTargetScoreVectors(const QString &outputFilePath);

private:

    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> m_pepStrWModsVsFrameIndexScoreResultOfTargets;
    QMap<PeptideStringWithMods, QVector<MS2Ion>> m_fragPreds;
    QMap<PeptideStringWithMods, bool> m_fragPredsIsDecoy;
    MsFrame m_msFrame;

    PythiaParameters m_params;
    QString m_msDataFilePath;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    QPair<double, double> m_mzTargetStartStop;

};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
