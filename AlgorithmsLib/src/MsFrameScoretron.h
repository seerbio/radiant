//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRON_H
#define PYTHIADIACPP_MSFRAMESCORETRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FeatureFinderHillBuilder.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "PeakIntegratomatic.h"
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
    friend class MsFrameScoretronTests;

    MsFrameScoretron();
    ~MsFrameScoretron() = default;

    Err init(
            const PythiaParameters &params,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
            );

    Err extractHillsForCandidtates(QString *frameHillsFilePath);


private:

    Err buildFragIonLibForTargetMz(const QString &fragLibUri);

    Err groupHillsForFrameCandidates(
            QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
            );

    Err getCandidateHills(
            const MS2IonsSeparated &ms2IonsTandemPred,
            QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> *featureFinderHills
            );

    Err getHillsIonType(
            const QMap<IonIndex, MS2Ion> &ions,
            QMap<IonIndex, QVector<FeatureFinderHill>> *featureFinderHills
    );


private:

    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> m_pepStrWModsVsFrameIndexScoreResultOfTargets;
    QMap<PeptideStringWithMods, MS2IonsSeparated> m_fragPreds;
    QMap<PeptideStringWithMods, bool> m_fragPredsIsDecoy;
    MsFrame m_msFrame;
    FeatureFinderHillBuilder m_featureFinderHillBuilder;
    FrameIndex m_frameIndexMin;
    FrameIndex m_frameIndexMax;

    PythiaParameters m_params;
    QString m_msDataFilePath;
    UniqueMsInfoScanKey m_uniqueMsInfoScanKey;
    QPair<double, double> m_mzTargetStartStop;

};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
