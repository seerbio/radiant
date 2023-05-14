//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
#define PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "ExtractedScansReader.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "MsFrameScoreVectorReader.h"
#include "MsUtils.h"
#include "PythiaParameterReader.h"
#include "TandemSpectraDeconvolvotron.h"


using namespace Error;

class PSMsReaderRow;

struct FrameStats {

    double scoreMean = -1.0;
    double scoreMedian = -1.0;
    double scoreStDev = -1.0;
    double scoreMin = -1.0;
    double scoreMax = -1.0;

    double discScoreMean = -1.0;
    double discScoreMedian = -1.0;
    double discScoreStDev = -1.0;
    double discScoreMin = -1.0;
    double discScoreMax = -1.0;

    int frameCandidateCount = -1;

};

struct ReScoreVals {
    ReScore reScore = -1.0;
    Score score = -1.0;
    double cosineSim = -1.0;
    double klDiv = -1.0;
    double fractionFound = -1.0;
    int ionsFound = -1;
};


class ALGORITHMSLIB_EXPORTS MsFrameScoretronProcessormatic {

public:

    friend class MsFrameScoretronProcessormaticTests;

    MsFrameScoretronProcessormatic() = default;
    ~MsFrameScoretronProcessormatic() = default;

    Err init(
            const QString &frameScoreVecFilePath,
            const QString &frameExtractedScansFilePath,
            const PythiaParameters &pythiaParameters,
            const QString &msDataFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            QPair<double, double> mzTargetStartStop
            );

    Err processFrameScoreVectors(QVector<PSMsReaderRow> *psmReaderRows);


private:

    Err buildPepStrWModsVsExtractedScanRow();

    Err buildFrameIndexVsApexScorePeptideStringWithMods();

    Err buildPepStringWithModsVsMS2Ions();

    Err rescoreCandidatesWithFullPrediction();

    Err processorLogic(
            const QVector<PeptideStringWithMods> &peptideStringWithMods,
            FrameIndex frameIndex
            );

    Err collateMS2IonsByPepStrWithModsForFrameIndex(
            const QVector<PeptideStringWithMods> &peptideStringWithMods,
            QMap<PeptideStringWithMods, QVector<MS2Ion>> *peptideByExtractedPoints
            );

    Err calculateDiscriminateScoreForFrameIndex(
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &peptideByExtractedPoints,
            const ScanPoints &scanPoints,
            FrameIndex frameIndex
    );

    Err buildFrameIndexVsFrameStats();


    Err compileScores(QVector<PSMsReaderRow> *psmReaderRows);




private:

    PythiaParameters m_params;
    QString m_scoreVectorsFilePath;
    QString m_frameExtractedScansFilePath;
    TandemSpectraDeconvolvotron m_deconvolvotron;

    MsFrame m_msFrame;
    QMap<PeptideStringWithMods, ExtractedScansReaderRow> m_pepStrWModsVsExtractedScanRow;
    QMap<PeptideStringWithMods, QVector<MS2Ion>> m_pepStrWModsVsMS2Ions;
    QMap<PeptideStringWithMods, Score> m_pepStrWModsVsOgScore;
    QMap<PeptideStringWithMods, Score> m_pepStrWModsVsCharge;
    QMap<PeptideStringWithMods, bool> m_pepStrWModsVsIsDecoy;
    QMap<FrameIndex, QVector<PeptideStringWithMods>> m_frameIndexVsApexScorePeptideStringWithMods;
    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, TandemDeconvolverResult>>> m_frameIndexVsPeptideDeconResult;
    QMap<PeptideStringWithMods, ReScoreVals> m_pepStringModsVsRescore;
    QMap<FrameIndex, FrameStats> m_frameIndexVsPeptideFrameStats;

};


#endif //PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
