//
// Created by anichols on 4/16/23.
//

#ifndef PYTHIADIACPP_MSCALIBRATOMATIC_H
#define PYTHIADIACPP_MSCALIBRATOMATIC_H


#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FragLibraryTronDIA.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "MsFrameScoreVectorReader.h"
#include "MsUtils.h"
#include "NearestNeighborsSearch.h"
#include "ProteinDigestomatic.h"
#include "PythiaParameterReader.h"

using namespace Error;

class PeptideSequence;

struct PeptideStringWithModsScoreResult {
    PeptideStringWithMods peptideStringWithMods;
    Score score = -1.0;
    FrameIndex frameIndex = -1;
    Charge charge = -1;
};

class ALGORITHMSLIB_EXPORTS MsCalibratomatic {

public:

    MsCalibratomatic();
    ~MsCalibratomatic() = default;

    Err init(
            const QMap<QString, QString> &scoreVectorsVsScanFrameFilePaths,
            const PythiaParameters &pythiaParameters,
            int calPointK,
            FragLibraryTronDIA *fragLibraryTronDia
            );

private:

    Err buildCalibrator(const QMap<QString, QString> &scoreVectorsVsScanFrameFilePaths);

    Err processLogicForFrameScores(
            const QString &scoreVectorsFilePath,
            const QString &msFrameScansFilePath
    );

    Err buildFrameIndexVsScanPoints(const QVector<MsFrameScanPointRows> &msFrameScanPointRows);

    Err getTopNCandidatesPerFrameIndex(
            int topN,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
    );

    Err getScoredPSMsUntilFirstDecoyIsFound(QVector<PeptideStringWithModsScoreResult> *scoresNoDecoys);

    Err buildCalibrationPoints(const QVector<PeptideStringWithModsScoreResult> &scoresNoDecoys);

    Err buildPeptideSequenceWithModsVsCharge();

    Err loadCalibrationPointsToKDTree();

private:

    //Need to be cleared w/ each new ScoreFrame iteration
    QVector<MsFrameScoreVectorReaderRow> m_scoreVectors;
    QMap<FrameIndex, ScanPoints> m_frameIndexVsScanPoints;
    QMap<PeptideStringWithMods, Charge> m_peptideWithModsVsCharge;
    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> m_topCandidatesInFrameIndex;

    //Never cleared
    PythiaParameters m_params;
    FragLibraryTronDIA *m_fragLibraryTronDia;
    QMap<FrameIndex, QVector<ExtractPoints>> m_calibrationPoints;
    NearestNeighborsSearch m_nnSearch;
    int m_calPointK;

};


#endif //PYTHIADIACPP_MSCALIBRATOMATIC_H
