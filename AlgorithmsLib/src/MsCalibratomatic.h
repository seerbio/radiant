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

using DiffPPM = double;
using Coors = QVector<double>;

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

    Err init(
            const QString &calFilePath,
            const QString &matFilePath
            );

    Err writeCalibratomatic(
            const QString &msDataFilePath,
            QString *calibrationMatFilePath,
            QString *calibarationCalFilePath
            );

    // either FrameIndex, or ScanNumber can be key as they are both ints.
    Err recalibratePoints(
            const QMap<FrameIndex, ScanPoints> &indexVsScanPoints,
            QMap<FrameIndex, ScanPoints> *recalIndexVsScanPoints
            );

    static Err recalibratePoints(
            const QMap<ScanNumber, ScanPoints> &scanPoints,
            const QString &calibrationMatFilePath,
            const QString &calibarationCalFilePath,
            QMap<ScanNumber, ScanPoints> *recalScanPoints
    );

    [[nodiscard]] int newStDev();

private:

    Err buildCalibrator(const QMap<QString, QString> &scoreVectorsVsScanFrameFilePaths);

    void filterNNInput(
            const QVector<double> &ppmDiffVals,
            QVector<QPair<DiffPPM, Coors>> *valuesVsTreePoints,
            QVector<QPair<DiffPPM, Coors>> *valuesVsTreePointsRemoved
    );

    Err processLogicForFrameScores(
            const QString &scoreVectorsFilePath,
            const QString &msFrameScansFilePath
    );

    Err getTopNCandidatesPerFrameIndex(
            int topN,
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
    );

    Err getScoredPSMsUntilFirstDecoyIsFound(QVector<PeptideStringWithModsScoreResult> *scoresNoDecoys);

    Err buildCalibrationPoints(const QVector<PeptideStringWithModsScoreResult> &scoresNoDecoys);

    Err buildPeptideSequenceWithModsVsCharge();

    Err loadCalibrationPointsToKDTree();

    Err calculateNewAccuracyMetrics();

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
    double m_stDevNew;

};


#endif //PYTHIADIACPP_MSCALIBRATOMATIC_H
