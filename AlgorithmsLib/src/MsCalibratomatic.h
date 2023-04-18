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
#include "ProteinDigestomatic.h"
#include "PythiaParameterReader.h"

using namespace Error;

class PeptideSequence;

struct PeptideStringWithModsScoreResult {
    PeptideStringWithMods peptideStringWithMods;
    Score score = -1.0;
    FrameIndex frameIndex = -1;
};

class ALGORITHMSLIB_EXPORTS MsCalibratomatic {

public:

    MsCalibratomatic() = default;
    ~MsCalibratomatic() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            FragLibraryTronDIA *fragLibraryTronDia
            );

    Err exec(const QMap<QString, QString> &scoreVectorsVsScanFrameFilePaths);

private:

    Err processLogic(
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

private:

    PythiaParameters m_params;
    QVector<MsFrameScoreVectorReaderRow> m_scoreVectors;
    QMap<FrameIndex, ScanPoints> m_frameIndexVsScanPoints;
    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> m_topCandidatesInFrameIndex;

    FragLibraryTronDIA *m_fragLibraryTronDia;

    QMap<FrameIndex, QVector<ExtractPoints>> m_calibrationPoints;

};


#endif //PYTHIADIACPP_MSCALIBRATOMATIC_H
