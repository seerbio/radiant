//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoreVectorReader.h"
#include "PeakIntegratomatic.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>


namespace {

    struct ScoringMatrices {
        Eigen::MatrixX<double> scoringMatrix;
        Eigen::MatrixX<double> theoMatrixMatrix;
    };

    Eigen::MatrixX<double> buildTheoMat(
            const QVector<MS2Ion> &ms2Ions,
            int frameScanCount,
            int topNMs2Ions
    ) {

        Eigen::MatrixX<double> theoMat(frameScanCount, topNMs2Ions);
        theoMat.setZero();

        Eigen::VectorX<double> theoIntensities(topNMs2Ions);
        for (int i = 0; i < topNMs2Ions; i++) {
            theoIntensities(i) = ms2Ions.at(i).intensity;
        }
        for (int i = 0; i < frameScanCount; i++) {
            theoMat.row(i) = theoIntensities;
        }

        return theoMat;
    }

    ScoringMatrices buildTargetScoringMatrices(
            const QVector<MS2Ion> &ms2Ions,
            int frameScanCount,
            int topNMs2Ions,
            double featureFinderTolerancePPM,
            TurboXIC *turboXic
    ) {

        Eigen::MatrixX<double> scoringMat(frameScanCount, topNMs2Ions);
        scoringMat.setZero();

        const Eigen::MatrixX<double> theoMat = buildTheoMat(
                ms2Ions,
                frameScanCount,
                topNMs2Ions
        );

        for (int colIdx = 0; colIdx < ms2Ions.size(); colIdx++) {

            const double mz = ms2Ions.at(colIdx).mz;
            const Intensity intensityTheo = ms2Ions.at(colIdx).intensity;

            const double massTol = MathUtils::calculatePPM(
                    mz,
                    featureFinderTolerancePPM
            );

            const double mzMin = mz - massTol;
            const double mzMax = mz + massTol;
            const FrameIndex frameIndexMin = 0;
            const FrameIndex frameIndexMax = frameScanCount;

            const XICPoints xic = turboXic->extractPoints(
                    mzMin,
                    mzMax,
                    frameIndexMin,
                    frameIndexMax
            );

            if (xic.isEmpty()) {
                continue;
            }

            for (auto it = xic.begin(); it != xic.end(); it++) {
                const FrameIndex frameIndex = it.key();
                const Intensity intensity = it.value();

                scoringMat(frameIndex, colIdx) = intensity;

            }
        }

        ScoringMatrices sm;
        sm.scoringMatrix = scoringMat;
        sm.theoMatrixMatrix = theoMat;

        return sm;
    }

    struct FrameIndexScoreResultOfTarget {
        QVector<double> cosineSimPerFrameIndexOfTargetVec;
        QVector<double> intensityPerFrameIndexOfTargetVec;
        QVector<int> foundIonsPerFrameIndexOfTargetVec;
        QVector<double> scorePerFrameIndexOfTargetVec;
        PeakIntegrationIndexes bestScorePeakLimits = {-1, -1};
    };

    Err buildPerFrameIndexScoreVectors(
            const ScoringMatrices &scoringMatrices,
            FrameIndexScoreResultOfTarget *frameIndexScoreResultOfTarget
    ) {

        ERR_INIT

        const Eigen::VectorX<double> cosineSimPerFrameIndexOfTarget = EigenUtils::rowWiseCosineSimilarOfMatrices(
                scoringMatrices.scoringMatrix,
                scoringMatrices.theoMatrixMatrix
        );

        const Eigen::VectorX<double> intensityPerFrameIndexOfTarget = scoringMatrices.scoringMatrix.rowwise().sum();

        const Eigen::VectorX<int> foundIonsPerFrameIndexOfTarget
                = (scoringMatrices.scoringMatrix.array() > 0.0).rowwise().count().cast<int>();

        const Eigen::VectorX<double> foundIonsPerFrameIndexOfTargetFactorial
                = foundIonsPerFrameIndexOfTarget.unaryExpr([](int x) { return MathUtils::factorial(x); }).cast<double>();

        const Eigen::VectorX<double> scorePerFrameIndexOfTarget = foundIonsPerFrameIndexOfTargetFactorial.array()
                                                                  * cosineSimPerFrameIndexOfTarget.array();

        frameIndexScoreResultOfTarget->cosineSimPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(cosineSimPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget->intensityPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(intensityPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget->foundIonsPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(foundIonsPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget->scorePerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(scorePerFrameIndexOfTarget);

        const double stopThresholdFraction = 0.1;
        const int filterLength = 7;
        const int smoothCount = 1;
        e = PeakIntegratomatic::simpleIntegrator(
                frameIndexScoreResultOfTarget->scorePerFrameIndexOfTargetVec,
                stopThresholdFraction,
                filterLength,
                smoothCount,
                &frameIndexScoreResultOfTarget->bestScorePeakLimits,
                &frameIndexScoreResultOfTarget->scorePerFrameIndexOfTargetVec
        ); ree;

        ERR_RETURN
    }

    Err writeFrameTargetScoreVectors(
            const QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> &pepStrWModsVsFrameIndexScoreResultOfTargets,
            const QString &outputFilePath,
            int minFoundMzCount
    ) {

        ERR_INIT

        QVector<MsFrameScoreVectorReaderRow> msFrameScoreVectorReaderRows;
        for (auto it = pepStrWModsVsFrameIndexScoreResultOfTargets.begin(); it != pepStrWModsVsFrameIndexScoreResultOfTargets.end(); it++) {

            const PeptideStringWithMods &peptideStringWithMods = it.key();
            const FrameIndexScoreResultOfTarget &frameIndexScoreResultOfTarget = it.value();

            const QVector<int> &founds = frameIndexScoreResultOfTarget.foundIonsPerFrameIndexOfTargetVec;
            const int maxFoundCount = *std::max_element(founds.begin(), founds.end());

            if (maxFoundCount < minFoundMzCount) {
                continue;
            }

            MsFrameScoreVectorReaderRow row;
            row.peptideStringWithMods = peptideStringWithMods;
            row.cosineSimPerFrameIndexOfTargetVec = frameIndexScoreResultOfTarget.cosineSimPerFrameIndexOfTargetVec;
            row.foundIonsPerFrameIndexOfTargetVec = frameIndexScoreResultOfTarget.foundIonsPerFrameIndexOfTargetVec;
            row.scorePerFrameIndexOfTargetVec = frameIndexScoreResultOfTarget.scorePerFrameIndexOfTargetVec;
            row.intensityPerFrameIndexOfTargetVec = frameIndexScoreResultOfTarget.intensityPerFrameIndexOfTargetVec;
            row.scorePeakStart = frameIndexScoreResultOfTarget.bestScorePeakLimits.first;
            row.scorePeakEnd = frameIndexScoreResultOfTarget.bestScorePeakLimits.second;

            msFrameScoreVectorReaderRows.push_back(row);
        }

        e = ParquetReader::write(
                msFrameScoreVectorReaderRows,
                outputFilePath
        ); ree;

        ERR_RETURN
    }

    QPair<Err, QString> scoreFrameTargets(
            const MsFrame &frame,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPreds,
            const PythiaParameters &params,
            const QString &msDataFilePath
    ) {

        ERR_INIT

        qDebug() << "Building matrix" << frame.uniqueMsInfoScanKey();
        qDebug() << "Frame size" << frame.scanCount();

        const QMap<FrameIndex, ScanPoints> scanPoints = frame.frameIndexVsScanPoints();

        TurboXIC turboXic;
        e = turboXic.init(scanPoints); rree;

        QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> pepStrWModsVsFrameIndexScoreResultOfTargets;
        for (auto it = tandemPreds.begin(); it != tandemPreds.end(); it++) {

            const PeptideStringWithMods &peptideStringWithMods = it.key();
            const QVector<MS2Ion> &ms2Ions = it.value();

//#define DEBUG_FRAME
#ifdef DEBUG_FRAME
            if (peptideStringWithMods != "AAANAVXQDK") {
                continue;
            }
            for (const MS2Ion &i : ms2Ions) {
                qDebug() << i.mz << i.intensity << i.ionLabel;
            }
#endif
            const ScoringMatrices targetScoringMatrices = buildTargetScoringMatrices(
                    ms2Ions,
                    frame.scanCount(),
                    params.topNMs2Ions,
                    params.featureFinderTolerancePPM,
                    &turboXic
            );

            FrameIndexScoreResultOfTarget frameIndexScoreResultOfTarget;
            e = buildPerFrameIndexScoreVectors(
                    targetScoringMatrices,
                    &frameIndexScoreResultOfTarget
            ); rree;

            pepStrWModsVsFrameIndexScoreResultOfTargets.insert(peptideStringWithMods, frameIndexScoreResultOfTarget);
        }

        const QString outputFilePath = msDataFilePath + "." + frame.uniqueMsInfoScanKey() + ".frame";

        e = writeFrameTargetScoreVectors(
                pepStrWModsVsFrameIndexScoreResultOfTargets,
                outputFilePath,
                params.minFoundMzPeaks
        ); rree;

        qDebug() << "Found Targets:" << pepStrWModsVsFrameIndexScoreResultOfTargets.size();
        qDebug() << "Frame written to:" << outputFilePath;
        return {e, outputFilePath};
    }

    QPair<Err, QString> processFrameLogic(
            const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk,
            const PythiaParameters &params,
            const QString &msDataFilePath
    ) {

        ERR_INIT

        const MsFrame &frame = chunk.first;
        const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPreds = chunk.second;

        const QPair<double, double> &mzTargetStartEnd = frame.precursorMzTargetStartEnd();
        qDebug() << "Processing window" << mzTargetStartEnd.first << mzTargetStartEnd.second;

        QPair<Err, QString> scoreFrameTargetsResult = scoreFrameTargets(
                frame,
                tandemPreds,
                params,
                msDataFilePath
        );

        e = scoreFrameTargetsResult.first; rree;

        return scoreFrameTargetsResult;
    }

}//namespace
QPair<Err, QString> MsFrameScoretron::scoreCandidatesFrame(
        const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk,
        const PythiaParameters &params,
        const QString &msDataFilePath
        ) {

    ERR_INIT

    QPair<Err, QString> scoreFrameTargetsResult = processFrameLogic(
            chunk,
            params,
            msDataFilePath
    );

    e = scoreFrameTargetsResult.first; rree;

    return scoreFrameTargetsResult;
}


