//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoreVectorReader.h"
#include "MsReaderPointerFactory.h"
#include "TurboXIC.h"


#include <QtConcurrent/QtConcurrent>

#include <iostream>


Err PythiaDIAWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &pepLibUri
        ) {

    ERR_INIT

    pythiaParameters.print();

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(fragLibUri); ree;
    e = ErrorUtils::isNotEmpty(pepLibUri); ree;

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;
    m_pepLibUri = pepLibUri;

    e = m_fragLibraryTronDia.init(
            pythiaParameters,
            fragLibUri
    ); ree;

    ERR_RETURN
}


Err PythiaDIAWorkflow::processFile(const QString &msDatalFilePath) {

    ERR_INIT

    QPair<Err, MsReaderPointer> msReaderPointerResult
            = MsReaderPointerFactory::createInstance(msDatalFilePath);
    e = msReaderPointerResult.first; ree;
    MsReaderPointer msReaderPointer = msReaderPointerResult.second;

    QVector<MsFrame> msFrames;
    e = preprocessDIAFramesParallel(
            msReaderPointer,
            &msFrames
            ); ree;

    e = scoreCandidatesPerFrameParallel(
            msFrames,
            msDatalFilePath
            ); ree;

    ERR_RETURN
}

namespace {

    struct ProcessFileLogicInput {
        MsScanInfo msScanInfo;
        QMap<ScanNumber, ScanPoints> frame;
        PythiaParameters pythiaParameters;
    };

    QPair<Err, MsFrame> processFileParallelLogic(const ProcessFileLogicInput &input) {

        ERR_INIT

        const MsScanInfo &msScanInfo = input.msScanInfo;
        const PythiaParameters &pythiaParameters = input.pythiaParameters;

        MsFrame frame;
        e = frame.init(
                pythiaParameters,
                input.frame,
                msScanInfo.targetScanKey(),
                msScanInfo.collisionEnergy,
                msScanInfo.precursorTargetMz,
                msScanInfo.isoWindowLower,
                msScanInfo.isoWindowUpper
        ); rree;

        e = frame.preprocessMsFrame(
                pythiaParameters.denoise,
                pythiaParameters.deisotope,
                pythiaParameters.smooth
        ); rree

        return {e, frame};
    }

    Err buildProcessFileLogicInputs(
            const MsReaderPointer &msReaderPointer,
            const PythiaParameters &pythiaParameters,
            QVector<ProcessFileLogicInput> *processFileLogicInputs
    ) {

        ERR_INIT

        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaTargetFrames;
        e = msReaderPointer->collateTandemPrecursorTargetsDIA(&diaTargetFrames); ree;

        for (const QMap<ScanNumber, ScanPoints> &frame : diaTargetFrames) {

            const ScanNumber firstScanNumber = frame.firstKey();

            ProcessFileLogicInput processFileLogicInput;

            e = msReaderPointer->getMsScanInfo(
                    firstScanNumber,
                    &processFileLogicInput.msScanInfo
            ); ree;

            processFileLogicInput.frame = frame;
            processFileLogicInput.pythiaParameters = pythiaParameters;

            processFileLogicInputs->push_back(processFileLogicInput);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::preprocessDIAFramesParallel(
        const MsReaderPointer &msReaderPointer,
        QVector<MsFrame> *msFrames
        ) {

    ERR_INIT

    QVector<ProcessFileLogicInput> processFileLogicInputs;
    e = buildProcessFileLogicInputs(
            msReaderPointer,
            m_pythiaParameters,
            &processFileLogicInputs
    ); ree;

    QFuture<QPair<Err, MsFrame>> futures = QtConcurrent::mapped(
            processFileLogicInputs,
            processFileParallelLogic
    );
    futures.waitForFinished();

    for (const QPair<Err, MsFrame> &pr : futures) {
        qDebug() << "FrameSize" << pr.second.scanCount();
        e = pr.first; ree;
        msFrames->push_back(pr.second);
    }

    ERR_RETURN
}

namespace {

    Err groupFramesWithTandemPredictions(
            const QVector<MsFrame> &msFrames,
            const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &framePredictions,
            QVector<QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>>> *processingChunks
            ) {

        ERR_INIT

        for (const MsFrame &frame : msFrames) {
            const UniqueMsInfoScanKey uniqueMsInfoScanKey = frame.uniqueMsInfoScanKey();
            e = ErrorUtils::isTrue(framePredictions.contains(uniqueMsInfoScanKey)); ree;
            processingChunks->push_back({frame, framePredictions.value(uniqueMsInfoScanKey)});
        }

        ERR_RETURN
    }

    struct ScoringMatrices {
        Eigen::MatrixX<double> scoringMatrix;
        Eigen::MatrixX<double> theoMatrixMatrix;
    };

    ScoringMatrices buildTargetScoringMatrices(
            const QVector<MS2Ion> &ms2Ions,
            int frameScanCount,
            int topNMs2Ions,
            double featureFinderTolerancePPM,
            TurboXIC *turboXic
            ) {

        Eigen::MatrixX<double> scoringMat(frameScanCount, topNMs2Ions);
        scoringMat.setZero();

        Eigen::MatrixX<double> theoMat(frameScanCount, topNMs2Ions);
        theoMat.setZero();

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
                theoMat(frameIndex, colIdx) = intensityTheo;
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
    };

    FrameIndexScoreResultOfTarget buildPerFrameIndexScoreVectors(const ScoringMatrices &scoringMatrices) {

        const Eigen::VectorX<double> cosineSimPerFrameIndexOfTarget = EigenUtils::rowWiseCosineSimilarOfMatrices(
                scoringMatrices.scoringMatrix,
                scoringMatrices.theoMatrixMatrix
                );

        const Eigen::VectorX<double> intensityPerFrameIndexOfTarget = scoringMatrices.scoringMatrix.rowwise().sum();

        const Eigen::VectorX<int> foundIonsPerFrameIndexOfTarget
                = (scoringMatrices.scoringMatrix.array() > 0.0).rowwise().count().cast<int>();

        const Eigen::VectorX<double> foundIonsPerFrameIndexOfTargetCubed = foundIonsPerFrameIndexOfTarget
                .cwiseProduct(foundIonsPerFrameIndexOfTarget
                .cwiseProduct(foundIonsPerFrameIndexOfTarget)
                ).cast<double>();

        const Eigen::VectorX<double> scorePerFrameIndexOfTarget = foundIonsPerFrameIndexOfTargetCubed.array()
                                                                * intensityPerFrameIndexOfTarget.array()
                                                                * cosineSimPerFrameIndexOfTarget.array();

        FrameIndexScoreResultOfTarget frameIndexScoreResultOfTarget;

        frameIndexScoreResultOfTarget.cosineSimPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(cosineSimPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget.intensityPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(intensityPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget.foundIonsPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(foundIonsPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget.scorePerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(scorePerFrameIndexOfTarget);

        return frameIndexScoreResultOfTarget;
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

            msFrameScoreVectorReaderRows.push_back(row);
        }

        e = ParquetReader::write(
                msFrameScoreVectorReaderRows,
                outputFilePath
                ); ree;

        ERR_RETURN
    }

    Err scoreFrameTargets(
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
        e = turboXic.init(scanPoints); ree;

        QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> pepStrWModsVsFrameIndexScoreResultOfTargets;
        for (auto it = tandemPreds.begin(); it != tandemPreds.end(); it++) {

            const PeptideStringWithMods &peptideStringWithMods = it.key();
            const QVector<MS2Ion> &ms2Ions = it.value();

            const ScoringMatrices targetScoringMatrices = buildTargetScoringMatrices(
                    ms2Ions,
                    frame.scanCount(),
                    params.topNMs2Ions,
                    params.featureFinderTolerancePPM,
                    &turboXic
                    );

            const FrameIndexScoreResultOfTarget frameIndexScoreResultOfTarget
                            = buildPerFrameIndexScoreVectors(targetScoringMatrices);

            pepStrWModsVsFrameIndexScoreResultOfTargets.insert(peptideStringWithMods, frameIndexScoreResultOfTarget);
        }

        const QString outputFilePath = msDataFilePath + "." + frame.uniqueMsInfoScanKey() + ".frame";

        e = writeFrameTargetScoreVectors(
                pepStrWModsVsFrameIndexScoreResultOfTargets,
                outputFilePath,
                params.minFoundMzPeaks
                ); ree;

        qDebug() << "Found Targets:" << pepStrWModsVsFrameIndexScoreResultOfTargets.size();
        qDebug() << "Frame written to:" << outputFilePath;
        return e;
    }

    Err processFrameLogic(
            const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk,
            const PythiaParameters &params,
            const QString &msDataFilePath
            ) {

        ERR_INIT

        const MsFrame &frame = chunk.first;
        const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPreds = chunk.second;

        const QPair<double, double> &mzTargetStartEnd = frame.precursorMzTargetStartEnd();
        qDebug() << "Processing window" << mzTargetStartEnd.first << mzTargetStartEnd.second;

        e = scoreFrameTargets(
                frame,
                tandemPreds,
                params,
                msDataFilePath
                ); ree;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::scoreCandidatesPerFrameParallel(
        const QVector<MsFrame> &msFrames,
        const QString &msDataFilePath
        ) {

    ERR_INIT

    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, QVector<MS2Ion>>> framePredictions;
    e = buildTargetCandidatesForFrame(
            msFrames,
            &framePredictions
    ); ree;

    QVector<QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>>> processingChunks;
    e = groupFramesWithTandemPredictions(
            msFrames,
            framePredictions,
            &processingChunks
            ); ree;

#define PARALLEL_DIA_SCORE
#ifdef PARALLEL_DIA_SCORE
    qDebug() << "Running scoreCandidatesPerFrame parallel";

    const auto scoringMatrixLogicBinder = std::bind(
            processFrameLogic,
            std::placeholders::_1,
            m_pythiaParameters,
            msDataFilePath
            );

    QFuture<Err> futures = QtConcurrent::mapped(
            processingChunks,
            scoringMatrixLogicBinder
    );
    futures.waitForFinished();

    for (const Err err : futures) {
        e = err; ree;
    }
#else
    qDebug() << "Running scoreCandidatesPerFrame serial";
    for (const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk : processingChunks) {

        e = processFrameLogic(
                chunk,
                m_pythiaParameters,
                msDataFilePath
                ); ree;

    }
#endif

    ERR_RETURN
}

Err PythiaDIAWorkflow::buildTargetCandidatesForFrame(
        const QVector<MsFrame> &msFrames,
        QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, QVector<MS2Ion>>> *framePredictions
        ) {

    ERR_INIT

    for (const MsFrame &frame : msFrames) {

        const QPair<double, double> mzTargetStartEnd = frame.precursorMzTargetStartEnd();
        qDebug() << "Collecting Targets for:" << mzTargetStartEnd.first << mzTargetStartEnd.second;

        QMap<PeptideStringWithMods, QVector<MS2Ion>> peptideStringWithModsVsMS2Ions;

        e = m_fragLibraryTronDia.getMS2Ions(
                mzTargetStartEnd.first,
                mzTargetStartEnd.second,
                m_pythiaParameters.topNMs2Ions,
                &peptideStringWithModsVsMS2Ions
        ); ree;

        framePredictions->insert(frame.uniqueMsInfoScanKey(), peptideStringWithModsVsMS2Ions);
    }

    ERR_RETURN
}
