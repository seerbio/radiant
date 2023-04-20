//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoretron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "MsReaderPointerFactory.h"
#include "PeakIntegratomatic.h"

#include <QtConcurrent/QtConcurrent>


Err PythiaDIAWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri
        ) {

    ERR_INIT

    pythiaParameters.print();

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(fragLibUri); ree;

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;

    e = m_fragLibraryTronDia.init(
            pythiaParameters,
            fragLibUri
    ); ree;

    ERR_RETURN
}

namespace {

    Err buildScoreVectorsVsScanFrameFilePaths(
            const QMap<UniqueMsInfoScanKey, QString> &uniqueMsInfoScanKeyVsScoreVecsFrameFilePaths,
            const QMap<UniqueMsInfoScanKey, QString> &uniqueMsInfoScanKeyVsMsFrameFilePath,
            QMap<QString, QString> *scoreVectorsVsScanFrameFilePaths
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKeyVsMsFrameFilePath); ree;
        e = ErrorUtils::isEqualString(
                uniqueMsInfoScanKeyVsMsFrameFilePath.keys(),
                uniqueMsInfoScanKeyVsScoreVecsFrameFilePaths.keys()
                ); ree;

        const QList<UniqueMsInfoScanKey> &uniqueMsInfoScanKeys = uniqueMsInfoScanKeyVsMsFrameFilePath.keys();
        for (const UniqueMsInfoScanKey &uniqueMsInfoScanKey : uniqueMsInfoScanKeys) {

            scoreVectorsVsScanFrameFilePaths->insert(
                    uniqueMsInfoScanKeyVsScoreVecsFrameFilePaths.value(uniqueMsInfoScanKey),
                    uniqueMsInfoScanKeyVsMsFrameFilePath.value(uniqueMsInfoScanKey)
                    );
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::processFile(const QString &msDatalFilePath) {

    ERR_INIT

    QPair<Err, MsReaderPointer> msReaderPointerResult
            = MsReaderPointerFactory::createInstance(msDatalFilePath);
    e = msReaderPointerResult.first; ree;
    MsReaderPointer msReaderPointer = msReaderPointerResult.second;

    QMap<UniqueMsInfoScanKey, QString> uniqueMsInfoScanKeyVsScoredFrameFilePathsCalibration;
    QMap<UniqueMsInfoScanKey, QString> uniqueMsInfoScanKeyVsMsFrameFilePathCalibration;

    const int numberOfFramesToProcessForCalibration = -1;
    e = buildCandidateScoreVectors(
            msReaderPointer,
            numberOfFramesToProcessForCalibration,
            &uniqueMsInfoScanKeyVsScoredFrameFilePathsCalibration,
            &uniqueMsInfoScanKeyVsMsFrameFilePathCalibration
            ); ree;

    QMap<QString, QString> scoreVectorsVsScanFrameFilePaths;
    e = buildScoreVectorsVsScanFrameFilePaths(
            uniqueMsInfoScanKeyVsScoredFrameFilePathsCalibration,
            uniqueMsInfoScanKeyVsMsFrameFilePathCalibration,
            &scoreVectorsVsScanFrameFilePaths
            ); ree;

    MsCalibratomatic msCalibratomatic;
    const int calPointK = 25;
    e = msCalibratomatic.init(
            scoreVectorsVsScanFrameFilePaths,
            m_pythiaParameters,
            calPointK,
            &m_fragLibraryTronDia
    );

    ERR_RETURN
}

namespace {

    Err writeMsFrames(
            const QVector<MsFrame> &msFrames,
            const QString &msDataFilePath,
            QMap<UniqueMsInfoScanKey, QString> *uniqueMsInfoScanKeyVsMsFrameFilePath
            ) {

        ERR_INIT

        for (const MsFrame &frame : msFrames) {

            const QString outputFilePath
                    = msDataFilePath + "." + frame.uniqueMsInfoScanKey() + ".frameScans";

            qDebug() << "Writing frame to" << outputFilePath;
            e = frame.writeFramScans(outputFilePath); ree;
            uniqueMsInfoScanKeyVsMsFrameFilePath->insert(
                    frame.uniqueMsInfoScanKey(),
                    outputFilePath
            );
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::buildCandidateScoreVectors(
        const MsReaderPointer &msReaderPointer,
        int numberOfFramesToProcess,
        QMap<UniqueMsInfoScanKey, QString> *uniqueMsInfoScanKeyVsScoredFrameFilePaths,
        QMap<UniqueMsInfoScanKey, QString> *uniqueMsInfoScanKeyVsMsFrameFilePath
        ) {

    ERR_INIT

    uniqueMsInfoScanKeyVsScoredFrameFilePaths->clear();

    QVector<MsFrame> msFrames;
    e = preprocessDIAFramesParallel(
            msReaderPointer,
            numberOfFramesToProcess,
            &msFrames
    ); ree;

    const QString msDataFilePath = msReaderPointer->filePath();
    e = ErrorUtils::isNotEmpty(msDataFilePath); ree;

    e = scoreCandidatesPerFrameParallelWrite(
            msFrames,
            msDataFilePath,
            uniqueMsInfoScanKeyVsScoredFrameFilePaths
    ); ree;

    e = writeMsFrames(
            msFrames,
            msDataFilePath,
            uniqueMsInfoScanKeyVsMsFrameFilePath
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
            int numberOfFramesToProcess,
            QVector<ProcessFileLogicInput> *processFileLogicInputs
    ) {

        ERR_INIT

        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaTargetFrames;
        e = msReaderPointer->collateTandemPrecursorTargetsDIA(&diaTargetFrames); ree;

        int counter = 0;
        for (const QMap<ScanNumber, ScanPoints> &frame : diaTargetFrames) {

            if (numberOfFramesToProcess > 0 && counter >= numberOfFramesToProcess) {
                break;
            }

            const ScanNumber firstScanNumber = frame.firstKey();

            ProcessFileLogicInput processFileLogicInput;

            e = msReaderPointer->getMsScanInfo(
                    firstScanNumber,
                    &processFileLogicInput.msScanInfo
            ); ree;

            processFileLogicInput.frame = frame;
            processFileLogicInput.pythiaParameters = pythiaParameters;

            processFileLogicInputs->push_back(processFileLogicInput);
            counter++;
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::preprocessDIAFramesParallel(
        const MsReaderPointer &msReaderPointer,
        int numberOfFramesToProcess,
        QVector<MsFrame> *msFrames
        ) {

    ERR_INIT

    QVector<ProcessFileLogicInput> processFileLogicInputs;
    e = buildProcessFileLogicInputs(
            msReaderPointer,
            m_pythiaParameters,
            numberOfFramesToProcess,
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

}//namespace
Err PythiaDIAWorkflow::scoreCandidatesPerFrameParallelWrite(
        const QVector<MsFrame> &msFrames,
        const QString &msDataFilePath,
        QMap<UniqueMsInfoScanKey, QString> *uniqueMsInfoScanKeyVsScoredFrameFilePaths
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
    qDebug() << "Running scoreCandidatesPerFrameWrite parallel";

    const auto scoringMatrixLogicBinder = std::bind(
            MsFrameScoretron::scoreCandidatesFrameWrite,
            std::placeholders::_1,
            m_pythiaParameters,
            msDataFilePath
    );

    QFuture<QPair<Err, QPair<UniqueMsInfoScanKey, QString>>> futures = QtConcurrent::mapped(
            processingChunks,
            scoringMatrixLogicBinder
    );
    futures.waitForFinished();

    for (const QPair<Err, QPair<UniqueMsInfoScanKey, QString>> &result : futures) {
        e = result.first; ree;
        uniqueMsInfoScanKeyVsScoredFrameFilePaths->insert(result.second.first, result.second.second);
    }
#else
    qDebug() << "Running scoreCandidatesPerFrame serial";
    for (const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk : processingChunks) {

        const QString frameID = chunk.first.uniqueMsInfoScanKey();
        if (frameID != "504979") {
            continue;
        }

        QPair<Err, QPair<UniqueMsInfoScanKey, QString>> scoreFrameTargetsResult = MsFrameScoretron::scoreCandidatesFrame(
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
