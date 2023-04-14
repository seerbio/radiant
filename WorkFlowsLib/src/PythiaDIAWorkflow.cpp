//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoretron.h"
#include "MsFrameScoreVectorReader.h"
#include "MsReaderPointerFactory.h"
#include "PeakIntegratomatic.h"
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

    MsFrameScoretron msFrameScoretron;
    e = msFrameScoretron.init(m_pythiaParameters); ree;

    e = msFrameScoretron.scoreCandidatesPerFrameParallel(
            msFrames,
            msDatalFilePath,
            &m_fragLibraryTronDia
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

