//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "ErrorUtils.h"
#include "FragLibraryTronDIA.h"
#include "MsReaderPointerFactory.h"
#include "ParallelUtils.h"

#include <QtConcurrent/QtConcurrent>


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

    ERR_RETURN

}


Err PythiaDIAWorkflow::processFile(const QString &msDatalFilePath) {

    ERR_INIT

    QPair<Err, MsReaderPointer> msReaderPointerResult
            = MsReaderPointerFactory::createInstance(msDatalFilePath);
    e = msReaderPointerResult.first; ree;
    MsReaderPointer msReaderPointer = msReaderPointerResult.second;

    QVector<MsFrame> msFrames;
    e = preprocessDIAFrames(
            msReaderPointer,
            &msFrames
            ); ree;

    e = scoreCandidatesPerFrame(msFrames); ree;

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
Err PythiaDIAWorkflow::preprocessDIAFrames(
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

        e = pr.first; ree;
        msFrames->push_back(pr.second);
    }

    ERR_RETURN
}


namespace {

    Err processFrame(
            const MsFrame &frame
            ) {

        ERR_INIT

        qDebug() << frame.precursorMzTargetStartEnd();

        ERR_RETURN
    }
}
Err PythiaDIAWorkflow::scoreCandidatesPerFrame(const QVector<MsFrame> &msFrames) {
    ERR_INIT

    FragLibraryTronDIA fragLibraryTronDia;
    e = fragLibraryTronDia.init(
            m_pythiaParameters,
            m_fragLibUri
            ); ree;

    for (const MsFrame &frame : msFrames) {

        const QPair<double, double> mzTargetStartEnd = frame.precursorMzTargetStartEnd();

        qDebug() << "Processing window" << mzTargetStartEnd.first << mzTargetStartEnd.second;

        QHash<PeptideStringWithMods, QVector<MS2Ion>> peptideStringWithModsVsMS2Ions;

        const int topNMs2Ions = 12;

        e = fragLibraryTronDia.getMS2Ions(
                mzTargetStartEnd.first,
                mzTargetStartEnd.second,
                topNMs2Ions,
                &peptideStringWithModsVsMS2Ions
                ); ree;

        processFrame(frame);
    }

    ERR_RETURN
}

