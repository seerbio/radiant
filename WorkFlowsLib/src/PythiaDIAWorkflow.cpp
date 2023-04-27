//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoretron.h"
#include "MsReaderParquet.h"
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


    ERR_RETURN
}

struct FrameParallelInput {
    PythiaParameters params;
    QString msDataFilePath;
    QString fragLibFilePath;
    UniqueMsInfoScanKey uniqueMsInfoScanKey;
    QPair<double, double> mzTargetStartStop;
};

namespace {

    Err buildParallelInput(
            const PythiaParameters &pythiaParameters,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            QVector<FrameParallelInput> *frameParallelInputs
    ) {

        ERR_INIT
        MsReaderParquet msReaderParquet;
        e = msReaderParquet.openFile(
                msDataFilePath,
                MsParquetReaderNamespace::TARGET_KEY
        ); ree;

        const QVector<MsScanInfo> uniqueScanInfos = msReaderParquet.getUniqueTandemMsScanInfos();

        for (const MsScanInfo &si : uniqueScanInfos) {

            FrameParallelInput fpi;
            fpi.msDataFilePath = msDataFilePath;
            fpi.params = pythiaParameters;
            fpi.fragLibFilePath = fragLibFilePath;
            fpi.uniqueMsInfoScanKey = si.targetScanKey();
            fpi.mzTargetStartStop
                    = {si.precursorTargetMz - si.isoWindowLower, si.precursorTargetMz + si.isoWindowUpper};

            e = ErrorUtils::isBelowThreshold(
                    fpi.mzTargetStartStop.first,
                    fpi.mzTargetStartStop.second,
                    ErrorUtilsParam::ExcludeThreshold
            ); ree;

            frameParallelInputs->push_back(fpi);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    QVector<FrameParallelInput> frameParallelInputs;
    e = buildParallelInput(
            m_pythiaParameters,
            msDataFilePath,
            m_fragLibUri,
            &frameParallelInputs
    ); ree;
    e = ErrorUtils::isNotEmpty(frameParallelInputs); ree;

    e = buildCalibrationFiles(frameParallelInputs); ree;

//    e = processDIAFramesParallel(frameParallelInputs); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::buildCalibrationFiles(const QVector<FrameParallelInput> &frameParallelInputs) {

    ERR_INIT

    const double calibrationFraction = 0.25;
    const int calibrationResize = static_cast<int>(std::round(frameParallelInputs.size() * calibrationFraction));

    QVector<FrameParallelInput> frameParallelInputsCalibration = frameParallelInputs;
    frameParallelInputsCalibration.resize(calibrationResize);

    e = processDIAFramesParallel(frameParallelInputsCalibration); ree;

    ERR_RETURN
}

namespace {

    QPair<Err, QPair<UniqueMsInfoScanKey, QString>> parallelFrameProcossingLogic(
            const FrameParallelInput &fpi
            ) {

       ERR_INIT

        MsFrameScoretron msFrameScoretron;

        QPair<Err, QPair<UniqueMsInfoScanKey, QString>> result = msFrameScoretron.scoreCandidates(
                fpi.params,
                fpi.msDataFilePath,
                fpi.fragLibFilePath,
                fpi.uniqueMsInfoScanKey,
                fpi.mzTargetStartStop
                ); rree;

        return result;
   }

}//namespace
Err PythiaDIAWorkflow::processDIAFramesParallel(const QVector<FrameParallelInput> &frameParallelInputs) {

    ERR_INIT

#define PARALLEL_RUN_SCORE_VEC
#ifdef PARALLEL_RUN_SCORE_VEC
    QFuture<QPair<Err, QPair<UniqueMsInfoScanKey, QString>>> futures = QtConcurrent::mapped(
            frameParallelInputs,
            parallelFrameProcossingLogic
            );
    futures.waitForFinished();
#else
    for (const FrameParallelInput &fpi : frameParallelInputs) {

        if (fpi.uniqueMsInfoScanKey != "474966") {
            continue;
        }

        QPair<Err, QPair<UniqueMsInfoScanKey, QString>> result
                = parallelFrameProcossingLogic(fpi); ree;

    }
#endif


    ERR_RETURN
}
