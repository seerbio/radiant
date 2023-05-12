//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoretron.h"
#include "MsReaderParquet.h"
#include "PeakIntegratomatic.h"
#include "PSMsReader.h"

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
    bool applySmooth2D = false;
};

namespace {

    Err buildParallelInput(
            const PythiaParameters &pythiaParameters,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            bool applySmooth2D,
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
            fpi.applySmooth2D = applySmooth2D;
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

    Err buildRecalibratedMsDataFile(
            const QString &msDataFilePath,
            const QString &firstPassResultsFilePath,
            PythiaParameters *pythiaParameters,
            QString *outputFilePath
            ) {

        ERR_INIT

        const int calibrationPoints = 3; //TODO add this to params.
        MsCalibratomatic msCalibratomatic;
        e = msCalibratomatic.init(
                *pythiaParameters,
                firstPassResultsFilePath,
                calibrationPoints
                ); ree;

        MsReaderParquet msReaderParquet;
        e = msReaderParquet.openFile(msDataFilePath); ree;

        const QMap<ScanNumber, ScanPoints> scanPoints = msReaderParquet.getScanPoints();

        QMap<ScanNumber, ScanPoints> scanPointsRecalibrated;
        e = msCalibratomatic.recalibratePoints(
                scanPoints,
                &scanPointsRecalibrated
                ); ree;

        msReaderParquet.setScanPoints(scanPointsRecalibrated);

        *outputFilePath = msDataFilePath + ".reCal";

        e = MsReaderParquet::writeMsReaderToParquet(
                *outputFilePath,
                QSharedPointer<MsReaderBase>(new MsReaderBase(msReaderParquet))
        ); ree;

//        const double ppmMultilplier = 4.0;
//        const double newPrecisionPPM = msCalibratomatic.newStDev() * ppmMultilplier;
//        const double oldPrecisionPPM = pythiaParameters->ms2ExtractionWidthPPM;
//        (*pythiaParameters).ms2ExtractionWidthPPM = newPrecisionPPM;
//        (*pythiaParameters).featureFinderTolerancePPM = newPrecisionPPM;
//        qDebug() << "PPM tolerance adjusted from" << oldPrecisionPPM << "->" << newPrecisionPPM;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    QString msDataFilePathRecalibrated = msDataFilePath + ".reCal"; //drewholio remove this path but keep var

//    const bool applySmooth2DCalibration = false;
//    QVector<FrameParallelInput> frameParallelInputs;
//    e = buildParallelInput(
//            m_pythiaParameters,
//            msDataFilePath,
//            m_fragLibUri,
//            applySmooth2DCalibration,
//            &frameParallelInputs
//    ); ree;
//    e = ErrorUtils::isNotEmpty(frameParallelInputs); ree;
//
//    const QString firstPassResultsFilePath = msDataFilePath + ".firstPass.pythiaDIA";
//
//    e = buildPSMResultsForCalibrationFile(
//            frameParallelInputs,
//            firstPassResultsFilePath
//            ); ree;
//
//    e = buildRecalibratedMsDataFile(
//            msDataFilePath,
//            firstPassResultsFilePath,
//            &m_pythiaParameters,
//            &msDataFilePathRecalibrated
//            ); ree;

    const bool applySmooth2D = true;
    QVector<FrameParallelInput> frameParallelInputsRecal;
    e = buildParallelInput(
            m_pythiaParameters,
            msDataFilePathRecalibrated,
            m_fragLibUri,
            applySmooth2D,
            &frameParallelInputsRecal
    ); ree;
    e = ErrorUtils::isNotEmpty(frameParallelInputsRecal); ree;

    QStringList frameScoreVectorsFilePaths;
    e = processDIAFramesParallel(
            frameParallelInputsRecal,
            &frameScoreVectorsFilePaths
            ); ree;

//    const QString resultsFilePath = msDataFilePath + ".pythiaDIA";
//    e = ParquetReader::write(
//            psmReaderRowsRecal,
//            resultsFilePath
//            ); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::buildPSMResultsForCalibrationFile(
        const QVector<FrameParallelInput> &frameParallelInputs,
        const QString &firstPassResultsFilePath
        ) {

    ERR_INIT

    const double calibrationFraction = 0.5; //TODO fix this. Find a better way.
    const int calibrationResize = static_cast<int>(std::round(frameParallelInputs.size() * calibrationFraction));

    QVector<FrameParallelInput> frameParallelInputsCalibration = frameParallelInputs;
    frameParallelInputsCalibration.resize(calibrationResize);

    QStringList frameScoreVectorsFilePaths;
    e = processDIAFramesParallel(
            frameParallelInputsCalibration,
            &frameScoreVectorsFilePaths
            ); ree;

//    qDebug() << "Rows to write:" << psmReaderRows.size(); //drewholio
//    e = ParquetReader::write(
//            psmReaderRows,
//            firstPassResultsFilePath
//            ); ree;

    ERR_RETURN
}

namespace {

    QPair<Err, QString> parallelFrameProcossingLogic(const FrameParallelInput &fpi) {

       ERR_INIT

       MsFrameScoretron msFrameScoretron;
       e = msFrameScoretron.init(
               fpi.params,
               fpi.msDataFilePath,
               fpi.fragLibFilePath,
               fpi.uniqueMsInfoScanKey,
               fpi.mzTargetStartStop
               ); rree;

       QString frameScoreVectorsFilePath;
       e = msFrameScoretron.buildFrameScoreVectors(&frameScoreVectorsFilePath); rree;

       return {e, frameScoreVectorsFilePath};
   }

}//namespace
Err PythiaDIAWorkflow::processDIAFramesParallel(
        const QVector<FrameParallelInput> &frameParallelInputs,
        QStringList *psmReaderRows
        ) {

    ERR_INIT

#define PARALLEL_RUN_SCORE_VEC
#ifdef PARALLEL_RUN_SCORE_VEC
    QFuture<QPair<Err, QString>> futures = QtConcurrent::mapped(
            frameParallelInputs, //.mid(20,8), //drewholio undu this
            parallelFrameProcossingLogic
            );
    futures.waitForFinished();

    for (const QPair<Err, QString> &result : futures) {
        e = result.first; ree;
        psmReaderRows->push_back(result.second);
    }
#else
    for (const FrameParallelInput &fpi : frameParallelInputs) {

//        if (fpi.uniqueMsInfoScanKey != "645043") {
//            continue;
//        }

        QPair<Err, QVector<PSMsReaderRow>> result
                = parallelFrameProcossingLogic(fpi); ree;

        e = result.first; ree;
        psmReaderRows->append(result.second);
    }
#endif

    ERR_RETURN
}
