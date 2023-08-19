//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MsFrameScoretron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
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
    e = ErrorUtils::fileExists(fragLibUri); ree;

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;

    ERR_RETURN
}

Err PythiaDIAWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &iRTReCalFilePath
) {

    ERR_INIT

    e = init(
            pythiaParameters,
            fragLibUri
            ); ree;

    e = ErrorUtils::fileExists(iRTReCalFilePath); ree;

    m_iRTReCalFilePath = iRTReCalFilePath;

    ERR_RETURN
}

struct FrameParallelInput {
    PythiaParameters params;
    QString msDataFilePath;
    QString fragLibFilePath;
    QString iRTReCalFilePath;
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

    Err buildParallelInput(
            const PythiaParameters &pythiaParameters,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const QString &iRTReCalFilePath,
            QVector<FrameParallelInput> *frameParallelInputs
    ) {

        ERR_INIT

        e = buildParallelInput(
                pythiaParameters,
                msDataFilePath,
                fragLibFilePath,
                frameParallelInputs
                ); ree;

        for (int i = 0; i < frameParallelInputs->size(); i++) {
            FrameParallelInput &fpi = (*frameParallelInputs)[i];
            fpi.iRTReCalFilePath = iRTReCalFilePath;
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

//        const int calibrationPoints = 3; //TODO add this to params.
//        MsCalibratomatic msCalibratomatic;
//        e = msCalibratomatic.init(
//                *pythiaParameters,
//                firstPassResultsFilePath,
//                calibrationPoints
//                ); ree;
//
//        MsReaderParquet msReaderParquet;
//        e = msReaderParquet.openFile(msDataFilePath); ree;
//
//        const QMap<ScanNumber, ScanPoints> scanPoints = msReaderParquet.getScanPoints();
//
//        QMap<ScanNumber, ScanPoints> scanPointsRecalibrated;
//        e = msCalibratomatic.recalibratePoints(
//                scanPoints,
//                &scanPointsRecalibrated
//                ); ree;
//
//        msReaderParquet.setScanPoints(scanPointsRecalibrated);
//
//        *outputFilePath = msDataFilePath + ".reCal";
//
//        e = MsReaderParquet::writeMsReaderToParquet(
//                *outputFilePath,
//                QSharedPointer<MsReaderBase>(new MsReaderBase(msReaderParquet))
//        ); ree;
//
////        const double ppmMultilplier = 4.0;
////        const double newPrecisionPPM = msCalibratomatic.newStDev() * ppmMultilplier;
////        const double oldPrecisionPPM = pythiaParameters->ms2ExtractionWidthPPM;
////        (*pythiaParameters).ms2ExtractionWidthPPM = newPrecisionPPM;
////        (*pythiaParameters).featureFinderTolerancePPM = newPrecisionPPM;
////        qDebug() << "PPM tolerance adjusted from" << oldPrecisionPPM << "->" << newPrecisionPPM;

        ERR_RETURN
    }

    QPair<Err, QPair<ScanNumber, ScanPoints>>  deisotopeLogic(
            const QPair<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
            double ppmTol
            ) {

        ERR_INIT

        const QString &chargeModelFilePath
                = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");

        const QString &monoModelFilePath
                = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");

        MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
        e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);

        ScanPoints scanPointsIterDeisotoped;
        e = ms2ChargeDeconvolvotron.deisotopeScanPoints(
                scanNumberVsScanPoints.second,
                &scanPointsIterDeisotoped
                ); rree

        return {e, {scanNumberVsScanPoints.first, scanPointsIterDeisotoped}};
    }

    Err parallelDeisotopeMsFile(
            const QString &msFilePath,
            double ppmTol,
            QString *outputFilePath
            ) {

        ERR_INIT
        e = ErrorUtils::fileExists(msFilePath); ree;

        *outputFilePath = msFilePath;
        *outputFilePath = outputFilePath->replace(".prq", ".deiso.prq");
        e = ErrorUtils::fileExists(*outputFilePath);
        if (e == eNoError) {
            ERR_RETURN
        }

        MsReaderParquet msReaderParquet;
        e = msReaderParquet.openFile(msFilePath); ree;

        const QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints = msReaderParquet.getScanPoints();
        e = ErrorUtils::isNotEmpty(scanNumberVsScanPoints); ree;

        const QVector<QPair<ScanNumber, ScanPoints>> scanNumberVsScanPointsVector
                = ParallelUtils::convertMapToVectorPairs(scanNumberVsScanPoints);

        const auto digestLogicBinder = std::bind(
                deisotopeLogic,
                std::placeholders::_1,
                ppmTol
        );

        QFuture<QPair<Err, QPair<ScanNumber, ScanPoints>>> futures = QtConcurrent::mapped(
                scanNumberVsScanPointsVector,
                digestLogicBinder
        );
        futures.waitForFinished();

        QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsDeisotoped;
        for (const QPair<Err, QPair<ScanNumber, ScanPoints>> &result : futures) {
            e = result.first;
            ree;

            const QPair<ScanNumber, ScanPoints> &res = result.second;
            scanNumberVsScanPointsDeisotoped.insert(res.first, res.second);
        }

        e = ErrorUtils::isNotEmpty(scanNumberVsScanPointsDeisotoped); ree;
        msReaderParquet.setScanPoints(scanNumberVsScanPointsDeisotoped);

        MsReaderBase msReaderBase;
        msReaderBase.setScanPoints(msReaderParquet.getScanPoints());
        msReaderBase.setMsScanInfo(msReaderParquet.getMsScanInfos());

        e = MsReaderParquet::writeMsReaderToParquet(
                *outputFilePath,
                QSharedPointer<MsReaderBase>(new MsReaderBase(msReaderBase))
                ); ree;

        ERR_RETURN
    }

    Err separateTargetWindows(
            const QString &msFilePath,
            QVector<QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>>> *frames
            ) {

        ERR_INIT

        MsReaderParquet msReaderParquet;
        e = msReaderParquet.openFile(msFilePath); ree;

        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaTargetFrame;
        e = msReaderParquet.collateTandemPrecursorTargetsDIA(&diaTargetFrame); ree

        *frames = ParallelUtils::convertMapToVectorPairs(diaTargetFrame);

        ERR_RETURN
    }


    Err parallelFeatureFind(
            const QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> &frame,
            const PythiaParameters &pythiaParameters
            ) {

        ERR_INIT

        FeatureFinderParameters ffParams;
        ffParams.tolerancePPM = pythiaParameters.ms2ExtractionWidthPPM;
        ffParams.skipScanCount = pythiaParameters.skipScanCount;
        ffParams.minScanCount = pythiaParameters.minScanCount;
        ffParams.filterLength = pythiaParameters.filterLength;
        ffParams.smoothCount = pythiaParameters.smoothCount;
        ffParams.sigma = pythiaParameters.sigma;
        ffParams.signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;

        FeatureFinderHillBuilder featureFinderHillBuilder;
        e = featureFinderHillBuilder.init(ffParams); ree;
        e = featureFinderHillBuilder.buildHills(frame.second); ree;
        e = featureFinderHillBuilder.refineHills(true); ree;

        ERR_RETURN
    }


}//namespace
Err PythiaDIAWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    const QString &msDataFilePathRecalibrated = msDataFilePath; //drewholio remove this path but keep var
    qDebug() << msDataFilePathRecalibrated;

    QVector<QPair<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>>> diaTargetFrame;
    e = separateTargetWindows(
            msDataFilePathRecalibrated,
            &diaTargetFrame
            ); ree;

    const auto featureFinderLogicBinder = std::bind(
            parallelFeatureFind,
            std::placeholders::_1,
            m_pythiaParameters
    );

    QFuture<Err> futures = QtConcurrent::mapped(
            diaTargetFrame,
            featureFinderLogicBinder
    );
    futures.waitForFinished();



//    e = parallelDeisotopeMsFile(
//            msDataFilePath,
//            m_pythiaParameters.ms2ExtractionWidthPPM,
//            &msDataFilePathRecalibrated
//            ); ree;

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

//    QVector<FrameParallelInput> frameParallelInputsRecal;
//
//    if (m_iRTReCalFilePath.isEmpty()) {
//
//        e = buildParallelInput(
//                m_pythiaParameters,
//                msDataFilePathRecalibrated,
//                m_fragLibUri,
//                &frameParallelInputsRecal
//        ); ree;
//
//    }
//
//    else {
//
//        e = buildParallelInput(
//                m_pythiaParameters,
//                msDataFilePathRecalibrated,
//                m_fragLibUri,
//                m_iRTReCalFilePath,
//                &frameParallelInputsRecal
//        ); ree;
//
//    }
//
//    e = ErrorUtils::isNotEmpty(frameParallelInputsRecal); ree;
//
//    QVector<ScoredCandidate> scoredCandidates;
//    e = buildFrameScoreVectors(
//            frameParallelInputsRecal,
//            &scoredCandidates
//            ); ree;
//
//    const QString resultsFilePath = msDataFilePath + ".pythiaDIA";
//    e = ParquetReader::write(
//            scoredCandidates,
//            resultsFilePath
//            ); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::buildPSMResultsForCalibrationFile(
        const QVector<FrameParallelInput> &frameParallelInputs,
        QVector<ScoredCandidate> *frameScoreVectorsAndExtractFilePaths
        ) {

    ERR_INIT

    const double calibrationFraction = 0.5; //TODO fix this. Find a better way.
    const int calibrationResize = static_cast<int>(std::round(frameParallelInputs.size() * calibrationFraction));

    QVector<FrameParallelInput> frameParallelInputsCalibration = frameParallelInputs;
    frameParallelInputsCalibration.resize(calibrationResize);

    e = buildFrameScoreVectors(
            frameParallelInputsCalibration,
            frameScoreVectorsAndExtractFilePaths
            ); ree;

//    qDebug() << "Rows to write:" << psmReaderRows.size(); //drewholio
//    e = ParquetReader::write(
//            psmReaderRows,
//            firstPassResultsFilePath
//            ); ree;

    ERR_RETURN
}

namespace {

    QPair<Err, QVector<ScoredCandidate>> parallelBuildFrameScoreVectorLogic(const FrameParallelInput &fpi) {

       ERR_INIT

       MsFrameScoretron msFrameScoretron;

       if (fpi.iRTReCalFilePath.isEmpty()) {
           e = msFrameScoretron.init(
                   fpi.params,
                   fpi.msDataFilePath,
                   fpi.fragLibFilePath,
                   fpi.uniqueMsInfoScanKey,
                   fpi.mzTargetStartStop
           ); rree;
       }

       else {
           e = msFrameScoretron.init(
                   fpi.params,
                   fpi.msDataFilePath,
                   fpi.fragLibFilePath,
                   fpi.iRTReCalFilePath,
                   fpi.uniqueMsInfoScanKey,
                   fpi.mzTargetStartStop
           ); rree;
       }

       QVector<ScoredCandidate> scoredCandidates;
       e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates); rree;
        
       return {e, scoredCandidates};
   }

}//namespace
Err PythiaDIAWorkflow::buildFrameScoreVectors(
        const QVector<FrameParallelInput> &frameParallelInputs,
        QVector<ScoredCandidate> *scoredCandidates
        ) {

    ERR_INIT

#define PARALLEL_RUN_SCORE_VEC
#ifdef PARALLEL_RUN_SCORE_VEC
    QFuture<QPair<Err, QVector<ScoredCandidate>>> futures = QtConcurrent::mapped(
            frameParallelInputs, //.mid(20,8),
            parallelBuildFrameScoreVectorLogic
            );
    futures.waitForFinished();

    for (const QPair<Err, QVector<ScoredCandidate>> &result : futures) {
        e = result.first; ree;
        scoredCandidates->append(result.second);
    }
#else
    for (const FrameParallelInput &fpi : frameParallelInputs) {

        if (fpi.uniqueMsInfoScanKey != "645043") {
            continue;
        }

        const QPair<Err, ScoreVectorsOutput> result = parallelBuildFrameScoreVectorLogic(fpi); ree;
        e = result.first; ree;
        scoreVectorsOutputs->push_back(result.second);
    }
#endif

    ERR_RETURN
}
