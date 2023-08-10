//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoretron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "PeakIntegratomatic.h"
#include "PSMsReader.h"

#include <QtConcurrent/QtConcurrent>

struct ScoreVectorsOutput {
    QVector<ScoredCandidate> scoredCandidates;
    UniqueMsInfoScanKey uniqueMsInfoScanKey;
    QPair<double, double> mzTargetStartStop = {-1.0, -1.0};
};


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

}//namespace
Err PythiaDIAWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    QString msDataFilePathRecalibrated; //drewholio remove this path but keep var

    e = parallelDeisotopeMsFile(
            msDataFilePath,
            m_pythiaParameters.ms2ExtractionWidthPPM,
            &msDataFilePathRecalibrated
            ); ree;

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

    QVector<FrameParallelInput> frameParallelInputsRecal;

    if (m_iRTReCalFilePath.isEmpty()) {

        e = buildParallelInput(
                m_pythiaParameters,
                msDataFilePathRecalibrated,
                m_fragLibUri,
                &frameParallelInputsRecal
        ); ree;

    }

    else {

        e = buildParallelInput(
                m_pythiaParameters,
                msDataFilePathRecalibrated,
                m_fragLibUri,
                m_iRTReCalFilePath,
                &frameParallelInputsRecal
        ); ree;

    }

    e = ErrorUtils::isNotEmpty(frameParallelInputsRecal); ree;

    QVector<ScoreVectorsOutput> frameScoreVectorOutput;
    e = buildFrameScoreVectors(
            frameParallelInputsRecal,
            &frameScoreVectorOutput
            ); ree;

//    QVector<PSMsReaderRow> psmReaderRowsRecal;
//
//    for (const ScoreVectorsOutput &svo  : frameScoreVectorOutput) {
//
//        const QMap<FrameIndex , QVector<ScoredCandidate>> &fivsc = svo.frameIndexVsScoredCandidates;
//
//        for (auto it = fivsc.begin(); it != fivsc.end(); it++) {
//
//            const FrameIndex frameIndex = it.key();
//            const QVector<ScoredCandidate> &scs = it.value();
//
//            for (const ScoredCandidate &sc : scs) {
//                PSMsReaderRow psMsReaderRow;
//
////                std::cout << sc.scanTime << " "
////                        << sc.scanNumber << " "
////                        << sc.peptideStringWithMods.toStdString() << " "
////                        << sc.isDecoy << " "
////                        << sc.frequencyPercentSum << " "
////                        << sc.charge << " " << std::endl;
//
//                psMsReaderRow.scanTime = sc.scanTime;
//                psMsReaderRow.scanNumber = sc.scanNumber;
//                psMsReaderRow.peptideStringWithMods = sc.peptideStringWithMods;
//                psMsReaderRow.isDecoy = sc.isDecoy;
//                psMsReaderRow.score = sc.frequencyPercentSum;
//                psMsReaderRow.charge = sc.charge;
//
//                psmReaderRowsRecal.push_back(psMsReaderRow);
//
//
//            }
//
//        }
//    }
//
////    e = processFrameScoreVectors(
////            frameScoreVectorOutput,
////            msDataFilePath,
////            m_pythiaParameters,
////            &psmReaderRowsRecal
////            ); ree;
//
//    const QString resultsFilePath = msDataFilePath + ".pythiaDIA";
//    e = ParquetReader::write(
//            psmReaderRowsRecal,
//            resultsFilePath
//            ); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::buildPSMResultsForCalibrationFile(
        const QVector<FrameParallelInput> &frameParallelInputs,
        QVector<ScoreVectorsOutput> *frameScoreVectorsAndExtractFilePaths
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

    QPair<Err, ScoreVectorsOutput> parallelBuildFrameScoreVectorLogic(const FrameParallelInput &fpi) {

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

       ScoreVectorsOutput output;
//       output.scoredCandidates = scoredCandidates;
//       output.uniqueMsInfoScanKey = fpi.uniqueMsInfoScanKey;
//       output.mzTargetStartStop = fpi.mzTargetStartStop;

       const QString resultsFilePath = fpi.msDataFilePath + "." +fpi.uniqueMsInfoScanKey + ".pythiaDIA";
       e = ParquetReader::write(
               scoredCandidates,
               resultsFilePath
       ); rree;

       return {e, output};
   }

}//namespace
Err PythiaDIAWorkflow::buildFrameScoreVectors(
        const QVector<FrameParallelInput> &frameParallelInputs,
        QVector<ScoreVectorsOutput> *scoreVectorsOutputs
        ) {

    ERR_INIT

#define PARALLEL_RUN_SCORE_VEC
#ifdef PARALLEL_RUN_SCORE_VEC
    QFuture<QPair<Err, ScoreVectorsOutput>> futures = QtConcurrent::mapped(
            frameParallelInputs, //.mid(20,8),
            parallelBuildFrameScoreVectorLogic
            );
    futures.waitForFinished();

    for (const QPair<Err, ScoreVectorsOutput> &result : futures) {
        e = result.first; ree;
//        scoreVectorsOutputs->push_back(result.second);
//
//        QVector<PSMsReaderRow> psmReaderRowsRecal;
//
//        const QVector<ScoredCandidate> &fivsc = result.second.scoredCandidates;
//
//        for (const ScoredCandidate &sc : fivsc) {
//
//            PSMsReaderRow psMsReaderRow;
//
//            psMsReaderRow.scanTime = sc.scanTime;
//            psMsReaderRow.scanNumber = sc.scanNumber;
//            psMsReaderRow.peptideStringWithMods = sc.peptideStringWithMods;
//            psMsReaderRow.isDecoy = sc.isDecoy;
//            psMsReaderRow.score = sc.frequencyPercentSum;
//            psMsReaderRow.charge = sc.charge;
//
//            psmReaderRowsRecal.push_back(psMsReaderRow);
//        }
//
//        const QString resultsFilePath = result.second.uniqueMsInfoScanKey + ".pythiaDIA";
//        e = ParquetReader::write(
//                psmReaderRowsRecal,
//                resultsFilePath
//        ); ree;

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

namespace {

    QPair<Err, QVector<PSMsReaderRow>> processScoreVectorsParallelLogic(
            const ScoreVectorsOutput &scoreVectorsOutput,
            const QString &msDataFilePath,
            const PythiaParameters &pythiaParameters
            ) {

        ERR_INIT

        MsFrameScoretronProcessormatic msFrameScoretronProcessormatic;
        e = msFrameScoretronProcessormatic.init(pythiaParameters); rree;

        QVector<PSMsReaderRow> psmsReaderRows;
        e = msFrameScoretronProcessormatic.processFrameScoreVectors(&psmsReaderRows); rree;

        return {e, psmsReaderRows};
    }

}//namespace
Err PythiaDIAWorkflow::processFrameScoreVectors(
        const QVector<ScoreVectorsOutput> &scoreVectorOutputs,
        const QString &msDataFilePath,
        const PythiaParameters &pythiaParameters,
        QVector<PSMsReaderRow> *psmsPreaderRows
        ) {

    ERR_INIT

#ifdef PARALLEL_RUN_SCORE_VEC
    const auto scoreVecProcessorLogicBinder = std::bind(
            processScoreVectorsParallelLogic,
            std::placeholders::_1,
            msDataFilePath,
            pythiaParameters
    );

    QFuture<QPair<Err, QVector<PSMsReaderRow>>> futures = QtConcurrent::mapped(
            scoreVectorOutputs, //.mid(20,8),
            scoreVecProcessorLogicBinder
    );
    futures.waitForFinished();

    for (const QPair<Err, QVector<PSMsReaderRow>> &result : futures) {
        e = result.first; ree;
        psmsPreaderRows->append(result.second);
    }
#else
    for (const ScoreVectorsOutput &svo : scoreVectorOutputs) {

        const QPair<Err, QVector<PSMsReaderRow>> result = processScoreVectorsParallelLogic(
                svo,
                msDataFilePath,
                pythiaParameters
                ); ree;

        e = result.first; ree;
        psmsPreaderRows->append(result.second);
    }
#endif

    ERR_RETURN
}

