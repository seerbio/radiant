//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "CandidateClassifier.h"
#include "ClassifierWeightsManager.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FastaFileToPeptidesListWorkFlow.h"
#include "FastaReader.h"
#include "FDRCLassifierNeuralNet.h"
#include "FeatureFinderHillBuilder.h"
#include "MS2ChargeDeconvolvotron.h"
#include "MS2DataExtractomatic.h"
#include "MsCalibrationReader.h"
#include "MsFrameScoretron.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "PeakIntegratomatic.h"
#include "TandemSpectraDeconvolvotron.h"

#include <QtConcurrent/QtConcurrent>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <iostream>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
using rTreeBox = bg::model::box<rTreeCoor>;
using rTreePoint = std::pair<rTreeBox, QString> ;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;


PythiaDIAWorkflow::PythiaDIAWorkflow() : m_minTopNMs2Ions(6) {}

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

    e = m_fragLibReader.init(m_fragLibUri); ree;

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

namespace{

    Err filterScoreCandidatesByFDR(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            double qValueThreshold,
            bool filterDecoys,
            QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;
        e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

        scoredCandidatesTargetsFDRThresholded->clear();

        for (const ScoredCandidate &sc : scoredCandidatesAll) {

            if (filterDecoys && sc.isDecoy) {
                continue;
            }

            else if (sc.qValue > qValueThreshold) {
                continue;
            }

            scoredCandidatesTargetsFDRThresholded->push_back(sc);
        }

        ERR_RETURN
    }

    Err outputFDRResults(const QVector<ScoredCandidate> &scoredCandidatesAll) {

        ERR_INIT

        const QVector<double> fdrFractions = {0.5, 0.2, 0.1, 0.01, 0.005};
        for (double fdrThresh : fdrFractions) {
            int foundAtThreshold;
            e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                    scoredCandidatesAll,
                    fdrThresh,
                    &foundAtThreshold
            ); ree;

            qDebug() << foundAtThreshold << "PSMs found at" << fdrThresh * 100 << "% FDR";
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::processFile(const QString &_msDataFilePath) {

    ERR_INIT

    m_msScanInfos.clear();

    QString msDataFilePath = _msDataFilePath;

#define USE_FILE_CACHING
#ifdef USE_FILE_CACHING
    {
        const QString msDataFilePathCached = msDataFilePath + S_GLOBAL_SETTINGS.DOT_CACHED_FILE_EXTENSION;
        const bool cacheFileExists = QFileInfo::exists(msDataFilePathCached);
        qDebug() << "Using cached msdatafile" << cacheFileExists;

        if (cacheFileExists) {
            msDataFilePath = msDataFilePathCached;
        }
        else {
            MsReaderParquet msReaderParquet;
            e = msReaderParquet.openFile(msDataFilePath); ree;
            e = deisotopeScans(&msReaderParquet); ree;
            e = MsReaderParquet::writeMsReaderToParquet(
                    msDataFilePathCached,
                    QSharedPointer<MsReaderBase>(new MsReaderBase(msReaderParquet))
            ); ree;
        }
    }
#endif
    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(msDataFilePath); ree;
    m_msScanInfos = msReaderParquet.getUniqueTandemMsScanInfos();
#ifndef USE_FILE_CACHING
    if (m_pythiaParameters.deisotopeScans) {
        e = deisotopeScans(&msReaderParquet); ree;
    }
#endif

    e = buildCalibration(&msReaderParquet); ree;

//#define BYPASS_OPTI
//#ifndef BYPASS_OPTI
//    e = optimizeParameters(&msReaderParquet); ree;
//#else
//    //Pythia Main Library
////    m_pythiaParameters.ms2ExtractionWidthPPM = 12.2715;
////    m_pythiaParameters.scanTimeWindowMinutes = 1.79397;
////    m_pythiaParameters.cosineSimToAnchorThreshold = 0.9;
//
//    //Entrapment libarary
//    m_pythiaParameters.ms2ExtractionWidthPPM = 11.945;
//    m_pythiaParameters.scanTimeWindowMinutes = 5.25372;
//    m_pythiaParameters.cosineSimToAnchorThreshold = 0.9;
//#endif
//
//    QVector<ScoredCandidate> scoredCandidatesTargetsFDRThresholded;
//    QVector<ScoredCandidate> scoredCandidatesAll;
//    e = mainAnalysis(
//            &msReaderParquet,
//            &scoredCandidatesTargetsFDRThresholded,
//            &scoredCandidatesAll
//            ); ree;
//
////#define WRITE_UNFILTERED_PSM
//#ifdef WRITE_UNFILTERED_PSM
//    qDebug() << scoredCandidatesAll.size() << "Candidates written";
//    e = ParquetReader::write(scoredCandidatesAll, "scoredCandidatesAll.parquet");
//#endif
//
//    QVector<ScoredCandidate> scoredCandidatesAllUpdated;
//    e = removeInterferingCandidates(
//            &msReaderParquet,
//            scoredCandidatesTargetsFDRThresholded,
//            scoredCandidatesAll,
//            &scoredCandidatesAllUpdated
//            ); ree;
//
////#define WRITE_FILTERED_PSM
//#ifdef WRITE_FILTERED_PSM
//    qDebug() << scoredCandidatesAllUpdated.size() << "Candidates written";
//    e = ParquetReader::write(scoredCandidatesAllUpdated, "scoredCandidatesAllUpdated.parquet");
//#endif
//
////#define USE_NEURAL_NET_CLASSIFIER
//#ifdef USE_NEURAL_NET_CLASSIFIER
//    QVector<ScoredCandidate> scoredCandidatesClassifierUpdated;
//    e = applyNeuralNetClassifier(
//            scoredCandidatesAllUpdated,
//            &msReaderParquet,
//            &scoredCandidatesClassifierUpdated
//            ); ree;
//
//    QVector<ScoredCandidate> scoredCandidatesAllThresholded;
//    filterScoreCandidatesByFDR(
//            scoredCandidatesClassifierUpdated,
//            0.01,
//            true,
//            &scoredCandidatesAllThresholded
//    );
//    scoredCandidatesClassifierUpdated = scoredCandidatesAllThresholded;
//
//    e = updateProteinGroupAnnotation(
//            "/home/anichols/Downloads/human_plasma_arath_entrapment.fasta", //TODO make this proper input
//            &scoredCandidatesClassifierUpdated
//            ); ree;
//
//    const QString resultsFilePath = msReaderParquet.filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
//    e = ParquetReader::write(scoredCandidatesClassifierUpdated, resultsFilePath); ree;
//#else
//    QVector<ScoredCandidate> scoredCandidatesAllThresholded;
//    filterScoreCandidatesByFDR(
//            scoredCandidatesAllUpdated,
//            0.01,
//            true,
//            &scoredCandidatesAllThresholded
//    );
//    scoredCandidatesAllUpdated = scoredCandidatesAllThresholded;
//
//    e = updateProteinGroupAnnotation(
//            "/home/anichols/Downloads/human_plasma_arath_entrapment.fasta", //TODO make this proper input
//            &scoredCandidatesAllUpdated
//            ); ree;
//
//    const QString resultsFilePath = msReaderParquet.filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
//    e = ParquetReader::write(scoredCandidatesAllUpdated, resultsFilePath); ree;
//#endif

    ERR_RETURN
}

namespace {

    QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> deisotopeLogic(const QVector<QPair<ScanNumber, ScanPoints>> &scanPointPairs, double ppmTol) {

        ERR_INIT

        const QString &chargeModelFilePath
                = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");

        const QString &monoModelFilePath
                = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");

        MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
        e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);

        QVector<QPair<ScanNumber, ScanPoints>> deisotopedScanPoints;
        for (const QPair<ScanNumber, ScanPoints> &pr : scanPointPairs) {

            ScanPoints scanPointsIterDeisotoped;
            e = ms2ChargeDeconvolvotron.deisotopeScanPoints(pr.second, &scanPointsIterDeisotoped); rree;

            deisotopedScanPoints.push_back({pr.first, scanPointsIterDeisotoped});
        }

        return {e, deisotopedScanPoints};
    }

}//namespace
Err PythiaDIAWorkflow::deisotopeScans(MsReaderParquet *msReaderParquet) {

    ERR_INIT

    qDebug() << "Deisotoping start";

    QElapsedTimer et;
    et.start();

    const auto deisotopeLogicBinder = std::bind(
            deisotopeLogic,
            std::placeholders::_1,
            m_pythiaParameters.ms2ExtractionWidthPPM
    );

    QVector<QVector<QPair<ScanNumber, ScanPoints>>> scanPointsTranched;
    ParallelUtils::trancheMapForParallelization(
            msReaderParquet->getScanPoints(),
            ParallelUtils::numberOfAvailableSystemProcessors(),
            &scanPointsTranched
    );

    QFuture<QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>>> futures
            = QtConcurrent::mapped(scanPointsTranched, deisotopeLogicBinder);
    futures.waitForFinished();

    QMap<ScanNumber, ScanPoints> deisotopedScanPoints;
    for (const QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> &res : futures) {

        e = res.first; ree;
        for (const QPair<ScanNumber, ScanPoints> &r : res.second) {
            const ScanNumber scanNumber = r.first;
            const ScanPoints &scanPoints = r.second;
            deisotopedScanPoints.insert(scanNumber, scanPoints);
        }
    }

    e = msReaderParquet->setScanPoints(deisotopedScanPoints); ree;

    qDebug() << "Scans deisotoped in" << et.elapsed() << "mSec";

    ERR_RETURN
}

namespace {

    Err buildMsCalibrationReaderRows(
            const QVector<ScoredCandidate> &scoredCandidatesFDRThresholded,
            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesFDRThresholded); ree;

        qDebug() << scoredCandidatesFDRThresholded.size() << "Found for recalibartion";

        const auto msCalibrationReaderRowsInsertLogic = [](const ScoredCandidate &sc){
            MsCalibarationReaderRow row;
            row.peptideStringWithMods = sc.peptideStringWithMods;
            row.iRTPredicted = static_cast<float>(sc.iRTPredicted);
            row.scanTime = sc.scanTime;
            row.mzSearchedVec = sc.mzSearchedVec;
            row.mzFoundMeanVec = sc.mzFoundMeanVec;
            row.mzFoundStDevVec = sc.mzFoundStDevVec;
            return row;
        };

        ;
        std::transform(
                scoredCandidatesFDRThresholded.begin(),
                scoredCandidatesFDRThresholded.end(),
                std::back_inserter(*msCalibrationReaderRows),
                msCalibrationReaderRowsInsertLogic
        );

        ERR_RETURN
    }


}//namespace
Err PythiaDIAWorkflow::buildCalibration(MsReaderParquet *msReaderParquet) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_fragLibReader.libarySize() > 0); ree;
    const double fdrThreshold = 0.01;

    const double calibrationSelectionFractionIncrement = 0.2;
    double calibrationSelectionFraction = calibrationSelectionFractionIncrement;
    const int topNMs2IonsCalibration = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
    );

    qDebug() << "Using top:" << topNMs2IonsCalibration << "fragments for calibration";

    const bool useExtendedScores = false;
    const bool useNeuralNetworkScores = false;
    QVector<ScoredCandidate> scoredCandidatesAll;
    QVector<ScoredCandidate> scoredCandidatesTargetsFDRThresholded;

    int onePercentFDRCount = 0;
    const double maxTrainingFraction = 1.1;
    const int minTrainingCount = 100;

    MS2DataExtractomatic ms2DataExtractomatic;
    e = ms2DataExtractomatic.init(
            m_pythiaParameters,
            topNMs2IonsCalibration,
            useExtendedScores,
            useNeuralNetworkScores,
            msReaderParquet
            ); ree;

    while (calibrationSelectionFraction < maxTrainingFraction && onePercentFDRCount < minTrainingCount) {

        QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptideCalibration;
        e = buildCandidates(
                topNMs2IonsCalibration,
                calibrationSelectionFraction,
                &uniqueInfoScanKeyVsCandidatePeptideCalibration
        ); ree;

        e = ms2DataExtractomatic.extractMS2ForCandidates(
                uniqueInfoScanKeyVsCandidatePeptideCalibration,
                fdrThreshold,
                &scoredCandidatesAll,
                &scoredCandidatesTargetsFDRThresholded
        ); ree;

        FDRCLassifierNeuralNet::countScoreCandidatesByFDR(scoredCandidatesAll, 0.01, &onePercentFDRCount);

        calibrationSelectionFraction += calibrationSelectionFractionIncrement;
    }

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = buildMsCalibrationReaderRows(
            scoredCandidatesTargetsFDRThresholded,
            &msCalibrationReaderRows
            ); ree;

//#define WRITE_CALIBRATION
#ifdef WRITE_CALIBRATION
    const QString resultsFilePath = msReaderParquet->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_CAL_FILE_EXTENSION;
    e = ParquetReader::write(msCalibrationReaderRows, resultsFilePath); ree;
    e = m_msCalibratomatic.init(resultsFilePath); ree;
#else
    e = m_msCalibratomatic.init(msCalibrationReaderRows); ree;
#endif

    ERR_RETURN
}

namespace {

    RTree loadScanInfoToRTree(const QVector<MsScanInfo> &msScanInfos) {

        ERR_INIT

        std::vector<rTreePoint> cloudLoader;

        std::transform(
                msScanInfos.begin(),
                msScanInfos.end(),
                std::back_inserter(cloudLoader),
                [](const MsScanInfo &msScanInfo){

                    rTreeBox box(
                            {msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower, 0.0},
                            {msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper, 0.0}
                    );

                    return rTreePoint(box, msScanInfo.targetScanKey());
                }
        );

        const int maxElements = 16;
        return RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    }

    int calculateuniqueInfoScanKeyVsCandidatePeptideSize(
            const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptide
    ){

        const auto accumulateLogic = [](int sum, const QMap<PeptideStringWithMods, CandidatePeptide> &c){return sum + c.size();};

        return std::accumulate(
                uniqueInfoScanKeyVsCandidatePeptide.begin(),
                uniqueInfoScanKeyVsCandidatePeptide.end(),
                0,
                accumulateLogic
        );
    }

} //namespace
Err PythiaDIAWorkflow::buildCandidates(
        int topNMs2Ions,
        double selectionListFraction,
        QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide = nullptr
) {

    ERR_INIT

    uniqueInfoScanKeyVsCandidatePeptide->clear();

    QMap<PeptideSequenceChargeKey, CandidatePeptide> peptideSequenceChargeKeyVsCandidatePeptide;

    QMap<Index, bool> selectionList;
    if (selectionListFraction > 0) {
        const int selectionCount
            = static_cast<int>(std::round(m_fragLibReader.libarySize() * selectionListFraction));

        selectionList = MathUtils::generateRandomSelectionList(m_fragLibReader.libarySize(), selectionCount);

        e = m_fragLibReader.getMS2IonsTopN(
                selectionList,
                topNMs2Ions,
                m_pythiaParameters.mzMinDataStructure,
                m_pythiaParameters.mzMaxDataStructure,
                &peptideSequenceChargeKeyVsCandidatePeptide
        ); ree;
    }

    else {
        e = m_fragLibReader.getMS2IonsTopN(
                topNMs2Ions,
                m_pythiaParameters.mzMinDataStructure,
                m_pythiaParameters.mzMaxDataStructure,
                &peptideSequenceChargeKeyVsCandidatePeptide
        ); ree;

    }

    const RTree rtree = loadScanInfoToRTree(m_msScanInfos);

    for (const CandidatePeptide &candidatePeptide : peptideSequenceChargeKeyVsCandidatePeptide) {

        const double mz = BiophysicalCalcs::calculateThomsonFromMass(candidatePeptide.mass, candidatePeptide.charge);

        const rTreeBox queryBox(
                rTreeCoor(mz, 0.0),
                rTreeCoor(mz, 0.0)
        );

        std::vector<rTreePoint> rTreeSearchResult;
        rtree.query(bgi::intersects(queryBox), std::back_inserter(rTreeSearchResult));

        for (const rTreePoint &rtp : rTreeSearchResult) {
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey = rtp.second;

            (*uniqueInfoScanKeyVsCandidatePeptide)[uniqueMsInfoScanKey].insert(candidatePeptide.peptideStringWithMods, candidatePeptide);
        }

    }

    qDebug() << "Candidates size" << calculateuniqueInfoScanKeyVsCandidatePeptideSize(*uniqueInfoScanKeyVsCandidatePeptide);

    ERR_RETURN
}

Err PythiaDIAWorkflow::buildCandidates(
        const QVector<PeptideStringWithMods> &inclusionList,
        int topNMs2Ions,
        QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide = nullptr
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(inclusionList); ree;

    uniqueInfoScanKeyVsCandidatePeptide->clear();

    QMap<PeptideStringWithMods, bool> inclusionMap;
    for (const PeptideStringWithMods &peptideStringWithMods : inclusionList) {
        inclusionMap.insert(peptideStringWithMods, true);
    }

    QMap<PeptideSequenceChargeKey, CandidatePeptide> peptideSequenceChargeKeyVsCandidatePeptide;

    e = m_fragLibReader.getMS2IonsTopN(
            topNMs2Ions,
            m_pythiaParameters.mzMinDataStructure,
            m_pythiaParameters.mzMaxDataStructure,
            &peptideSequenceChargeKeyVsCandidatePeptide
    ); ree;

    e = ErrorUtils::isNotEmpty(m_msScanInfos); ree;
    const RTree rtree = loadScanInfoToRTree(m_msScanInfos);

    for (const CandidatePeptide &candidatePeptide : peptideSequenceChargeKeyVsCandidatePeptide) {

        const double mz = BiophysicalCalcs::calculateThomsonFromMass(candidatePeptide.mass, candidatePeptide.charge);

        const rTreeBox queryBox(
                rTreeCoor(mz, 0.0),
                rTreeCoor(mz, 0.0)
        );

        std::vector<rTreePoint> rTreeSearchResult;
        rtree.query(bgi::intersects(queryBox), std::back_inserter(rTreeSearchResult));

        for (const rTreePoint &rtp : rTreeSearchResult) {
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey = rtp.second;

            if (!inclusionMap.value(candidatePeptide.peptideStringWithMods)) {
                continue;
            }

            (*uniqueInfoScanKeyVsCandidatePeptide)[uniqueMsInfoScanKey].insert(candidatePeptide.peptideStringWithMods, candidatePeptide);
        }

    }

    qDebug() << "Candidates size" << calculateuniqueInfoScanKeyVsCandidatePeptideSize(*uniqueInfoScanKeyVsCandidatePeptide);

    ERR_RETURN
}

namespace {

    Err buildDOE(
            const PythiaParameters &pythiaParameters,
            double mzPPMStDev,
            double scanTimeStDev,
            QVector<PythiaParameters> *pythiaParametersExperiments
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(mzPPMStDev > 0.0); ree;
        e = ErrorUtils::isTrue(scanTimeStDev > 0.0); ree;
        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

        const QVector<QVector<double>> experiments = {
                {1.5, 1.0,  2.0},
                {3.5, 1.0,  2.0},
                {1.5,  3.0,  2.0},
                {3.5,  3.0,  2.0},
                {1.5,  2.0, 1.0},
                {3.5,  2.0, 1.0},
                {1.5,  2.0,  3.0},
                {3.5,  2.0,  3.0},
                {2.5, 1.0, 1.0},
                {2.5,  3.0, 1.0},
                {2.5, 1.0,  3.0},
                {2.5,  3.0,  3.0},
                {2.5,  2.0,  2.0}
        };

        for (const QVector<double> &exp : experiments) {

            PythiaParameters params = pythiaParameters;
            params.ms2ExtractionWidthPPM = mzPPMStDev * exp.at(0);
            params.scanTimeWindowMinutes = scanTimeStDev * exp.at(1);

            switch (static_cast<int>(exp.at(2))) {
                case 1:
                    params.cosineSimToAnchorThreshold = 0.9;
                    break;
                case 2:
                    params.cosineSimToAnchorThreshold = 0.935;
                    break;
                case 3:
                    params.cosineSimToAnchorThreshold = 0.97;
                    break;
                default:
                    params.cosineSimToAnchorThreshold = pythiaParameters.cosineSimToAnchorThreshold;
            }

            pythiaParametersExperiments->push_back(params);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::optimizeParameters(MsReaderParquet *msReaderParquet) {

    ERR_INIT

//    const int topNMs2IonsOptimization = std::max(
//            m_minTopNMs2Ions,
//            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
//    );
//
//    qDebug() << "Using top:" << topNMs2IonsOptimization << "fragments for optimization";
//
//    const double selectionFractionValue = 0.1;
//    const double fdrThreshold = 0.01;
//
//    QVector<PythiaParameters> pythiaParametersExperiments;
//    e = buildDOE(
//            m_pythiaParameters,
//            m_msCalibratomatic.mzStDev(),
//            m_msCalibratomatic.scanTimeStDev(),
//            &pythiaParametersExperiments
//            ); ree;
//
//    struct DOEResult {
//        double mzStDev = -1.0;
//        double scanTimeStDev = -1.0;
//        double cosineSimAnchor = -1.0;
//        int fdrCount = -1;
//    };
//
//    const bool useExtendedScores = true;
//    const bool useNeuralNetworkScores = false;
//
//    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptideCalibration;
//    e = buildCandidates(
//            topNMs2IonsOptimization,
//            selectionFractionValue,
//            &uniqueInfoScanKeyVsCandidatePeptideCalibration
//    ); ree;
//
//    QVector<DOEResult> results;
//    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {
//
//        QVector<ScoredCandidate> scoredCandidatesOptimization;
//        e = extractTargetDecoyData(
//                uniqueInfoScanKeyVsCandidatePeptideCalibration,
//                pythiaParams,
//                topNMs2IonsOptimization,
//                msReaderParquet,
//                &scoredCandidatesOptimization
//        ); ree;
//
//        QVector<ScoredCandidate> scoredCandidatesAll;
//        e = buildScoredCandidatesFDR(
//                scoredCandidatesOptimization,
//                useExtendedScores,
//                useNeuralNetworkScores,
//                topNMs2IonsOptimization,
//                msReaderParquet->scanTimeMinMax(),
//                &scoredCandidatesAll
//        ); ree;
//
//       int targetCountAboveFDRQValueThreshold;
//        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
//                scoredCandidatesAll,
//                fdrThreshold,
//                &targetCountAboveFDRQValueThreshold
//                ); ree;
//
//        e = outputFDRResults(scoredCandidatesAll); ree;
//
//        DOEResult res;
//        res.mzStDev = pythiaParams.ms2ExtractionWidthPPM;
//        res.scanTimeStDev = pythiaParams.scanTimeWindowMinutes;
//        res.cosineSimAnchor = pythiaParams.cosineSimToAnchorThreshold;
//        res.fdrCount = targetCountAboveFDRQValueThreshold;
//        results.push_back(res);
//    }
//
//    std::sort(
//            results.rbegin(),
//            results.rend(),
//            [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
//            );
//
//    for (const DOEResult &r : results) {
//        qDebug() << r.mzStDev << r.scanTimeStDev << r.cosineSimAnchor << r.fdrCount;
//    }
//
//    //TODO replace this with response surface derived DOE parameters when you figure out how to do it.
//    const DOEResult bestParametersFDR = *std::max_element(
//            results.rbegin(),
//            results.rend(),
//            [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
//    );
//
//    m_pythiaParameters.ms2ExtractionWidthPPM = bestParametersFDR.mzStDev;
//    m_pythiaParameters.scanTimeWindowMinutes = bestParametersFDR.scanTimeStDev;
//    m_pythiaParameters.cosineSimToAnchorThreshold = bestParametersFDR.cosineSimAnchor;
//    ///////////////////////////////////////////////////////////////////////////////////////////////////

    ERR_RETURN
}

Err PythiaDIAWorkflow::mainAnalysis(
        MsReaderParquet *msReaderParquet,
        QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded,
        QVector<ScoredCandidate> *scoredCandidatesAll
        ) {

    ERR_INIT

//    m_pythiaParameters.print();
//
//    const double selectionFractionBypassValue = -1.0;
//    const double fdrThreshold = 0.2;
//
//    const int topNMs2IonsMainAnalysis = std::max(
//            m_minTopNMs2Ions,
//            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
//    );
//
//    qDebug() << "Using top:" << topNMs2IonsMainAnalysis << "fragments for main analysis";
//
//    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptides;
//    e = buildCandidates(
//            topNMs2IonsMainAnalysis,
//            selectionFractionBypassValue,
//            &uniqueInfoScanKeyVsCandidatePeptides
//    ); ree;
//
//    const bool useExtendedScores = true;
//    const bool useNeuralNetworkScores = false;
//
//    e = extractionLoopLogic(
//            uniqueInfoScanKeyVsCandidatePeptides,
//            fdrThreshold,
//            useExtendedScores,
//            useNeuralNetworkScores,
//            topNMs2IonsMainAnalysis,
//            msReaderParquet,
//            scoredCandidatesAll,
//            scoredCandidatesTargetsFDRThresholded
//    ); ree;
//
//    int psmCountTenPercentFDR;
//    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
//            *scoredCandidatesTargetsFDRThresholded,
//            0.1,
//            &psmCountTenPercentFDR
//            ); ree;

    ERR_RETURN
}

namespace {

    Err buildPeptideStringWithModsVsScoreCandidatesDecoys(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            QMap<PeptideStringWithMods, ScoredCandidate> *peptideStringWithModsVsScoreCandidatesDecoys
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;

        for (const ScoredCandidate &sc : scoredCandidatesAll) {

            if (!sc.isDecoy) {
                continue;
            }

            const QStringList splitPepStrWithMods = sc.peptideStringWithMods.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP);
            e = ErrorUtils::isEqual(splitPepStrWithMods.size(), 2); ree;

            const PeptideStringWithMods &peptideStringWithMods = splitPepStrWithMods.front();
            peptideStringWithModsVsScoreCandidatesDecoys->insert(peptideStringWithMods, sc);

        }

        e = ErrorUtils::isNotEmpty(*peptideStringWithModsVsScoreCandidatesDecoys); ree

        ERR_RETURN
    }

    Err buildMs2IonsFromScoredCandidate(
            const ScoredCandidate &scoredCandidate,
            QVector<MS2Ion> *ms2Ions
    ) {

        ERR_INIT

        if (scoredCandidate.theoIntensityVec.isEmpty()) {
            qDebug() << scoredCandidate.peptideStringWithMods;
            qDebug() << scoredCandidate.mzSearchedVec;
            qDebug() << scoredCandidate.intensityFoundMaxVec;
        }

        e = ErrorUtils::isNotEmpty(scoredCandidate.theoIntensityVec); ree;
        e = ErrorUtils::isEqual(scoredCandidate.theoIntensityVec.size(), scoredCandidate.mzSearchedVec.size()); ree;

        for (int i = 0; i < scoredCandidate.theoIntensityVec.size(); i++) {
            const double intensity = scoredCandidate.theoIntensityVec.at(i);
            const double mz = scoredCandidate.mzSearchedVec.at(i);

            if (MathUtils::tZero(mz) || MathUtils::tZero(intensity)) {
                continue;
            }

            MS2Ion ms2Ion;
            ms2Ion.mz = mz;
            ms2Ion.intensity = intensity;

            ms2Ions->push_back(ms2Ion);
        }

        ERR_RETURN
    }

    Err buildTandemDeconvolutionInput(
            const QVector<ScoredCandidate> &scoredCandidatesTargetsFDRThresholded,
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            QMap<ScanNumber, QMap<PeptideStringWithMods, QVector<MS2Ion>>> *scanNumberVsTandemPredictions
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesTargetsFDRThresholded); ree;
        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;

        QMap<PeptideStringWithMods, ScoredCandidate> peptideStringWithModsVsScoreCandidatesDecoys;
        e = buildPeptideStringWithModsVsScoreCandidatesDecoys(
                scoredCandidatesAll,
                &peptideStringWithModsVsScoreCandidatesDecoys
        ); ree;

        for (const ScoredCandidate &sc : scoredCandidatesTargetsFDRThresholded) {
            e = ErrorUtils::isTrue(peptideStringWithModsVsScoreCandidatesDecoys.contains(sc.peptideStringWithMods)); ree;

            const ScoredCandidate &scoredCandidateDecoy
                    = peptideStringWithModsVsScoreCandidatesDecoys.value(sc.peptideStringWithMods);

            QVector<MS2Ion> ms2IonsTarget;
            e = buildMs2IonsFromScoredCandidate(sc, &ms2IonsTarget); ree;
            (*scanNumberVsTandemPredictions)[sc.scanNumber].insert(sc.peptideStringWithMods, ms2IonsTarget);

            if(scoredCandidateDecoy.mzSearchedVec.isEmpty()) {
                continue;
            }

            QVector<MS2Ion> ms2IonsDecoy;
            e = buildMs2IonsFromScoredCandidate(scoredCandidateDecoy, &ms2IonsDecoy); ree;
            (*scanNumberVsTandemPredictions)[sc.scanNumber].insert(scoredCandidateDecoy.peptideStringWithMods, ms2IonsDecoy);
        }

        ERR_RETURN
    }

    struct DeconVol {
        ScanNumber scanNumber = -1;
        ScanPoints scanPoints;
        QMap<PeptideStringWithMods, QVector<MS2Ion>> tandemPredictions;
    };

    struct DeconResult {
        Err e = eNoError;
        ScanNumber scanNumber = -1;
        QMap<PeptideStringWithMods, TandemDeconvolverResult> tandemDeconvolverResult;
    };

    DeconResult tandemDeconvolutionLogic(
            const DeconVol &deconVol,
            const PythiaParameters &params
    ) {

        ERR_INIT

        const int maxIteration= 20;
        const double stopTol = 0.000000001;
        const double pValThreshold = 0.05;

        DeconResult deconResult;

        TandemSpectraDeconvolvotron tandemSpectraDeconvolvotron;
        e = tandemSpectraDeconvolvotron.init(
                S_GLOBAL_SETTINGS.HASHING_PRECISION,
                params.mzMaxDataStructure,
                params.ms2ExtractionWidthPPM,
                maxIteration,
                stopTol,
                pValThreshold
        );
        if (e != eNoError){
            return deconResult;
        }

        QMap<PeptideStringWithMods, TandemDeconvolverResult> pepSeqVsWeight;
        e = tandemSpectraDeconvolvotron.deconvolveTandemSpectra(
                deconVol.scanPoints,
                deconVol.tandemPredictions,
                &pepSeqVsWeight
        );
        if (e != eNoError){
            return deconResult;
        }

        deconResult.e = e;
        deconResult.scanNumber = deconVol.scanNumber;
        deconResult.tandemDeconvolverResult = pepSeqVsWeight;

        return deconResult;
    }

    Err deconvolveTandemSpectra(
            const QMap<ScanNumber, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &scanNumberVsTandemPredictions,
            const PythiaParameters &params,
            MsReaderParquet *msReaderParquet,
            QMap<ScanNumber, QMap<PeptideStringWithMods, TandemDeconvolverResult>> *scanNumberVsTandemDeconvolverResult
    ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(scanNumberVsTandemPredictions); ree;

        QElapsedTimer et;
        et.start();

        const auto insertLogic = [msReaderParquet, scanNumberVsTandemPredictions](const ScanNumber sn){

            DeconVol deconVol;
            deconVol.scanNumber = sn;
            deconVol.tandemPredictions = scanNumberVsTandemPredictions.value(sn);
            msReaderParquet->getScanPoints(sn, &deconVol.scanPoints);

            return deconVol;
        };

        QVector<DeconVol> scanPointsVsTandemPredictions;
        const QList<ScanNumber> &scanNumbers = scanNumberVsTandemPredictions.keys();
        std::transform(
                scanNumbers.begin(),
                scanNumbers.end(),
                std::back_inserter(scanPointsVsTandemPredictions),
                insertLogic
        );

#define PARALLEL_DECONVOLVE
#ifdef PARALLEL_DECONVOLVE
        const auto deconvolutionLogicBinder = std::bind(
            tandemDeconvolutionLogic,
            std::placeholders::_1,
            params
        );

        QFuture<DeconResult> futures = QtConcurrent::mapped(
                scanPointsVsTandemPredictions,
                deconvolutionLogicBinder
                );
        futures.waitForFinished();

        for (const DeconResult &result : futures) {
            e = result.e; ree;
            (*scanNumberVsTandemDeconvolverResult)[result.scanNumber] = result.tandemDeconvolverResult;
        }
#else
        for (const DeconVol &deconVol : scanPointsVsTandemPredictions) {

            const DeconResult result = tandemDeconvolutionLogic(
                    deconVol,
                    params
            );

            e = result.e; ree;
            (*scanNumberVsTandemDeconvolverResult)[result.scanNumber] = result.tandemDeconvolverResult;
        }
#endif

        qDebug() << "Deconvolution finished in" << et.elapsed() << "mSec";

        ERR_RETURN
    }

    void filterScoredCandidatesByWeightAndPVal(
            QVector<ScoredCandidate> *scoredCandidatesAllUpdated,
            double pValThreshold) {

        const double weightThreshold = 0.0;
        const auto terminatorLogic = [weightThreshold, pValThreshold](const ScoredCandidate &s){
            return s.matrixWeight < weightThreshold || s.matrixPValue > pValThreshold;
        };

        const auto terminator = std::remove_if(
                scoredCandidatesAllUpdated->begin(),
                scoredCandidatesAllUpdated->end(),
                terminatorLogic
        );
        scoredCandidatesAllUpdated->erase(terminator, scoredCandidatesAllUpdated->end());
    }

}//namespace
Err PythiaDIAWorkflow::removeInterferingCandidates(
        MsReaderParquet *msReaderParquet,
        const QVector<ScoredCandidate> &scoredCandidatesTargetsFDRThresholded,
        const QVector<ScoredCandidate> &scoredCandidatesAll,
        QVector<ScoredCandidate> *scoredCandidatesAllUpdated
        ) {

    ERR_INIT

    QMap<ScanNumber, QMap<PeptideStringWithMods, QVector<MS2Ion>>> scanNumberVsTandemPredictions;
    e = buildTandemDeconvolutionInput(
            scoredCandidatesTargetsFDRThresholded,
            scoredCandidatesAll,
            &scanNumberVsTandemPredictions
    ); ree;

    QMap<ScanNumber, QMap<PeptideStringWithMods, TandemDeconvolverResult>> scanNumberVsTandemDeconvolverResult;
    e = deconvolveTandemSpectra(
            scanNumberVsTandemPredictions,
            m_pythiaParameters,
            msReaderParquet,
            &scanNumberVsTandemDeconvolverResult
    ); ree;

    for (const ScoredCandidate &sc : scoredCandidatesAll) {

        const QMap<PeptideStringWithMods, TandemDeconvolverResult> &tandemResult
                = scanNumberVsTandemDeconvolverResult.value(sc.scanNumber);

        const TandemDeconvolverResult &tandemDeconvolverResult = tandemResult.value(sc.peptideStringWithMods);

        ScoredCandidate scNew = sc;
        scNew.matrixWeight = tandemDeconvolverResult.discScore;
        scNew.matrixPValue = tandemDeconvolverResult.pVal;
        scNew.matrixError = tandemDeconvolverResult.frameError;
        scNew.scanNumberCandidateCount = tandemDeconvolverResult.scanNumberCandidateCount;

        scoredCandidatesAllUpdated->push_back(scNew);
    }

    filterScoredCandidatesByWeightAndPVal(
            scoredCandidatesAllUpdated,
            m_pythiaParameters.pValThreshold
            );

    e = outputFDRResults(*scoredCandidatesAllUpdated); ree;

    ERR_RETURN
}

namespace {

    Err buildKeyVsScoredCandidateCulled(
            const QVector<ScoredCandidate> &scoredCandidatesCulled,
            QMap<QString, ScoredCandidate> *keyVsScoredCandidateCulled
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesCulled); ree;

        for (const ScoredCandidate &sc : scoredCandidatesCulled) {
            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    sc.peptideStringWithMods,
                    sc.targetKey,
                    sc.charge
            );
            keyVsScoredCandidateCulled->insert(key, sc);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::applyNeuralNetClassifier(
        const QVector<ScoredCandidate> &scoredCandidatesCulled,
        MsReaderParquet *msReaderParquet,
        QVector<ScoredCandidate> *scoredCandidatesClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scoredCandidatesCulled); ree;

    QVector<ScoredCandidate> scoredCandidatesAllFullFragIons;
    e = returnAllCandidatesScoredFullFragIons(
            msReaderParquet,
            &scoredCandidatesAllFullFragIons
            ); ree;

    QMap<QString, ScoredCandidate> keyVsScoredCandidateCulled;
    e = buildKeyVsScoredCandidateCulled(
            scoredCandidatesCulled,
            &keyVsScoredCandidateCulled
            ); ree;

#define WRITE_NEURAL_NET_INPUT
#ifdef WRITE_NEURAL_NET_INPUT
      e = ParquetReader::write(
              scoredCandidatesAllFullFragIons,
              "scoredCandidatesAllFullFragIons.parquet"
              ); ree;

    e = ParquetReader::write(
            keyVsScoredCandidateCulled.values().toVector(),
            "keyVsScoredCandidateCulled.parquet"
    ); ree;
#endif

    const int epochs = 20;
    const int baggingSize = 6;
    const double batchFraction = 0.01;
    const double learningRate = 0.001;
    FDRCLassifierNeuralNet fdrClassifierNeuralNet;
    e = fdrClassifierNeuralNet.init(
            m_pythiaParameters.topNMs2Ions,
            epochs,
            baggingSize,
            batchFraction,
            learningRate,
            msReaderParquet->scanTimeMinMax()
            ); ree;

    e = fdrClassifierNeuralNet.exec(
            keyVsScoredCandidateCulled,
            scoredCandidatesAllFullFragIons,
            scoredCandidatesClassifier
            ); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::returnAllCandidatesScoredFullFragIons(
        MsReaderParquet *msReaderParquet,
        QVector<ScoredCandidate> *scoredCandidatesAllTemp
        ) {

    ERR_INIT

//    e = ErrorUtils::isTrue(m_fragLibReader.libarySize() > 0); ree;
//
//    const double selectionFractionBypassValue = -1.0;
//    const int topNMs2IonsFull = std::max(
//            m_minTopNMs2Ions,
//            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions))
//    );
//    qDebug() << "Using top:" << topNMs2IonsFull << "fragments for calibration";
//
//    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptides;
//    e = buildCandidates(
//            topNMs2IonsFull,
//            selectionFractionBypassValue,
//            &uniqueInfoScanKeyVsCandidatePeptides
//    ); ree;
//
//    const bool useExtendedScores = true;
//    const bool useNeuralNetworkScores = true;
//    const double fdrThreshold = 1.0;
//
//    QVector<ScoredCandidate> scoredCandidatesTargetsFDRThresholdedUnused;
//    e = extractionLoopLogic(
//            uniqueInfoScanKeyVsCandidatePeptides,
//            fdrThreshold,
//            useExtendedScores,
//            useNeuralNetworkScores,
//            topNMs2IonsFull,
//            msReaderParquet,
//            scoredCandidatesAllTemp,
//            &scoredCandidatesTargetsFDRThresholdedUnused
//    ); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::updateProteinGroupAnnotation(
        const QString &fastaFilePath,
        QVector<ScoredCandidate> *scoredCandidates
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(*scoredCandidates); ree;

    FastaReader fastaReader;
    e = fastaReader.parseFastaFile(fastaFilePath); ree;

    QVector<PeptideSequence> peptideSequences;
    QMap<PeptideStringWithMods, QVector<FastaEntry>> peptideStringWithModsVsFastaEntries;
    e = FastaFileToPeptidesListWorkFlow::digestFastaEntries(
            m_pythiaParameters,
            fastaReader.fastaEntries(),
            &peptideSequences,
            &peptideStringWithModsVsFastaEntries
            ); ree;

    for (int i = 0; i < scoredCandidates->size(); i++) {

        ScoredCandidate &sc = (*scoredCandidates)[i];

        QString peptideSeqReplacedLeucines = sc.peptideStringWithMods;
        peptideSeqReplacedLeucines = peptideSeqReplacedLeucines.replace('L', 'J').replace('I', 'J');

        const QVector<FastaEntry> fastaEntries = peptideStringWithModsVsFastaEntries.value(peptideSeqReplacedLeucines);

        QStringList fastaDescriptions;
        std::transform(
                fastaEntries.begin(),
                fastaEntries.end(),
                std::back_inserter(fastaDescriptions),
                [](const FastaEntry &fe){return fe.fastaDescription;}
                );

        sc.proteinGroup = fastaDescriptions.join(';');
    }

    ERR_RETURN
}
