//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MS2ChargeDeconvolvotron.h"
#include "MsFrameScoretron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "PeakIntegratomatic.h"
#include "PSMsReader.h"
#include "ClassifierWeightsManager.h"

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

Err PythiaDIAWorkflow::processFile(const QString &_msDataFilePath) {

    ERR_INIT

    m_msScanInfos.clear();

    QString msDataFilePath = _msDataFilePath;

#define USE_FILE_CACHING
#ifdef USE_FILE_CACHING
    {
        const QString msDataFilePathCached = msDataFilePath + ".cached";
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

#ifndef USE_FILE_CACHING
    e = deisotopeScans(&msReaderParquet); ree;
#endif

    e = buildCalibration(&msReaderParquet); ree;










//    const QString resultsFilePath = msDataFilePath + ".pythiaDIA";
//    e = ParquetReader::write(scoredCandidatesCalibration, resultsFilePath); ree;

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

    QElapsedTimer et;
    et.start();

    const auto deisotopeLogicBinder = std::bind(
            deisotopeLogic,
            std::placeholders::_1,
            m_pythiaParameters.ms2ExtractionWidthPPM
    );

    QVector<QVector<QPair<ScanNumber, ScanPoints>>> scanPointsTranched;
    ParallelUtils::tranchMapForParallelization(
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

    QString buildTargetDecoyKey(
            const PeptideStringWithMods &peptideStringWithMods,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            Charge charge
            ) {

        const QString &sep = S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP;

        return peptideStringWithMods + sep + QString::number(charge) + sep + uniqueMsInfoScanKey;
    }

    Err buildClassifierInput(
            const QVector<ScoredCandidate> &scoredCandidates,
            QVector<QVector<double>> *targetScoresVector,
            QVector<QVector<double>> *decoyScoresVector
            ) {

        ERR_INIT

        QMap<PeptideStringWithMods, ScoredCandidate> targets;
        QMap<PeptideStringWithMods, ScoredCandidate> decoys;
        for (const ScoredCandidate &sc : scoredCandidates) {

            const QString key = buildTargetDecoyKey(sc.peptideStringWithMods, sc.targetKey, sc.charge);
            if(sc.isDecoy) {
                decoys.insert(key, sc);
            }
            else {
                targets.insert(key, sc);
            }
        }

        e = ErrorUtils::isEqual(targets.size(), decoys.size()); ree;

        for (auto it = targets.begin(); it != targets.end(); it++) {

            const ScoredCandidate &scTarget = it.value();
            const QString key = buildTargetDecoyKey(scTarget.peptideStringWithMods, scTarget.targetKey, scTarget.charge);

            e = ErrorUtils::isTrue(decoys.contains(key)); ree;

            const ScoredCandidate scDecoy = decoys.value(key);

            const QVector<double> targetScores = {
                    scTarget.cosineSimSum,
                    scTarget.cosineSimMS1,
                    std::pow(scTarget.cosineSimSpectrum, 3)
            };

            const QVector<double> decoyScores = {
                    std::max(scDecoy.cosineSimSum, 0.0),
                    std::max(scDecoy.cosineSimMS1, 0.0),
                    std::max(std::pow(scDecoy.cosineSimSpectrum, 3), 0.0)
            };

            targetScoresVector->push_back(targetScores);
            decoyScoresVector->push_back(decoyScores);
        }

        e = ErrorUtils::isEqual(targetScoresVector->size(), decoyScoresVector->size()); ree;

        ERR_RETURN
    }


}//namespace
Err PythiaDIAWorkflow::buildCalibration(MsReaderParquet *msReaderParquet) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_fragLibReader.libarySize() > 0); ree;

    const double calibrationSelectionFraction = -0.5;
    const int minTopNMs2Ions = 6;
    const int topNMs2IonsCalibration = std::max(
            minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
    );

    qDebug() << "Using top:" << topNMs2IonsCalibration << "fragments";

    QVector<ScoredCandidate> scoredCandidatesCalibration;
    e = extractTargetDecoyData(
            topNMs2IonsCalibration,
            calibrationSelectionFraction,
            msReaderParquet,
            &scoredCandidatesCalibration
    ); ree;

//    QVector<QVector<double>> targetScoresVector;
//    QVector<QVector<double>> decoyScoresVector;
//    e = buildClassifierInput(
//            scoredCandidatesCalibration,
//            &targetScoresVector,
//            &decoyScoresVector
//            ); ree;
//
//
//    QVector<QVector<double>> A;
//    QVector<double> b;
//    e = ClassifierWeightsManager::buildDataClassifier1(
//            targetScoresVector,
//            decoyScoresVector,
//            &A,
//            &b
//            ); ree;
//
//    QVector<double> weights;
//    e = ClassifierWeightsManager::fitWeights(A, b, &weights); ree;
//
//    QVector<double> results;
//
//    for (const ScoredCandidate &sc : scoredCandidatesCalibration) {
//        e = ClassifierWeightsManager::applyWeights({{
//            std::max(sc.cosineSimSum, 0.0),
//            std::max(sc.cosineSimMS1, 0.0),
//            std::pow(std::max(0.0, sc.cosineSimSpectrum), 3)
//        }}, weights, &results); ree;
//
//        qDebug() << sc.cosineSimSum << results.front();
//    }
//
//    qDebug() << weights;


//    qDebug() << weights;
//    for (int i = 0; i < results.size(); i++) {
//        qDebug() << b.at(i) << results.at(i);
//    }

    ClassifierWeightsManager classifierWeightsManager;


#define WRITE_CALIBRATION
#ifdef WRITE_CALIBRATION
    const QString resultsFilePath = msReaderParquet->filePath() + ".pythiaDIA";
    e = ParquetReader::write(scoredCandidatesCalibration, resultsFilePath); ree;
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

    struct ParallelProcessingInput {
        UniqueMsInfoScanKey uniqueMsInfoScanKey;
        PythiaParameters pythiaParameters;
        QString iRTReCalFilePath;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
        QMap<ScanNumber , ScanPoints> *scanPoints = nullptr;
        QMap<PeptideStringWithMods, CandidatePeptide> *peptideStringWithModsVsCandidatePeptide = nullptr;
    };

    QPair<Err, QVector<ScoredCandidate>> parallelProciessingLogic(const ParallelProcessingInput &ppi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

//        qDebug() << "Processing" << ppi.uniqueMsInfoScanKey;

        MsFrameScoretron msFrameScoretron;

        e = msFrameScoretron.init(
                ppi.uniqueMsInfoScanKey,
                ppi.pythiaParameters,
                *ppi.scanPoints,
                ppi.scanNumberVsScanTimeMS1,
                *ppi.peptideStringWithModsVsCandidatePeptide,
                ppi.scanNumberVsScanTime
        ); rree;

//        e = msFrameScoretron.init(
//                ppi.uniqueMsInfoScanKey,
//                ppi.pythiaParameters,
//                *ppi.scanPoints,
//                ppi.scanNumberVsScanTimeMS1,
//                *ppi.peptideStringWithModsVsCandidatePeptide,
//                ppi.scanNumberVsScanTime,
//                ppi.iRTReCalFilePath
//                ); rree;

        QVector<ScoredCandidate> scoredCandidates;
        e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates); rree;

        qDebug() << "Processed" << ppi.uniqueMsInfoScanKey << "Count" << scoredCandidates.size() << et.elapsed() << "mSec";

        return {e, scoredCandidates};
    }

}//namespace
Err PythiaDIAWorkflow::extractTargetDecoyData(
        int topNMs2Ions,
        double selectionFraction,
        MsReaderParquet *msReaderParquet,
        QVector<ScoredCandidate> *combinedResults
        ) {

    ERR_INIT

    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptideCalibration;
    e = buildCandidates(
            topNMs2Ions,
            selectionFraction,
            msReaderParquet,
            &uniqueInfoScanKeyVsCandidatePeptideCalibration
    ); ree;

    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> uniqueInfoScanKeyVsScanPoints;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
    e = buildUniqueMsInfoScanKeyVsScanPoints(
            msReaderParquet,
            &uniqueInfoScanKeyVsScanPoints,
            &scanNumberVsScanTime,
            &scanNumberVsScanTimeMS1
    ); ree;

    QVector<ParallelProcessingInput> parallelProcessingInputs;
    for (const UniqueMsInfoScanKey &uniqueMsInfoScanKey : uniqueInfoScanKeyVsCandidatePeptideCalibration.keys()) {
        ParallelProcessingInput ppi;
        ppi.uniqueMsInfoScanKey = uniqueMsInfoScanKey;
        ppi.scanNumberVsScanTime = scanNumberVsScanTime;
        ppi.scanNumberVsScanTimeMS1 = scanNumberVsScanTimeMS1;
        ppi.iRTReCalFilePath = m_iRTReCalFilePath;
        ppi.pythiaParameters = m_pythiaParameters;
        ppi.scanPoints = &uniqueInfoScanKeyVsScanPoints[uniqueMsInfoScanKey];
        ppi.peptideStringWithModsVsCandidatePeptide = &uniqueInfoScanKeyVsCandidatePeptideCalibration[uniqueMsInfoScanKey];

        parallelProcessingInputs.push_back(ppi);
    }

#define PROCESS_PARALLEL
#ifdef PROCESS_PARALLEL
    QFuture<QPair<Err, QVector<ScoredCandidate>>> futures = QtConcurrent::mapped(
            parallelProcessingInputs,
            parallelProciessingLogic
    );
    futures.waitForFinished();

    for (const QPair<Err, QVector<ScoredCandidate>> &res : futures) {

        e = res.first; ree;
        combinedResults->append(res.second);
    }
#else
    for (const ParallelProcessingInput &inp : parallelProcessingInputs) {

        if (inp.uniqueMsInfoScanKey != "695066") {
            continue;
        }

        const QPair<Err, QVector<ScoredCandidate>> res = parallelProciessingLogic(inp);
        e = res.first; ree;
        combinedResults->append(res.second);
    }
#endif

    ERR_RETURN
}

Err PythiaDIAWorkflow::buildCandidates(
        int topNMs2Ions,
        double selectionListFraction,
        MsReaderParquet *msReaderParquet,
        QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> *uniqueInfoScanKeyVsCandidatePeptide = nullptr
) {

    ERR_INIT

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
                100.0,
                2000.0,
                &peptideSequenceChargeKeyVsCandidatePeptide
        ); ree;

    }

    if (m_msScanInfos.isEmpty()) {
        m_msScanInfos = msReaderParquet->getUniqueTandemMsScanInfos();
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

Err PythiaDIAWorkflow::buildUniqueMsInfoScanKeyVsScanPoints(
        MsReaderParquet *msReaderParquet,
        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
        QMap<ScanNumber, ScanTime> *scanNumberVsScanTime,
        QMap<ScanNumber, ScanPoints > *scanNumberVsScanTimeMS1
        ) {

    ERR_INIT

    e = msReaderParquet->collateTandemPrecursorTargetsDIA(diaTargetFrames); ree

    const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderParquet->getMsScanInfos();
    for (auto it = msScanInfos.begin(); it != msScanInfos.end(); it++) {
        scanNumberVsScanTime->insert(it.key(), it.value().scanTime);
    }

    const int msLevel = 1;
    e = msReaderParquet->getScanPoints(msLevel, scanNumberVsScanTimeMS1); ree;

    ERR_RETURN
}
