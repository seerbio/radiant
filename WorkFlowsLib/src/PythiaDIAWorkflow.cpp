//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MS2ChargeDeconvolvotron.h"
#include "MsCalibrationReader.h"
#include "MsFrameScoretron.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "PeakIntegratomatic.h"
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

#ifndef USE_FILE_CACHING
    e = deisotopeScans(&msReaderParquet); ree;
#endif

    e = buildCalibration(&msReaderParquet); ree;

    e = optimizeParameters(&msReaderParquet); ree;

//    const QString resultsFilePath = msDataFilePath + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
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

    template <typename T>
    QVector<double> extractScoresFromVecFeatures(
            const QVector<T> &featureVec,
            int theoMzIonsSize
            ) {

        QVector<double> vec(theoMzIonsSize, 0.0);

        for (int i = 0; i < theoMzIonsSize; i++) {

            if (i >= featureVec.size()) {
                break;
            }

            vec[i] = static_cast<double>(featureVec[i]);
        }

        return vec;
    }

    QVector<double> buildScoreVector(
            const ScoredCandidate &scoreCandidate,
            bool useExtendedScores,
            int theoMzIonsSize
            ) {

        QVector<double> scores = {
                std::max(scoreCandidate.cosineSimSum, 0.0),
                std::max(scoreCandidate.cosineSimMS1, 0.0),
                std::pow(std::max(0.0, scoreCandidate.cosineSimSpectrum), 3)
        };

        if (useExtendedScores) {
            scores.push_back(std::abs(scoreCandidate.scanTime - scoreCandidate.scanTimePredicted));
            scores.push_back(scoreCandidate.theoFragmentCount);
            scores.push_back(scoreCandidate.charge);
            scores.push_back(scoreCandidate.peptideStringWithMods.size());

            scores.push_back(scoreCandidate.mass);

//            const QVector<double> foundIntensities
//                    = extractScoresFromVecFeatures(scoreCandidate.intensityFoundMaxVec, theoMzIonsSize);
//            const double intensitySum = std::accumulate(foundIntensities.begin(), foundIntensities.end(), 0.0);
//            scores.push_back(std::log(std::max(1.0, intensitySum)));
//            scores.push_back(std::log(std::max(1.0, scoreCandidate.xCorr)));

            const QVector<double> cosineSimToAnchors
                    = extractScoresFromVecFeatures(scoreCandidate.cosineSimToAnchorVec, theoMzIonsSize);
            scores.append(cosineSimToAnchors);

            const QVector<double> topHalfCosineSimScores = cosineSimToAnchors.mid(0, theoMzIonsSize / 2);
            scores.push_back(std::accumulate(topHalfCosineSimScores.begin(), topHalfCosineSimScores.end(), 0.0));

//            else if (theoMzIonsSize % 3 == 0) {
//
//                for (int i = 0; i < theoMzIonsSize; i += 3) {
//                    const QVector<double> cos1 = cosineSimToAnchors.mid(i, 3);
//                    scores.push_back(std::accumulate(cos1.begin(), cos1.end(), 0.0));
//                }
//            }

//            const QVector<double> individualPeakPointCount
//                    = extractScoresFromVecFeatures(scoreCandidate.peakPointCountFoundVec, theoMzIonsSize);
//            scores.append(individualPeakPointCount);

//            const QVector<double> frameIndexMaxDiffFromAnchorVec
//                    = extractScoresFromVecFeatures(scoreCandidate.frameIndexMaxDiffFromAnchorVec, theoMzIonsSize);
//            scores.append(frameIndexMaxDiffFromAnchorVec);

            const QVector<double> theoApexIntensity
                    = extractScoresFromVecFeatures(scoreCandidate.theoIntensityVec, theoMzIonsSize);
            scores.append(theoApexIntensity);

            const QVector<double> intensityFoundMaxVec
                    = extractScoresFromVecFeatures(scoreCandidate.intensityFoundMaxVec, theoMzIonsSize);
            const double maxIntensity = std::max(
                    *std::max(intensityFoundMaxVec.begin(), intensityFoundMaxVec.end()),
                    1.0
                    );
            QVector<double> intensityFoundMaxVecNorm;
            std::transform(
                    intensityFoundMaxVec.begin(),
                    intensityFoundMaxVec.end(),
                    std::back_inserter(intensityFoundMaxVecNorm),
                    [maxIntensity](double d){return d / maxIntensity;}
                    );
            scores.append(intensityFoundMaxVecNorm);


//            QVector<double> ppmVec;
//            for (int i = 0; i < scoreCandidate.mzSearchedVec.size(); i++) {
//
//                const double mzSearched = scoreCandidate.mzSearchedVec.at(i);
//                if (i >= scoreCandidate.mzFoundMeanVec.size()) {
//                    break;
//                }
//
//                const double mzFound = scoreCandidate.mzFoundMeanVec[i];
//
//                const double ppm = 1e6 * (mzFound - mzSearched) / mzSearched;
//                ppmVec.push_back(ppm);
//            }
//            const QVector<double> ppmMz
//                    = extractScoresFromVecFeatures(ppmVec, theoMzIonsSize);
//            scores.append(ppmMz);

            const QVector<double> mzStDev
                    = extractScoresFromVecFeatures(scoreCandidate.mzFoundStDevVec, theoMzIonsSize);
            scores.append(mzStDev);

        }

        return scores;
    }

    Err buildClassifierInput(
            const QVector<ScoredCandidate> &scoredCandidates,
            bool useExtendedScores,
            int theoMzIonsSize,
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

            const QVector<double> targetScores = buildScoreVector(scTarget, useExtendedScores, theoMzIonsSize);
            const QVector<double> decoyScores = buildScoreVector(scDecoy, useExtendedScores,theoMzIonsSize);

            targetScoresVector->push_back(targetScores);
            decoyScoresVector->push_back(decoyScores);
        }

        e = ErrorUtils::isEqual(targetScoresVector->size(), decoyScoresVector->size()); ree;

        ERR_RETURN
    }

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

    Err separateScoredTargetDecoys(
            const QVector<ScoredCandidate> &scoredCandidatesCalibration,
            const QVector<double> &weights,
            bool useExtendedScores,
            const int theoMzIonsSize,
            QMap<PeptideStringWithMods, ScoredCandidate> *peptideStringWithModsVsScoredCandidateTargets,
            QMap<PeptideStringWithMods, double> *peptideStringWithModsVsDiscScoreTargets,
            QMap<PeptideStringWithMods, double> *peptideStringWithModsVsDiscScoreDecoys
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesCalibration); ree;
        e = ErrorUtils::isNotEmpty(weights); ree;

        for (const ScoredCandidate &sc : scoredCandidatesCalibration) {

            QVector<double> scores = buildScoreVector(sc, useExtendedScores, theoMzIonsSize);

            QVector<double> results;
            e = ClassifierWeightsManager::applyWeights({scores}, weights, &results); ree;

            const QString key = buildTargetDecoyKey(sc.peptideStringWithMods, sc.targetKey, sc.charge);

            if(sc.isDecoy) {
                peptideStringWithModsVsDiscScoreDecoys->insert(key, results.front());
                continue;
            }

            peptideStringWithModsVsDiscScoreTargets->insert(key, results.front());
            peptideStringWithModsVsScoredCandidateTargets->insert(key, sc);
        }

        ERR_RETURN
    }

    Err thresholdQValues(
            const QMap<PeptideStringWithMods, ScoredCandidate> &peptideStringWithModsVsScoredCandidateTargets,
            const QMap<QString, double> &identifierVsQValue,
            double fdrThreshold,
            QVector<ScoredCandidate> *scoredCandidatesFDRThresholded
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithModsVsScoredCandidateTargets); ree;
        e = ErrorUtils::isNotEmpty(identifierVsQValue); ree;
        e = ErrorUtils::isTrue(fdrThreshold > 0.0); ree;

        for (auto it = identifierVsQValue.begin(); it != identifierVsQValue.end(); it++) {

            const QString &key = it.key();
            const double qVal = it.value();

            if (qVal > fdrThreshold) {
                continue;
            }

            const ScoredCandidate &sc = peptideStringWithModsVsScoredCandidateTargets.value(key);
            scoredCandidatesFDRThresholded->push_back(sc);
        }

        ERR_RETURN
    }

    Err buildScoredCandidatesFDRThresholded(
            const QVector<ScoredCandidate> &scoredCandidatesCalibration,
            double fdrThreshold,
            bool useExtendedScores,
            int theoMzIonsSize,
            QVector<ScoredCandidate> *scoredCandidatesFDRThresholded
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesCalibration); ree;
        e = ErrorUtils::isTrue(fdrThreshold > 0.0); ree;

        QVector<QVector<double>> targetScoresVector;
        QVector<QVector<double>> decoyScoresVector;
        e = buildClassifierInput(
                scoredCandidatesCalibration,
                useExtendedScores,
                theoMzIonsSize,
                &targetScoresVector,
                &decoyScoresVector
        ); ree;

        QVector<QVector<double>> A;
        QVector<double> b;
        e = ClassifierWeightsManager::buildDataClassifier1(
                targetScoresVector,
                decoyScoresVector,
                &A,
                &b
        ); ree;

        QVector<double> weights;
        e = ClassifierWeightsManager::fitWeights(A, b, &weights); ree;

        QMap<PeptideStringWithMods, ScoredCandidate> peptideStringWithModsVsScoredCandidateTargets;
        QMap<PeptideStringWithMods, double> peptideStringWithModsVsDiscScoreTargets;
        QMap<PeptideStringWithMods, double> peptideStringWithModsVsDiscScoreDecoys;
        e = separateScoredTargetDecoys(
                scoredCandidatesCalibration,
                weights,
                useExtendedScores,
                theoMzIonsSize,
                &peptideStringWithModsVsScoredCandidateTargets,
                &peptideStringWithModsVsDiscScoreTargets,
                &peptideStringWithModsVsDiscScoreDecoys
        ); ree;

        QMap<QString, double> identifierVsQValue;
        QMap<QString, double> identifierVsDecoyRatio;
        e = MathUtils::calculateQValue(
                peptideStringWithModsVsDiscScoreTargets,
                peptideStringWithModsVsDiscScoreDecoys,
                &identifierVsQValue,
                &identifierVsDecoyRatio
        ); ree;

        e = thresholdQValues(
                peptideStringWithModsVsScoredCandidateTargets,
                identifierVsQValue,
                fdrThreshold,
                scoredCandidatesFDRThresholded
        ); ree;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::buildCalibration(MsReaderParquet *msReaderParquet) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_fragLibReader.libarySize() > 0); ree;
    const double fdrThreshold = 0.01;

    const double calibrationSelectionFraction = 0.2;
    const int minTopNMs2Ions = 6;
    const int topNMs2IonsCalibration = std::max(
            minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
    );

    qDebug() << "Using top:" << topNMs2IonsCalibration << "fragments for calibration";

    QVector<ScoredCandidate> scoredCandidatesCalibration;
    e = extractTargetDecoyData(
            m_pythiaParameters,
            topNMs2IonsCalibration,
            calibrationSelectionFraction,
            msReaderParquet,
            &scoredCandidatesCalibration
    ); ree;

    const bool useExtendedScores = false;

    QVector<ScoredCandidate> scoredCandidatesFDRThresholded;
    e = buildScoredCandidatesFDRThresholded(
            scoredCandidatesCalibration,
            fdrThreshold,
            useExtendedScores,
            topNMs2IonsCalibration,
            &scoredCandidatesFDRThresholded
            ); ree;

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = buildMsCalibrationReaderRows(
            scoredCandidatesFDRThresholded,
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

    struct ParallelProcessingInput {
        UniqueMsInfoScanKey uniqueMsInfoScanKey;
        PythiaParameters pythiaParameters;
        MsCalibratomatic msCalibratomatic;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
        QMap<ScanNumber , ScanPoints> *scanPoints = nullptr;
        QMap<PeptideStringWithMods, CandidatePeptide> *peptideStringWithModsVsCandidatePeptide = nullptr;
    };

    QPair<Err, QVector<ScoredCandidate>> parallelProciessingLogic(const ParallelProcessingInput &ppi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        MsFrameScoretron msFrameScoretron;

        e = msFrameScoretron.init(
                ppi.uniqueMsInfoScanKey,
                ppi.pythiaParameters,
                *ppi.scanPoints,
                ppi.scanNumberVsScanTimeMS1,
                *ppi.peptideStringWithModsVsCandidatePeptide,
                ppi.scanNumberVsScanTime,
                ppi.msCalibratomatic
                ); rree;

        QVector<ScoredCandidate> scoredCandidates;
        e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates); rree;

        qDebug() << "Processed" << ppi.uniqueMsInfoScanKey << "Count" << scoredCandidates.size() << et.elapsed() << "mSec";

        return {e, scoredCandidates};
    }

}//namespace
Err PythiaDIAWorkflow::extractTargetDecoyData(
        const PythiaParameters &pythiaParameters,
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
        ppi.msCalibratomatic = m_msCalibratomatic;
        ppi.pythiaParameters = pythiaParameters;
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
                m_pythiaParameters.mzMinDataStructure,
                m_pythiaParameters.mzMaxDataStructure,
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

    const int minTopNMs2Ions = 6;
    const int topNMs2IonsOptimization = std::max(
            minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
    );

    qDebug() << "Using top:" << topNMs2IonsOptimization << "fragments for optimization";

    const double selectionFractionBypassValue = 0.1;
    const double fdrThreshold = 0.01;

    QVector<PythiaParameters> pythiaParametersExperiments;
    e = buildDOE(
            m_pythiaParameters,
            m_msCalibratomatic.mzStDev(),
            m_msCalibratomatic.scanTimeStDev(),
            &pythiaParametersExperiments
            ); ree;

    struct DOEResult {
        double mzStDev = -1.0;
        double scanTimeStDev = -1.0;
        double cosineSimAnchor = -1.0;
        int fdrCount = -1;
    };

    const bool useExtendedScores = true;

    QVector<DOEResult> results;
    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {

        QVector<ScoredCandidate> scoredCandidatesCalibration;
        e = extractTargetDecoyData(
                pythiaParams,
                topNMs2IonsOptimization,
                selectionFractionBypassValue,
                msReaderParquet,
                &scoredCandidatesCalibration
        ); ree;

        QVector<ScoredCandidate> scoredCandidatesFDRThresholded;
        e = buildScoredCandidatesFDRThresholded(
                scoredCandidatesCalibration,
                fdrThreshold,
                useExtendedScores,
                topNMs2IonsOptimization,
                &scoredCandidatesFDRThresholded
        ); ree;

        DOEResult res;
        res.mzStDev = pythiaParams.ms2ExtractionWidthPPM;
        res.scanTimeStDev = pythiaParams.scanTimeWindowMinutes;
        res.cosineSimAnchor = pythiaParams.cosineSimToAnchorThreshold;
        res.fdrCount = scoredCandidatesFDRThresholded.size();
        results.push_back(res);
    }

    std::sort(
            results.rbegin(),
            results.rend(),
            [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
            );

    for (const DOEResult &r : results) {
        qDebug() << r.mzStDev << r.scanTimeStDev << r.cosineSimAnchor << r.fdrCount;
    }

    //TODO replace this with response surface derived DOE parameters when you figure out how to do it.
    const DOEResult bestParametersFDR = *std::max_element(
            results.rbegin(),
            results.rend(),
            [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
    );

    m_pythiaParameters.ms2ExtractionWidthPPM = bestParametersFDR.mzStDev;
    m_pythiaParameters.scanTimeWindowMinutes = bestParametersFDR.scanTimeStDev;
    m_pythiaParameters.cosineSimToAnchorThreshold = bestParametersFDR.cosineSimAnchor;
    m_pythiaParameters.print();
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    ERR_RETURN
}
