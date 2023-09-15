//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "ClassifierWeightsManager.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MS2ChargeDeconvolvotron.h"
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
    e = deisotopeScans(&msReaderParquet); ree;
#endif

    e = buildCalibration(&msReaderParquet); ree;

//#define BYPASS_OPTI
#ifndef BYPASS_OPTI
    e = optimizeParameters(&msReaderParquet); ree;
#else
    m_pythiaParameters.ms2ExtractionWidthPPM = 12.6052;
    m_pythiaParameters.scanTimeWindowMinutes = 1.64017;
    m_pythiaParameters.cosineSimToAnchorThreshold = 0.9;
#endif

    QVector<ScoredCandidate> scoredCandidatesTargetsFDRThresholded;
    QVector<ScoredCandidate> scoredCandidatesAll;
    e = mainAnalysis(
            &msReaderParquet,
            &scoredCandidatesTargetsFDRThresholded,
            &scoredCandidatesAll
            ); ree;

    QVector<ScoredCandidate> scoredCandidatesAllUpdated;
    e = removeInterferingCandidates(
            &msReaderParquet,
            scoredCandidatesTargetsFDRThresholded,
            scoredCandidatesAll,
            &scoredCandidatesAllUpdated
            ); ree;


//#define WRITE_FDR_PMSS
#ifdef WRITE_FDR_PMSS
    const QString resultsFilePath = msReaderParquet.filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
    e = ParquetReader::write(scoredCandidatesTargetsFDRThresholded, resultsFilePath); ree;
#endif

#define WRITE_ALL_PMSS
#ifdef WRITE_ALL_PMSS
    const QString resultsFilePath = msReaderParquet.filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
    e = ParquetReader::write(scoredCandidatesAllUpdated, resultsFilePath); ree;
#endif

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
            bool useNeuralNetworkScores,
            int theoMzIonsSize
            ) {

        QVector<double> scores = {
                std::max(scoreCandidate.cosineSimSum, 0.0),
                std::max(scoreCandidate.cosineSimMS1, 0.0),
                //TODO make sure you explore this again if FMR is too big during entrapment experiments
//                std::pow(std::max(0.0, scoreCandidate.cosineSimSpectrum), 3),
                std::pow(std::max(0.0, scoreCandidate.klDivSpectrum), 3)
        };

        if (useExtendedScores || useNeuralNetworkScores) {
            scores.push_back(std::abs(scoreCandidate.scanTime - scoreCandidate.scanTimePredicted));
            scores.push_back(scoreCandidate.theoFragmentCount);
            scores.push_back(scoreCandidate.charge);
            scores.push_back(scoreCandidate.peptideStringWithMods.size());
            scores.push_back(scoreCandidate.mass);

//            scores.push_back(std::log(std::max(1.0, scoreCandidate.xCorr)));

            const QVector<double> cosineSimToAnchors
                    = extractScoresFromVecFeatures(scoreCandidate.cosineSimToAnchorVec, theoMzIonsSize);
            scores.append(cosineSimToAnchors);

            const QVector<double> topHalfCosineSimScores = cosineSimToAnchors.mid(0, theoMzIonsSize / 2);
            scores.push_back(std::accumulate(topHalfCosineSimScores.begin(), topHalfCosineSimScores.end(), 0.0));

            const QVector<double> bottomHalfCosineSimScores = cosineSimToAnchors.mid(theoMzIonsSize / 2, theoMzIonsSize / 2);
            scores.push_back(std::accumulate(bottomHalfCosineSimScores.begin(), bottomHalfCosineSimScores.end(), 0.0));

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

            const QVector<double> mzStDev
                    = extractScoresFromVecFeatures(scoreCandidate.mzFoundStDevVec, theoMzIonsSize);
            scores.append(mzStDev);

        }

        if (useNeuralNetworkScores) {
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
//                ppmVec.push_back(std::min(ppm, 100.0));
//            }
//            const QVector<double> ppmMz
//                    = extractScoresFromVecFeatures(ppmVec, theoMzIonsSize);
//            scores.append(ppmMz);
        }

        return scores;
    }

    Err buildClassifierInput(
            const QVector<ScoredCandidate> &scoredCandidates,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
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

            const QVector<double> targetScores = buildScoreVector(
                    scTarget,
                    useExtendedScores,
                    useNeuralNetworkScores,
                    theoMzIonsSize
                    );

            const QVector<double> decoyScores = buildScoreVector(
                    scDecoy,
                    useExtendedScores,
                    useNeuralNetworkScores,
                    theoMzIonsSize
                    );

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
            bool useNeuralNetworkScores,
            const int theoMzIonsSize,
            QMap<PeptideStringWithMods, ScoredCandidate> *peptideStringWithModsVsScoredCandidateTargets,
            QMap<PeptideStringWithMods, ScoredCandidate> *peptideStringWithModsVsScoredCandidateDecoys,
            QMap<PeptideStringWithMods, double> *peptideStringWithModsVsDiscScoreTargets,
            QMap<PeptideStringWithMods, double> *peptideStringWithModsVsDiscScoreDecoys
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesCalibration); ree;
        e = ErrorUtils::isNotEmpty(weights); ree;

        for (const ScoredCandidate &sc : scoredCandidatesCalibration) {

            QVector<double> scores = buildScoreVector(
                    sc,
                    useExtendedScores,
                    useNeuralNetworkScores,
                    theoMzIonsSize
                    );

            QVector<double> results;
            e = ClassifierWeightsManager::applyWeights({scores}, weights, &results); ree;

            const QString key = buildTargetDecoyKey(sc.peptideStringWithMods, sc.targetKey, sc.charge);

            ScoredCandidate scoredCandidate = sc;
            scoredCandidate.discriminateScore = results.front();

            if(sc.isDecoy) {
                peptideStringWithModsVsDiscScoreDecoys->insert(key, scoredCandidate.discriminateScore);
                peptideStringWithModsVsScoredCandidateDecoys->insert(key, scoredCandidate);
                continue;
            }

            peptideStringWithModsVsDiscScoreTargets->insert(key, scoredCandidate.discriminateScore);
            peptideStringWithModsVsScoredCandidateTargets->insert(key, scoredCandidate);
        }

        ERR_RETURN
    }

    Err assignQValuesAndDecoyRatios(
            const QMap<QString, double> &identifierVsQValue,
            const QMap<QString, double> &identifierVsDecoyRatio,
            QVector<ScoredCandidate> *scoredCandidatesQValsUpdated
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(identifierVsDecoyRatio); ree;
        e = ErrorUtils::isNotEmpty(identifierVsQValue); ree;
        e = ErrorUtils::isNotEmpty(*scoredCandidatesQValsUpdated); ree;

        for (int i = 0; i < scoredCandidatesQValsUpdated->size(); i++) {

            ScoredCandidate &sc = (*scoredCandidatesQValsUpdated)[i];

            if (sc.isDecoy) {
                continue;
            }

            const QString key = buildTargetDecoyKey(sc.peptideStringWithMods, sc.targetKey, sc.charge);

            e = ErrorUtils::isTrue(identifierVsQValue.contains(key)); ree;
            e = ErrorUtils::isTrue(identifierVsDecoyRatio.contains(key)); ree;

            sc.qValue = identifierVsQValue.value(key);
            sc.decoyRatio = identifierVsQValue.value(key);
        }

        ERR_RETURN
    }

    Err buildScoredCandidatesAll(
            const QMap<PeptideStringWithMods, ScoredCandidate> &peptideStringWithModsVsScoredCandidateTargets,
            const QMap<PeptideStringWithMods, ScoredCandidate> &peptideStringWithModsVsScoredCandidateDecoys,
            QVector<ScoredCandidate> *scoredCandidatesAll
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithModsVsScoredCandidateTargets); ree;
        e = ErrorUtils::isNotEmpty(peptideStringWithModsVsScoredCandidateDecoys); ree;

        for (const ScoredCandidate &sc : peptideStringWithModsVsScoredCandidateTargets) {
            scoredCandidatesAll->push_back(sc);
        }

        for (const ScoredCandidate &sc : peptideStringWithModsVsScoredCandidateDecoys) {
            ScoredCandidate scNew = sc;
            scNew.peptideStringWithMods += (static_cast<QString>(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP) + "DECOY");
            scoredCandidatesAll->push_back(scNew);
        }

        ERR_RETURN
    }

    Err buildScoredCandidatesFDR(
            const QVector<ScoredCandidate> &scoredCandidatesCalibration,
            double fdrThreshold,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            int theoMzIonsSize,
            QVector<ScoredCandidate> *scoredCandidatesAll
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesCalibration); ree;
        e = ErrorUtils::isTrue(fdrThreshold > 0.0); ree;

        scoredCandidatesAll->clear();

        QVector<QVector<double>> targetScoresVector;
        QVector<QVector<double>> decoyScoresVector;
        e = buildClassifierInput(
                scoredCandidatesCalibration,
                useExtendedScores,
                useNeuralNetworkScores,
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
        QMap<PeptideStringWithMods, ScoredCandidate> peptideStringWithModsVsScoredCandidateDecoys;
        QMap<PeptideStringWithMods, double> peptideStringWithModsVsDiscScoreTargets;
        QMap<PeptideStringWithMods, double> peptideStringWithModsVsDiscScoreDecoys;
        e = separateScoredTargetDecoys(
                scoredCandidatesCalibration,
                weights,
                useExtendedScores,
                useNeuralNetworkScores,
                theoMzIonsSize,
                &peptideStringWithModsVsScoredCandidateTargets,
                &peptideStringWithModsVsScoredCandidateDecoys,
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

        e = buildScoredCandidatesAll(
                peptideStringWithModsVsScoredCandidateTargets,
                peptideStringWithModsVsScoredCandidateDecoys,
                scoredCandidatesAll
                ); ree;

        e = assignQValuesAndDecoyRatios(
                identifierVsQValue,
                identifierVsDecoyRatio,
                scoredCandidatesAll
        ); ree;

        ERR_RETURN
    }

    Err filterScoreCandidatesByFDR(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            double qValueThreshold,
            QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;
        e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

        scoredCandidatesTargetsFDRThresholded->clear();

        for (const ScoredCandidate &sc : scoredCandidatesAll) {

            if (sc.isDecoy || sc.qValue > qValueThreshold) {
                continue;
            }

            scoredCandidatesTargetsFDRThresholded->push_back(sc);
        }

        ERR_RETURN
    }

    Err countScoreCandidatesByFDR(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            double qValueThreshold,
            int *targetCountBelowFDRThreshold
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;
        e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

        const auto countLogic = [qValueThreshold](const ScoredCandidate &sc){
            return !sc.isDecoy && sc.qValue < qValueThreshold;
        };

        *targetCountBelowFDRThreshold
            = static_cast<int>(std::count_if(scoredCandidatesAll.begin(), scoredCandidatesAll.end(), countLogic));

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::buildCalibration(MsReaderParquet *msReaderParquet) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_fragLibReader.libarySize() > 0); ree;
    const double fdrThreshold = 0.01;

    const double calibrationSelectionFraction = 0.2;
    const int topNMs2IonsCalibration = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
    );

    qDebug() << "Using top:" << topNMs2IonsCalibration << "fragments for calibration";

    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptideCalibration;
    e = buildCandidates(
            topNMs2IonsCalibration,
            calibrationSelectionFraction,
            &uniqueInfoScanKeyVsCandidatePeptideCalibration
    ); ree;

    const bool useExtendedScores = false;
    const bool useNeuralNetworkScores = false;
    QVector<ScoredCandidate> scoredCandidatesAll;
    QVector<ScoredCandidate> scoredCandidatesTargetsFDRThresholded;
    e = extractionLoopLogic(
            uniqueInfoScanKeyVsCandidatePeptideCalibration,
            fdrThreshold,
            useExtendedScores,
            useNeuralNetworkScores,
            topNMs2IonsCalibration,
            msReaderParquet,
            &scoredCandidatesAll,
            &scoredCandidatesTargetsFDRThresholded
            );

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

Err PythiaDIAWorkflow::extractionLoopLogic(
        const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptide,
        double fdrThreshold,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int topNMs2IonsMainAnalysis,
        MsReaderParquet *msReaderParquet,
        QVector<ScoredCandidate> *scoredCandidatesAll,
        QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded
        ) {
    ERR_INIT

    QVector<ScoredCandidate> scoredCandidatesOptimization;
    e = extractTargetDecoyData(
            uniqueInfoScanKeyVsCandidatePeptide,
            m_pythiaParameters,
            msReaderParquet,
            &scoredCandidatesOptimization
    ); ree;

    e = buildScoredCandidatesFDR(
            scoredCandidatesOptimization,
            fdrThreshold,
            useExtendedScores,
            useNeuralNetworkScores,
            topNMs2IonsMainAnalysis,
            scoredCandidatesAll
    ); ree;

    e = filterScoreCandidatesByFDR(
            *scoredCandidatesAll,
            fdrThreshold,
            scoredCandidatesTargetsFDRThresholded
    ); ree;

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
        QMap<PeptideStringWithMods, CandidatePeptide> peptideStringWithModsVsCandidatePeptide;
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
                ppi.peptideStringWithModsVsCandidatePeptide,
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
        const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptideCalibration,
        const PythiaParameters &pythiaParameters,
        MsReaderParquet *msReaderParquet,
        QVector<ScoredCandidate> *combinedResults
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(uniqueInfoScanKeyVsCandidatePeptideCalibration); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

    combinedResults->clear();

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
        ppi.peptideStringWithModsVsCandidatePeptide = uniqueInfoScanKeyVsCandidatePeptideCalibration[uniqueMsInfoScanKey];

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

    const int topNMs2IonsOptimization = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
    );

    qDebug() << "Using top:" << topNMs2IonsOptimization << "fragments for optimization";

    const double selectionFractionValue = 0.1;
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
    const bool useNeuralNetworkScores = false;

    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptideCalibration;
    e = buildCandidates(
            topNMs2IonsOptimization,
            selectionFractionValue,
            &uniqueInfoScanKeyVsCandidatePeptideCalibration
    ); ree;

    QVector<DOEResult> results;
    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {

        QVector<ScoredCandidate> scoredCandidatesOptimization;
        e = extractTargetDecoyData(
                uniqueInfoScanKeyVsCandidatePeptideCalibration,
                pythiaParams,
                msReaderParquet,
                &scoredCandidatesOptimization
        ); ree;

        QVector<ScoredCandidate> scoredCandidatesAll;
        e = buildScoredCandidatesFDR(
                scoredCandidatesOptimization,
                fdrThreshold,
                useExtendedScores,
                useNeuralNetworkScores,
                topNMs2IonsOptimization,
                &scoredCandidatesAll
        ); ree;

       int targetCountAboveFDRQValueThreshold;
        e = countScoreCandidatesByFDR(
                scoredCandidatesAll,
                fdrThreshold,
                &targetCountAboveFDRQValueThreshold
                ); ree;

        DOEResult res;
        res.mzStDev = pythiaParams.ms2ExtractionWidthPPM;
        res.scanTimeStDev = pythiaParams.scanTimeWindowMinutes;
        res.cosineSimAnchor = pythiaParams.cosineSimToAnchorThreshold;
        res.fdrCount = targetCountAboveFDRQValueThreshold;
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
    ///////////////////////////////////////////////////////////////////////////////////////////////////

    ERR_RETURN
}

Err PythiaDIAWorkflow::mainAnalysis(
        MsReaderParquet *msReaderParquet,
        QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded,
        QVector<ScoredCandidate> *scoredCandidatesAll
        ) {

    ERR_INIT

    m_pythiaParameters.print();

    const double selectionFractionBypassValue = -1.0;
    const double fdrThreshold = 0.1;

    const int topNMs2IonsMainAnalysis = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
    );

    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptideCalibration;
    e = buildCandidates(
            topNMs2IonsMainAnalysis,
            selectionFractionBypassValue,
            &uniqueInfoScanKeyVsCandidatePeptideCalibration
    ); ree;

    //TODO determine if this should be false w/ Entrapment experiements
    const bool useExtendedScores = true;
    const bool useNeuralNetworkScores = false;

    e = extractionLoopLogic(
            uniqueInfoScanKeyVsCandidatePeptideCalibration,
            fdrThreshold,
            useExtendedScores,
            useNeuralNetworkScores,
            topNMs2IonsMainAnalysis,
            msReaderParquet,
            scoredCandidatesAll,
            scoredCandidatesTargetsFDRThresholded
    );

    const double fdrOnePercent = 0.01;
    int foundAtOnePercentFDR;
    e = countScoreCandidatesByFDR(
            *scoredCandidatesAll,
            fdrOnePercent,
            &foundAtOnePercentFDR
            ); ree;

    qDebug() <<  foundAtOnePercentFDR << "PSMs found at 1 % FDR";
    qDebug() << scoredCandidatesTargetsFDRThresholded->size() << "PSMs found at " << fdrThreshold * 100 << "% FDR";

    double lastMinQVal = 1.0;
    while (true) {

        QVector<PeptideStringWithMods> scoredCandidatesTargetsFDRThresholdedSequences;
        std::transform(
                scoredCandidatesTargetsFDRThresholded->begin(),
                scoredCandidatesTargetsFDRThresholded->end(),
                std::back_inserter(scoredCandidatesTargetsFDRThresholdedSequences),
                [](const ScoredCandidate &sc){return sc.peptideStringWithMods;}
        );

        e = buildCandidates(
                scoredCandidatesTargetsFDRThresholdedSequences,
                topNMs2IonsMainAnalysis,
                &uniqueInfoScanKeyVsCandidatePeptideCalibration
        ); ree;

        e = extractionLoopLogic(
                uniqueInfoScanKeyVsCandidatePeptideCalibration,
                fdrThreshold,
                useExtendedScores,
                useNeuralNetworkScores,
                topNMs2IonsMainAnalysis,
                msReaderParquet,
                scoredCandidatesAll,
                scoredCandidatesTargetsFDRThresholded
        );

        e = countScoreCandidatesByFDR(
                *scoredCandidatesAll,
                fdrOnePercent,
                &foundAtOnePercentFDR
                ); ree;


        const double minQVal = std::min_element(
                scoredCandidatesTargetsFDRThresholded->begin(),
                scoredCandidatesTargetsFDRThresholded->end(),
                [](const ScoredCandidate &l, const ScoredCandidate &r){return l.qValue < r.qValue;}
                )->qValue;

        qDebug() << "Min QValue" << minQVal;
        qDebug() <<  foundAtOnePercentFDR << "PSMs found at 1 % FDR";
        qDebug() << scoredCandidatesTargetsFDRThresholded->size() << "PSMs found at " << fdrThreshold * 100 << "% FDR";

        if (minQVal > lastMinQVal || MathUtils::tSame(minQVal, lastMinQVal)) {
            break;
        }

        lastMinQVal = minQVal;
    }

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

    void filterScoredCandidatesByWeightAndPVal(QVector<ScoredCandidate> *scoredCandidatesAllUpdated) {

        const double pValThreshold = 0.05;
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

    filterScoredCandidatesByWeightAndPVal(scoredCandidatesAllUpdated);

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
