//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "BiophysicalCalcs.h"
#include "CalibrationReader.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FrameExtractsReader.h"
#include "IRTPredictron.h"
#include "MsFrameScoretronProcessormatic.h"
#include "NearestNeighborsSearch.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"
#include "TandemSpectraDeconvolvotron.h"
#include "TurboXIC.h"

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
using rTreePoint = std::pair<rTreeBox, int> ;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;


Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QVector<FeatureFinderHill> &featureFinderHills,
        const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(featureFinderHills); ree;
    e = ErrorUtils::isNotEmpty(peptideStringWithModsVsCandidatePeptide); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;

    m_params = params;
    m_featureFinderHills = featureFinderHills;
    m_fragPredsTopN = peptideStringWithModsVsCandidatePeptide;
    m_scanNumberVsScanTime = scanNumberVsScanTime;

    e = FragLibReader::generateFragmentFrequencies(
            m_fragPredsTopN,
            m_params.ms2ExtractionWidthPPM,
            &m_fragmentFrequencies
    ); ree

    ERR_RETURN
}

Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QVector<FeatureFinderHill> &featureFinderHills,
        const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const QString &iRTRecalibrationFilePath
        ) {

    ERR_INIT

    e = init(
            params,
            featureFinderHills,
            peptideStringWithModsVsCandidatePeptide,
            scanNumberVsScanTime
            ); ree;

    qDebug() << "updating rt vals from iRT";

    QVector<QPair<double, Coors>> nnInputData;
    e = IRTPredictron::buildNearestNeighborsIRTData(
            iRTRecalibrationFilePath,
            &nnInputData
    );

    NearestNeighborsSearch nearestNeighborsSearch;
    e = nearestNeighborsSearch.init(nnInputData); ree

    const int kNearestPoints = 10;

    for (auto it = peptideStringWithModsVsCandidatePeptide.begin(); it != peptideStringWithModsVsCandidatePeptide.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const double iRT = it.value().iRt;

        QVector<NNSearchResult> nnSearchResult;
        nearestNeighborsSearch.kNearestNeighborsSearch(
                {{iRT, 0.0}},
                kNearestPoints,
                &nnSearchResult
        );

        m_fragPredsPredictedScanTime.insert(peptideStringWithMods, nnSearchResult.front().meanValues);
    }

    ERR_RETURN
}


namespace {

    QMap<PeptideStringWithMods, QVector<MS2IonPeak>> buildScanNumberVsMS2IonPeaks(const QVector<MS2IonPeak> &ms2IonPeaks) {

        QMap<PeptideStringWithMods , QVector<MS2IonPeak>> scanNumberVsMS2IonPeaks;
        for (const MS2IonPeak &ms2IonPeak : ms2IonPeaks) {
            scanNumberVsMS2IonPeaks[ms2IonPeak.peptideStringWithMods].push_back(ms2IonPeak);
        }

        return scanNumberVsMS2IonPeaks;
    }

    RTree loadRTree(
            const QVector<MS2IonPeak> &ms2IonPeaks,
            QMap<Id, MS2IonPeak> *idVsMs2IonPeak,
            QSet<int> *frameIndexMaxesUnique
            ) {

        *idVsMs2IonPeak =  ParallelUtils::convertVectorToMap(ms2IonPeaks);

        std::vector<rTreePoint> cloudLoader;

        for (auto it = idVsMs2IonPeak->begin(); it != idVsMs2IonPeak->end(); it++) {

            const Id id = it.key();
            const MS2IonPeak &ms2ip = it.value();

            frameIndexMaxesUnique->insert(ms2ip.frameIndexMax);

            rTreeBox box(
                    {static_cast<double>(ms2ip.frameIndexStart), ms2ip.mzSearched},
                    {static_cast<double>(ms2ip.frameIndexEnd), ms2ip.mzSearched}
            );

            cloudLoader.push_back(rTreePoint({box, id}));
        }

        const int maxElements = 16;
        return RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    }

    QSet<QVector<Id>> buildUniqueIDs(
            const RTree &rtree,
            int mzPeaksFoundMin,
            const QSet<int> &frameIndexMaxesUnique
            ) {

        const double mzMin = 0.0;
        const double mzMax = 2000.0;

        QSet<QVector<Id>> idsUnique;
        for (int frameIndexMax : frameIndexMaxesUnique) {

            const rTreeBox queryBox(
                    rTreeCoor(frameIndexMax, mzMin),
                    rTreeCoor(frameIndexMax, mzMax)
            );

            std::vector<rTreePoint> rTreeSearchResult;
            rtree.query(
                    bgi::intersects(queryBox),
                    std::back_inserter(rTreeSearchResult)
            );

            QVector<Id> ids;
            std::transform(
                    rTreeSearchResult.begin(),
                    rTreeSearchResult.end(),
                    std::back_inserter(ids),
                    [&](const rTreePoint &p){return p.second;}
            );

            if (ids.size() < mzPeaksFoundMin) {
                continue;
            }

            std::sort(ids.begin(), ids.end());
            idsUnique.insert(ids);
        }

        return idsUnique;
    }

    Err buildComparisonMatrix(
            const MS2IonPeak &p1,
            const MS2IonPeak &p2,
            Eigen::MatrixX<double> *comparisonMatrix
    ) {

        ERR_INIT

        const QPair<int,int> p1ScanNumberIndexMinMax = {p1.frameIndexStart, p1.frameIndexEnd};
        const QPair<int,int> p2ScanNumberIndexMinMax = {p2.frameIndexStart, p2.frameIndexEnd};

        const int lowerScanNumberIndex = std::min(
                p1ScanNumberIndexMinMax.first,
                p2ScanNumberIndexMinMax.first
        );

        const int upperScanNumberIndex = std::max(
                p1ScanNumberIndexMinMax.second,
                p2ScanNumberIndexMinMax.second
        );

        const int colCount = 2;
        const int rowBuffer = 1;
        Eigen::MatrixX<double> mat(upperScanNumberIndex + rowBuffer, colCount);
        mat.setZero();

        const QVector<double> p1Intensities = p1.intensityVals;
        const QVector<double> p2Intensities = p2.intensityVals;
        QVector<int> p1ScanNumberIndexes(p1Intensities.size());
        std::iota(p1ScanNumberIndexes.begin(), p1ScanNumberIndexes.end(), p1ScanNumberIndexMinMax.first);

        QVector<int> p2ScanNumberIndexes(p2Intensities.size());
        std::iota(p2ScanNumberIndexes.begin(), p2ScanNumberIndexes.end(), p2ScanNumberIndexMinMax.first);

        e = ErrorUtils::isEqual(p1Intensities.size(), p1ScanNumberIndexes.size()); ree
        e = ErrorUtils::isEqual(p2Intensities.size(), p2ScanNumberIndexes.size()); ree

        const int insertColumn1 = 0;
        for (int i = 0; i < p1Intensities.size(); i++) {

            const int scanNumberIndex = p1ScanNumberIndexes.at(i);
            const double intensity = p1Intensities.at(i);

            mat.coeffRef(scanNumberIndex, insertColumn1) = intensity;
        }

        const int insertColumn2 = 1;
        for (int i = 0; i < p2Intensities.size(); i++) {

            const int scanNumberIndex = p2ScanNumberIndexes.at(i);
            const double intensity = p2Intensities.at(i);

            mat.coeffRef(scanNumberIndex, insertColumn2) = intensity;
        }

        const int filterLength = 3;
        const double sigma = 1.0;
        const Eigen::VectorX<double> gaussKernel = EigenKernelUtils::buildGaussianFilter1D(
                filterLength,
                sigma,
                false
        );

        mat = EigenKernelUtils::applyKernelRowWiseToMatrix(mat, gaussKernel);

        *comparisonMatrix = mat.block(
                lowerScanNumberIndex,
                insertColumn1,
                (upperScanNumberIndex - lowerScanNumberIndex + rowBuffer),
                colCount
        );

        ERR_RETURN
    }

    Err calculateCosineSimBetweenPeaks(
            const MS2IonPeak &p1,
            const MS2IonPeak &p2,
            double *cosineSim
    ) {

        ERR_INIT

        Eigen::MatrixX<double> comparisonMatrix;
        e = buildComparisonMatrix(
                p1,
                p2,
                &comparisonMatrix
        ); ree

        *cosineSim = EigenUtils::cosineSimilarity(comparisonMatrix.col(0), comparisonMatrix.col(1));

        ERR_RETURN
    }

    Err clusterLogic(
            const QVector<MS2IonPeak> &ms2ip,
            int mzPeaksFoundMin,
            double cosineSimMinThreshold,
            QVector<MS2IonPeak> *bestMS2IonPeaksClusteringOutput,
            double *bestCosineSimSumOutput
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2ip); ree

        QMap<Id, MS2IonPeak> idVsMs2IonPeak;
        QSet<int> frameIndexMaxesUnique;
        const RTree rtree = loadRTree(ms2ip, &idVsMs2IonPeak, &frameIndexMaxesUnique);

        const QSet<QVector<Id>> uniqueIDs = buildUniqueIDs(
                rtree,
                mzPeaksFoundMin,
                frameIndexMaxesUnique
                );

        if (uniqueIDs.isEmpty()) {
            ERR_RETURN
        }

        QVector<QVector<MS2IonPeak>> clusteredMs2IonPeaksForCosineSimCalc;
        const auto transformLogic = [&](const QVector<Id> &ids){
            QVector<MS2IonPeak> peaks;
            for (Id id : ids) {
                peaks.push_back(idVsMs2IonPeak.value(id));
            }
            return peaks;
        };
        std::transform(
                uniqueIDs.begin(),
                uniqueIDs.end(),
                std::back_inserter(clusteredMs2IonPeaksForCosineSimCalc),
                transformLogic
                );

        const int maxRank = 5;

        QVector<MS2IonPeak> bestMS2IonPeaksClustering;
        double bestCosineSimSum = 0.0;

        for (const QVector<MS2IonPeak> &ms2ips : clusteredMs2IonPeaksForCosineSimCalc) {

            for (const MS2IonPeak &ms2IonPeakAnchor : ms2ips) {

                if (ms2IonPeakAnchor.rank > maxRank) {
                    continue;
                }

                double cosineSimSum = 0.0;
                QVector<MS2IonPeak> clustering;

                for (const MS2IonPeak &ms2IonPeak : ms2ips) {

                    double cosineSim;
                    e = calculateCosineSimBetweenPeaks(ms2IonPeakAnchor, ms2IonPeak, &cosineSim); ree

                    if (cosineSim < cosineSimMinThreshold) {
                        continue;
                    }

                    MS2IonPeak ms2IonPeakNew = ms2IonPeak;
                    ms2IonPeakNew.cosineSimToAnchor = cosineSim;
                    ms2IonPeakNew.frameIndexMaxDiffFromAnchor = ms2IonPeakAnchor.frameIndexMax - ms2IonPeak.frameIndexMax;

                    clustering.push_back(ms2IonPeakNew);
                    cosineSimSum += ms2IonPeakNew.cosineSimToAnchor;
                }

                if (cosineSimSum > bestCosineSimSum) {
                    bestCosineSimSum = cosineSimSum;
                    bestMS2IonPeaksClustering = clustering;
                }
            }
        }

        *bestMS2IonPeaksClusteringOutput = bestMS2IonPeaksClustering;
        *bestCosineSimSumOutput = bestCosineSimSum;

        ERR_RETURN
    }

    Err findBestClusterGroupingPerPeptideStringWithMods(
            const QVector<MS2IonPeak> &ms2IonPeaks,
            int mzPeaksFoundMin,
            double cosineSimMinThreshold,
            QMap<PeptideStringWithMods, QPair<double, QVector<MS2IonPeak>>> *bestClusters
            ) {

        ERR_INIT

        const QMap<PeptideStringWithMods, QVector<MS2IonPeak>> scanNumberVsMS2IonPeaks
            = buildScanNumberVsMS2IonPeaks(ms2IonPeaks);

        for (auto it = scanNumberVsMS2IonPeaks.begin(); it != scanNumberVsMS2IonPeaks.end(); it++) {

            const PeptideStringWithMods &peptideStringWithMods = it.key();
            const QVector<MS2IonPeak> &ms2ip = it.value();

            QVector<MS2IonPeak> bestMS2IonPeaksClustering;
            double bestCosineSimSum;
            e = clusterLogic(
                    ms2ip,
                    mzPeaksFoundMin,
                    cosineSimMinThreshold,
                    &bestMS2IonPeaksClustering,
                    &bestCosineSimSum
                    ); ree

            if (bestMS2IonPeaksClustering.size() < mzPeaksFoundMin) {
                continue;
            }

            bestClusters->insert(peptideStringWithMods, {bestCosineSimSum, bestMS2IonPeaksClustering});
        }


        ERR_RETURN
    }

    Err fillUnFoundMS2IonPeaks(
            const QVector<MS2IonPeak> &bestClusterPeaks,
            const QVector<MS2Ion> &fragPredTopNOfCluster,
            const QMap<MzHashed, double> &fragFrequencies,
            QVector<MS2IonPeak> *bestClusterPeaksComplete
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(fragPredTopNOfCluster); ree
        e = ErrorUtils::isNotEmpty(bestClusterPeaks); ree

        QMap<MzHashed, MS2IonPeak> completeCluster;

        for (const MS2Ion &ms2Ion : fragPredTopNOfCluster) {
            const MzHashed mzHashed
                    = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

            MS2IonPeak ms2IonPeakUnfound;
            ms2IonPeakUnfound.mzSearched = ms2Ion.mz;
            ms2IonPeakUnfound.theoIntensity = ms2Ion.intensity;
            ms2IonPeakUnfound.pointCountFound = 0;
            ms2IonPeakUnfound.fragmentFrequency
                = fragFrequencies.value(MathUtils::hashDecimal(ms2IonPeakUnfound.mzSearched, S_GLOBAL_SETTINGS.HASHING_PRECISION));

            completeCluster.insert(mzHashed, ms2IonPeakUnfound);
        }

        for (const MS2IonPeak &ms2IonPeak : bestClusterPeaks) {

            const MzHashed &mzHashed
                    = MathUtils::hashDecimal(ms2IonPeak.mzSearched, S_GLOBAL_SETTINGS.HASHING_PRECISION);

            e = ErrorUtils::isTrue(completeCluster.contains(mzHashed)); ree

            completeCluster[mzHashed] = ms2IonPeak;
        }

        *bestClusterPeaksComplete = completeCluster.values().toVector();

        const auto sortLogicIntensityAsc
                = [](const MS2IonPeak &l, const MS2IonPeak &r){return l.theoIntensity < r.theoIntensity;};

        std::sort(bestClusterPeaksComplete->rbegin(), bestClusterPeaksComplete->rend(), sortLogicIntensityAsc);

        ERR_RETURN
    }

    MS2IonPeak findBestMS2IonPeak(const QVector<MS2IonPeak> &ms2IonPeaksBestCluster) {

        const auto sortLogic = [](const MS2IonPeak &l, const MS2IonPeak &r){
            return l.cosineSimToAnchor < r.cosineSimToAnchor;
        };

        const MS2IonPeak frameIndexMaxMax = *std::max_element(
                ms2IonPeaksBestCluster.rbegin(),
                ms2IonPeaksBestCluster.rend(),
                sortLogic
        );

        return frameIndexMaxMax;
    }

    QPair<double, double> calculateFreqPercentFoundVsBest(const QVector<MS2IonPeak> &ms2IonPeaksBestCluster) {

        double frequencyPercentFound = 0.0;
        double frequencyPercentBestPossible = 0.0;

        for (const MS2IonPeak &ms2IonPeak : ms2IonPeaksBestCluster) {

            if (ms2IonPeak.pointCountFound > 0) {
                frequencyPercentFound += ms2IonPeak.fragmentFrequency;
            }

            frequencyPercentBestPossible += ms2IonPeak.fragmentFrequency;
        }

        return {frequencyPercentFound, frequencyPercentBestPossible};
    }

    Err extractMs2PeakToVectors(
            const QVector<MS2IonPeak> &ms2IonPeaksBestCluster,
            QVector<double> *mzSearchedVec,
            QVector<double> *theoIntensityVec,
            QVector<double> *mzFoundMeanVec,
            QVector<double> *mzFoundStDevVec,
            QVector<double> *intensityFoundMaxVec,
            QVector<int> *frameIndexMaxDiffFromAnchorVec,
            QVector<double> *cosineSimToAnchorVec,
            QVector<int> *peakPointCountFoundVec,
            QVector<double> *fragmentFrequencyVec,
            QVector<int> *rankVec
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2IonPeaksBestCluster); ree

        for (const MS2IonPeak &mip : ms2IonPeaksBestCluster) {
            mzSearchedVec->push_back(mip.mzSearched);
            theoIntensityVec->push_back(mip.theoIntensity);
            mzFoundMeanVec->push_back(mip.mzFoundMean);
            mzFoundStDevVec->push_back(mip.mzFoundStDev);
            intensityFoundMaxVec->push_back(*std::max_element(mip.intensityVals.begin(), mip.intensityVals.end()));
            frameIndexMaxDiffFromAnchorVec->push_back(mip.frameIndexMaxDiffFromAnchor);
            cosineSimToAnchorVec->push_back(mip.cosineSimToAnchor);
            peakPointCountFoundVec->push_back(mip.pointCountFound);
            fragmentFrequencyVec->push_back(mip.fragmentFrequency);
            rankVec->push_back(mip.rank);
        }

        ERR_RETURN
    }



}//namespace
Err MsFrameScoretron::scoreFrameCandidates(QVector<ScoredCandidate> *scoredCandidates) {

    ERR_INIT

    QVector<MS2IonPeak> ms2IonPeaks;
    e = buildMS2Peaks(&ms2IonPeaks); ree;

    QMap<PeptideStringWithMods, QPair<CosineSimSum, QVector<MS2IonPeak>>> bestCosineSimSumVsClusters;
    e = findBestClusterGroupingPerPeptideStringWithMods(
            ms2IonPeaks,
            m_params.minFoundMzPeaks,
            m_params.cosineSimThreshold,
            &bestCosineSimSumVsClusters
            ); ree

    for (auto it = bestCosineSimSumVsClusters.begin(); it != bestCosineSimSumVsClusters.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const QPair<CosineSimSum, QVector<MS2IonPeak>> &bestCluster = it.value();
        const CosineSimSum cosineSimSum = bestCluster.first;
        QVector<MS2IonPeak> ms2IonPeaksBestCluster = bestCluster.second;

        const MS2IonPeak anchorMS2IonPeak = findBestMS2IonPeak(ms2IonPeaksBestCluster);
        const CandidatePeptide &candidatePeptide = m_fragPredsTopN.value(peptideStringWithMods);

        //NOTE: this is placed here so that unfound peaks are not included in frameIndexmaxMean calculation
        e = fillUnFoundMS2IonPeaks(
                ms2IonPeaksBestCluster,
                candidatePeptide.ms2Ions,
                m_fragmentFrequencies,
                &ms2IonPeaksBestCluster
                ); ree

        const QPair<double, double> freqScores = calculateFreqPercentFoundVsBest(ms2IonPeaksBestCluster);

        ScoredCandidate sc;
        sc.peptideStringWithMods = peptideStringWithMods;
        sc.frequencyPercentSum = freqScores.first;
        sc.frequencyPercentSumBestPossible = freqScores.second;
        sc.cosineSimSum = cosineSimSum;
        sc.isDecoy = candidatePeptide.isDecoy;
        sc.mass = candidatePeptide.mass;
        sc.charge = candidatePeptide.charge;
        sc.scanNumber = anchorMS2IonPeak.scanNumberMax;
        sc.scanTime = m_scanNumberVsScanTime.value(sc.scanNumber);
        sc.theoreticalFragmentCount = candidatePeptide.totalFragmentCount;
        sc.iRTPredicted = candidatePeptide.iRt;
        sc.scanTimePredicted = m_fragPredsPredictedScanTime.value(peptideStringWithMods);
        sc.targetKey = m_uniqueMsInfoScanKey;

        e = extractMs2PeakToVectors(
                ms2IonPeaksBestCluster,
                &sc.mzSearchedVec,
                &sc.theoIntensityVec,
                &sc.mzFoundMeanVec,
                &sc.mzFoundStDevVec,
                &sc.intensityFoundMaxVec,
                &sc.frameIndexMaxDiffFromAnchorVec,
                &sc.cosineSimToAnchorVec,
                &sc.peakPointCountFoundVec,
                &sc.fragmentFrequencyVec,
                &sc.rankVec
        ); ree

        scoredCandidates->push_back(sc);
    }

//#define WRITE_SCAN_FRAME
#ifdef WRITE_SCAN_FRAME
    e = MsFrame::writeFrameScans(m_msFrame.scanNumberVsScanPoints(), "frame.prq"); ree
#endif

    ERR_RETURN
}

namespace {

    RTree loadHillsToRTree(
            const QMap<Id, FeatureFinderHill> &idVsFeatureFinderHills,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
            ) {

        std::vector<rTreePoint> cloudLoader;

        for (auto it = idVsFeatureFinderHills.begin(); it != idVsFeatureFinderHills.end(); it++) {

            const Id id = it.key();
            const FeatureFinderHill &ffh = it.value();

            const QPair<int, int> scanNumberMinMax = ffh.scanNumberMinMax();
            const double mz = ffh.mzMean();

            rTreeBox box(
                    {mz, scanNumberVsScanTime.value(scanNumberMinMax.first)},
                    {mz, scanNumberVsScanTime.value(scanNumberMinMax.second)}
            );

            cloudLoader.push_back(rTreePoint({box, id}));
        }

        const int maxElements = 16;
        return RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

    }

}//namespace
Err MsFrameScoretron::buildMS2Peaks(QVector<MS2IonPeak> *ms2IonPeaks) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragPredsTopN); ree;

    const QMap<Id, FeatureFinderHill> idVsFeatureFinderHills
        = ParallelUtils::convertVectorToMap(m_featureFinderHills);

    e = ErrorUtils::isNotEmpty(idVsFeatureFinderHills); ree;
    const RTree rTreeHills = loadHillsToRTree(
            idVsFeatureFinderHills,
            m_scanNumberVsScanTime
            );

    const double scanTimeBuffer = 5.0; //TODO make this settable.

    for (auto it = m_fragPredsTopN.begin(); it != m_fragPredsTopN.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();

        const QVector<MS2Ion> &ms2IonsTopN = it.value().ms2Ions;
        const ScanTime scanTime = m_fragPredsPredictedScanTime.value(peptideStringWithMods);

        for (const MS2Ion &ms2Ion : ms2IonsTopN) {

            const double mzTol = MathUtils::calculatePPM(ms2Ion.mz, m_params.ms2ExtractionWidthPPM);
            const double mzMin = ms2Ion.mz - mzTol;
            const double mzMax = ms2Ion.mz + mzTol;

            const rTreeBox queryBox(
                    rTreeCoor(mzMin, scanTime - scanTimeBuffer),
                    rTreeCoor(mzMax, scanTime + scanTimeBuffer)
            );

            std::vector<rTreePoint> rTreeSearchResult;
            rTreeHills.query(
                    bgi::intersects(queryBox),
                    std::back_inserter(rTreeSearchResult)
            );

            for (const rTreePoint &rtp : rTreeSearchResult) {

                const FeatureFinderHill &ffh = idVsFeatureFinderHills.value(rtp.second);

                MS2IonPeak ms2IonPeak;
                ms2IonPeak.peptideStringWithMods = peptideStringWithMods;
                ms2IonPeak.mzSearched = ms2Ion.mz;
                ms2IonPeak.theoIntensity = ms2Ion.intensity;
                ms2IonPeak.rank = ms2Ion.rank;
                ms2IonPeak.frameIndexMax = ffh.maxIntensityScanNumberIndex();
                ms2IonPeak.frameIndexStart = ffh.scanNumberIndexMinMax().first;
                ms2IonPeak.frameIndexEnd = ffh.scanNumberIndexMinMax().second;
                ms2IonPeak.scanNumberMax = ffh.maxIntensityScanNumber();
                ms2IonPeak.scanNumberStart = ffh.scanNumberMinMax().first;
                ms2IonPeak.scanNumberEnd = ffh.scanNumberMinMax().second;
                ms2IonPeak.intensityVals = ffh.intensities();

                ms2IonPeak.mzFoundMean = ffh.mzMean();
                ms2IonPeak.mzFoundStDev = ffh.mzStDev();
                ms2IonPeak.pointCountFound = ffh.scanCount();
                ms2IonPeak.fragmentFrequency = m_fragmentFrequencies.value(MathUtils::hashDecimal(
                        ms2IonPeak.mzSearched,
                        S_GLOBAL_SETTINGS.HASHING_PRECISION
                        ));

                ms2IonPeaks->push_back(ms2IonPeak);
            }

        }
    }

    ERR_RETURN
}
