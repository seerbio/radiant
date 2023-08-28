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
#include "ScanTimeFromIRtMapper.h"
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


namespace{

    PeakIntegratomaticParameters buildPeakIntegratomaticParams(const PythiaParameters  &pythiaParameters) {

        PeakIntegratomaticParameters params;
        params.smoothCount = pythiaParameters.smoothCount;
        params.filterLength = pythiaParameters.filterLength;
        params.sigma = pythiaParameters.sigma;
        params.signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;

        return params;
    }

}//namespace
Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
        const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanPoints); ree;
    e = ErrorUtils::isNotEmpty(peptideStringWithModsVsCandidatePeptide); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;

    m_params = params;
    m_fragPredsTopN = peptideStringWithModsVsCandidatePeptide;
    m_scanNumberVsScanTime = scanNumberVsScanTime;

    e = FragLibReader::generateFragmentFrequencies(
            m_fragPredsTopN,
            m_params.ms2ExtractionWidthPPM,
            &m_fragmentFrequencies
    ); ree

    m_msFrame.init(
            scanNumberVsScanPoints,
            scanNumberVsScanTime
            ); ree;

    const PeakIntegratomaticParameters ffParams = buildPeakIntegratomaticParams(m_params);
    e = m_peakIntegratomatic.init(ffParams); ree;

//    e = msFrame.deisotopeMsFrame(ppmTol); ree;
//    e = msFrame.smoothFrame(
//            ffParams.filterLength,
//            ffParams.sigma,
//            ffParams.smoothCount,
//            mzMax
//            ); ree;


    ERR_RETURN
}

Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
        const QMap<PeptideStringWithMods, CandidatePeptide> &peptideStringWithModsVsCandidatePeptide,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const QString &iRTRecalibrationFilePath
        ) {

    ERR_INIT

    e = init(
            params,
            scanNumberVsScanPoints,
            peptideStringWithModsVsCandidatePeptide,
            scanNumberVsScanTime
            ); ree;

    qDebug() << "updating rt vals from iRT";

    ScanTimeFromIRtMapper scanTimeFromIRtMapper;
    e = scanTimeFromIRtMapper.init(iRTRecalibrationFilePath); ree;

    for (auto it = peptideStringWithModsVsCandidatePeptide.begin(); it != peptideStringWithModsVsCandidatePeptide.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const double iRT = it.value().iRt;

        double predictedScanTime;
        e = scanTimeFromIRtMapper.predictScanTime(iRT, &predictedScanTime); ree;

        m_fragPredsPredictedScanTime.insert(peptideStringWithMods, predictedScanTime);
    }

    ERR_RETURN
}

namespace {

    QMap<PeptideStringWithMods, QVector<MS2IonPeak>> buildScanNumberVsMS2IonPeaks(const QVector<MS2IonPeak> &ms2IonPeaks) {

        QMap<PeptideStringWithMods , QVector<MS2IonPeak>> scanNumberVsMS2IonPeaks;
//        for (const MS2IonPeak &ms2IonPeak : ms2IonPeaks) {
//            scanNumberVsMS2IonPeaks[ms2IonPeak.peptideStringWithMods].push_back(ms2IonPeak);
//        }

        return scanNumberVsMS2IonPeaks;
    }

    RTree loadRTree(
            const QVector<MS2IonPeak> &ms2IonPeaks,
            QMap<Id, MS2IonPeak> *idVsMs2IonPeak,
            QSet<int> *frameIndexMaxesUnique
            ) {

        *idVsMs2IonPeak =  ParallelUtils::convertVectorToMap(ms2IonPeaks);

        std::vector<rTreePoint> cloudLoader;

//        for (auto it = idVsMs2IonPeak->begin(); it != idVsMs2IonPeak->end(); it++) {
//
//            const Id id = it.key();
//            const MS2IonPeak &ms2ip = it.value();
//
//            frameIndexMaxesUnique->insert(ms2ip.frameIndexMax);
//
//            rTreeBox box(
//                    {static_cast<double>(ms2ip.frameIndexStart), ms2ip.mzSearched},
//                    {static_cast<double>(ms2ip.frameIndexEnd), ms2ip.mzSearched}
//            );
//
//            cloudLoader.push_back(rTreePoint({box, id}));
//        }

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

        const QVector<double> p1Intensities = p1.intensityValsSmoothed;
        const QVector<double> p2Intensities = p2.intensityValsSmoothed;
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

//        QMap<Id, MS2IonPeak> idVsMs2IonPeak;
//        QSet<int> frameIndexMaxesUnique;
//        const RTree rtree = loadRTree(ms2ip, &idVsMs2IonPeak, &frameIndexMaxesUnique);
//
//        const QSet<QVector<Id>> uniqueIDs = buildUniqueIDs(
//                rtree,
//                mzPeaksFoundMin,
//                frameIndexMaxesUnique
//                );
//
//        if (uniqueIDs.isEmpty()) {
//            ERR_RETURN
//        }
//
//        QVector<QVector<MS2IonPeak>> clusteredMs2IonPeaksForCosineSimCalc;
//        const auto transformLogic = [&](const QVector<Id> &ids){
//            QVector<MS2IonPeak> peaks;
//            for (Id id : ids) {
//                peaks.push_back(idVsMs2IonPeak.value(id));
//            }
//            return peaks;
//        };
//        std::transform(
//                uniqueIDs.begin(),
//                uniqueIDs.end(),
//                std::back_inserter(clusteredMs2IonPeaksForCosineSimCalc),
//                transformLogic
//                );
//
//        const int maxRank = 5;
//
//        QVector<MS2IonPeak> bestMS2IonPeaksClustering;
//        double bestCosineSimSum = 0.0;
//
//        for (const QVector<MS2IonPeak> &ms2ips : clusteredMs2IonPeaksForCosineSimCalc) {
//
//            for (const MS2IonPeak &ms2IonPeakAnchor : ms2ips) {
//
//                if (ms2IonPeakAnchor.rank > maxRank) {
//                    continue;
//                }
//
//                double cosineSimSum = 0.0;
//                QVector<MS2IonPeak> clustering;
//
//                for (const MS2IonPeak &ms2IonPeak : ms2ips) {
//
//                    double cosineSim;
//                    e = calculateCosineSimBetweenPeaks(ms2IonPeakAnchor, ms2IonPeak, &cosineSim); ree
//
//                    if (cosineSim < cosineSimMinThreshold) {
//                        continue;
//                    }
//
//                    MS2IonPeak ms2IonPeakNew = ms2IonPeak;
//                    ms2IonPeakNew.cosineSimToAnchor = cosineSim;
//                    ms2IonPeakNew.frameIndexMaxDiffFromAnchor = ms2IonPeakAnchor.frameIndexMax - ms2IonPeak.frameIndexMax;
//
//                    clustering.push_back(ms2IonPeakNew);
//                    cosineSimSum += ms2IonPeakNew.cosineSimToAnchor;
//                }
//
//                if (cosineSimSum > bestCosineSimSum) {
//                    bestCosineSimSum = cosineSimSum;
//                    bestMS2IonPeaksClustering = clustering;
//                }
//            }
//        }
//
//        *bestMS2IonPeaksClusteringOutput = bestMS2IonPeaksClustering;
//        *bestCosineSimSumOutput = bestCosineSimSum;

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

//        e = ErrorUtils::isNotEmpty(fragPredTopNOfCluster); ree
//        e = ErrorUtils::isNotEmpty(bestClusterPeaks); ree
//
//        QMap<MzHashed, MS2IonPeak> completeCluster;
//
//        for (const MS2Ion &ms2Ion : fragPredTopNOfCluster) {
//            const MzHashed mzHashed
//                    = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
//
//            MS2IonPeak ms2IonPeakUnfound;
//            ms2IonPeakUnfound.mzSearched = ms2Ion.mz;
//            ms2IonPeakUnfound.theoIntensity = ms2Ion.intensity;
//            ms2IonPeakUnfound.pointCountFound = 0;
//            ms2IonPeakUnfound.fragmentFrequency
//                = fragFrequencies.value(MathUtils::hashDecimal(ms2IonPeakUnfound.mzSearched, S_GLOBAL_SETTINGS.HASHING_PRECISION));
//
//            completeCluster.insert(mzHashed, ms2IonPeakUnfound);
//        }
//
//        for (const MS2IonPeak &ms2IonPeak : bestClusterPeaks) {
//
//            const MzHashed &mzHashed
//                    = MathUtils::hashDecimal(ms2IonPeak.mzSearched, S_GLOBAL_SETTINGS.HASHING_PRECISION);
//
//            e = ErrorUtils::isTrue(completeCluster.contains(mzHashed)); ree
//
//            completeCluster[mzHashed] = ms2IonPeak;
//        }
//
//        *bestClusterPeaksComplete = completeCluster.values().toVector();
//
//        const auto sortLogicIntensityAsc
//                = [](const MS2IonPeak &l, const MS2IonPeak &r){return l.theoIntensity < r.theoIntensity;};
//
//        std::sort(bestClusterPeaksComplete->rbegin(), bestClusterPeaksComplete->rend(), sortLogicIntensityAsc);

        ERR_RETURN
    }

//    MS2IonPeak findBestMS2IonPeak(const QVector<MS2IonPeak> &ms2IonPeaksBestCluster) {
//
//        const auto sortLogic = [](const MS2IonPeak &l, const MS2IonPeak &r){
//            return l.cosineSimToAnchor < r.cosineSimToAnchor;
//        };
//
//        const MS2IonPeak frameIndexMaxMax = *std::max_element(
//                ms2IonPeaksBestCluster.rbegin(),
//                ms2IonPeaksBestCluster.rend(),
//                sortLogic
//        );
//
//        return frameIndexMaxMax;
//    }

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

//    Err extractMs2PeakToVectors(
//            const QVector<MS2IonPeak> &ms2IonPeaksBestCluster,
//            QVector<double> *mzSearchedVec,
//            QVector<double> *theoIntensityVec,
//            QVector<double> *mzFoundMeanVec,
//            QVector<double> *mzFoundStDevVec,
//            QVector<double> *intensityFoundMaxVec,
//            QVector<int> *frameIndexMaxDiffFromAnchorVec,
//            QVector<double> *cosineSimToAnchorVec,
//            QVector<int> *peakPointCountFoundVec,
//            QVector<double> *fragmentFrequencyVec,
//            QVector<int> *rankVec
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(ms2IonPeaksBestCluster); ree
//
//        for (const MS2IonPeak &mip : ms2IonPeaksBestCluster) {
//            mzSearchedVec->push_back(mip.mzSearched);
//            theoIntensityVec->push_back(mip.theoIntensity);
//            mzFoundMeanVec->push_back(mip.mzFoundMean);
//            mzFoundStDevVec->push_back(mip.mzFoundStDev);
//            intensityFoundMaxVec->push_back(*std::max_element(mip.intensityVals.begin(), mip.intensityVals.end()));
//            frameIndexMaxDiffFromAnchorVec->push_back(mip.frameIndexMaxDiffFromAnchor);
//            cosineSimToAnchorVec->push_back(mip.cosineSimToAnchor);
//            peakPointCountFoundVec->push_back(mip.pointCountFound);
//            fragmentFrequencyVec->push_back(mip.fragmentFrequency);
//            rankVec->push_back(mip.rank);
//        }
//
//        ERR_RETURN
//    }

}//namespace
Err MsFrameScoretron::scoreFrameCandidates(QVector<ScoredCandidate> *scoredCandidates) {

    ERR_INIT

    QMap<MzHashed, XICPoints> mzHashedVsXICPoints;
    QMap<MzHashed, QVector<double>> mzHashedVsIonPresence;
    e = buildMS2Peaks(&mzHashedVsXICPoints, &mzHashedVsIonPresence); ree;

    for (auto it = m_fragPredsTopN.begin(); it != m_fragPredsTopN.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const CandidatePeptide &candidatePeptide = it.value();

//#define TROUBLESHOOT_CANDIDATE
#ifdef TROUBLESHOOT_CANDIDATE
        qDebug() << "troubleshoot cand" << candidatePeptide.peptideStringWithMods;
        const PeptideStringWithMods peptideStringWithModsTroubleshoot = "SLDFKDSJL";
        if (peptideStringWithModsTroubleshoot != candidatePeptide.peptideStringWithMods) {
            continue;
        }
#endif

        e = processCandidate(
                candidatePeptide,
                mzHashedVsXICPoints,
                mzHashedVsIonPresence
                ); ree;

    }

//    QMap<PeptideStringWithMods, QPair<CosineSimSum, QVector<MS2IonPeak>>> bestCosineSimSumVsClusters;
//    e = findBestClusterGroupingPerPeptideStringWithMods(
//            ms2IonPeaks,
//            m_params.minFoundMzPeaks,
//            m_params.cosineSimThreshold,
//            &bestCosineSimSumVsClusters
//            ); ree
//
//    for (auto it = bestCosineSimSumVsClusters.begin(); it != bestCosineSimSumVsClusters.end(); it++) {
//
//        const PeptideStringWithMods &peptideStringWithMods = it.key();
//        const QPair<CosineSimSum, QVector<MS2IonPeak>> &bestCluster = it.value();
//        const CosineSimSum cosineSimSum = bestCluster.first;
//        QVector<MS2IonPeak> ms2IonPeaksBestCluster = bestCluster.second;
//
//        const MS2IonPeak anchorMS2IonPeak = findBestMS2IonPeak(ms2IonPeaksBestCluster);
//        const CandidatePeptide &candidatePeptide = m_fragPredsTopN.value(peptideStringWithMods);
//
//        //NOTE: this is placed here so that unfound peaks are not included in frameIndexmaxMean calculation
//        e = fillUnFoundMS2IonPeaks(
//                ms2IonPeaksBestCluster,
//                candidatePeptide.ms2Ions,
//                m_fragmentFrequencies,
//                &ms2IonPeaksBestCluster
//                ); ree
//
//        const QPair<double, double> freqScores = calculateFreqPercentFoundVsBest(ms2IonPeaksBestCluster);
//
//        ScoredCandidate sc;
//        sc.peptideStringWithMods = peptideStringWithMods;
//        sc.frequencyPercentSum = freqScores.first;
//        sc.frequencyPercentSumBestPossible = freqScores.second;
//        sc.cosineSimSum = cosineSimSum;
//        sc.isDecoy = candidatePeptide.isDecoy;
//        sc.mass = candidatePeptide.mass;
//        sc.charge = candidatePeptide.charge;
//        sc.scanNumber = anchorMS2IonPeak.scanNumberMax;
//        sc.scanTime = m_scanNumberVsScanTime.value(sc.scanNumber);
//        sc.theoreticalFragmentCount = candidatePeptide.totalFragmentCount;
//        sc.iRTPredicted = candidatePeptide.iRt;
//        sc.scanTimePredicted = m_fragPredsPredictedScanTime.value(peptideStringWithMods);
//        sc.targetKey = m_uniqueMsInfoScanKey;
//
//        e = extractMs2PeakToVectors(
//                ms2IonPeaksBestCluster,
//                &sc.mzSearchedVec,
//                &sc.theoIntensityVec,
//                &sc.mzFoundMeanVec,
//                &sc.mzFoundStDevVec,
//                &sc.intensityFoundMaxVec,
//                &sc.frameIndexMaxDiffFromAnchorVec,
//                &sc.cosineSimToAnchorVec,
//                &sc.peakPointCountFoundVec,
//                &sc.fragmentFrequencyVec,
//                &sc.rankVec
//        ); ree
//
//        scoredCandidates->push_back(sc);
//    }

//#define WRITE_SCAN_FRAME
#ifdef WRITE_SCAN_FRAME
    e = MsFrame::writeFrameScans(m_msFrame.scanNumberVsScanPoints(), "frame.prq"); ree
#endif

    ERR_RETURN
}

namespace {

    Err buildMzHashedVsMzIon(
            const QMap<PeptideStringWithMods, CandidatePeptide> &fragPredsTopN,
            QMap<MzHashed, MZION> *mzHashedVsMzIon
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(fragPredsTopN); ree;
        mzHashedVsMzIon->clear();

        for (auto it = fragPredsTopN.begin(); it != fragPredsTopN.end(); it++) {

            const QVector<MS2Ion> &ms2IonsTopN = it.value().ms2Ions;

            for (const MS2Ion &ms2Ion: ms2IonsTopN) {
                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                mzHashedVsMzIon->insert(mzHashed, ms2Ion.mz);
            }
        }

        ERR_RETURN
    }

    Err buildMzHashedVsXICPoints(
            const QMap<MzHashed, MZION> &mzHashedVsMzIon,
            const MsFrame &msFrame,
            double ppmTol,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints
            ) {

        ERR_INIT

        TurboXIC turboXic;
        e = turboXic.init(msFrame.frameIndexVsScanPoints()); ree;

        double scanNumberMin;
        double scanNumberMax;
        double mzMinRTree;
        double mzMaxRTree;
        e = turboXic.getRTreeLimits(
                &scanNumberMin,
                &scanNumberMax,
                &mzMinRTree,
                &mzMaxRTree
        ); ree;

        for (auto it = mzHashedVsMzIon.begin(); it != mzHashedVsMzIon.end(); it++) {

            const MzHashed mzHashed = it.key();
            const MZION mz = it.value();
            const double mzTol = MathUtils::calculatePPM(mz, ppmTol);

            const double mzMin = mz - mzTol;
            const double mzMax = mz + mzTol;

            const XICPoints xicPoints = turboXic.extractPointsXIC(
                    mzMin,
                    mzMax,
                    static_cast<int>(scanNumberMin),
                    static_cast<int>(scanNumberMax)
            );

            mzHashedVsXICPoints->insert(mzHashed, xicPoints);

        }

        ERR_RETURN
    }

    Err buildMzHashedVsXICPointsNormalized(
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            int scanCount,
            QMap<MzHashed, Eigen::VectorX<double>> *mzHashedVsXICPointsNormalized
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;

        for (auto it = mzHashedVsXICPoints.begin(); it != mzHashedVsXICPoints.end(); it++) {

            const MzHashed mzHashed = it.key();
            const XICPoints &xicPoints = it.value();

            if (xicPoints.scanNumbersVsIntensityVals.isEmpty()) {
                mzHashedVsXICPointsNormalized->insert(mzHashed, {});
                continue;
            }

            Eigen::VectorX<double> vecNormalized = EigenUtils::convertQMapToEigenVector(
                    xicPoints.scanNumbersVsIntensityVals,
                    scanCount
            );
            const double denom = vecNormalized.maxCoeff();
            vecNormalized = vecNormalized.array() / denom;

            mzHashedVsXICPointsNormalized->insert(mzHashed, vecNormalized);
        }

        ERR_RETURN
    }

    Err buildMzHashedVsIonPresence(
            const QMap<MzHashed, Eigen::VectorX<double>> &mzHashedVsXICPointsNormalized,
            QMap<MzHashed, QVector<double>> *mzHashedVsIonPresence
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzHashedVsXICPointsNormalized); ree;

        for (auto it = mzHashedVsXICPointsNormalized.begin(); it != mzHashedVsXICPointsNormalized.end(); it++) {

            const MzHashed mzHashed = it.key();
            const Eigen::VectorX<double> &normVec = it.value();

            Eigen::VectorX<double> vecPresence = normVec.array() / normVec.array();
            vecPresence = EigenUtils::setNANToZero(vecPresence);
            mzHashedVsIonPresence->insert(mzHashed, EigenUtils::convertEigenVectorToQVector(vecPresence));
        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::buildMS2Peaks(
        QMap<MzHashed, XICPoints> *mzHashedVsXICPoints,
        QMap<MzHashed, QVector<double>> *mzHashedVsIonPresence
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragPredsTopN); ree;

    QMap<MzHashed, MZION> mzHashedVsMzIon;
    e = buildMzHashedVsMzIon(m_fragPredsTopN, &mzHashedVsMzIon); ree;

    qDebug() << "Extracting" << mzHashedVsMzIon.size() << "XICs";

    e = buildMzHashedVsXICPoints(
            mzHashedVsMzIon,
            m_msFrame,
            m_params.ms2ExtractionWidthPPM,
            mzHashedVsXICPoints
            ); ree;

    QMap<MzHashed, Eigen::VectorX<double>> mzHashedVsXICPointsNormalized;
    e = buildMzHashedVsXICPointsNormalized(
            *mzHashedVsXICPoints,
            m_msFrame.scanCount(),
            &mzHashedVsXICPointsNormalized
            ); ree;

    e = buildMzHashedVsIonPresence(
            mzHashedVsXICPointsNormalized,
            mzHashedVsIonPresence
            );ree


    ERR_RETURN
}

namespace {

    int getMaxRowCount(
            const CandidatePeptide &candidatePeptide,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence
            ) {

        int maxRowCount = 0;
        for (const MS2Ion &ms2Ion : candidatePeptide.ms2Ions) {

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const QVector<double> &ionPresenceVec = mzHashedVsIonPresence.value(mzHashed);

            maxRowCount = std::max(maxRowCount, ionPresenceVec.size());
        }

        return maxRowCount;
    }

    Eigen::MatrixX<double> buildSummingMatrix(
            const CandidatePeptide &candidatePeptide,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
            int topNMs2Ions
            ) {

        const int cols = topNMs2Ions;
        const int rows = getMaxRowCount(
                candidatePeptide,
                mzHashedVsIonPresence
                );

        if (rows == 0) {
            return {};
        }

        Eigen::MatrixX<double> mat(rows, cols);
        mat.setZero();

        for (int i = 0; i < cols; i++) {

            if (i >= candidatePeptide.ms2Ions.size()) {
                break;
            }

            const MS2Ion &ms2Ion = candidatePeptide.ms2Ions.at(i);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const QVector<double> &ionPresenceVec = mzHashedVsIonPresence.value(mzHashed);

            if (ionPresenceVec.isEmpty()) {
                continue;
            }

            const Eigen::VectorX<double> eVec = EigenUtils::convertQVectorToEigenVector(ionPresenceVec);
            mat.col(i) = eVec;
        }

        return mat;
    }

    void filterSummedVecPeakIntegrations(
            const QVector<double> &summedMatToVec,
            int summedMzPresenceMin,
            int peakWidthMin,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
            ) {

        const auto terminatorLogic = [summedMatToVec, summedMzPresenceMin, peakWidthMin](const PeakIntegrationIndexes &pii){

            const int peakWidth = pii.second - pii.first;
            if (peakWidth < peakWidthMin) {
                return true;
            }

            const QVector<double> summedMatVecMax = summedMatToVec.mid(pii.first, peakWidth);
            const double summedPresenceIntegrationMax = *std::max_element(summedMatVecMax.begin(), summedMatVecMax.end());

            return summedPresenceIntegrationMax < summedMzPresenceMin;
        };

        const auto terminator = std::remove_if(peakIntegrationIndexes->begin(), peakIntegrationIndexes->end(), terminatorLogic);
        peakIntegrationIndexes->erase(terminator, peakIntegrationIndexes->end());
    }

    void sortPeakIntegrationsDescMaxSumFound(
            const QVector<double> &summedMatToVec,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
            ) {

        const auto sortLogic = [summedMatToVec](const PeakIntegrationIndexes &l, const PeakIntegrationIndexes &r){

            const int peakWidthL = l.second - l.first;
            const int peakWidthR = r.second - r.first;

            const QVector<double> summedMatVecMaxL = summedMatToVec.mid(l.first, peakWidthL);
            const QVector<double> summedMatVecMaxR = summedMatToVec.mid(r.first, peakWidthR);

            const double summedPresenceIntegrationMaxL = *std::max_element(summedMatVecMaxL.begin(), summedMatVecMaxL.end());
            const double summedPresenceIntegrationMaxR = *std::max_element(summedMatVecMaxR.begin(), summedMatVecMaxR.end());

            if (MathUtils::tSame(summedPresenceIntegrationMaxL, summedPresenceIntegrationMaxR)) {
                return std::accumulate(summedMatVecMaxL.begin(), summedMatVecMaxL.end(), 0.0)
                     < std::accumulate(summedMatVecMaxR.begin(), summedMatVecMaxR.end(), 0.0);
            }

            return summedPresenceIntegrationMaxL < summedPresenceIntegrationMaxR;
        };

        std::sort(peakIntegrationIndexes->rbegin(), peakIntegrationIndexes->rend(), sortLogic);
    }

    Eigen::MatrixX<double> buildIntensityVecMatrix(
            const CandidatePeptide &candidatePeptide,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            int vecSize,
            int topNMzFrags
    ) {

        const int cols = topNMzFrags;
        const int rows = vecSize;

        Eigen::MatrixX<double> mat(rows, cols);
        mat.setZero();

        for (int i = 0; i < cols; i++) {

            if (i >= candidatePeptide.ms2Ions.size()) {
                break;
            }

            const MS2Ion &ms2Ion = candidatePeptide.ms2Ions.at(i);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

            QVector<double> intensityVals = ParallelUtils::convertMapToVector(
                    mzHashedVsXICPoints.value(mzHashed).scanNumbersVsIntensityVals,
                    rows
                    );

            if (intensityVals.isEmpty()) {
                intensityVals = QVector<double>(vecSize, 0.0);
            }

            const Eigen::VectorX<double> eVec = EigenUtils::convertQVectorToEigenVector(intensityVals);
            mat.col(i) = eVec;
        }

        return mat;
    }

    Err calculateCandidateAllignementMetrics(
            const Eigen::MatrixX<double> &intensityMatrix,
            const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
            int topNMs2FragPeaks
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrix.nonZeros() > 1); ree;
        e = ErrorUtils::isNotEmpty(peakIntegrationIndexes); ree;

        for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {

            const int rowStart = pii.first;
            const int rowCount = pii.second - pii.first + 1;

            const int colStart = 0;
            const int colCount = topNMs2FragPeaks;

            Eigen::MatrixXd intensityMatrixIntegratedLimits = intensityMatrix.block(
                    rowStart,
                    colStart,
                    rowCount,
                    colCount
                    );

//            std::cout << intensityMatrixIntegratedLimits << std::endl;
//            std::cout << std::endl;
        }
//        std::cout << "*********" << std::endl;

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::processCandidate(
        const CandidatePeptide &candidatePeptide,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
        const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(mzHashedVsIonPresence); ree;
    e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;

    const Eigen::MatrixX<double> presenceMatrix = buildSummingMatrix(
            candidatePeptide,
            mzHashedVsIonPresence,
            m_params.topNMs2Ions
            );

    if (presenceMatrix.rows() == 0) {
        // TODO when you figue out what to return here, add it here for null result.
        ERR_RETURN
    }

    const Eigen::VectorX<double> summedMat = presenceMatrix.rowwise().sum();
    const QVector<double> summedMatToVec = EigenUtils::convertEigenVectorToQVector(summedMat);

    QVector<PeakIntegrationIndexes> peakIntegrationIndexes;
    QVector<double> summedPresenceVecSmoothed;
    e = m_peakIntegratomatic.findAllPeaksLimitsInXIC(
            summedMatToVec,
            &peakIntegrationIndexes,
            &summedPresenceVecSmoothed
    ); ree;

    const int peakWidthMin = 3;
    filterSummedVecPeakIntegrations(
            summedMatToVec,
            m_params.minFoundMzPeaks,
            peakWidthMin,
            &peakIntegrationIndexes
            );

    if (peakIntegrationIndexes.isEmpty()) {
        ERR_RETURN
    }

    sortPeakIntegrationsDescMaxSumFound(
            summedMatToVec,
            &peakIntegrationIndexes
            );

//#define TROUBLESHOOT_MATRICIES
#ifdef TROUBLESHOOT_MATRICIES
    for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {
        qDebug() << pii << summedMatToVec.mid(pii.first, pii.second - pii.first + 1);
    }
    qDebug() << "****";
#endif

    const Eigen::MatrixX<double> intensityMatrix = buildIntensityVecMatrix(
            candidatePeptide,
            mzHashedVsXICPoints,
            static_cast<int>(summedMat.size()),
            m_params.topNMs2Ions
    );

    e = calculateCandidateAllignementMetrics(
            intensityMatrix,
            peakIntegrationIndexes,
            m_params.topNMs2Ions
            ); ree;


    ERR_RETURN
}
