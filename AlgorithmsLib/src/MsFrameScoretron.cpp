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
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
) {

    ERR_INIT

    params.print();

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    e = ErrorUtils::isAboveThreshold(
            mzTargetStartStop.second,
            mzTargetStartStop.first,
            ErrorUtilsParam::ExcludeThreshold
    ); ree;

    m_params = params;
    m_msDataFilePath = msDataFilePath;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;
    m_mzTargetStartStop = mzTargetStartStop;

    e = FragLibReader::buildFragIonLibForCandidates(
            fragLibFilePath,
            m_params.chargeStateMin,
            m_params.chargeStateMax,
            m_mzTargetStartStop.first,
            m_mzTargetStartStop.second,
            &m_fragPreds,
            &m_fragPredsIsDecoy,
            &m_fragPredsMass,
            &m_fragPredsIRT
    ); ree

    e = loadFragPredsTopN(m_params.topNMs2Ions); ree

    e = initMsFrame(
        msDataFilePath,
        uniqueMsInfoScanKey,
        mzTargetStartStop
        ); ree

    e = m_msFrame.deisotopeMsFrame(m_params.ms2ExtractionWidthPPM);

    m_mzStartStopMean = (mzTargetStartStop.second + mzTargetStartStop.first) / 2.0;

    QMap<MzHashed, FrequencyPercent> fragmentFrequencies;
    e = FragLibReader::generateFragmentFrequencies(
            m_fragPredsTopN,
            m_params.ms2ExtractionWidthPPM,
            &fragmentFrequencies
            ); ree

    qDebug() << "TargetKey" << uniqueMsInfoScanKey;
    qDebug() << "Scan Count" << m_msFrame.scanCount();
    qDebug() << "Candidate Count:" << m_fragPreds.size();

    ERR_RETURN
}

namespace {

    Err buildNearestNeighborsIRTData(
            const QString &iRTRecalibrationFilePath,
            QVector<QPair<double, Coors>> *nnInputData
    ) {

        ERR_INIT

        e = ErrorUtils::fileExists(iRTRecalibrationFilePath); ree

        QVector<IRTReCalibrationRow> iRTReCalibrationReaderRows;
        e  = CSVReader::read(
                iRTRecalibrationFilePath,
                &iRTReCalibrationReaderRows
        ); ree;

        for (const IRTReCalibrationRow &row : iRTReCalibrationReaderRows) {
            nnInputData->push_back({row.scanTime, {row.iRT, 0.0}});
        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const QString &iRTRecalibrationFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
        ) {
    ERR_INIT

    e = init(
            params,
            msDataFilePath,
            fragLibFilePath,
            uniqueMsInfoScanKey,
            mzTargetStartStop
            ); ree;

    qDebug() << "updating rt vals from iRT";

    QVector<QPair<double, Coors>> nnInputData;
    e = buildNearestNeighborsIRTData(
            iRTRecalibrationFilePath,
            &nnInputData
    );

    NearestNeighborsSearch nearestNeighborsSearch;
    e = nearestNeighborsSearch.init(nnInputData); ree

    const int kNearestPoints = 10;

    for (auto it = m_fragPredsIRT.begin(); it != m_fragPredsIRT.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const double iRT = it.value();

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


Err MsFrameScoretron::initMsFrame(
        const QString &msDataFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
) {

    ERR_INIT

    e = MsFrame::buildMsFrame(
            msDataFilePath,
            uniqueMsInfoScanKey,
            mzTargetStartStop,
            &m_msFrame
    ); ree;

    e = m_msFrame.filterFrameByMz(
            m_params.mzMinDataStructure,
            m_params.mzMaxDataStructure
    ); ree;

    e = ErrorUtils::isAboveThreshold(
            m_msFrame.scanCount(),
            0,
            ErrorUtilsParam::ExcludeThreshold
    );ree;

    ERR_RETURN
}

Err MsFrameScoretron::loadFragPredsTopN(int n) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragPreds); ree;

    const double mzMin = 300.0;
    const double mzMax = 2000.0;

    for (auto it = m_fragPreds.begin(); it != m_fragPreds.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const MS2IonsSeparated &ms2IonsSeparated = it.value();

        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.yIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.bIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.y2Ions.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.b2Ions.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.aIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.yNH3Ions.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.yH2OIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.bNH3Ions.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.bH2OIons.values().toVector());
        m_fragPredsTopN[peptideStringWithMods].append(ms2IonsSeparated.precursorIons.values().toVector());

        QVector<MS2Ion> *ms2Ions = &m_fragPredsTopN[peptideStringWithMods];
        FragLibReader::getTopNMostIntenseMs2Ions(n, mzMin, mzMax, ms2Ions);
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

        QVector<MS2IonPeak> bestMS2IonPeaksClustering;
        double bestCosineSimSum = 0.0;

        for (const QVector<MS2IonPeak> &ms2ips : clusteredMs2IonPeaksForCosineSimCalc) {

            for (const MS2IonPeak &ms2IonPeakAnchor : ms2ips) {

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

                    clustering.push_back(ms2IonPeakNew);
                    cosineSimSum += cosineSim;
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
        QVector<MS2IonPeak> ms2IonPeaksMap = bestCluster.second;

        const double frameIndexMaxMean = std::accumulate(
                ms2IonPeaksMap.rbegin(),
                ms2IonPeaksMap.rend(),
                0.0,
                [](double sum, const MS2IonPeak &m){return sum + m.frameIndexMax;}) / ms2IonPeaksMap.size();

        const ScanNumber scanNumber
                = m_msFrame.scanNumberFromFrameIndex(static_cast<int>(std::round(frameIndexMaxMean)));

        const ScanTime scanTime = m_msFrame.scanTimeFromScanNumber(scanNumber);

        ScoredCandidate sc;
        sc.cosineSim = cosineSimSum;
        sc.peptideStringWithMods = peptideStringWithMods;
        sc.scanNumber = scanNumber;
        sc.scanTime = scanTime;
        sc.isDecoy = m_fragPredsIsDecoy.value(sc.peptideStringWithMods);

        scoredCandidates->push_back(sc);
    }


//    int falsers = 0;
//    int trueers = 0;
//    for (auto it = bestClusters.begin(); it != bestClusters.end(); it++) {
//
//        if (it.value().first < 6) {
//            continue;
//        }
//
//        qDebug() << "drewholio" << it.key() << it.value().first << m_fragPredsIsDecoy.value(it.key())
//                << it.value().second.size() << it.value().first / it.value().second.size();
//
//        if (m_fragPredsIsDecoy.value(it.key())) {
//            trueers++;
//        }
//        else{
//            falsers++;
//        }
//    }
//    qDebug() << "drewholio" << falsers << trueers;

//#define WRITE_MS2_ION_PEAKS
#ifdef WRITE_MS2_ION_PEAKS
    e = ParquetReader::write(ms2IonPeaks, "ms2IonPeaks.prq"); ree
#endif

//#define WRITE_SCAN_FRAME
#ifdef WRITE_SCAN_FRAME
    e = MsFrame::writeFrameScans(m_msFrame.scanNumberVsScanPoints(), "frame.prq"); ree
#endif

    ERR_RETURN
}

namespace {

    struct ConversionXICVecs {
        QVector<double> intensityVals;
        QVector<double> mzValsMeans;
        QVector<double> mzValsStDevs;
    };

    ConversionXICVecs convertXicLimitsToVec(
            const XICPoints &xicPoints,
            MsFrame *msFrame
            ) {

        const int buffer = 1;
        const int frameIndexVecSize = msFrame->frameIndexFromScanNumber(xicPoints.scanNumbersVsIntensityVals.lastKey()) + buffer;
        QVector<double> intensityVals(frameIndexVecSize, 0.0);
        QVector<double> mzValsMeans(frameIndexVecSize, 0.0);
        QVector<double> mzValsStDevs(frameIndexVecSize, 0.0);

        for (auto it = xicPoints.scanNumbersVsIntensityVals.begin(); it != xicPoints.scanNumbersVsIntensityVals.end(); it++) {

            const ScanNumber scanNumber = it.key();
            const FrameIndex frameIndex = msFrame->frameIndexFromScanNumber(scanNumber);
            const double intensity = it.value();

            intensityVals[frameIndex] = intensity;
        }

        for (auto it = xicPoints.scanNumberVsMzVals.begin(); it != xicPoints.scanNumberVsMzVals.end(); it++) {

            const ScanNumber scanNumber = it.key();
            const FrameIndex frameIndex = msFrame->frameIndexFromScanNumber(scanNumber);
            const QVector<double> &mzVals = it.value();

            mzValsMeans[frameIndex] = MathUtils::mean(mzVals);
            mzValsStDevs[frameIndex] = MathUtils::stDev(mzVals);
        }

        ConversionXICVecs cv;
        cv.intensityVals = intensityVals;
        cv.mzValsMeans = mzValsMeans;
        cv.mzValsStDevs = mzValsStDevs;

        return cv;
    }

}//namespace
Err MsFrameScoretron::buildMS2Peaks(QVector<MS2IonPeak> *ms2IonPeaks) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_fragPredsTopN); ree;

    TurboXIC turboXic;
    e = turboXic.init(m_msFrame.scanNumberVsScanPoints()); ree

    PeakIntegratomaticParameters params;
    params.smoothCount = 2;
    params.filterLength = 5;
    params.sigma = 1;
    params.signalToNoiseRatio = 1;

    PeakIntegratomatic peakIntegratomatic;
    e = peakIntegratomatic.init(params); ree;

    double scanNumberMinRTree;
    double scanNumberMaxRTree;
    double mzMinRTree;
    double mzMaxRTree;
    e = turboXic.getRTreeLimits(
            &scanNumberMinRTree,
            &scanNumberMaxRTree,
            &mzMinRTree,
            &mzMaxRTree
    );

    const double scanTimeBuffer = 5.0; //TODO make this settable.

    for (auto it = m_fragPredsTopN.begin(); it != m_fragPredsTopN.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const QVector<MS2Ion> &ms2IonsTopN = it.value();
        const ScanTime scanTime = m_fragPredsPredictedScanTime.value(peptideStringWithMods);

        const ScanNumber scanNumberMin = m_msFrame.scanNumberFromScanTime(scanTime - scanTimeBuffer);
        const ScanNumber scanNumberMax = m_msFrame.scanNumberFromScanTime(scanTime + scanTimeBuffer);

        for (const MS2Ion &ms2Ion : ms2IonsTopN) {

            const double mzTol = MathUtils::calculatePPM(ms2Ion.mz, m_params.ms2ExtractionWidthPPM);
            const double mzMin = ms2Ion.mz - mzTol;
            const double mzMax = ms2Ion.mz + mzTol;

            const XICPoints xicPoints = turboXic.extractPointsXIC(
                    mzMin,
                    mzMax,
                    scanNumberMin,
                    scanNumberMax
            );

            if (xicPoints.scanNumbersVsIntensityVals.isEmpty()) {
                continue;
            }

            const ConversionXICVecs cxv = convertXicLimitsToVec(
                    xicPoints,
                    &m_msFrame
                    );

            QVector<PeakIntegrationIndexes> peakLimits;
            QVector<double> intensityVecSmoothed;
            e = peakIntegratomatic.findAllPeaksLimitsInXIC(
                    cxv.intensityVals,
                    &peakLimits,
                    &intensityVecSmoothed
            ); ree;

            for (const PeakIntegrationIndexes &pii : peakLimits) {

                if (pii.second <= 0) {
                    continue;
                }

                e = ErrorUtils::isTrue(pii.second > pii.first); ree;

                const int pointsLength = pii.second - pii.first + 1;

                MS2IonPeak ms2IonPeak;
                ms2IonPeak.peptideStringWithMods = peptideStringWithMods;
                ms2IonPeak.mzSearched = ms2Ion.mz;
                ms2IonPeak.theoIntensity = ms2Ion.intensity;
                ms2IonPeak.frameIndexStart = pii.first;
                ms2IonPeak.frameIndexEnd = pii.second;
                ms2IonPeak.intensityVals = cxv.intensityVals.mid(pii.first, pointsLength);
                ms2IonPeak.frameIndexMax = MathUtils::findMaxIndexInVector(ms2IonPeak.intensityVals) + pii.first;
                ms2IonPeak.mzFoundMean = MathUtils::mean(cxv.mzValsMeans.mid(pii.first, pointsLength));
                ms2IonPeak.mzFoundStDev = MathUtils::stDev(cxv.mzValsMeans.mid(pii.first, pointsLength));

                if (ms2IonPeak.intensityVals.isEmpty()) {
                    continue;
                }

                ms2IonPeaks->push_back(ms2IonPeak);
            }
        }
    }

    ERR_RETURN
}



//namespace{
//
//    Err buildScanFragPreds(
//            const QVector<ScoredCandidate> &scoredCandidatesTargets,
//            QMap<PeptideStringWithMods, QVector<MS2Ion>> *fragPredsFlattened,
//            QMap<PeptideStringWithMods, QVector<MS2Ion>> *scanFragPredsFlattened
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(*fragPredsFlattened); ree;
//
//        for (const ScoredCandidate &sc : scoredCandidatesTargets) {
//
//            e = ErrorUtils::isTrue(fragPredsFlattened->contains(sc.peptideStringWithMods)); ree;
//
//            const QVector<MS2Ion> &pred = fragPredsFlattened->value(sc.peptideStringWithMods);
//
//            scanFragPredsFlattened->insert(sc.peptideStringWithMods, pred);
//        }
//
//        ERR_RETURN
//    }
//
//    void removeNegativeDiscScores(QMap<PeptideStringWithMods, TandemDeconvolverResult> *pepSeqVsWeight) {
//
//        QMap<PeptideStringWithMods, TandemDeconvolverResult> pepSeqVsWeightFiltered;
//        for (auto it = pepSeqVsWeight->begin(); it != pepSeqVsWeight->end(); it++) {
//
//            const PeptideStringWithMods &peptideStringWithMods = it.key();
//            const TandemDeconvolverResult &res = it.value();
//
//            if (res.discScore < 0 || res.pVal > 0.05) {
//                continue;
//            }
//
//            pepSeqVsWeightFiltered.insert(peptideStringWithMods, res);
//        }
//
//        *pepSeqVsWeight = pepSeqVsWeightFiltered;
//    }
//
//}//namespace
//Err MsFrameScoretron::deconvolveScan(
//        const QVector<ScoredCandidate> &scoredCandidatesTargets,
//        const ScanPoints &scanPointsDeisotoped,
//        QVector<ScoredCandidate> *scoredCandidatesTargetsFiltered
//        ) {
//
//    ERR_INIT
//
//    QMap<PeptideStringWithMods, QVector<MS2Ion>> scanFragPredsTopN;
//    e = buildScanFragPreds(
//            scoredCandidatesTargets,
//            &m_fragPredsTopN,
//            &scanFragPredsTopN
//            ); ree
//
//    const int precision = 2;
//    const double mzMax = 1400.0;
//    const int iters = 4;
//    const double stopTol = 1e-8;
//    const double pVal = 0.05;
//
//    if (scanFragPredsTopN.isEmpty()) {
//        ERR_RETURN
//    }
//
//    TandemSpectraDeconvolvotron decon;
//    e = decon.init(
//            precision,
//            mzMax,
//            iters,
//            stopTol,
//            pVal
//            ); ree
//
//    QMap<PeptideStringWithMods, TandemDeconvolverResult> pepSeqVsWeight;
//    e = decon.deconvolveTandemSpectra(scanPointsDeisotoped, scanFragPredsTopN, &pepSeqVsWeight); ree
//
//    removeNegativeDiscScores(&pepSeqVsWeight);
//
//    for (const ScoredCandidate &sc : scoredCandidatesTargets) {
//
//        if (!pepSeqVsWeight.contains(sc.peptideStringWithMods)) {
//            continue;
//        }
//
//        const TandemDeconvolverResult &tdr = pepSeqVsWeight.value(sc.peptideStringWithMods);
//
//        ScoredCandidate scNew = sc;
//
//        scNew.discScore = tdr.discScore;
//        scNew.pVal = tdr.pVal;
//        scNew.frameError = tdr.frameError;
//        scNew.tTestVal = tdr.tTestVal;
//
//        scoredCandidatesTargetsFiltered->push_back(scNew);
//    }
//
//    ERR_RETURN
//}
