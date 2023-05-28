//
// Created by anichols on 5/26/23.
//

#include "FeatureFinderHillClusterTron.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsUtils.h"
#include "ParallelUtils.h"

#include <Eigen/Dense>


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

Err FeatureFinderHillClusterTron::init(const FeatureFinderParameters &params) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid()); ree;

    m_params = params;

    ERR_RETURN
}

namespace {

    RTree buildHillsRTree(const QMap<int, FeatureFinderHill> &featureFinderHills) {

        std::vector<rTreePoint> cloudLoader;

        for (auto it = featureFinderHills.begin(); it != featureFinderHills.end(); it++) {

            const int id = it.key();
            const FeatureFinderHill &ffh = it.value();

            const QPair<int, int> frameIndexMinMax = ffh.scanNumberIndexMinMax();
            const double mz = ffh.mzMean();

            rTreeBox box(
                    {static_cast<double>(frameIndexMinMax.first), mz},
                    {static_cast<double>(frameIndexMinMax.second), mz}
            );

            cloudLoader.push_back(rTreePoint({box, id}));
        }

        const int maxElements = 16;
        RTree rTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

        return rTree;
    }

     Err buildComparisonMatrix(
            const FeatureFinderHill &h1,
            const FeatureFinderHill &h2,
            Eigen::MatrixX<double> *comparisonMatrix
            ) {

        ERR_INIT

        const QPair<int,int> h1ScanNumberIndexMinMax = h1.scanNumberIndexMinMax();
        const QPair<int,int> h2ScanNumberIndexMinMax = h2.scanNumberIndexMinMax();

        const int lowerScanNumberIndex = std::min(
             h1ScanNumberIndexMinMax.first,
             h2ScanNumberIndexMinMax.first
        );

        const int upperScanNumberIndex = std::max(
             h1ScanNumberIndexMinMax.second,
             h2ScanNumberIndexMinMax.second
        );

        const int colCount = 2;
        const int rowBuffer = 1;
        Eigen::MatrixX<double> mat(upperScanNumberIndex + rowBuffer, colCount);
        mat.setZero();

        const QVector<double> h1Intensities = h1.intensities();
        const QVector<double> h2Intensities = h2.intensities();
        const QVector<int> h1ScanNumberIndexes = h1.scanNumberIndexes();
        const QVector<int> h2ScanNumberIndexes = h2.scanNumberIndexes();

        e = ErrorUtils::isEqual(h1Intensities.size(), h1ScanNumberIndexes.size()); ree
        e = ErrorUtils::isEqual(h2Intensities.size(), h2ScanNumberIndexes.size()); ree

        const int insertColumn1 = 0;
        for (int i = 0; i < h1Intensities.size(); i++) {

            const int scanNumberIndex = h1ScanNumberIndexes.at(i);
            const double intensity = h1Intensities.at(i);

            mat.coeffRef(scanNumberIndex, insertColumn1) = intensity;
        }

        const int insertColumn2 = 1;
        for (int i = 0; i < h2Intensities.size(); i++) {

             const int scanNumberIndex = h2ScanNumberIndexes.at(i);
             const double intensity = h2Intensities.at(i);

             mat.coeffRef(scanNumberIndex, insertColumn2) = intensity;
        }

         *comparisonMatrix = mat.block(
                lowerScanNumberIndex,
                insertColumn1,
                (upperScanNumberIndex - lowerScanNumberIndex + rowBuffer),
                2
                );

        ERR_RETURN
    }

    Err calculateCosineSimBetweenHills(
            const FeatureFinderHill &h1,
            const FeatureFinderHill &h2,
            double *cosineSim
            ) {

        ERR_INIT

        Eigen::MatrixX<double> comparisonMatrix;
        e = buildComparisonMatrix(
                h1,
                h2,
                &comparisonMatrix
                ); ree

        *cosineSim = EigenUtils::cosineSimilarity(comparisonMatrix.col(0), comparisonMatrix.col(1));

        ERR_RETURN
    }

    Err removePointsBelowCosineSimCorrThreshold(
            const FeatureFinderHill &featureFinderHillAnchor,
            const QMap<int, FeatureFinderHill> &featureFinderHillMap,
            const std::vector<rTreePoint> &foundRTreeFeatureFinderHills,
            double cosineSimThreshold,
            QVector<FeatureFinderHill> *correlatedHills
            ) {

        ERR_INIT

        for (const rTreePoint &rtp : foundRTreeFeatureFinderHills) {

            const int &ffhId = rtp.second;

            double cosineSim;
            e = calculateCosineSimBetweenHills(
                    featureFinderHillAnchor,
                    featureFinderHillMap.value(ffhId),
                    &cosineSim
                    ); ree

            if (cosineSim < cosineSimThreshold) {
                continue;
            }

            correlatedHills->push_back(featureFinderHillMap.value(ffhId));
        }

        ERR_RETURN
    }

    ScanPoints convertHillApexesToScanPoints(const QVector<FeatureFinderHill> &featureFinderHills) {

        ScanPoints hillApexes;

        const auto transformLogic = [](const FeatureFinderHill &ffh){
            return ScanPoint (ffh.mzMean(), ffh.intensityValueMax());
        };

        std::transform(
                featureFinderHills.begin(),
                featureFinderHills.end(),
                std::back_inserter(hillApexes),
                transformLogic
                );

        return hillApexes;
    }

    rTreePoint findNearestRTreePoint(
            const std::vector<rTreePoint> &rTreeResults,
            const QMap<int, FeatureFinderHill> &featureFinderHillMap,
            double mz
            ) {

        for (const rTreePoint &rtp : rTreeResults){

            if (MathUtils::tZero(featureFinderHillMap.value(rtp.second).mzMean() - mz)) {
                return rtp;
            }
        }

        return {};
    }


}//namespace
Err FeatureFinderHillClusterTron::clusterHills(
        const QVector<FeatureFinderHill> &featureFinderHills,
        bool isMs2Hills,
        QVector<HillsClustering> *hillClusterings,
        QMap<int, FeatureFinderHill> *featureFinderHillsMap
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(featureFinderHills); ree;

    QVector<FeatureFinderHill> hillsSortedHiLoIntensity = featureFinderHills;
    FeatureFinderHillUtils::sortFeatureFinderHillsIntensityDesc(&hillsSortedHiLoIntensity);

    *featureFinderHillsMap = ParallelUtils::convertVectorToMap(featureFinderHills);

    RTree featureFinderHillRTree = buildHillsRTree(*featureFinderHillsMap);

    const double leftDistanceThomsons = 2.1;
    const double rightDistanceThomsons = 3.1;

    for (const FeatureFinderHill &ffh : featureFinderHills) {

        const int maxIntensityScanNumberIndex = ffh.maxIntensityScanNumberIndex();
        const double mzMean = ffh.mzMean();

        const rTreeBox queryBox(
                rTreeCoor(maxIntensityScanNumberIndex, mzMean - leftDistanceThomsons),
                rTreeCoor(maxIntensityScanNumberIndex, mzMean + rightDistanceThomsons)
        );

        std::vector<rTreePoint> rTreeSearchResult;
        featureFinderHillRTree.query(
                bgi::intersects(queryBox),
                std::back_inserter(rTreeSearchResult)
                );

        const rTreePoint rtpSearched = findNearestRTreePoint(
                rTreeSearchResult,
                *featureFinderHillsMap,
                ffh.mzMean()
                );

        if (!MathUtils::tZero(ffh.mzMean() - featureFinderHillsMap->value(rtpSearched.second).mzMean())) {
            continue;
        }

        if (rTreeSearchResult.size() == 1) {
            featureFinderHillRTree.remove(rtpSearched);

            if (isMs2Hills) {

                HillsClustering hillsClustering;
                hillsClustering.hillsMapIndexes.push_back(rtpSearched.second);
                hillClusterings->push_back(hillsClustering);
            }

            continue;
        }

        QVector<FeatureFinderHill> correlatedHills;
        e = removePointsBelowCosineSimCorrThreshold(
                ffh,
                *featureFinderHillsMap,
                rTreeSearchResult,
                m_params.cosineSimThreshold,
                &correlatedHills
        ); ree;

        const ScanPoints hillApexPoints = convertHillApexesToScanPoints(correlatedHills);

        const ScanPoint centerPoint = {ffh.mzMean(), ffh.intensityValueMax()};

        int charge;
        e = MsUtils::chargeDeterminator(
                centerPoint,
                hillApexPoints,
                m_params.tolerancePPM,
                m_params.chargeMin,
                m_params.chargeMax,
                &charge
                ); ree

        int monoOffset;
        ScanPoints subtractPoints;
        double bestCosineSim;
        e = MsUtils::monoIsotopeDeterminator(
                centerPoint,
                hillApexPoints,
                m_params.tolerancePPM,
                charge,
                &monoOffset,
                &subtractPoints,
                &bestCosineSim
                ); ree;

        HillsClustering hillCluster;
        hillCluster.monoOffset = monoOffset;
        hillCluster.cosineSim = bestCosineSim;
        hillCluster.charge = charge;
        for (const ScanPoint &sp : subtractPoints) {

            if (sp.x() < 1) {
                continue;
            }

            const rTreePoint rtp = findNearestRTreePoint(
                    rTreeSearchResult,
                    *featureFinderHillsMap,
                    sp.x()
                    );

            hillCluster.hillsMapIndexes.push_back(rtp.second);
            featureFinderHillRTree.remove(rtp);
        }

        hillClusterings->push_back(hillCluster);
    }

    ERR_RETURN
}
