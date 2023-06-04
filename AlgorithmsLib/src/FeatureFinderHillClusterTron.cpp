//
// Created by anichols on 5/26/23.
//

#include "FeatureFinderHillClusterTron.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FragLibReader.h"
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


FeatureFinderHillClusterTron::FeatureFinderHillClusterTron()
: m_chargeMin(1)
, m_chargeMax(5)
{}

Err FeatureFinderHillClusterTron::init(
        const FeatureFinderParameters &params,
        int chargeMin,
        int chargeMax
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid()); ree

    const double minCosineSim = 0.0;
    const double maxCosineSim = 1.0;
    e = ErrorUtils::isWithinRange(
            params.cosineSimThreshold,
            minCosineSim,
            maxCosineSim,
            ErrorUtilsParam::IncludeThreshold
            ); ree

    e = ErrorUtils::isTrue(chargeMin < chargeMax); ree
    e = ErrorUtils::isTrue(chargeMin > 0); ree
    e = ErrorUtils::isTrue(chargeMax > 0); ree

    m_params = params;

    ERR_RETURN
}

Err FeatureFinderHillClusterTron::init(const FeatureFinderParameters &params) {

    ERR_INIT

    e = init(params, 1, 666); ree

    ERR_RETURN
}


namespace {

    QVector<FeatureFinderHillPlus> buildFeatureFinderHillPlussesSortedHiLo(
            const QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> &featureFinderHills
    ) {

        QVector<FeatureFinderHillPlus> featureFinderHillsPlusses;

        for (auto it = featureFinderHills.begin(); it != featureFinderHills.end(); it++) {

            const IonType &ionType = it.key();
            const QMap<IonIndex, QVector<FeatureFinderHill>> &hills = it.value();

            for (auto itt = hills.begin(); itt != hills.end(); itt++) {

                const IonIndex ionIndex = itt.key();
                const QVector<FeatureFinderHill> &ffhs = itt.value();

                for (const FeatureFinderHill &ffh : ffhs) {

                    FeatureFinderHillPlus ffhp;
                    ffhp.featureFinderHill = ffh;
                    ffhp.ionType = ionType;
                    ffhp.ionIndex = ionIndex;

                    featureFinderHillsPlusses.push_back(ffhp);
                }

            }

        }

        FeatureFinderHillUtils::sortFeatureFinderHillsPlussesIntensityDesc(&featureFinderHillsPlusses);

        return featureFinderHillsPlusses;
    }

    RTree buildHillsRTree(
            const QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> &featureFinderHills,
            QMap<Id, FeatureFinderHillPlus> *featureFinderHillsPlussesMap
            ) {

        QVector<FeatureFinderHillPlus> hillsSortedHiLoIntensity = buildFeatureFinderHillPlussesSortedHiLo(featureFinderHills);

        *featureFinderHillsPlussesMap = ParallelUtils::convertVectorToMap(hillsSortedHiLoIntensity);

        std::vector<rTreePoint> cloudLoader;

        for (auto it = featureFinderHillsPlussesMap->begin(); it != featureFinderHillsPlussesMap->end(); it++) {

            const int id = it.key();
            const FeatureFinderHillPlus &ffhp = it.value();

            const QPair<int, int> frameIndexMinMax = ffhp.featureFinderHill.scanNumberIndexMinMax();
            const double mz = ffhp.featureFinderHill.mzMean();

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
                colCount
                );

        ERR_RETURN
    }

//    Err removePointsBelowCosineSimCorrThreshold(
//            const FeatureFinderHill &featureFinderHillAnchor,
//            const QMap<int, FeatureFinderHill> &featureFinderHillMap,
//            const std::vector<rTreePoint> &foundRTreeFeatureFinderHills,
//            double cosineSimThreshold,
//            QVector<FeatureFinderHill> *correlatedHills
//            ) {
//
//        ERR_INIT
//
//        for (const rTreePoint &rtp : foundRTreeFeatureFinderHills) {
//
//            const int &ffhId = rtp.second;
//
//            double cosineSim;
//            e = FeatureFinderHillClusterTron::calculateCosineSimBetweenHills(
//                    featureFinderHillAnchor,
//                    featureFinderHillMap.value(ffhId),
//                    &cosineSim
//                    ); ree
//
//            if (cosineSim < cosineSimThreshold) {
//                continue;
//            }
//
//            correlatedHills->push_back(featureFinderHillMap.value(ffhId));
//        }
//
//        ERR_RETURN
//    }

//    ScanPoints convertHillApexesToScanPoints(const QVector<FeatureFinderHillPlus> &featureFinderHills) {
//
//        ScanPoints hillApexes;
//
//        const auto transformLogic = [](const FeatureFinderHillPlus &ffh){
//            return ScanPoint (ffh.featureFinderHill.mzMean(), ffh.featureFinderHill.intensityValueMax());
//        };
//
//        std::transform(
//                featureFinderHills.begin(),
//                featureFinderHills.end(),
//                std::back_inserter(hillApexes),
//                transformLogic
//                );
//
//        return hillApexes;
//    }
//
//    rTreePoint findNearestRTreePoint(
//            const std::vector<rTreePoint> &rTreeResults,
//            const QMap<int, FeatureFinderHillPlus> &featureFinderHillPlussesMap,
//            double mz
//            ) {
//
//        for (const rTreePoint &rtp : rTreeResults){
//
//            if (MathUtils::tZero(featureFinderHillPlussesMap.value(rtp.second).featureFinderHill.mzMean() - mz)) {
//                return rtp;
//            }
//        }
//
//        return {};
//    }


}//namespace
Err FeatureFinderHillClusterTron::clusterHillsMS1(
        const QVector<FeatureFinderHill> &featureFinderHills,
        bool isMs2Hills,
        QVector<HillsClustering> *hillClusterings,
        QMap<int, FeatureFinderHill> *featureFinderHillsMap
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(featureFinderHills); ree;

    QVector<FeatureFinderHill> hillsSortedHiLoIntensity = featureFinderHills;
//    FeatureFinderHillUtils::sortFeatureFinderHillsIntensityDesc(&hillsSortedHiLoIntensity);
//
//    *featureFinderHillsMap = ParallelUtils::convertVectorToMap(featureFinderHills);
//
//    RTree featureFinderHillRTree = buildHillsRTree(*featureFinderHillsMap);
//
//    const double leftDistanceThomsons = 2.1;
//    const double rightDistanceThomsons = 3.1;
//
//    for (const FeatureFinderHill &ffh : hillsSortedHiLoIntensity) {
//
//        const int maxIntensityScanNumberIndex = ffh.maxIntensityScanNumberIndex();
//        const double mzMean = ffh.mzMean();
//
//        const rTreeBox queryBox(
//                rTreeCoor(maxIntensityScanNumberIndex, mzMean - leftDistanceThomsons),
//                rTreeCoor(maxIntensityScanNumberIndex, mzMean + rightDistanceThomsons)
//        );
//
//        std::vector<rTreePoint> rTreeSearchResult;
//        featureFinderHillRTree.query(
//                bgi::intersects(queryBox),
//                std::back_inserter(rTreeSearchResult)
//                );
//
//        const rTreePoint rtpSearched = findNearestRTreePoint(
//                rTreeSearchResult,
//                *featureFinderHillsMap,
//                ffh.mzMean()
//                );
//
//        if (!MathUtils::tZero(ffh.mzMean() - featureFinderHillsMap->value(rtpSearched.second).mzMean())) {
//            continue;
//        }
//
//        if (rTreeSearchResult.size() == 1) {
//            featureFinderHillRTree.remove(rtpSearched);
//
//            if (isMs2Hills) {
//
//                HillsClustering hillsClustering;
//                hillsClustering.hillsMapIndexes.push_back(rtpSearched.second);
//                hillClusterings->push_back(hillsClustering);
//            }
//
//            continue;
//        }
//
//        QVector<FeatureFinderHill> correlatedHills;
//        e = removePointsBelowCosineSimCorrThreshold(
//                ffh,
//                *featureFinderHillsMap,
//                rTreeSearchResult,
//                m_cosineSimThreshold,
//                &correlatedHills
//        ); ree;
//
//        const ScanPoints hillApexPoints = convertHillApexesToScanPoints(correlatedHills);
//
//        const ScanPoint centerPoint = {ffh.mzMean(), ffh.intensityValueMax()};
//
//        int charge;
//        e = MsUtils::chargeDeterminator(
//                centerPoint,
//                hillApexPoints,
//                m_params.tolerancePPM,
//                m_chargeMin,
//                m_chargeMax,
//                &charge
//                ); ree
//
//        int monoOffset;
//        ScanPoints subtractPoints;
//        double bestCosineSim;
//        e = MsUtils::monoIsotopeDeterminator(
//                centerPoint,
//                hillApexPoints,
//                m_params.tolerancePPM,
//                charge,
//                &monoOffset,
//                &subtractPoints,
//                &bestCosineSim
//                ); ree;
//
//        HillsClustering hillCluster;
//        hillCluster.monoOffset = monoOffset;
//        hillCluster.cosineSim = bestCosineSim;
//        hillCluster.charge = charge;
//        for (const ScanPoint &sp : subtractPoints) {
//
//            if (sp.x() < 1) {
//                continue;
//            }
//
//            const rTreePoint rtp = findNearestRTreePoint(
//                    rTreeSearchResult,
//                    *featureFinderHillsMap,
//                    sp.x()
//                    );
//
//            hillCluster.hillsMapIndexes.push_back(rtp.second);
//            featureFinderHillRTree.remove(rtp);
//        }
//
//        hillClusterings->push_back(hillCluster);
//    }

    ERR_RETURN
}

Err FeatureFinderHillClusterTron::writeClustersToMzRt(
        const QMap<ScanNumber, double> &scanNumberVsScanTime,
        const QVector<HillsClustering> &hillClustersByIndexs,
        const QMap<int, FeatureFinderHill> &featureFinderHillsMap,
        const QString &destinationFilePath
        ) {

    ERR_INIT

    struct MzRtRow : public CSVReaderInputBase {
        double mzLo = -1.0;
        double mzHi = -1.0;
        double rtLo = -1.0;
        double rtHi = -1.0;

        QMap<QString, QVariant> map() override {
            return {
                    {"mzLo", mzLo},
                    {"mzHi", mzHi},
                    {"rtLo", rtLo},
                    {"rtHi", rtHi}
            };
        }
    };

    QVector<MzRtRow> mzRtRows;
    for (const HillsClustering &hc : hillClustersByIndexs) {

        const QVector<int> &hillIndexes = hc.hillsMapIndexes;

        const auto insertLogic = [featureFinderHillsMap](int ind){
            return featureFinderHillsMap.value(ind);
        };

        QVector<FeatureFinderHill> ffhs;
        std::transform(
                hillIndexes.begin(),
                hillIndexes.end(),
                std::back_inserter(ffhs),
                insertLogic
                );

        const FeatureFinderHill maxIntensityHill = *std::max_element(
                ffhs.rbegin(),
                ffhs.rend(),
                [](const FeatureFinderHill &l, const FeatureFinderHill &r){
                    return l.intensityValueMax() < r.intensityValueMax();
                }
        );

        const auto minMaxMzHills = std::minmax_element(
                ffhs.begin(),
                ffhs.end(),
                [](const FeatureFinderHill &l, const FeatureFinderHill &r){
                    return l.mzMean() < r.mzMean();
                }
        );

        MzRtRow mzRtRow;
        mzRtRow.mzLo = minMaxMzHills.first->mzMean() - 0.025;
        mzRtRow.mzHi = minMaxMzHills.second->mzMean() + 0.025;
        mzRtRow.rtLo = scanNumberVsScanTime.value(maxIntensityHill.scanNumberMinMax().first);
        mzRtRow.rtHi = scanNumberVsScanTime.value(maxIntensityHill.scanNumberMinMax().second);

        mzRtRows.push_back(mzRtRow);
    }

    e = CSVReader::write(
            mzRtRows,
            destinationFilePath
    ); ree

    ERR_RETURN

}


namespace {

    std::vector<rTreePoint> searchRTree(
            double minScanNumberIndex,
            double maxScanNumberIndex,
            double mzSearchMin,
            double mzSearchMax,
            const RTree &featureFinderHillRTree,
            const QMap<Id, FeatureFinderHillPlus> &featureFinderHillsPlussesMap
            ) {

        const rTreeBox queryBox(
                rTreeCoor(minScanNumberIndex, mzSearchMin),
                rTreeCoor(maxScanNumberIndex, mzSearchMax)
        );

        std::vector<rTreePoint> rTreeSearchResult;
        featureFinderHillRTree.query(
                bgi::intersects(queryBox),
                std::back_inserter(rTreeSearchResult)
        );

        return rTreeSearchResult;
    }

    Err removePointsBelowCosineSimCorrThresholdMS2(
            const FeatureFinderHillPlus &featureFinderHillAnchor,
            const QMap<int, FeatureFinderHillPlus> &featureFinderHillMap,
            const std::vector<rTreePoint> &foundRTreeFeatureFinderHills,
            double cosineSimThreshold,
            QVector<FeatureFinderHillPlus> *correlatedHills
    ) {

        ERR_INIT

        using MZKEY = int;
        using BEST_COSINE_SIM = double;

        QMap<MZKEY, BEST_COSINE_SIM> mzKeyVsBestCosineSim;
        QMap<MZKEY, FeatureFinderHillPlus> mzKeyVsBestFeatureFinderHillPlus;

        for (const rTreePoint &rtp : foundRTreeFeatureFinderHills) {

            const int &ffhpId = rtp.second;
            const FeatureFinderHillPlus &ffhp = featureFinderHillMap.value(ffhpId);

            double cosineSim;
            e = FeatureFinderHillClusterTron::calculateCosineSimBetweenHills(
                    featureFinderHillAnchor.featureFinderHill,
                    ffhp.featureFinderHill,
                    &cosineSim
            ); ree

            if (cosineSim < cosineSimThreshold) {
                continue;
            }

            const int mzKey = MathUtils::hashDecimal(ffhp.featureFinderHill.mzMean(), S_GLOBAL_SETTINGS.HASHING_PRECISION);

            if (cosineSim > mzKeyVsBestCosineSim.value(mzKey)) {
                mzKeyVsBestCosineSim[mzKey] = cosineSim;
                mzKeyVsBestFeatureFinderHillPlus[mzKey] = ffhp;
            }

        }

        const QList<double> &bestCosineSims = mzKeyVsBestCosineSim.values();
        const QList<MZKEY> &mzKeys = mzKeyVsBestCosineSim.keys();

        const QList<FeatureFinderHillPlus> &bestFeatureFinderHillsPlus = mzKeyVsBestFeatureFinderHillPlus.values();

        const auto transformLogic = [&](const FeatureFinderHillPlus &ffhp){

            FeatureFinderHillPlus ffhpNew = ffhp;
            const int mzKey = MathUtils::hashDecimal(
                    ffhp.featureFinderHill.mzMean(),
                    S_GLOBAL_SETTINGS.HASHING_PRECISION
            );

            ffhpNew.cosineSimToAnchor = mzKeyVsBestCosineSim.value(mzKey);
            return ffhpNew;
        };

        std::transform(
                bestFeatureFinderHillsPlus.begin(),
                bestFeatureFinderHillsPlus.end(),
                std::back_inserter(*correlatedHills),
                transformLogic
                );

        ERR_RETURN
    }

    Err findBestMS2IonAnchor(
            const QVector<MS2Ion> &ms2IonsAnchors,
            const QMap<Id, FeatureFinderHillPlus> &featureFinderHillsPlussesMap,
            const RTree &featureFinderHillRTree,
            double ppmTolerance,
            QVector<FeatureFinderHillPlus> *ffhp
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2IonsAnchors); ree
        e = ErrorUtils::isFalse(featureFinderHillRTree.empty()); ree

        const double minScanNumberIndex = 0.0;
        const double maxScanNumberIndex = featureFinderHillRTree.bounds().max_corner().get<0>();

        const auto transformLogic
            = [featureFinderHillsPlussesMap](const rTreePoint &p){return featureFinderHillsPlussesMap.value(p.second);};

        for (const MS2Ion &potentialMS2Anchor : ms2IonsAnchors) {

            const double mzPotentialAnchor = potentialMS2Anchor.mz;

            const double massTolerance = MathUtils::calculatePPM(
                    mzPotentialAnchor,
                    ppmTolerance
                    );

            const double mzSearchMin = mzPotentialAnchor - massTolerance;
            const double mzSearchMax = mzPotentialAnchor + massTolerance;

            std::vector<rTreePoint> searchResultPoints = searchRTree(
                    minScanNumberIndex,
                    maxScanNumberIndex,
                    mzSearchMin,
                    mzSearchMax,
                    featureFinderHillRTree,
                    featureFinderHillsPlussesMap
                    );

            if (searchResultPoints.empty()) {
                continue;
            }

            std::transform(
                    searchResultPoints.begin(),
                    searchResultPoints.end(),
                    std::back_inserter(*ffhp),
                    transformLogic
                    );

            ERR_RETURN
        }

        ERR_RETURN
    }

    Err getClosestMS2Ion(
            const QVector<MS2Ion> &ms2Ions,
            const FeatureFinderHillPlus &ffhp,
            MS2Ion *foundIon
            ) {

        ERR_INIT

        auto it = std::min_element(ms2Ions.begin(), ms2Ions.end(), [ffhp] (const MS2Ion &a, const MS2Ion &b) {
            return std::abs(a.mz - ffhp.featureFinderHill.mzMean()) <= std::abs(b.mz - ffhp.featureFinderHill.mzMean());
        });

        if(it == ms2Ions.end()) {
            rrr(eValueError);
        }

        *foundIon = *it;

        ERR_RETURN
    }

    Err removeDuplicateCorrelatedHills(
            const MS2IonsSeparated &ms2IonsTandemPred,
            QVector<FeatureFinderHillPlus> *correlatedHills
            ) {

        ERR_INIT

        const int ionIndexThreshold = 0;
        const int topX = 500;
        const QVector<MS2Ion> ms2Ions = ms2IonsTandemPred.getTopXMS2Ions(topX, ionIndexThreshold);

        using MZIONKEY = int;
        QMap<MZIONKEY, QVector<FeatureFinderHillPlus>> mzKeyVsFeatureFinderHillPlusses;

        for (const FeatureFinderHillPlus &ffhp : *correlatedHills) {

            MS2Ion closestMS2Ion;
            e = getClosestMS2Ion(
                    ms2Ions,
                    ffhp,
                    &closestMS2Ion
            ); ree

            const MZIONKEY mzKey = MathUtils::hashDecimal(closestMS2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

            mzKeyVsFeatureFinderHillPlusses[mzKey].push_back(ffhp);
        }

        QVector<FeatureFinderHillPlus> correlatedHillsFiltered;

        const auto cosineSimSortLogic = [](const FeatureFinderHillPlus &l, const FeatureFinderHillPlus &r){
            return l.cosineSimToAnchor < r.cosineSimToAnchor;
        };

        for (const QVector<FeatureFinderHillPlus> &ffhps : mzKeyVsFeatureFinderHillPlusses) {

            if (ffhps.size() == 1) {
                correlatedHillsFiltered.push_back(ffhps.front());
                continue;
            }

            const FeatureFinderHillPlus ffhpMaxCosineSim = *std::max_element(ffhps.begin(), ffhps.end(), cosineSimSortLogic);
            correlatedHillsFiltered.push_back(ffhpMaxCosineSim);
        }

        *correlatedHills = correlatedHillsFiltered;

        ERR_RETURN
    }

    double calculateWeightedHillCosineSimSum(const QVector<FeatureFinderHillPlus> &correlatedHills) {

        const auto accLogic = [](double weightedSum, const FeatureFinderHillPlus &ffhp){
            return weightedSum + (ffhp.cosineSimToAnchor * ffhp.featureFinderHill.mzMean());
        };

        return std::accumulate(correlatedHills.begin(), correlatedHills.end(), 0.0, accLogic);
    }


}//namespace
Err FeatureFinderHillClusterTron::clusterHillsByBestMS2IonAnchor(
        const QMap<IonType, QMap<IonIndex, QVector<FeatureFinderHill>>> &featureFinderHills,
        const QVector<MS2Ion> &ms2IonsAnchors,
        const MS2IonsSeparated &ms2IonsTandemPred,
        HillsClusteringMS2 *bestHillsClusteringMS2
        ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(featureFinderHills); ree
    e = ErrorUtils::isNotEmpty(ms2IonsAnchors); ree

    bestHillsClusteringMS2->cosineSimSumWeighted = 0.0;

    QMap<Id, FeatureFinderHillPlus> featureFinderHillsPlussesMap;
    RTree featureFinderHillRTree = buildHillsRTree(
            featureFinderHills,
            &featureFinderHillsPlussesMap
            );

    const double mzMin = featureFinderHillRTree.bounds().min_corner().get<1>();
    const double mzMax = featureFinderHillRTree.bounds().max_corner().get<1>();

    QVector<FeatureFinderHillPlus> anchorHillsAtBestMz;
    e = findBestMS2IonAnchor(
            ms2IonsAnchors,
            featureFinderHillsPlussesMap,
            featureFinderHillRTree,
            m_params.tolerancePPM,
            &anchorHillsAtBestMz
            ); ree

    if (anchorHillsAtBestMz.isEmpty()) {
        ERR_RETURN
    }

    for (const FeatureFinderHillPlus &bestMS2IonAnchor : anchorHillsAtBestMz) {

        const double maxIntensityScanNumberIndex
                = bestMS2IonAnchor.featureFinderHill.maxIntensityScanNumberIndex();

        std::vector<rTreePoint> rTreeSearchResult = searchRTree(
                maxIntensityScanNumberIndex,
                maxIntensityScanNumberIndex,
                mzMin,
                mzMax,
                featureFinderHillRTree,
                featureFinderHillsPlussesMap
                );

        if (rTreeSearchResult.size() == 1) {
            continue;
        }

        QVector<FeatureFinderHillPlus> correlatedHills;
        e = removePointsBelowCosineSimCorrThresholdMS2(
                bestMS2IonAnchor,
                featureFinderHillsPlussesMap,
                rTreeSearchResult,
                m_params.cosineSimThreshold,
                &correlatedHills
        ); ree;

        e = removeDuplicateCorrelatedHills(
                ms2IonsTandemPred,
                &correlatedHills
                ); ree

        const double cosineSimSumWeighted = calculateWeightedHillCosineSimSum(correlatedHills);

        if (cosineSimSumWeighted > bestHillsClusteringMS2->cosineSimSumWeighted) {

            FeatureFinderHillUtils::sortFeatureFinderHillsPlussesMzAsc(&correlatedHills);

            bestHillsClusteringMS2->apexFeatureFinderHillPlus = bestMS2IonAnchor;
            bestHillsClusteringMS2->cosineSimSumWeighted = cosineSimSumWeighted;
            bestHillsClusteringMS2->correlatedHills = correlatedHills;
        }

    }

    ERR_RETURN
}

Err FeatureFinderHillClusterTron::calculateCosineSimBetweenHills(
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
