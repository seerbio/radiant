//
// Created by anichols on 4/29/23.
//

#include "Ms1FeatureFinder.h"

#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MsFrame.h"
#include "MsReaderParquet.h"
#include "MsUtils.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <QElapsedTimer>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
using rTreeSearchBox = bg::model::box<rTreeCoor>;
using rTreePoint = std::pair<rTreeCoor, double> ;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

Err Ms1FeatureFinder::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    ERR_RETURN
}

namespace {

    Err buildMsFrame(
            const QString &msDataFilePath,
            const PythiaParameters &params,
            MsFrame *msFrame,
            QMap<ScanNumber, ScanTime> *scanNumberVsScanTime
            ) {

        ERR_INIT

        MsReaderParquet msReader;
        e = msReader.openFile(
                msDataFilePath,
                "msLevel",
                {1, 1}
        ); ree;

        const int msLevel = 1;
        QMap<ScanNumber, ScanPoints> ms1Scans;
        e = msReader.getScanPoints(
                msLevel,
                &ms1Scans
        ); ree;

        e = msFrame->init(
                params,
                "-1000",
                ms1Scans,
                {-1.0, -1.0}
        ); ree;


        const bool denoise = true;
        const bool deisotope = false;
        const bool smooth = true;
        e = msFrame->preprocessMsFrame(
                denoise,
                deisotope,
                smooth
        ); ree;

        *scanNumberVsScanTime = msReader.getScanNumberVsScanTime();

        e = ErrorUtils::isNotEmpty(*scanNumberVsScanTime); ree;

        ERR_RETURN
    }

    void sortHillsIntensityHiLo(QVector<FeatureFinderHill> *hills) {

        const auto sortLogic = [](const FeatureFinderHill &l, const FeatureFinderHill &r){
            return l.maxIntensityValue() < r.maxIntensityValue();
        };

        std::sort(hills->rbegin(), hills->rend(), sortLogic);
    }

    RTree buildHillMaximaRtree(const QVector<FeatureFinderHill> &featureFinderHillsSortedHiLo) {

        QElapsedTimer et;
        et.start();

        std::vector<rTreePoint> cloudLoader;

        const auto loadLogic = [&](const FeatureFinderHill &ffh) {
            rTreeCoor coor(static_cast<double>(ffh.maxIntensityScanNumber()), ffh.mzMean());
            return std::make_pair(coor, ffh.maxIntensityValue());
        };

        std::transform(
                featureFinderHillsSortedHiLo.begin(),
                featureFinderHillsSortedHiLo.end(),
                back_inserter(cloudLoader),
                loadLogic
        );

        const int maxElements = 16;
        RTree rTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

        qDebug() << featureFinderHillsSortedHiLo.size()
                 << "Hills loaded into rtree in"
                 << et.restart()
                 << "mSec";

        return rTree;
    }

    std::vector<rTreePoint> getTreePoints(
            ScanNumber scanNumber,
            int scanBuffer,
            double mz,
            double mzBuffer,
            RTree *rTree
            ) {

        const double scanNumberMin = scanNumber - scanBuffer;
        const double scanNumberMax = scanNumber + scanBuffer;

        const double mzMin = mz - mzBuffer;
        const double mzMax = mz + mzBuffer;

        const rTreeSearchBox queryBox(
                rTreeCoor(scanNumberMin, mzMin),
                rTreeCoor(scanNumberMax, mzMax)
        );

        std::vector<rTreePoint> result;
        rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

        return result;
    }

    ScanPoints convertTreePointsToScanPoints(const std::vector<rTreePoint> &rTreePoints) {

        ScanPoints scanPoints;
        for (const rTreePoint &rtp : rTreePoints) {
            scanPoints.push_back({rtp.first.get<1>(), rtp.second});
        }

        return scanPoints;
    }

    void removePointsFromRTree(
            const QVector<QPointF> &subtractionPoints,
            ScanNumber scanNumber,
            RTree *rTree
            ) {

        for (const QPointF &p : subtractionPoints) {

            rTreeCoor deletePoint(static_cast<double>(scanNumber), p.x());

            auto first = bgi::qbegin(*rTree, bgi::nearest(deletePoint, 1));
            auto last  = bgi::qend(*rTree);

            if (first != last) {
                rTree->remove(*first); // assuming single result
            }

        }

    }

    Err buildHillCluster(
            const FeatureFinderHill &hill,
            const FeatureFinderParameters &ffParams,
            RTree *rTree
            ) {

        ERR_INIT

        const std::vector<rTreePoint> treePoints = getTreePoints(
                hill.maxIntensityScanNumber(),
                ffParams.scanBuffer,
                hill.mzMean(),
                ffParams.mzBuffer,
                rTree
        );

        const ScanPoints scanPoints = convertTreePointsToScanPoints(treePoints);
        const QPointF mzCenterPoint = {hill.mzMean(), hill.maxIntensityValue()};

        int charge;
        e = MsUtils::chargeDeterminator(
                mzCenterPoint,
                scanPoints,
                ffParams.tolerancePPM,
                &charge
        ); ree;

        if (charge < 1) {
//            qDebug() << charge;
            ERR_RETURN
        }

        int monoIsoOffset;
        QVector<QPointF> subtractionPoints;
        e = MsUtils::monoIsotopeDeterminator(
                mzCenterPoint,
                scanPoints,
                ffParams.tolerancePPM,
                charge,
                &monoIsoOffset,
                &subtractionPoints
                ); ree;

        removePointsFromRTree(
                subtractionPoints,
                hill.maxIntensityScanNumber(),
                rTree
                );

//        qDebug() << charge * mzCenterPoint.x()
//                 << hill.maxIntensityScanNumber() << mzCenterPoint
//                 << charge << monoIsoOffset;

        ERR_RETURN
    }

    Err clusterHills(
            const FeatureFinderParameters &ffParams,
            const QVector<FeatureFinderHill> &featureFinderHills
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(featureFinderHills); ree;

        RTree rTree = buildHillMaximaRtree(featureFinderHills);

        for (const FeatureFinderHill &h : featureFinderHills) {

//            if(h.maxIntensityScanNumber() !=  6742 /*|| !MathUtils::tSame(h.mzMean(), 605.997)*/) {
//                continue;
//            }

            rTreeCoor queryPoint(static_cast<double>(h.maxIntensityScanNumber()), h.mzMean());
            auto firstFoundQueryPoint = bgi::qbegin(rTree, bgi::nearest(queryPoint, 1));

            if(!MathUtils::tSame(firstFoundQueryPoint->first.get<1>(), h.mzMean())) {
                continue;
            }

            e = buildHillCluster(
                    h,
                    ffParams,
                    &rTree
                    ); ree;
//            break;

        }

        ERR_RETURN
    }


}//namespace
Err Ms1FeatureFinder::exec(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;

    MsFrame msFrame;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    e = buildMsFrame(
            msDataFilePath,
            m_params,
            &msFrame,
            &scanNumberVsScanTime
    ); ree;

    //TODO set these accordingly in PythiaParams file.
    FeatureFinderParameters ffParams;
    ffParams.tolerancePPM = m_params.featureFinderTolerancePPM;
    ffParams.skipScanCount = 0;
    ffParams.minScanCount = 2;
    ffParams.useMeanMz = true;

    FeatureFinderHillBuilder hillBuilder;
    e = hillBuilder.init(ffParams); ree;
    hillBuilder.setRunParallel(true);

    QVector<FeatureFinderHill> featureFinderHills;
    e = hillBuilder.buildHills(
            msFrame.scanNumberVsScanPoints(),
            &featureFinderHills
            ); ree;

    e = hillBuilder.refineHills(&featureFinderHills); ree;
    sortHillsIntensityHiLo(&featureFinderHills);

//#define WRITE_HILLS
#ifdef WRITE_HILLS
    const QString batmassHillsFilePath = msDataFilePath + ".newRefine" + ".mzrt.csv";
    e = FeatureFinderHillBuilder::writeHillsToBatmassMzMrtFile(
            scanNumberVsScanTime,
            featureFinderHills,
            batmassHillsFilePath
            ); ree;
#endif

    e = clusterHills(
            ffParams,
            featureFinderHills
            ); ree;

    ERR_RETURN
}
