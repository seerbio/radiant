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

#include <iostream>

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

    Err buildMs1Frame(
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
                "-1000",
                ms1Scans,
                {-1.0, -1.0}
        ); ree;

        *scanNumberVsScanTime = msReader.getScanNumberVsScanTime();

        e = ErrorUtils::isNotEmpty(*scanNumberVsScanTime); ree;

        ERR_RETURN
    }

    void sortHillsIntensityHiLo(QVector<FeatureFinderHill> *hills) {

        const auto sortLogic = [](const FeatureFinderHill &l, const FeatureFinderHill &r){
            return l.intensityValueMax() < r.intensityValueMax();
        };

        std::sort(hills->rbegin(), hills->rend(), sortLogic);
    }

    RTree buildHillMaximaRtree(const QVector<FeatureFinderHill> &featureFinderHillsSortedHiLo) {

        QElapsedTimer et;
        et.start();

        std::vector<rTreePoint> cloudLoader;

        const auto loadLogic = [&](const FeatureFinderHill &ffh) {
            rTreeCoor coor(static_cast<double>(ffh.maxIntensityScanNumber()), ffh.mzMean());
            return std::make_pair(coor, ffh.intensityValueMax());
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
            RTree *rTree,
            FeatureFinderMS1Feature *ms1Feature
            ) {

        ERR_INIT

        const std::vector<rTreePoint> treePoints = getTreePoints(
                hill.maxIntensityScanNumber(),
                ffParams.scanBuffer,
                hill.mzMean(),
                ffParams.mzBuffer,
                rTree
        );

        if (treePoints.size() < 2) {
            ERR_RETURN
        }

        const ScanPoints scanPoints = convertTreePointsToScanPoints(treePoints);
        const QPointF mzCenterPoint = {hill.mzMean(), hill.intensityValueMax()};

        e = MsUtils::chargeDeterminator(
                mzCenterPoint,
                scanPoints,
                ffParams.tolerancePPM,
                1,
                5,
                &ms1Feature->charge
        ); ree;

        if (ms1Feature->charge < 1) {
            ERR_RETURN
        }

        e = MsUtils::monoIsotopeDeterminator(
                mzCenterPoint,
                scanPoints,
                ffParams.tolerancePPM,
                ms1Feature->charge,
                &ms1Feature->monoIsotopeOffset,
                &ms1Feature->foundIsotopes,
                &ms1Feature->cosineSim
                ); ree;

        ms1Feature->mzFound = hill.mzMean();
        ms1Feature->scanNumberMinMax = hill.scanNumberMinMax();

        removePointsFromRTree(
                ms1Feature->foundIsotopes,
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
            const QVector<FeatureFinderHill> &featureFinderHills,
            QVector<FeatureFinderMS1Feature> *ms1Features
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(featureFinderHills); ree;

        RTree rTree = buildHillMaximaRtree(featureFinderHills);

        const double sameNessTol = 0.001;
        for (const FeatureFinderHill &h : featureFinderHills) {

//            if(h.maxIntensityScanNumber() != 6734
//                /*|| !MathUtils::tSame(h.mzMean(), 605.997)*/) {
//                continue;
//            }

            rTreeCoor queryPoint(static_cast<double>(h.maxIntensityScanNumber()), h.mzMean());
            auto firstFoundQueryPoint = bgi::qbegin(rTree, bgi::nearest(queryPoint, 1));

            if(!MathUtils::tSame(firstFoundQueryPoint->first.get<1>(), h.mzMean(),sameNessTol)
               || !MathUtils::tSame(static_cast<int>(firstFoundQueryPoint->first.get<0>()), h.maxIntensityScanNumber())
            ) {
                continue;
            }
//            std::cout << "SDFSDFSD " << firstFoundQueryPoint->first.get<1>() << " " << h.mzMean()
//                      << firstFoundQueryPoint->first.get<0>() << " " <<  h.maxIntensityScanNumber() << std::endl;

            FeatureFinderMS1Feature ms1Feature;
            e = buildHillCluster(
                    h,
                    ffParams,
                    &rTree,
                    &ms1Feature
                    ); ree;
//            break;

            if (ms1Feature.charge < 1) {
                continue;
            }

//            qDebug() << ms1Feature.mzFound << ms1Feature.charge << ms1Feature.cosineSim;

            ms1Features->push_back(ms1Feature);
        }

        ERR_RETURN
    }

    Err loadFeatureFinderParams(
            const PythiaParameters &params,
            FeatureFinderParameters *featureFinderParameters
            ) {

        ERR_INIT

        featureFinderParameters->tolerancePPM = params.featureFinderTolerancePPM;
        featureFinderParameters->skipScanCount = 0; //set to zero because MS1 scans are smoothed.
        featureFinderParameters->minScanCount = params.minScanCount;
        featureFinderParameters->useMeanMz = params.useMeanMz;
        featureFinderParameters->filterLength = params.filterLength;
        featureFinderParameters->smoothCount = params.smoothCount;
        featureFinderParameters->sigma = params.sigma;
        featureFinderParameters->signalToNoiseRatio = params.signalToNoiseRatio;

        e = ErrorUtils::isTrue(featureFinderParameters->isValid()); ree;

        ERR_RETURN
    }

}//namespace
Err Ms1FeatureFinder::exec(
        const QString &msDataFilePath,
        QVector<FeatureFinderHill> *featureFinderHills
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;

    MsFrame msFrame;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    e = buildMs1Frame(
            msDataFilePath,
            m_params,
            &msFrame,
            &scanNumberVsScanTime
    ); ree;

    FeatureFinderParameters ffParams;
    e = loadFeatureFinderParams(
            m_params,
            &ffParams
            ); ree;

    FeatureFinderHillBuilder hillBuilder;
    e = hillBuilder.init(ffParams); ree;
    hillBuilder.setRunParallel(true);

    e = hillBuilder.buildHills(
            msFrame.scanNumberVsScanPoints(),
            featureFinderHills
            ); ree;

    e = hillBuilder.refineHills(featureFinderHills); ree;

    qDebug() << "Hill count" << featureFinderHills->size();

//    sortHillsIntensityHiLo(&featureFinderHills);
//
////#define WRITE_HILLS
//#ifdef WRITE_HILLS
//    const QString batmassHillsFilePath = msDataFilePath + ".newRefine" + ".mzrt.csv";
//    e = FeatureFinderHillBuilder::writeHillsToBatmassMzMrtFile(
//            scanNumberVsScanTime,
//            featureFinderHills,
//            batmassHillsFilePath
//            ); ree;
//#endif
//
//    QVector<FeatureFinderMS1Feature> ms1Features;
//    e = clusterHills(
//            ffParams,
//            featureFinderHills,
//            &ms1Features
//            ); ree;
//
//    qDebug() << ms1Features.size();

    ERR_RETURN
}
