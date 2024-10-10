//
// Created by andrewnichols on 10/8/24.
//

#include "MsReaderBrukerTims.h"

#include "ParallelUtils.h"

#include "CppSQLite3.h"
#include "timsdata_cpp.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <boost/math/constants/constants.hpp>

#include "EigenUtils.h"


#define OUTPUT_PRECISION 8
#define OUTPUT_WIDTH 9

MsReaderBrukerTims::MsReaderBrukerTims() {
}

MsReaderBrukerTims::~MsReaderBrukerTims() {
}

namespace {

    std::vector<TimsMS2WindowsInfo> buildTimsWindowsInfos(CppSQLite3DB *db) {

        CppSQLite3Query query = db->execQuery(
            "SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth, CollisionEnergy "
            "FROM DiaFrameMsMsWindows;"
            );

        std::vector<TimsMS2WindowsInfo> diaFrameInfos;

        while(!query.eof()) {

            TimsMS2WindowsInfo dfi;

            dfi.windowGroup = query.getIntField("WindowGroup");
            dfi.scanNumberBegin = query.getIntField("ScanNumBegin");
            dfi.scanNumberEnd = query.getIntField("ScanNumEnd");
            dfi.isolationMz = query.getFloatField("IsolationMz");
            dfi.isolationWidth = query.getFloatField("IsolationWidth");
            dfi.collisionEnergy = query.getFloatField("CollisionEnergy");

            diaFrameInfos.push_back(dfi);

            query.nextRow();
        }

        return diaFrameInfos;
    }

    std::vector<TimsFrameInfo> buildTimsFrameInfo(CppSQLite3DB *db) {

        CppSQLite3Query query = db->execQuery("SELECT Id, Time, MsMsType, NumScans FROM Frames;");

        std::vector<TimsFrameInfo> brukerFrameInfos;

        while(!query.eof()) {

            TimsFrameInfo frameInfo;
            frameInfo.frameId = query.getIntField("Id");
            frameInfo.scanTime = query.getFloatField("Time") / 60.0f;
            frameInfo.msmsType = query.getIntField("MsMsType");
            frameInfo.numScans = query.getIntField("NumScans");

            brukerFrameInfos.push_back(frameInfo);

            query.nextRow();
        }

        return brukerFrameInfos;
    }

    Err buildWindowGroupIndexes(
        const std::vector<TimsMS2WindowsInfo> &timeWindowsInfos,
        std::vector<int> *windowGroupIndexes
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(timeWindowsInfos); ree;

        const auto windowGroupMinMax = std::minmax_element(
            timeWindowsInfos.begin(),
            timeWindowsInfos.end(),
            [](const TimsMS2WindowsInfo &l, const TimsMS2WindowsInfo &r){return l.windowGroup < r.windowGroup;}
            );

        const int windowGroupCount = windowGroupMinMax.second->windowGroup - windowGroupMinMax.first->windowGroup + 1;
        std::vector<int> windowGroups(windowGroupCount);
        std::iota(windowGroups.begin(), windowGroups.end(), windowGroupMinMax.first->windowGroup);

        *windowGroupIndexes = windowGroups;
        e = ErrorUtils::isFalse(windowGroupIndexes->empty()); ree;

        ERR_RETURN
    }

    Err assignWindowGroupToTimsFrameInfos(
        timsdata::TimsData *data,
        std::vector<TimsMS2WindowsInfo> *timsMS2WindowsInfos,
        std::vector<TimsFrameInfo> *timsFramesInfos
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(timsMS2WindowsInfos->empty()); ree;
        e = ErrorUtils::isFalse(timsFramesInfos->empty()); ree;

        std::vector<int> windowGroupIndexes;
        e = buildWindowGroupIndexes(
            *timsMS2WindowsInfos,
            &windowGroupIndexes
            ); ree

        const size_t windowGroupCount = windowGroupIndexes.size();

        int ms2Counter = 0;
        for (int i = 0; i < timsFramesInfos->size(); i++) {

            TimsFrameInfo &tfi = (*timsFramesInfos)[i];

            if (tfi.msmsType < 1) {
                continue;
            }

            tfi.windowGroup = windowGroupIndexes.at(ms2Counter++ % windowGroupCount);
        }

        ERR_RETURN

    }

    Err weightedMean(const ScanPoints &scanPoints, float *weightedMean) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scanPoints); ree;

        const float numerator = std::accumulate(
            scanPoints.begin(),
            scanPoints.end(),
            0.0f,
            [](float sum, const ScanPoint &sp){return sum + (sp.x() * sp.y());}
            );

        const float denom = std::accumulate(
            scanPoints.begin(),
            scanPoints.end(),
            0.0f,
            [](float sum, const ScanPoint &sp){return sum + sp.y();}
            );

        e = ErrorUtils::isFalse(MathUtils::tZero(denom)); ree;

        *weightedMean = numerator / denom;

        ERR_RETURN
    }

    Err scanClustererer(
        const ScanPoints &scanPoints,
        ScanPoints *scanPointsClustered
        ) {

        ERR_INIT

        QVector<float> mzVals;
        std::transform(
            scanPoints.begin(),
            scanPoints.end(),
            std::back_inserter(mzVals),
            [](const ScanPoint &scanPoint){return scanPoint.x();}
            );

        const bool isSorted = std::is_sorted(mzVals.begin(), mzVals.end());
        e = ErrorUtils::isTrue(isSorted); ree;

        const Eigen::VectorX<float> v = EigenUtils::convertQVectorToEigenVector(mzVals);

        const Eigen::VectorX<float> v1 = v.head(v.size() - 1);
        const Eigen::VectorX<float> v2 = v.tail(v.size() - 1);

        const Eigen::VectorX<float> ppmDiffs = 1e6 * ((v2 - v1).array() / v1.array());

        QVector<ScanPoint> scanPointsCurrentCluster;
        for (int i = 0; i < scanPoints.size(); i++) {

            const ScanPoint &scanPoint = scanPoints.at(i);
            if (scanPointsCurrentCluster.isEmpty()) {
                scanPointsCurrentCluster.push_back(scanPoint);
                continue;
            }

            const float ppmDiff = ppmDiffs.coeffRef(i - 1);
            if (constexpr float ppmTol = 18.0; ppmDiff <= ppmTol) {
                scanPointsCurrentCluster.push_back(scanPoint);
                continue;
            }

            QVector<float> intensityValsCurrentCluster;
            std::transform(
                scanPointsCurrentCluster.begin(),
                scanPointsCurrentCluster.end(), std::back_inserter(intensityValsCurrentCluster),
                [](const ScanPoint &sp){return sp.y();}
                );

            float mzValsCurrentClusterMean;
            e = weightedMean(scanPointsCurrentCluster, &mzValsCurrentClusterMean); ree;

            const float intensitySum = std::accumulate(
                intensityValsCurrentCluster.begin(),
                intensityValsCurrentCluster.end(),
                0.0f
                );
            scanPointsClustered->push_back({mzValsCurrentClusterMean, intensitySum});

            scanPointsCurrentCluster.clear();
            scanPointsCurrentCluster.push_back(scanPoint);
        }

        if (!scanPointsCurrentCluster.isEmpty()) {

            float mzValsCurrentClusterMean;
            e = weightedMean(scanPointsCurrentCluster, &mzValsCurrentClusterMean); ree;

            QVector<float> intensityValsCurrentCluster;
            std::transform(
                scanPointsCurrentCluster.begin(),
                scanPointsCurrentCluster.end(), std::back_inserter(intensityValsCurrentCluster),
                [](const ScanPoint &sp){return sp.y();}
                );

            const float intensitySum = std::accumulate(
                intensityValsCurrentCluster.begin(),
                intensityValsCurrentCluster.end(),
                0.0f
                );
            scanPointsClustered->push_back({mzValsCurrentClusterMean, intensitySum});
        }

        ERR_RETURN
    }

    Err extractMS1ScanPointsFromFrame(
        const TimsFrameInfo &timsFrameInfo,
        timsdata::TimsData *data,
        ScanPoints *scanPointsMS1
        ) {

        ERR_INIT

        const timsdata::FrameProxy scans = data->readScans(timsFrameInfo.frameId, 0, timsFrameInfo.numScans);

        ScanPoints scanPointsFrame;
        scanPointsFrame.reserve(static_cast<int>(scans.getTotalNbrPeaks()));

        constexpr int nonZeroIndexArrayOffset = -1;
        for(int scanNumber = 0 + nonZeroIndexArrayOffset ; scanNumber < scans.getNbrScans(); scanNumber++) {

            if (scans.getNbrPeaks(scanNumber) < 1) {
                continue;
            }

            timsdata::FrameProxy::FrameIteratorRange xAxis = scans.getScanX(scanNumber);
            const timsdata::FrameProxy::FrameIteratorRange yAxis = scans.getScanY(scanNumber);

            std::vector<double> xAxisMasses;
            std::vector<double> indices(xAxis.first, xAxis.second);
            data->indexToMz(timsFrameInfo.frameId, indices, xAxisMasses);

            const size_t numberOfPeaks = scans.getNbrPeaks(scanNumber);

            ScanPoints scanPointsScan(static_cast<int>(numberOfPeaks));
            for(size_t pkNum = 0; pkNum < numberOfPeaks; ++pkNum) {
                scanPointsScan[static_cast<int>(pkNum)] = {static_cast<float>(xAxisMasses[pkNum]), static_cast<float>(yAxis.first[pkNum])};
            }

            scanPointsFrame.append(scanPointsScan);
        }

        std::sort(scanPointsFrame.begin(), scanPointsFrame.end(), [](const ScanPoint &l, const ScanPoint &r){return l.x() < r.x();});

        e = scanClustererer(scanPointsFrame, scanPointsMS1); ree;

        ERR_RETURN
    }

    Err extractMS2ScanPointsFromFrame(
        const TimsFrameInfo &timsFrameInfo,
        const QHash<int, QVector<TimsMS2WindowsInfo>> &windowGroupIndexVsTimsMs2WindowsInfoses,
        timsdata::TimsData *data,
        QMap<MzTargetKey, ScanPoints> *scanPointFrameClustered
        ) {

        ERR_INIT

        if (timsFrameInfo.msmsType < 1) {
            ERR_RETURN;
        }

        e = ErrorUtils::contains(timsFrameInfo.windowGroup, windowGroupIndexVsTimsMs2WindowsInfoses); ree;
        const QVector<TimsMS2WindowsInfo> &timsMs2WindowsInfos
                        = windowGroupIndexVsTimsMs2WindowsInfoses.value(timsFrameInfo.windowGroup);

        const timsdata::FrameProxy scans = data->readScans(timsFrameInfo.frameId, 0, timsFrameInfo.numScans);
        const size_t totalPeaks = scans.getTotalNbrPeaks();

        for (const TimsMS2WindowsInfo &tmwi : timsMs2WindowsInfos) {

            const MzTargetKey mzTargetKey = QString::number(MathUtils::hashDecimal(
                tmwi.isolationMz,
                S_GLOBAL_SETTINGS.HASHING_PRECISION)
                );

            ScanPoints scanPointsFrame;
            scanPointsFrame.reserve(static_cast<int>(totalPeaks));

            constexpr int nonZeroIndexArrayOffset = -1;
            for(int scanNumber = tmwi.scanNumberBegin + nonZeroIndexArrayOffset ; scanNumber < tmwi.scanNumberEnd; scanNumber++) {

                if (scans.getNbrPeaks(scanNumber) < 1) {
                    continue;
                }

                timsdata::FrameProxy::FrameIteratorRange xAxis = scans.getScanX(scanNumber);
                const timsdata::FrameProxy::FrameIteratorRange yAxis = scans.getScanY(scanNumber);

                std::vector<double> xAxisMasses;
                std::vector<double> indices(xAxis.first, xAxis.second);
                data->indexToMz(timsFrameInfo.frameId, indices, xAxisMasses);

                const size_t numberOfPeaks = scans.getNbrPeaks(scanNumber);

                ScanPoints scanPointsScan(static_cast<int>(numberOfPeaks));
                for(size_t pkNum = 0; pkNum < numberOfPeaks; ++pkNum) {
                    scanPointsScan[static_cast<int>(pkNum)] = {static_cast<float>(xAxisMasses[pkNum]), static_cast<float>(yAxis.first[pkNum])};
                }

                scanPointsFrame.append(scanPointsScan);
            }

            std::sort(scanPointsFrame.begin(), scanPointsFrame.end(), [](const ScanPoint &l, const ScanPoint &r){return l.x() < r.x();});

            ScanPoints scanPointsClustered;
            e = scanClustererer(scanPointsFrame, &scanPointsClustered); ree;
            scanPointFrameClustered->insert(mzTargetKey, scanPointsClustered);
        }

        ERR_RETURN
    }

    Err parallelReadTimsLogic(
        const std::string &tdfDirectory,
        const std::vector<TimsFrameInfo> &timsFramesInfos,
        const QHash<int, QVector<TimsMS2WindowsInfo>> &windowGroupIndexVsTimsMs2WindowsInfoses
        ) {

        ERR_INIT

        timsdata::TimsData data(tdfDirectory);

        for (const TimsFrameInfo &t : timsFramesInfos) {

            std::cout << t.frameId << " " << t.windowGroup << std::endl;

            if (t.msmsType < 1) {

                ScanPoints scanPointsMS1;
                e = extractMS1ScanPointsFromFrame(
                    t,
                    &data,
                    &scanPointsMS1
                    ); ree;
                continue;
            }

            QMap<MzTargetKey, ScanPoints> scanPointFrameClustered;
            e = extractMS2ScanPointsFromFrame(
                t,
                windowGroupIndexVsTimsMs2WindowsInfoses,
                &data,
                &scanPointFrameClustered
                ); ree;
        }

        ERR_RETURN
    }

}//namespace
Err MsReaderBrukerTims::openFile(const QString &filePath) {

    ERR_INIT

    std::string tdfDirectory(filePath.toStdString());
    std::string tdfFile = tdfDirectory + "/analysis.tdf";

    e = ErrorUtils::fileExists(QString::fromStdString(tdfFile)); ree;

    try {
        timsdata::TimsData data(tdfDirectory);  // NOTE: UTF-8 conversion needed here!

        CppSQLite3DB db;
        db.open(tdfFile.c_str());  // NOTE: UTF-8 conversion needed here!

        std::vector<TimsMS2WindowsInfo> timsMS2WindowsInfos = buildTimsWindowsInfos(&db);
        e = ErrorUtils::isNotEmpty(timsMS2WindowsInfos); ree;

        QHash<int, QVector<TimsMS2WindowsInfo>> windowGroupIndexVsTimsMs2WindowsInfoses;
        for (const TimsMS2WindowsInfo &w : timsMS2WindowsInfos) {
            windowGroupIndexVsTimsMs2WindowsInfoses[w.windowGroup].push_back(w);
        }

        std::vector<TimsFrameInfo> timsFramesInfos = buildTimsFrameInfo(&db);
        e = ErrorUtils::isNotEmpty(timsFramesInfos); ree;
        db.close();

        e = assignWindowGroupToTimsFrameInfos(
            &data,
            &timsMS2WindowsInfos,
            &timsFramesInfos
            ); ree;


#define TIMS_PARALLEL
#ifdef TIMS_PARALLEL

        const auto readerLogicBinder = std::bind(
            parallelReadTimsLogic,
            tdfDirectory,
            std::placeholders::_1,
            windowGroupIndexVsTimsMs2WindowsInfoses
            );

        QVector<std::vector<TimsFrameInfo>> timsFramesInfosTranched;
        ParallelUtils::trancheVectorForParallelizationInOrder(
            timsFramesInfos,
            ParallelUtils::numberOfAvailableSystemProcessors(),
            0,
            &timsFramesInfosTranched
            );

        QFuture<Err> future = QtConcurrent::mapped(
            timsFramesInfosTranched,
            readerLogicBinder
            );
        future.waitForFinished();
#else
        e = parallelReadTimsLogic(
            tdfDirectory,
            timsFramesInfos,
            windowGroupIndexVsTimsMs2WindowsInfoses
            ); ree;
#endif

    }
    catch(const std::exception& ee) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << ee.what();
        rrr(eFileError);
    }

    ERR_RETURN
}

Err MsReaderBrukerTims::closeFile() {
    ERR_INIT

    ERR_RETURN
}
