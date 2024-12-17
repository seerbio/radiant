//
// Created by andrewnichols on 10/8/24.
//

#include "MsReaderBrukerTims.h"

#include "ParallelUtils.h"
#include "MsReaderParquet.h"

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
        ScanPoints *scanPointsMS1,
        Ms1FrameTIMS *frameTims
        ) {

        ERR_INIT

        const timsdata::SpectrumData scansCollapsed = data->readScansSummed(timsFrameInfo.frameId, 0, timsFrameInfo.numScans);

        const int numberOfPeaks = static_cast<int>(scansCollapsed.area_values.size());
        scanPointsMS1->reserve(numberOfPeaks);

        for (int i = 0; i < numberOfPeaks; i++) {
            scanPointsMS1->push_back({static_cast<float>(scansCollapsed.mz_values.at(i)), scansCollapsed.area_values.at(i)});
        }

        const timsdata::FrameProxy scans = data->readScans(timsFrameInfo.frameId, 0, timsFrameInfo.numScans);



        for(IonMobilityIndex ionMobilityIndex = 0; ionMobilityIndex < timsFrameInfo.numScans; ionMobilityIndex++) {

            if (scans.getNbrPeaks(ionMobilityIndex) < 1) {
                continue;
            }

            timsdata::FrameProxy::FrameIteratorRange xAxis = scans.getScanX(ionMobilityIndex);
            const timsdata::FrameProxy::FrameIteratorRange yAxis = scans.getScanY(ionMobilityIndex);

            std::vector<double> xAxisMasses;
            std::vector<double> indices(xAxis.first, xAxis.second);
            data->indexToMz(timsFrameInfo.frameId, indices, xAxisMasses);

            const size_t numberOfPeaks = scans.getNbrPeaks(ionMobilityIndex);

            ScanPoints scanPointsScan(static_cast<int>(numberOfPeaks));
            for(size_t pkNum = 0; pkNum < numberOfPeaks; ++pkNum) {
                scanPointsScan[static_cast<int>(pkNum)] = {static_cast<float>(xAxisMasses[pkNum]), static_cast<float>(yAxis.first[pkNum])};
            }

            frameTims->insert(ionMobilityIndex, scanPointsScan);
        }

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

        for (const TimsMS2WindowsInfo &tmwi : timsMs2WindowsInfos) {

            const MzTargetKey mzTargetKey = QString::number(MathUtils::hashDecimal(
                tmwi.isolationMz,
                S_GLOBAL_SETTINGS.HASHING_PRECISION)
                );

            constexpr int buffer = 1;
            const timsdata::SpectrumData spectrumData = data->readScansSummed(
                timsFrameInfo.frameId,
                tmwi.scanNumberBegin + buffer,
                tmwi.scanNumberEnd - buffer
                );

            int scanPointsSize = static_cast<int>(spectrumData.mz_values.size());

            ScanPoints scanPointsFrame;
            scanPointsFrame.reserve(scanPointsSize);

            for (int i = 0; i < scanPointsSize; i++) {
                scanPointsFrame.push_back({static_cast<float>(spectrumData.mz_values.at(i)), spectrumData.area_values.at(i)});
            }

            std::sort(scanPointsFrame.begin(), scanPointsFrame.end(), [](const ScanPoint &l, const ScanPoint &r){return l.x() < r.x();});

            scanPointFrameClustered->insert(mzTargetKey, scanPointsFrame);
        }

        ERR_RETURN
    }

    Err buildMzTargetVsTimsMs2WindowsInfos(
        const QHash<int, QVector<TimsMS2WindowsInfo>> &windowGroupIndexVsTimsMs2WindowsInfoses,
        QHash<MzTargetKey, TimsMS2WindowsInfo> *mzTargetVsTimsMs2WindowsInfos
        ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(windowGroupIndexVsTimsMs2WindowsInfoses); ree;

        for (const QVector<TimsMS2WindowsInfo> &infosVector : windowGroupIndexVsTimsMs2WindowsInfoses) {
            for(const TimsMS2WindowsInfo &tmwi : infosVector) {
                const MzTargetKey mzTargetKey
                    = QString::number(MathUtils::hashDecimal(tmwi.isolationMz, S_GLOBAL_SETTINGS.HASHING_PRECISION));
                mzTargetVsTimsMs2WindowsInfos->insert(mzTargetKey, tmwi);
            }
        }

        ERR_RETURN
    }

    QPair<Err, QVector<std::tuple<MsScanInfo, ScanPoints, Ms1FrameTIMS>>> parallelReadTimsLogic(
        const std::string &tdfDirectory,
        const std::vector<TimsFrameInfo> &timsFramesInfos,
        const QMap<FrameIndex, double> &frameIndexVsDriftTime,
        const QHash<MzTargetKey, TimsMS2WindowsInfo> &mzTargetVsTimsMs2WindowsInfos,
        const QHash<int, QVector<TimsMS2WindowsInfo>> &windowGroupIndexVsTimsMs2WindowsInfoses
        ) {

        ERR_INIT

        timsdata::TimsData data(tdfDirectory);

        e = ErrorUtils::isNotEmpty(timsFramesInfos); rree;
        e = ErrorUtils::isNotEmpty(frameIndexVsDriftTime); rree;
        e = ErrorUtils::isNotEmpty(mzTargetVsTimsMs2WindowsInfos); rree;
        e = ErrorUtils::isNotEmpty(windowGroupIndexVsTimsMs2WindowsInfoses); rree;

        QVector<std::tuple<MsScanInfo, ScanPoints, Ms1FrameTIMS>> msScanInfoScanPointsPairs;

        for (const TimsFrameInfo &tfi : timsFramesInfos) {

            MsScanInfo msScanInfo;
            msScanInfo.scanNumber = tfi.frameId;
            msScanInfo.scanTime = tfi.scanTime;

            if (tfi.msmsType < 1) {

                ScanPoints scanPointsMS1;
                Ms1FrameTIMS ms1FrameTims;
                e = extractMS1ScanPointsFromFrame(
                    tfi,
                    &data,
                    &scanPointsMS1,
                    &ms1FrameTims
                    ); rree;

                msScanInfo.msLevel = 1;

                msScanInfoScanPointsPairs.push_back({msScanInfo, scanPointsMS1, ms1FrameTims});
                continue;
            }

            QMap<MzTargetKey, ScanPoints> mzTargetKeyVsScanPoints;
            e = extractMS2ScanPointsFromFrame(
                tfi,
                windowGroupIndexVsTimsMs2WindowsInfoses,
                &data,
                &mzTargetKeyVsScanPoints
                ); rree;

            for (auto it = mzTargetKeyVsScanPoints.begin(); it != mzTargetKeyVsScanPoints.end(); ++it) {
                const MzTargetKey &mzTargetKey = it.key();
                const ScanPoints &scanPoints = it.value();

                e = ErrorUtils::contains(mzTargetKey, mzTargetVsTimsMs2WindowsInfos); rree;
                const TimsMS2WindowsInfo& timsMs2WindowsInfo = mzTargetVsTimsMs2WindowsInfos.value(mzTargetKey);

                msScanInfo.msLevel = 2;
                msScanInfo.collisionEnergy = timsMs2WindowsInfo.collisionEnergy;
                msScanInfo.precursorTargetMz = timsMs2WindowsInfo.isolationMz;
                msScanInfo.isoWindowLower = timsMs2WindowsInfo.isolationWidth / 2.0f;
                msScanInfo.isoWindowUpper = msScanInfo.isoWindowLower;
                msScanInfo.ionMobilityIndex
                    = static_cast<int>(std::round((timsMs2WindowsInfo.scanNumberBegin + timsMs2WindowsInfo.scanNumberEnd) / 2.0f));

                e = ErrorUtils::contains(msScanInfo.ionMobilityIndex, frameIndexVsDriftTime); rree;
                msScanInfo.ionMobilityDriftTime = static_cast<float>(frameIndexVsDriftTime.value(msScanInfo.ionMobilityIndex));

                msScanInfoScanPointsPairs.push_back({msScanInfo, scanPoints, {}});
            }
        }

        return {e, msScanInfoScanPointsPairs};
    }

    Err buildFrameIndexVsDriftTime(
        const TimsFrameInfo &timsFrameInfo,
        timsdata::TimsData *data,
        QMap<FrameIndex, double> *frameIndexVsDriftTime
        ) {

        ERR_INIT

        e = ErrorUtils::isAboveThreshold(timsFrameInfo.numScans, 0, ErrorUtilsParam::ExcludeThreshold); ree;

        const timsdata::FrameProxy scans = data->readScans(timsFrameInfo.frameId, 0, timsFrameInfo.numScans);
        for(IonMobilityIndex ionMobilityIndex = 0; ionMobilityIndex < scans.getNbrScans(); ionMobilityIndex++) {
            std::vector<double> mobility;
            data->scanNumToOneOverK0(timsFrameInfo.frameId, { static_cast<double>(ionMobilityIndex) }, mobility);
            frameIndexVsDriftTime->insert(ionMobilityIndex, mobility[0]);
        }

        ERR_RETURN
    }

}//namespace
Err MsReaderBrukerTims::openFile(const QString &filePath) {

    ERR_INIT

    std::string tdfDirectory(filePath.toStdString());
    std::string tdfFile = tdfDirectory + "/analysis.tdf";

    e = ErrorUtils::fileExists(QString::fromStdString(tdfFile)); ree;
    m_filePath = filePath;
    m_isTIMS = true;

    try {
        timsdata::TimsData data(tdfDirectory);  // NOTE: UTF-8 conversion needed here!

        CppSQLite3DB db;
        db.open(tdfFile.c_str());  // NOTE: UTF-8 conversion needed here!

        std::vector<TimsMS2WindowsInfo> timsMS2WindowsInfos = buildTimsWindowsInfos(&db);
        e = ErrorUtils::isNotEmpty(timsMS2WindowsInfos); ree;


        for (const TimsMS2WindowsInfo &w : timsMS2WindowsInfos) {
            m_windowGroupIndexVsTimsMs2WindowsInfoses[w.windowGroup].push_back(w);
        }

        e = buildMzTargetVsTimsMs2WindowsInfos(
            m_windowGroupIndexVsTimsMs2WindowsInfoses,
            &m_mzTargetVsTimsMs2WindowsInfos
            ); ree;

        m_timsFramesInfos = buildTimsFrameInfo(&db);
        e = ErrorUtils::isNotEmpty(m_timsFramesInfos); ree;
        db.close();

        e = buildFrameIndexVsDriftTime(
            m_timsFramesInfos.front(),
            &data,
            &m_frameIndexVsDriftTime
            ); ree;

        e = assignWindowGroupToTimsFrameInfos(
            &data,
            &timsMS2WindowsInfos,
            &m_timsFramesInfos
            ); ree;

#define TIMS_PARALLEL
#ifdef TIMS_PARALLEL
        const auto readerLogicBinder = std::bind(
            parallelReadTimsLogic,
            tdfDirectory,
            std::placeholders::_1,
            m_frameIndexVsDriftTime,
            m_mzTargetVsTimsMs2WindowsInfos,
            m_windowGroupIndexVsTimsMs2WindowsInfoses
            );

        QVector<std::vector<TimsFrameInfo>> timsFramesInfosTranched;
        ParallelUtils::trancheVectorForParallelizationInOrder(
            m_timsFramesInfos,
            ParallelUtils::numberOfAvailableSystemProcessors(),
            0,
            &timsFramesInfosTranched
            );

        QFuture<QPair<Err, QVector<std::tuple<MsScanInfo, ScanPoints, Ms1FrameTIMS>>>> future = QtConcurrent::mapped(
            timsFramesInfosTranched,
            readerLogicBinder
            );
        future.waitForFinished();

        int scanNumber = 0;
        for (const QPair<Err, QVector<std::tuple<MsScanInfo, ScanPoints, Ms1FrameTIMS>>> &res : future) {
            e = res.first; ree;

            const QVector<std::tuple<MsScanInfo, ScanPoints, Ms1FrameTIMS>> &resVector = res.second;
            for (int i = 0; i < resVector.size(); i++) {
                const std::tuple<MsScanInfo, ScanPoints, Ms1FrameTIMS> &tup = resVector.at(i);
                MsScanInfo msi = std::get<0>(tup);
                msi.scanNumber = ++scanNumber;
                m_msScanInfo.insert(msi.scanNumber, msi);
                m_scanPoints.insert(msi.scanNumber, std::get<1>(tup));

                if (msi.msLevel > 1) {
                    continue;
                }
                m_frameNumberVsFrameTIMS.insert(msi.scanNumber, std::get<2>(tup));
            }
        }
#else
        QPair<Err, QVector<std::tuple<MsScanInfo, ScanPoints, Ms1FrameTIMS>>> result = parallelReadTimsLogic(
            tdfDirectory,
            timsFramesInfos,
            frameIndexVsDriftTime,
            mzTargetVsTimsMs2WindowsInfos,
            windowGroupIndexVsTimsMs2WindowsInfoses
            );
        e = result.first; ree;
#endif

    }
    catch(const std::exception& ee) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << ee.what();
        rrr(eFileError);
    }

    e = printFileInfo(); ree;

    ERR_RETURN
}

Err MsReaderBrukerTims::closeFile() {
    ERR_INIT

    ERR_RETURN
}

Err MsReaderBrukerTims::writeFrame(
    const QString &filePath,
    float scanTime,
    int msLevel
    ) {

    ERR_INIT

    timsdata::TimsData data(m_filePath.toStdString());

    TimsFrameInfo timsFrameInfo;
    for (const TimsFrameInfo &tfi : m_timsFramesInfos) {

        if (msLevel > 1) {
            if (tfi.msmsType == 0) {
                continue;
            }
        }
        else {
            if (tfi.msmsType != 0) {
                continue;
            }
        }

        timsFrameInfo = tfi;
        if (tfi.scanTime > scanTime || MathUtils::tSame(scanTime, tfi.scanTime, 0.00001)) {
            break;
        }
    }

    qDebug() << "Frame index range" << m_timsFramesInfos.front().frameId << m_timsFramesInfos.back().frameId << timsFrameInfo.frameId;

    const timsdata::FrameProxy scans = data.readScans(timsFrameInfo.frameId, 0, timsFrameInfo.numScans);

    e = ErrorUtils::contains(timsFrameInfo.windowGroup, m_windowGroupIndexVsTimsMs2WindowsInfoses); ree;
    const QVector<TimsMS2WindowsInfo> &timsMs2WindowsInfos
                    = m_windowGroupIndexVsTimsMs2WindowsInfoses.value(timsFrameInfo.windowGroup);

    QVector<MsParquetReaderRow> msParquetReaderRows;
    for (const TimsMS2WindowsInfo &tmwi : timsMs2WindowsInfos) {

        const MzTargetKey mzTargetKey = QString::number(MathUtils::hashDecimal(
            tmwi.isolationMz,
            S_GLOBAL_SETTINGS.HASHING_PRECISION)
            );

        constexpr int nonZeroIndexArrayOffset = -1;
        for(int scanNumber = tmwi.scanNumberBegin + nonZeroIndexArrayOffset ; scanNumber < tmwi.scanNumberEnd; scanNumber++) {

            if (scans.getNbrPeaks(scanNumber) < 1) {
                continue;
            }

            timsdata::FrameProxy::FrameIteratorRange xAxis = scans.getScanX(scanNumber);
            const timsdata::FrameProxy::FrameIteratorRange yAxis = scans.getScanY(scanNumber);

            std::vector<double> xAxisMasses;
            std::vector<double> indices(xAxis.first, xAxis.second);
            data.indexToMz(timsFrameInfo.frameId, indices, xAxisMasses);

            QVector<float> mzVals;
            QVector<float> intensityVals;
            const size_t numberOfPeaks = scans.getNbrPeaks(scanNumber);
            for(size_t pkNum = 0; pkNum < numberOfPeaks; ++pkNum) {
                 mzVals.push_back(static_cast<float>(xAxisMasses[pkNum]));
                 intensityVals.push_back(static_cast<float>(yAxis.first[pkNum]));
            }

            MsParquetReaderRow msParquetReaderRow;
            msParquetReaderRow.scanNumber = timsFrameInfo.frameId;
            msParquetReaderRow.scanTime = timsFrameInfo.scanTime;
            msParquetReaderRow.msLevel = 2;
            msParquetReaderRow.collisionEnergy = tmwi.collisionEnergy;
            msParquetReaderRow.precursorTargetMz = tmwi.isolationMz;
            msParquetReaderRow.isoWindowLower = tmwi.isolationWidth / 2.0f;
            msParquetReaderRow.isoWindowUpper = msParquetReaderRow.isoWindowLower;
            msParquetReaderRow.ionMobilityIndex = scanNumber;
            msParquetReaderRow.ionMobilityDriftTime = m_frameIndexVsDriftTime.value(msParquetReaderRow.ionMobilityIndex);
            msParquetReaderRow.mzVals = mzVals;
            msParquetReaderRow.intensityVals = intensityVals;
            msParquetReaderRow.targetKey = mzTargetKey;

            msParquetReaderRows.push_back(msParquetReaderRow);
        }
    }

    e = ParquetReader::write(msParquetReaderRows, filePath); ree;

    ERR_RETURN
}
