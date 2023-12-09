//
// Created by anichols on 4/8/23.
//

#include "MsFrame.h"

#include "EigenKernelUtils.h"
#include "EigenSparseUtils.h"
#include "EigenUtils.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"

#include <nanoflann.hpp>

#include <QElapsedTimer>


class Q_DECL_HIDDEN MsFrame::Private
{


public:

    using KDTree = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::MatrixXd>;

    Private();
    ~Private();

    Err init(const QMap<FrameIndex, ScanTime> &frameIndexVsScanTime);

    Err frameIndexFromScanTime(ScanTime scanTime, FrameIndex *frameIndex);

    bool isInit();

private:

    QMap<FrameIndex, ScanTime> m_frameIndexVsScanTime;
    QMap<Index, FrameIndex > m_indexVsFrameIndex;
    bool m_isInit;

    Eigen::MatrixX<double> *m_mat;
    KDTree *m_kdTree;

};

MsFrame::Private::Private() : m_isInit(false) {}

MsFrame::Private::~Private() {
    delete m_kdTree;
    delete m_mat;
}


Err MsFrame::Private::init(const QMap<FrameIndex, ScanTime> &frameIndexVsScanTime) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(frameIndexVsScanTime); ree;
    m_frameIndexVsScanTime = frameIndexVsScanTime;

    const int nDims = 2;
    m_mat = new Eigen::MatrixX<double>(frameIndexVsScanTime.size(), nDims);

    int index = 0;
    for (auto it = m_frameIndexVsScanTime.begin(); it != m_frameIndexVsScanTime.end(); it++) {

        const FrameIndex frameIndex = it.key();
        const ScanTime  scanTime = it.value();

        m_mat->coeffRef(index, 0) = scanTime;
        m_mat->coeffRef(index, 1) = 0.0;

        m_indexVsFrameIndex.insert(index++, frameIndex);
    }

    const int treeLeafSize = 15;
    m_kdTree = new KDTree(nDims, *m_mat, treeLeafSize);
    m_isInit = true;

    ERR_RETURN
}

bool MsFrame::Private::isInit() {
    return m_isInit;
}

Err MsFrame::Private::frameIndexFromScanTime(ScanTime scanTime, FrameIndex *frameIndex) {

    ERR_INIT

    e = ErrorUtils::isTrue(isInit());

    const size_t numResults = 1;
    std::vector<double> queryPt = {scanTime, 0.0};
    std::vector<long> retIndex(numResults);
    std::vector<double> outDistSqr(numResults);

    std::vector<std::pair<Eigen::Index, double>> matches;

    const size_t resultsSize = m_kdTree->index->knnSearch(
            queryPt.data(),
            numResults,
            retIndex.data(),
            outDistSqr.data()
    );

    *frameIndex = m_indexVsFrameIndex.value(static_cast<int>(retIndex.front()));

    ERR_RETURN
}



///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

MsFrame::MsFrame() : d_ptr(new Private()) {}

MsFrame::~MsFrame() {}

Err MsFrame::init(
        const QMap<ScanNumber, ScanPoints*> &scanPoints,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;

    m_frame = scanPoints;

    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree
    m_scanNumberVsScanTime = scanNumberVsScanTime;

    e = buildFrameIndexVsScanNumber(); ree
    e = initFrameIndexVsScanTimeKDTree(); ree;

    ERR_RETURN
}

Err MsFrame::buildFrameIndexVsScanNumber() {

    ERR_INIT

    e = ErrorUtils::isFalse(m_frame.isEmpty()); ree;

    for (const ScanNumber sn : m_frame.keys()) {
        m_frameIndexVsScanNumber.insert(m_frameIndexVsScanNumber.size(), sn);
    }

    ERR_RETURN
}

int MsFrame::scanCount() const {
    return m_frameIndexVsScanNumber.size();
}

QMap<FrameIndex, ScanPoints*> MsFrame::frameIndexVsScanPoints() const {

    QMap<FrameIndex, ScanPoints*> frameIndexVsScanPoints;

    for (ScanPoints *sp : m_frame) {
        frameIndexVsScanPoints.insert(frameIndexVsScanPoints.size(), sp);
    }

    return frameIndexVsScanPoints;
}

Err MsFrame::writeFrameScans(
        const QMap<FrameIndex, ScanPoints> &framesVsScanPoints,
        const QString &outputFilePath
        ) {

    ERR_INIT

    QVector<MsFrameScanPointRows> rowsToWrite;
    for (auto it = framesVsScanPoints.begin(); it != framesVsScanPoints.end(); it++) {
        MsFrameScanPointRows row;
        row.frameIndex = it.key();

        const ScanPoints &sp = it.value();

        if (sp.isEmpty()) {
            continue;
        }

        e = MsReaderBase::splitScanPoints(
                sp,
                &row.mzVals,
                &row.intensityVals
                ); ree;

        rowsToWrite.push_back(row);
    }

    e = ParquetReader::write(
            rowsToWrite,
            outputFilePath
            ); ree;

    ERR_RETURN
}

ScanNumber MsFrame::scanNumberFromFrameIndex(FrameIndex frameIndex) const {
    return m_frameIndexVsScanNumber.value(frameIndex);
}

ScanTime MsFrame::scanTimeFromScanNumber(ScanNumber scanNumber) const {
    return m_scanNumberVsScanTime.value(scanNumber);
}

ScanNumber MsFrame::frameIndexFromScanNumber(ScanNumber scanNumber) const {
    return m_frameIndexVsScanNumber.key(scanNumber);
}

ScanPoints* MsFrame::getScanPointsByScanNumber(ScanNumber scanNumber) const {
    return m_frame.value(scanNumber);
}

QMap<ScanNumber, ScanPoints*> MsFrame::scanNumberVsScanPoints() const {
    return m_frame;
}

bool MsFrame::isValid() const {
    const int minScanCount = 1;
    return scanCount() > minScanCount;
}

ScanNumber MsFrame::scanNumberFromScanTime(ScanTime scanTime) const {

    const QVector<ScanNumber> scanNumbers = m_scanNumberVsScanTime.keys().toVector();
    const QVector<ScanTime> scanTimes = m_scanNumberVsScanTime.values().toVector();

    const int closestIndex = MathUtils::closest(scanTimes, scanTime);

    return scanNumbers.at(closestIndex);
}

Err MsFrame::frameIndexFromScanTime(ScanTime scanTime, FrameIndex *frameIndex) const {
    ERR_INIT
    e = d_ptr->frameIndexFromScanTime(scanTime, frameIndex); ree;
    ERR_RETURN
}

namespace {

    Err buildFrameIndexVsScanTime(
            const QMap<FrameIndex, ScanNumber> &frameIndexVsScanNumber,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            QMap<FrameIndex, ScanTime> *frameIndexVsScanTime
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(frameIndexVsScanNumber); ree;
        e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;

        for (auto it = frameIndexVsScanNumber.begin(); it != frameIndexVsScanNumber.end(); it++) {
            const FrameIndex frameIndex = it.key();
            const ScanNumber scanNumber = it.value();
            const ScanTime scanTime = scanNumberVsScanTime.value(scanNumber);
            frameIndexVsScanTime->insert(frameIndex, scanTime);
        }

        ERR_RETURN
    }

}//namespace
Err MsFrame::initFrameIndexVsScanTimeKDTree() {

    ERR_INIT

    QMap<FrameIndex, ScanTime> frameIndexVsScanTime;
    e = buildFrameIndexVsScanTime(
            m_frameIndexVsScanNumber,
            m_scanNumberVsScanTime,
            &frameIndexVsScanTime
            ); ree;

    d_ptr->init(frameIndexVsScanTime); ree;

    ERR_RETURN
}
