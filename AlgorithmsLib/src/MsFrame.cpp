//
// Created by anichols on 4/8/23.
//

#include "MsFrame.h"

#include "EigenKernelUtils.h"
#include "EigenSparseUtils.h"
#include "EigenUtils.h"
#include "GlobalSettings.h"
#include "MS2ChargeDeconvolvotron.h"
#include "MsReaderBase.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"

#include <QElapsedTimer>

using SparseMatrixPoint = EigenSparseUtils::SparseMatrixPoint;

const int PRECISION = 3;


Err MsFrame::init(
        QMap<ScanNumber, ScanPoints> &scanPoints,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;

    m_frame = scanPoints;

    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree
    m_scanNumberVsScanTime = scanNumberVsScanTime;

    e = buildFrameIndexVsScanNumber(); ree

    ERR_RETURN
}

Err MsFrame::buildFrameIndexVsScanNumber() {

    ERR_INIT

    e = ErrorUtils::isFalse(m_frame.isEmpty()); ree;

    for (const ScanNumber &sn : m_frame.keys()) {
        m_frameIndexVsScanNumber.insert(m_frameIndexVsScanNumber.size(), sn);
    }

    ERR_RETURN
}

int MsFrame::scanCount() const {
    return m_frameIndexVsScanNumber.size();
}

QMap<FrameIndex, ScanPoints> MsFrame::frameIndexVsScanPoints() const {

    QMap<FrameIndex, ScanPoints> frameIndexVsScanPoints;

    for (const ScanPoints &sp : m_frame) {
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

ScanPoints MsFrame::getScanPointsByScanNumber(ScanNumber scanNumber) const {
    return m_frame.value(scanNumber);
}

QMap<ScanNumber, ScanPoints> MsFrame::scanNumberVsScanPoints() const {
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