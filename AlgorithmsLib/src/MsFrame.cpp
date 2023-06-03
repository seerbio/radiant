//
// Created by anichols on 4/8/23.
//

#include "MsFrame.h"

#include "DeisotoperTandem.h"
#include "EigenKernelUtils.h"
#include "EigenSparseUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "GlobalSettings.h"
#include "MsReaderBase.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include <QElapsedTimer>

using SparseMatrixPoint = EigenSparseUtils::SparseMatrixPoint;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using rTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
using rTreeSearchBox = bg::model::box<rTreeCoor>;
using rTreePoint = std::pair<rTreeCoor, double> ;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;

const int PRECISION = 3;


MsFrame::MsFrame()
: m_mzWindowLower(-1.0)
, m_mzWindowUpper(-1.0)

{}

Err MsFrame::init(
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QMap<ScanNumber, ScanPoints> &scanPoints,
        const QPair<double, double> &frameMzStartStop
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    m_frame = scanPoints;

    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;

    m_mzWindowLower = frameMzStartStop.first;
    m_mzWindowUpper = frameMzStartStop.second;

    e = buildFrameIndexVsScanNumber(); ree;

    ERR_RETURN
}

QPair<double, double> MsFrame::precursorMzTargetStartEnd() const {
    return {m_mzWindowLower, m_mzWindowUpper};
}

UniqueMsInfoScanKey MsFrame::uniqueMsInfoScanKey() const {
    return m_uniqueMsInfoScanKey;
}

Err MsFrame::buildFrameIndexVsScanNumber() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_frame); ree;

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

Err MsFrame::writeFrameScans(const QString &outputFilePath) const {

    ERR_INIT

    const QMap<FrameIndex, ScanPoints> framesVsScanPoints = frameIndexVsScanPoints();

    QVector<MsFrameScanPointRows> rowsToWrite;
    for (auto it = framesVsScanPoints.begin(); it != framesVsScanPoints.end(); it++) {
        MsFrameScanPointRows row;
        row.frameIndex = it.key();

        const ScanPoints &sp = it.value();

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

double MsFrame::meanPrecursorRange() const {
    return (m_mzWindowLower + m_mzWindowUpper) / 2.0;
}

ScanNumber MsFrame::scanNumberFromFrameIndex(FrameIndex frameIndex) const {
    return m_frameIndexVsScanNumber.value(frameIndex);
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

Err MsFrame::buildMsFrame(
        const QString &msDataFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        QPair<double, double> mzTargetStartStop,
        MsFrame *msFrame
        ) {

    ERR_INIT

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(
            msDataFilePath,
            MsParquetReaderNamespace::PERCURSOR_TARGET_MZ,
            {mzTargetStartStop.first, mzTargetStartStop.second}
    ); ree;

    const QMap<ScanNumber, ScanPoints> targetScanPoints = msReaderParquet.getScanPoints();
    e = ErrorUtils::isNotEmpty(targetScanPoints); ree;

    e = msFrame->init(
            uniqueMsInfoScanKey,
            targetScanPoints,
            mzTargetStartStop
    ); ree;

    ERR_RETURN
}

namespace {

    Err smoothEigenMatrix(
            int gaussFilterLength,
            double sigma,
            int smoothCount,
            Eigen::SparseMatrix<double, Eigen::ColMajor> *mat
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(mat->nonZeros() > 0); ree

        const Eigen::VectorX<double> gaussKernel = EigenKernelUtils::buildGaussianFilter1D(
                gaussFilterLength,
                sigma,
                false
                );

        for (int i = 0; i <= smoothCount; i++) {

            *mat = EigenKernelUtils::applyKernelColumnWiseToMatrix(
                    *mat,
                    gaussKernel,
                    false
                    );
        }

        ERR_RETURN
    }

}//namespace
Err MsFrame::smoothFrame(
        int gaussFilterLength,
        double sigma,
        int smoothCount,
        double mzMax
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(isValid()); ree

    const QMap<FrameIndex, ScanPoints> scanPoints = frameIndexVsScanPoints();

    Eigen::SparseMatrix<double, Eigen::ColMajor> mat = EigenSparseUtils::loadFrameToSparseMatrixColMajor(
            scanPoints,
            PRECISION,
            mzMax
            );

    e = smoothEigenMatrix(
            gaussFilterLength,
            sigma,
            smoothCount,
            &mat
            ); ree

    const QMap<int, ScanPoints> reFrame = EigenSparseUtils::loadSparseMatrixToFrame(
            mat,
            PRECISION
            );

    e = reloadMFrame(reFrame); ree;

    ERR_RETURN
}

Err MsFrame::reloadMFrame(const QMap<FrameIndex, ScanPoints> &scanPoints) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree

    m_frame.clear();

    for (auto it = scanPoints.begin(); it != scanPoints.end(); it++) {

        const FrameIndex frameIndex = it.key();
        const ScanNumber scanNumber = scanNumberFromFrameIndex(frameIndex);
        const ScanPoints &fiScanPoints = it.value();

        m_frame.insert(scanNumber, fiScanPoints);
    }

    ERR_RETURN
}
