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
        const QMap<ScanNumber, ScanPoints> &scanPoints,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree

    m_frame = scanPoints;

    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree
    m_scanNumberVsScanTime = scanNumberVsScanTime;

    e = buildFrameIndexVsScanNumber(); ree

    ERR_RETURN
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

Err MsFrame::filterFrameByMz(
        double mzMin,
        double mzMax
        ) {

    ERR_INIT

    const auto terminatorLogic = [&](const ScanPoint &p){return !(mzMin < p.x() && p.x() < mzMax);};

    for (ScanPoints &pnts : m_frame) {

        const auto terminator = std::remove_if(pnts.begin(), pnts.end(), terminatorLogic);

        pnts.erase(terminator, pnts.end());
    }

    ERR_RETURN
}

ScanNumber MsFrame::scanNumberFromScanTime(ScanTime scanTime) const {

    const QVector<ScanNumber> scanNumbers = m_scanNumberVsScanTime.keys().toVector();
    const QVector<ScanTime> scanTimes = m_scanNumberVsScanTime.values().toVector();

    const int closestIndex = MathUtils::closest(scanTimes, scanTime);

    return scanNumbers.at(closestIndex);
}

Err MsFrame::deisotopeMsFrame(double ppmTol) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_frame);

    const QString &chargeModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");

    const QString &monoModelFilePath
            = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");

    MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
    e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);

    QMap<ScanNumber, ScanPoints> scanPointsDeisotoped;
    for (auto it = m_frame.begin(); it != m_frame.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const ScanPoints &scanPointsIter = it.value();

        ScanPoints scanPointsIterDeisotoped;
        e = ms2ChargeDeconvolvotron.deisotopeScanPoints(scanPointsIter, &scanPointsIterDeisotoped);

        scanPointsDeisotoped.insert(scanNumber, scanPointsIterDeisotoped);
    }

    m_frame = scanPointsDeisotoped;
    e = buildFrameIndexVsScanNumber(); ree

    ERR_RETURN
}
