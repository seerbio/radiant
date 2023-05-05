//
// Created by anichols on 4/8/23.
//

#include "MsFrame.h"

#include "DeisotoperTandem.h"
#include "EigenKernelUtils.h"
#include "EigenSparseUtils.h"
#include "MsReaderBase.h"
#include "MsScansDenoiseTron.h"
#include "ParallelUtils.h"

#include <QElapsedTimer>


MsFrame::MsFrame()
: m_mzWindowLower(-1.0)
, m_mzWindowUpper(-1.0)

{}

Err MsFrame::init(
        const PythiaParameters &pythiaParameters,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QMap<ScanNumber, ScanPoints> &scanPoints,
        const QPair<double, double> &frameMzStartStop
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    m_frame = scanPoints;

    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;

    m_mzWindowLower = frameMzStartStop.first;
    m_mzWindowUpper = frameMzStartStop.second;

    e = buildFrameIndexVsScanNumber(); ree;

    ERR_RETURN
}

Err MsFrame::init(
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanPoints> &scanPoints,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        double collisionEnergy,
        double precursorTargetMz,
        double isoWindowLower,
        double isoWindowUpper
        ) {

    ERR_INIT

    e = init(
            pythiaParameters,
            uniqueMsInfoScanKey,
            scanPoints,
            {precursorTargetMz - isoWindowLower, precursorTargetMz + isoWindowUpper}
            ); ree;

    const double minVal = 0.0;

    for (const double val : {collisionEnergy, precursorTargetMz, isoWindowLower, isoWindowUpper}) {
        e = ErrorUtils::isAboveThreshold(
                val,
                minVal,
                ErrorUtilsParam::ExcludeThreshold
        ); ree;
    }

    ERR_RETURN
}

Err MsFrame::preprocessMsFrame(
        bool denoise,
        bool deisotope,
        bool smooth
        ) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    if (denoise) {
        e = denoiseFrame(); ree;
    }

    if (deisotope) {
        e = deisotopeFrame(); ree;
    }

    if (smooth) {
        e = smoothFrame(); ree;
    }

    qDebug() << "File preprocessed in" << et.elapsed() << "mSec";

    ERR_RETURN

}

Err MsFrame::denoiseFrame() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_frame); ree;

    MsScansDenoiseTron denoiser;
    e = denoiser.init(m_params); ree;

    QMap<ScanNumber, ScanPoints> denoisedFrame;
    e = denoiser.denoiseScansFrame(
            m_frame,
            &denoisedFrame
    ); ree;

    m_frame = denoisedFrame;

    ERR_RETURN
}

Err MsFrame::deisotopeFrame() {

    ERR_INIT

    const bool runParallel = false;

    e = ErrorUtils::isNotEmpty(m_frame); ree;

    QMap<ScanNumber, ScanPoints> deisotopedTandemScans;
    e = DeisotoperTandem::deisotopeTandemScansParallel(
            m_frame,
            m_params.ms2ExtractionWidthPPM,
            runParallel,
            &deisotopedTandemScans
    ); ree;

    m_frame = deisotopedTandemScans;

    ERR_RETURN
}

namespace {

    Eigen::SparseMatrix<double> rowWiseGaussianSmoothMatrix(
            const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat
            ) {

        //TODO find a way to autoset filters
        const int filterLen = 3;
        const double sigma = 2;
        Eigen::VectorX<double> gaussianFilter = EigenKernelUtils::buildGaussianFilter1D(
                filterLen,
                sigma
                );

        return EigenKernelUtils::applyKernelRowWiseToMatrix(
                mat,
                gaussianFilter,
                false
                );
    }

    Eigen::SparseMatrix<double> colWiseGaussianSmoothMatrix(
            const Eigen::SparseMatrix<double> &mat,
            double ppm
    ) {

       const auto ppmInt = static_cast<int>(std::round(ppm));

        const int filterLen =  ppmInt % 2 == 0 ? ppmInt - 1 : ppmInt;
        const double sigma = 2;
        Eigen::VectorX<double> gaussianFilter = EigenKernelUtils::buildGaussianFilter1D(
                filterLen,
                sigma
        );

        return EigenKernelUtils::applyKernelColumnWiseToMatrix(
                mat,
                gaussianFilter,
                false
        );
    }

}//namespace
Err MsFrame::smoothFrame() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_frame); ree;

    const int precision = 4;
    const double mzMax = 2000;

    const QMap<FrameIndex, ScanPoints> frame = frameIndexVsScanPoints();

    Eigen::SparseMatrix<double, Eigen::RowMajor> mat
        = EigenSparseUtils::loadFrameToSparseMatrixRowMajor(
            frame,
            precision,
            mzMax
            );

    //TODO add smoothing params to pythiaParams.
    mat = rowWiseGaussianSmoothMatrix(mat); ree;

    QMap<FrameIndex, ScanPoints> frameSmooth
        = EigenSparseUtils::loadSparseMatrixToFrame(mat, precision);

    m_frame.clear();

    for (auto it = frameSmooth.begin(); it != frameSmooth.end(); it++) {
        const FrameIndex frameIndex = it.key();
        const ScanPoints &scanPoints = it.value();
        const ScanNumber scanNumber = m_frameIndexVsScanNumber.value(frameIndex);
        m_frame.insert(scanNumber, scanPoints);
    }

    ERR_RETURN
}

namespace {

QMap<int, double> eigenUtilsApexWrapper(const Eigen::SparseVector<double> &vec) {
    return EigenSparseUtils::apexes(vec);
}

}//namespace
Err MsFrame::gaussianSmooth2D() {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isNotEmpty(m_frame); ree;
    e = ErrorUtils::isTrue(m_params.isValid()); ree;

    const int precision = 3;
    const double mzMax = 2000;

    const QMap<FrameIndex, ScanPoints> frame = frameIndexVsScanPoints();

    Eigen::SparseMatrix<double, Eigen::ColMajor> mat = EigenSparseUtils::loadFrameToSparseMatrixColMajor(
            frame,
            precision,
            mzMax
    );

    //TODO add smoothing params to pythiaParams.
    mat = colWiseGaussianSmoothMatrix(mat, m_params.featureFinderTolerancePPM); ree;
    mat = rowWiseGaussianSmoothMatrix(mat); ree;

    QVector<Eigen::SparseVector<double>> vecs;
    for (int i = 0; i < mat.cols(); ++i) {
        Eigen::SparseVector<double> scanVec = mat.col(i);
        vecs.push_back(scanVec);
    }

//#define PARALLEL_APEX_EXTRACT
#ifdef PARALLEL_APEX_EXTRACT
    QFuture<QMap<int, double>> apexVecs = QtConcurrent::mapped(
            vecs,
            eigenUtilsApexWrapper
            );
    apexVecs.waitForFinished();
#else
    QVector<QMap<int, double>> apexVecs;
    for (const Eigen::SparseVector<double> &vec : vecs) {
        apexVecs.push_back(eigenUtilsApexWrapper(vec));
    }
#endif

    m_frame.clear();

    FrameIndex frameIndex = 0;
    for (const QMap<int, double> &res : apexVecs) {

        for (auto it = res.begin(); it != res.end(); it++){
            const int hashedMz = it.key();
            const auto mz = MathUtils::unHashDecimal<double>(hashedMz, precision);
            const double intensityVal = it.value();

            m_frame[scanNumberFromFrameIndex(frameIndex)].push_back({mz, intensityVal});
        }

        frameIndex++;
    }

    qDebug() << "frame 2d smoothed in" << et.elapsed() << "mSec";

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

Err MsFrame::writeFramScans(const QString &outputFilePath) const {

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

Err MsFrame::buildFrameIndexVsScanPoints(
        const QVector<MsFrameScanPointRows> &msFrameScanPointRows,
        QMap<FrameIndex, ScanPoints> *frameIndexVsScanPoints
        ) {

    ERR_INIT

    frameIndexVsScanPoints->clear();
    e = ErrorUtils::isNotEmpty(msFrameScanPointRows); ree;

    for (const MsFrameScanPointRows &row : msFrameScanPointRows) {

        e = ErrorUtils::isEqual(
                row.mzVals.size(),
                row.intensityVals.size()
                ); ree;

        ScanPoints scanPoints;
        e = ParallelUtils::zip(
                row.mzVals,
                row.intensityVals,
                &scanPoints
        ); ree

        frameIndexVsScanPoints->insert(row.frameIndex, scanPoints);
    }

    ERR_RETURN
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
