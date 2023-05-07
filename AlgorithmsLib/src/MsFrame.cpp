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
#include "MsScansDenoiseTron.h"
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

namespace {

    Eigen::SparseMatrix<double> rowWiseGaussianSmoothMatrix(
            const Eigen::SparseMatrix<double, Eigen::RowMajor> &mat
    ) {

        //TODO find a way to autoset filters
        const int filterLen = 3;
        const double sigma = 2.0;
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

        const int filterLen =  ppmInt % 2 == 0 ? ppmInt + 1 : ppmInt;
        const double sigma = 1;
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

    Eigen::SparseMatrix<double> frameToSparseMatrixSmoothed(
            const QMap<FrameIndex, ScanPoints> &frame,
            double featureFinderTolerancePPM
            ) {

        const double mzMax = 2000;

        Eigen::SparseMatrix<double, Eigen::ColMajor> mat = EigenSparseUtils::loadFrameToSparseMatrixColMajor(
                frame,
                PRECISION,
                mzMax
        );

        //TODO add smoothing params to pythiaParams.
        mat = colWiseGaussianSmoothMatrix(mat, featureFinderTolerancePPM);
        mat = rowWiseGaussianSmoothMatrix(mat);

        return mat;
    }

    void sortHillsFrameIndexThenMzAsc(QVector<FeatureFinderHillPoint> *hillsPoints) {

        const auto sortLogic = [](const FeatureFinderHillPoint &l, const FeatureFinderHillPoint &r){

            if (l.frameIndex == r.frameIndex) {
                return l.mz < r.mz;
            }

            return l.frameIndex < r.frameIndex;
        };

        std::sort(hillsPoints->begin(), hillsPoints->end(), sortLogic);
    }

    RTree loadRTree(QVector<FeatureFinderHillPoint> &featureFinderHillPoints) {

        QElapsedTimer et;
        et.start();

        std::vector<rTreePoint> cloudLoader;

        for (const FeatureFinderHillPoint &h : featureFinderHillPoints) {

            rTreeCoor coor(h.frameIndex, h.mz);
            cloudLoader.emplace_back(coor, h.intensity);
        };

        const int maxElements = 16;
        RTree rTree(cloudLoader, bgi::dynamic_quadratic(maxElements));

        qDebug() << rTree.size()
                 << "Frame Apexes loaded into rtree in"
                 << et.restart()
                 << "mSec";

        return rTree;
    }

    const auto rTreeMaxPointLogic = [](const rTreePoint &l, const rTreePoint &r){
        return l.second < r.second;
    };

    rTreePoint getMaxRTreePointResultValue(const std::vector<rTreePoint> &result) {

        const rTreePoint searchedMaxIntensityPoint = *std::max_element(
                result.rbegin(),
                result.rend(),
                rTreeMaxPointLogic
        );

        return searchedMaxIntensityPoint;
    }

    Err iterativelyDeisotopeTree(
            double featureFinderPPM,
            RTree *rtree,
            QVector<FeatureFinderHillPoint> *featureFinderHillApexes
            ) {

        ERR_INIT

        const Charge chargeMax = 2;
        const int maxSearchDepth = 4;
        const int frameIndexBuffer = 1;

        while (!rtree->empty()) {

            const rTreePoint rTreeMaxPoint = *std::max_element(
                    rtree->begin(),
                    rtree->end(),
                    rTreeMaxPointLogic
                    );

            FeatureFinderHillPoint searchPoint(
                    static_cast<int>(rTreeMaxPoint.first.get<0>()),
                    rTreeMaxPoint.first.get<1>(),
                    rTreeMaxPoint.second
                    );
            featureFinderHillApexes->push_back(searchPoint);

            for (int charge = chargeMax; charge > 0; charge--) {

                const double chargeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / charge;
                const double mzTol = MathUtils::calculatePPM(searchPoint.mz, featureFinderPPM);

                double lastSearchedMaxIntensity = searchPoint.intensity;

                for (int depth = 1; depth < maxSearchDepth; depth++) {

                    const double mzMin = searchPoint.mz + (depth * chargeDistance) - mzTol;
                    const double mzMax = searchPoint.mz + (depth * chargeDistance) + mzTol;
                    const double frameIndexMin = searchPoint.frameIndex - frameIndexBuffer;
                    const double frameIndexMax = searchPoint.frameIndex + frameIndexBuffer;

                    const rTreeSearchBox queryBox(
                            rTreeCoor(frameIndexMin, mzMin),
                            rTreeCoor(frameIndexMax, mzMax)
                    );

                    std::vector<rTreePoint> result;
                    rtree->query(bgi::intersects(queryBox), std::back_inserter(result));

                    if (result.empty()) {
                        break;
                    }

                    const rTreePoint searchedMaxIntensityPoint = getMaxRTreePointResultValue(result);
                    const double searchedMaxIntensityPointValue = searchedMaxIntensityPoint.second;

                    if (searchedMaxIntensityPointValue < lastSearchedMaxIntensity) {
                        lastSearchedMaxIntensity = searchedMaxIntensityPointValue;
                        rtree->remove(searchedMaxIntensityPoint);
                        continue;
                    }

                    break;
                }
            }

            rtree->remove(rTreeMaxPoint);
        }

        ERR_RETURN
    }

}//namespace
Err MsFrame::deisotopeFrame() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_frame); ree;
    e = ErrorUtils::isTrue(m_params.isValid()); ree;

    FeatureFinderHillBuilder featureFinderHillBuilder;

    FeatureFinderParameters ffParams;
    ffParams.tolerancePPM = m_params.featureFinderTolerancePPM;
    ffParams.minScanCount = 1;
    ffParams.skipScanCount = 1;
    ffParams.useMeanMz = true;
    ffParams.filterLength = 3;
    ffParams.smoothCount = 1;

    e = featureFinderHillBuilder.init(ffParams); ree;
    featureFinderHillBuilder.setRunParallel(false);

    QVector<FeatureFinderHill> featureFinderHills;
    e = featureFinderHillBuilder.buildHills(
            scanNumberVsScanPoints(),
            &featureFinderHills
    ); ree;

//    e = featureFinderHillBuilder.refineHills(&featureFinderHills); ree;

    QVector<FeatureFinderHillPoint> featureFinderHillPoints;
    e = FeatureFinderHillBuilder::featureFinderHillPoints(
            featureFinderHills,
            &featureFinderHillPoints
            ); ree;
    sortHillsFrameIndexThenMzAsc(&featureFinderHillPoints);

    RTree rtree = loadRTree(featureFinderHillPoints);
    QVector<FeatureFinderHillPoint> featureFinderHillApexes;
    e = iterativelyDeisotopeTree(
            ffParams.tolerancePPM,
            &rtree,
            &featureFinderHillApexes
            ); ree;

    ERR_RETURN
}

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

    const Eigen::SparseMatrix<double> mat = frameToSparseMatrixSmoothed(
            frameIndexVsScanPoints(),
            m_params.featureFinderTolerancePPM
            );

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
            const auto mz = MathUtils::unHashDecimal<double>(hashedMz, PRECISION);
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
