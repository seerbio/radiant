//
// Created by anichols on 2/5/24.
//

#include "DeMultiplexScanerator.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MathUtils.h"
#include "MsReaderBase.h"

#include <Eigen/Dense>
#include <unsupported/Eigen/NNLS>

#include <iostream>

class Q_DECL_HIDDEN DeMultiplexScanerator::Private
{

public:

    Private(
            double ppmTol,
            double intensityFractionThreshold
            );

    ~Private();

    Err buildScanMaskMatrix(
            const QVector<MsScanInfo> &msScanInfos,
            Eigen::MatrixX<float> *scanMaskMatrix,
            QPair<float, float> *totalScanRange,
            double *amuStagger
            );

    Err buildTransitionsMatrix(
            const QVector<ScanPoints*> &scans,
            Eigen::MatrixX<float> *transitionsMatrix,
            QVector<float> *mzTransitions
            );

    Err deMultiplexScans(
            const QVector<ScanPoints*> &scanPoints,
            const QVector<MsScanInfo> &msScanInfos,
            QVector<ScanPoints> *demultiplexedScans,
            QVector<MzTargetKey> *mzTargetKeys,
            double *windowSize
            );

private:

    double m_ppmTol;
    double m_intensityFractionThreshold;

};

DeMultiplexScanerator::Private::Private(
        double ppmTol,
        double intensityFractionThreshold
        )
        : m_ppmTol(ppmTol)
        , m_intensityFractionThreshold(intensityFractionThreshold)
        {}

DeMultiplexScanerator::Private::~Private() {}

namespace {

    struct ScanPointsSorter {
        ScanPoint *scanPoint;
        int scansIndex = -1;
    };

}//namespace
Err DeMultiplexScanerator::Private::buildTransitionsMatrix(
        const QVector<ScanPoints*> &scans,
        Eigen::MatrixX<float> *transitionsMatrix,
        QVector<float> *mzTransitions
        ) {

    ERR_INIT

    QVector<ScanPointsSorter> scanPointsSorters;

    for (int scanIndex = 0; scanIndex < scans.size(); scanIndex++) {

        ScanPoints *scanPoints = scans.at(scanIndex);

        for (ScanPoint &sp : *scanPoints) {
            ScanPointsSorter sps;
            sps.scanPoint = &sp;
            sps.scansIndex = scanIndex;

            scanPointsSorters.push_back(sps);
        }
    }

    const auto sortLogic = [](const ScanPointsSorter &l, const ScanPointsSorter &r){
        return l.scanPoint->x() < r.scanPoint->x();
    };
    std::sort(scanPointsSorters.begin(), scanPointsSorters.end(), sortLogic);

    QVector<QVector<ScanPointsSorter>> transitions;

    QVector<ScanPointsSorter> currentScanPointsSorter = {scanPointsSorters.at(0)};
    for (int i = 1; i < scanPointsSorters.size() - 1; i++) {

        const ScanPointsSorter &spsL = currentScanPointsSorter.back();
        const ScanPointsSorter &spsR = scanPointsSorters.at(i);

        const double ppmDiff
            = 1e6 * std::abs(spsL.scanPoint->x() - spsR.scanPoint->x()) / spsL.scanPoint->x();

        if (ppmDiff <= m_ppmTol) {

            if (spsL.scansIndex == spsR.scansIndex) {
                continue;
            }

            currentScanPointsSorter.push_back(spsR);
            continue;
        }

        transitions.push_back(currentScanPointsSorter);
        currentScanPointsSorter.clear();
        currentScanPointsSorter.push_back(spsR);
    }

    transitionsMatrix->resize(scans.size(), transitions.size());
    transitionsMatrix->setZero();

    for (int col = 0; col < transitions.size(); col++) {

        const QVector<ScanPointsSorter> &sps = transitions.at(col);

        const double transitionMean = std::accumulate(
                sps.begin(),
                sps.end(),
                0.0,
                [](double sum, const ScanPointsSorter &sps){return sum + sps.scanPoint->x();}
                ) / sps.size();

        mzTransitions->push_back(transitionMean);

        for (const ScanPointsSorter &sp : sps) {
            transitionsMatrix->coeffRef(sp.scansIndex, col) = sp.scanPoint->y();
        }
    }

    ERR_RETURN
}

namespace {

    Err getTotalScanRange(
            const QVector<MsScanInfo> &msScanInfos,
            QPair<float, float> *totalScanRange
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(msScanInfos); ree;

        *totalScanRange = {std::numeric_limits<float>::max(), std::numeric_limits<float>::lowest()};

        for (const MsScanInfo &msScanInfo : msScanInfos) {
            const float mzLo = msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower;
            totalScanRange->first = std::min(totalScanRange->first, mzLo);

            const float mzHi = msScanInfo.precursorTargetMz + msScanInfo.isoWindowLower;
            totalScanRange->second = std::max(totalScanRange->second, mzHi);
        }

        ERR_RETURN
    }

    Err calculateWindowSize(const QVector<MsScanInfo> &msScanInfos, double *windowSize) {

        ERR_INIT

        QVector<double> mzStarts;
        std::transform(
                msScanInfos.begin(),
                msScanInfos.end(),
                std::back_inserter(mzStarts),
                [](const MsScanInfo &msScanInfo){return msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower;}
                );

        std::sort(mzStarts.begin(), mzStarts.end());

        QVector<double> distances;
        for (int i = 0; i < mzStarts.size() - 1; i++) {
            const double mzLeft = mzStarts.at(i);
            const double mzRight = mzStarts.at(i+1);
            distances.append(mzRight - mzLeft);
        }

        const double meanDistance = MathUtils::mean(distances);

        const bool allDistancesAreEqual = std::all_of(
                distances.begin(),
                distances.end(),
                [meanDistance](double d){return MathUtils::tSame(meanDistance, d);}
                );

        e = ErrorUtils::isTrue(allDistancesAreEqual); ree;

        *windowSize = meanDistance;

        ERR_RETURN
    }

}//namespace
Err DeMultiplexScanerator::Private::buildScanMaskMatrix(
        const QVector<MsScanInfo> &msScanInfos,
        Eigen::MatrixX<float> *scanMaskMatrix,
        QPair<float, float> *totalScanRange,
        double *windowSize
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msScanInfos); ree;

    e = getTotalScanRange(msScanInfos, totalScanRange); ree;
    e = ErrorUtils::isBelowThreshold(
            totalScanRange->first,
            totalScanRange->second,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;

    e = calculateWindowSize(msScanInfos, windowSize); ree;

    const int rows = msScanInfos.size();
    const int cols = static_cast<int>((totalScanRange->second - totalScanRange->first) / *windowSize);
    scanMaskMatrix->resize(rows, cols);
    scanMaskMatrix->setZero();

    for (int i = 0; i < msScanInfos.size(); i++) {

        const MsScanInfo &msi = msScanInfos.at(i);

        const int mzStart
            = static_cast<int>((msi.precursorTargetMz - msi.isoWindowLower) - totalScanRange->first);

        const int mzEnd
            = static_cast<int>((msi.precursorTargetMz + msi.isoWindowUpper) - totalScanRange->first);

        for (int j = mzStart; j < mzEnd; j++) {
            scanMaskMatrix->coeffRef(i, j) = 1.0;
        }
    }

    ERR_RETURN
}

namespace{

    Err buildDemultiplexedScans(
            const Eigen::MatrixX<float> &xMat,
            const QVector<float> &mzTransitions,
            QVector<ScanPoints> *demultiplexedScans
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzTransitions); ree;
        e = ErrorUtils::isEqual(
                static_cast<int>(xMat.cols()),
                mzTransitions.size()
        ); ree;

        demultiplexedScans->clear();

        const int xMatRows = static_cast<int>(xMat.rows());

        demultiplexedScans->resize(xMatRows);
        demultiplexedScans->reserve(xMatRows);

        for (int row = 0; row < xMatRows; row++) {

            const Eigen::VectorX<float> &mzRow = xMat.row(row);
            for (int col = 0; col < xMat.cols(); col++) {

                const float intensity = mzRow.coeffRef(col);

                if (MathUtils::tZero(intensity)) {
                    continue;
                }

                (*demultiplexedScans)[row].push_back({mzTransitions.at(col), intensity});
            }
        }

        ERR_RETURN
    }

}//namespace
Err DeMultiplexScanerator::Private::deMultiplexScans(
        const QVector<ScanPoints*> &scanPoints,
        const QVector<MsScanInfo> &msScanInfos,
        QVector<ScanPoints> *demultiplexedScans,
        QVector<MzTargetKey> *mzTargetKeys,
        double *windowSize
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    e = ErrorUtils::isNotEmpty(msScanInfos); ree;

    Eigen::MatrixX<float> scanMaskMatrix;
    QPair<float, float> totalScanRange;

    e = buildScanMaskMatrix(
            msScanInfos,
            &scanMaskMatrix,
            &totalScanRange,
            windowSize
            ); ree;

    Eigen::MatrixX<float> transitionsMatrix;
    QVector<float> mzTransitions;
    e = buildTransitionsMatrix(scanPoints, &transitionsMatrix, &mzTransitions); ree;

    Eigen::MatrixX<float> X(scanMaskMatrix.cols(), transitionsMatrix.cols());
    X.setZero();

    const int maxIters = 10;
    const float tolerance = 10.0;
    Eigen::NNLS<Eigen::MatrixX<float>> nnls(
            scanMaskMatrix,
            maxIters,
            tolerance
            );

    for(int i = 0; i < transitionsMatrix.cols(); i++) {

        Eigen::VectorX<float> b = transitionsMatrix.col(i);
        Eigen::VectorX<float> x = nnls.solve(b);

        e = ErrorUtils::isTrue(nnls.info() == Eigen::Success); ree;
        X.col(i) = x;
    }

    e = buildDemultiplexedScans(
            X,
            mzTransitions,
            demultiplexedScans
            ); ree;

    const int distance = static_cast<int>(totalScanRange.second - totalScanRange.first);

    QVector<float> demuxWindows;
    for (int i = 0; i < distance; i++) {
        const double window = totalScanRange.first + (*windowSize * i) + (*windowSize / 2.0);
        demuxWindows.push_back(window);
    }

    std::transform(
            demuxWindows.begin(),
            demuxWindows.end(),
            std::back_inserter(*mzTargetKeys),
            [](float f){return QString::number(MathUtils::hashDecimal(f, S_GLOBAL_SETTINGS.HASHING_PRECISION));}
            );

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


DeMultiplexScanerator::DeMultiplexScanerator(
        double ppmTol,
        double intensityFractionThreshold
        )
        : d_ptr(new Private(ppmTol, intensityFractionThreshold))
        {}

DeMultiplexScanerator::~DeMultiplexScanerator() {}

Err DeMultiplexScanerator::deMultiplexScans(
        const QVector<ScanPoints*> &scans,
        const QVector<MsScanInfo> &msScanInfos,
        QVector<ScanPoints> *demultiplexedScans,
        QVector<MzTargetKey> *mzTargetKeys,
        double *windowSize
        ) {
    ERR_INIT
    e = d_ptr->deMultiplexScans(
            scans,
            msScanInfos,
            demultiplexedScans,
            mzTargetKeys,
            windowSize
            ); ree;
    ERR_RETURN
}

Err DeMultiplexScanerator::_buildTransitionMatrixTestAccess(
        const QVector<ScanPoints*> &scans,
        QVector<QVector<float>> *transitionMatrixVecs
        ) {
    ERR_INIT

    Eigen::MatrixX<float> transitionsMatrix;
    QVector<float> mzTransitions;

    e = d_ptr->buildTransitionsMatrix(
            scans,
            &transitionsMatrix,
            &mzTransitions
            ); ree;

    transitionsMatrix.conservativeResize(4, 4);
    *transitionMatrixVecs = EigenUtils::convertEigenMatrixToQVectors(transitionsMatrix);

    ERR_RETURN
}

Err DeMultiplexScanerator::_buildScanMaskMatrixTestAccess(
        const QVector<MsScanInfo> &msScanInfos,
        QVector<QVector<float>> *scanMaskMatrixVecs
        ) {
    ERR_INIT

    Eigen::MatrixX<float> scanMaskMatrix;
    QPair<float, float> totalScanRange;
    double windowSize;
    e = d_ptr->buildScanMaskMatrix(
            msScanInfos,
            &scanMaskMatrix,
            &totalScanRange,
            &windowSize
            ); ree;

    *scanMaskMatrixVecs = EigenUtils::convertEigenMatrixToQVectors(scanMaskMatrix);

    ERR_RETURN
}
