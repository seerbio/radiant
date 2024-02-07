//
// Created by anichols on 2/5/24.
//

#include "DeMultiplexScanerator.h"

#include "ErrorUtils.h"
#include "MathUtils.h"

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
            Eigen::MatrixX<double> *scanMaskMatrix
            );

    Err buildTransitionsMatrix(
            const QVector<ScanPoints*> &scans,
            Eigen::MatrixX<double> *transitionsMatrix
            );

    Err deMultiplexScans(
            const QVector<ScanPoints*> &scanPoints,
            const QVector<MsScanInfo> &msScanInfos
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
        Eigen::MatrixX<double> *transitionsMatrix
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
        const ScanPointsSorter &spsR = scanPointsSorters.at(i + 1);

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

    transitions.resize(100);

    transitionsMatrix->resize(scans.size(), transitions.size());
    transitionsMatrix->setZero();

    for (int col = 0; col < transitions.size(); col++) {

        const QVector<ScanPointsSorter> &sps = transitions.at(col);

//        const double transitionMean = std::accumulate(
//                sps.begin(),
//                sps.end(),
//                0.0,
//                [](double sum, const ScanPointsSorter &sps){return sum + sps.scanPoint->x();}
//                ) / sps.size();
//        qDebug() << transitionMean;

        for (const ScanPointsSorter &sp : sps) {
            transitionsMatrix->coeffRef(sp.scansIndex, col) = sp.scanPoint->y();
        }
    }

    ERR_RETURN
}

namespace {

    Err getTotalScanRange(
            const QVector<MsScanInfo> &msScanInfos,
            QPair<double, double> *totalScanRange
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(msScanInfos); ree;

        *totalScanRange = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};

        for (const MsScanInfo &msScanInfo : msScanInfos) {
            const double mzLo = msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower;
            totalScanRange->first = std::min(totalScanRange->first, mzLo);

            const double mzHi = msScanInfo.precursorTargetMz + msScanInfo.isoWindowLower;
            totalScanRange->second = std::max(totalScanRange->second, mzHi);
        }

        ERR_RETURN
    }

    Err calculateWindowStaggering(const QVector<MsScanInfo> &msScanInfos, double *amuStagger) {

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

        *amuStagger = meanDistance;

        ERR_RETURN
    }

}//namespace
Err DeMultiplexScanerator::Private::buildScanMaskMatrix(
        const QVector<MsScanInfo> &msScanInfos,
        Eigen::MatrixX<double> *scanMaskMatrix
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msScanInfos); ree;

    QPair<double, double> totalScanRange;
    e = getTotalScanRange(msScanInfos, &totalScanRange); ree;
    e = ErrorUtils::isBelowThreshold(
            totalScanRange.first,
            totalScanRange.second,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;

    double amuStagger;
    e = calculateWindowStaggering(msScanInfos, &amuStagger); ree;

    qDebug() << totalScanRange << amuStagger;

    const int rows = msScanInfos.size();
    const int cols = static_cast<int>(totalScanRange.second - totalScanRange.first);
    scanMaskMatrix->resize(rows, cols);
    scanMaskMatrix->setZero();


    for (int i = 0; i < msScanInfos.size(); i++) {

        const MsScanInfo &msi = msScanInfos.at(i);

        const int mzStart
            = static_cast<int>((msi.precursorTargetMz - msi.isoWindowLower) - totalScanRange.first);

        const int mzEnd
            = static_cast<int>((msi.precursorTargetMz + msi.isoWindowUpper) - totalScanRange.first);

        for (int j = mzStart; j < mzEnd; j++) {
            scanMaskMatrix->coeffRef(i, j) = 1.0;
        }
    }

    ERR_RETURN
}

Err DeMultiplexScanerator::Private::deMultiplexScans(
        const QVector<ScanPoints*> &scanPoints,
        const QVector<MsScanInfo> &msScanInfos
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    e = ErrorUtils::isNotEmpty(msScanInfos); ree;

    Eigen::MatrixX<double> scanMaskMatrix;
    e = buildScanMaskMatrix(msScanInfos, &scanMaskMatrix); ree;

    Eigen::MatrixX<double> transitionsMatrix;
    e = buildTransitionsMatrix(scanPoints, &transitionsMatrix); ree;

    transitionsMatrix.coeffRef(1, 0) = 5295.0;

    std::cout << transitionsMatrix << std::endl;
    std::cout << " " << std::endl;
    std::cout << scanMaskMatrix << std::endl;

    Eigen::MatrixX<double> X(scanMaskMatrix.cols(), transitionsMatrix.cols());
    X.setZero();

    Eigen::NNLS<Eigen::MatrixX<double>> nnls(scanMaskMatrix, 10, 10.0);

    for(int i = 0; i < transitionsMatrix.cols(); i++) {

        Eigen::VectorX<double> b = transitionsMatrix.col(i);
        Eigen::VectorX<double> x = nnls.solve(b);

        e = ErrorUtils::isTrue(nnls.info() == Eigen::Success); ree;
        X.col(i) = x;
    }


    std::cout << scanMaskMatrix << std::endl;
    std::cout << "" << std::endl;
    std::cout << X << std::endl;
    std::cout << "" << std::endl;
    std::cout << transitionsMatrix << std::endl;
    std::cout << "" << std::endl;
    std::cout << scanMaskMatrix * X << std::endl;

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
        const QVector<MsScanInfo> &msScanInfos
        ) {
    ERR_INIT
    e = d_ptr->deMultiplexScans(scans, msScanInfos);
    ERR_RETURN
}

Err DeMultiplexScanerator::_buildTransitionMatrixTestAccess(
        const QVector<ScanPoints*> &scans
        ) {
    ERR_INIT
    Eigen::MatrixX<double> transitionsMatrix;
    e = d_ptr->buildTransitionsMatrix(scans, &transitionsMatrix); ree;
    ERR_RETURN
}

Err DeMultiplexScanerator::_buildScanMaskMatrixTestAccess(const QVector<MsScanInfo> &msScanInfos) {
    ERR_INIT
    Eigen::MatrixX<double> scanMaskMatrix;
    e = d_ptr->buildScanMaskMatrix(msScanInfos, &scanMaskMatrix); ree;
    ERR_RETURN
}
