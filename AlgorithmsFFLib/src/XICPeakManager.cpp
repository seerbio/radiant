//
// Created by anichols on 5/9/24.
//

#include "XICPeakManager.h"

#include "EigenKernelUtils.h"
#include "ErrorUtils.h"
#include "MsFrame.h"
#include "PeakIntegratomatic.h"
#include "TurboXIC.h"


XICPeakManager::XICPeakManager()
: m_filterLength(-1)
, m_polynomialOrder(-1)
, m_smoothCount(-1)
, m_stopThresholdValue(-1.0)
{}

Err XICPeakManager::init(
    int filterLength,
    int polynomialOrder,
    int smoothCount,
    double stopThresholdFraction
    ) {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(filterLength, 3, ErrorUtilsParam::IncludeThreshold); ree
    e = ErrorUtils::isAboveThreshold(polynomialOrder, 1, ErrorUtilsParam::IncludeThreshold); ree
    e = ErrorUtils::isAboveThreshold(smoothCount, 1, ErrorUtilsParam::IncludeThreshold); ree
    e = ErrorUtils::isAboveThreshold(stopThresholdFraction, 0.0, ErrorUtilsParam::ExcludeThreshold); ree

    m_filterLength = filterLength;
    m_polynomialOrder = polynomialOrder;
    m_smoothCount = smoothCount;
    m_stopThresholdValue = stopThresholdFraction;

    ERR_RETURN
}

namespace {

    PeakIntegratomaticParameters buildPeakIntegratomaticParameters(double stopThresholdFraction) {
        PeakIntegratomaticParameters params;
        params.sigma = 1.0;
        params.filterLength = 3;
        params.smoothCount = 0;
        params.signalToNoiseRatio = 1.0;
        params.stopThresholdFraction = stopThresholdFraction;

        return params;
    }

}//namespace
Err XICPeakManager::findPeaks(
    const MsFrame &msFrame,
    const QVector<float> &mzValsToExtract,
    float ppmTolerance
    ){

    ERR_INIT

    e = ErrorUtils::isTrue(msFrame.isValid()); ree;
    e = ErrorUtils::isNotEmpty(mzValsToExtract); ree;
    e = ErrorUtils::isAboveThreshold(ppmTolerance, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    TurboXIC turboXic;
    e = turboXic.init(msFrame.frameIndexVsScanPoints()); ree;

    PeakIntegratomatic peakIntegratomatic;
    e = peakIntegratomatic.init(buildPeakIntegratomaticParameters(m_stopThresholdValue)); ree;

    e = extractPeaks(
        mzValsToExtract,
        peakIntegratomatic,
        ppmTolerance,
        &turboXic
        ); ree;

    ERR_RETURN
}

namespace {

    Err buildSavGolayFilter(
        int filterLength,
        int polynomialOrder,
        Eigen::VectorX<float> *savitskyGolayKernel
        ) {

        ERR_INIT

        e = ErrorUtils::isAboveThreshold(filterLength, 3, ErrorUtilsParam::IncludeThreshold); ree
        e = ErrorUtils::isAboveThreshold(polynomialOrder, 1, ErrorUtilsParam::IncludeThreshold); ree

        Eigen::MatrixX<float> savitskyGolayKernelMat;
        e = EigenKernelUtils::buildSavitzkyGolayKernel(
            filterLength,
            polynomialOrder,
            0,
            1,
            &savitskyGolayKernelMat
            ); ree;

        *savitskyGolayKernel = savitskyGolayKernelMat;

        ERR_RETURN
    }

    Eigen::VectorX<float> loadXicPointsToVectorX(const XICPoints &xicPoints) {

        const FrameIndex frameIndexMax = std::max_element(
            xicPoints.begin(),
            xicPoints.end(),
            [](const XICPoint &l, const XICPoint &r){return l.scanNumber < r.scanNumber;}
            )->scanNumber + 1;

        Eigen::VectorX<float> vec(frameIndexMax);
        vec.setZero();

        for (const XICPoint &xicPoint : xicPoints) {
            vec.coeffRef(xicPoint.scanNumber) = xicPoint.intensity;
        }

        return vec;
    }

}//namespace
Err XICPeakManager::extractPeaks(
    const QVector<float> &mzValsToExtract,
    const PeakIntegratomatic &peakIntegratomatic,
    float ppmTolerance,
    TurboXIC *turboXic
) {

    ERR_INIT

    Eigen::VectorX<float> savitskyGolayKernel;
    e = buildSavGolayFilter(
        m_filterLength,
        m_polynomialOrder,
        &savitskyGolayKernel
        ); ree;

    for (float mzVal : mzValsToExtract) {

        const float massTol = MathUtils::calculatePPM(mzVal, ppmTolerance);
        const float mzMin = mzVal - massTol;
        const float mzMax = mzVal + massTol;

        const XICPoints xicPoints = turboXic->extractPointsXIC(mzMin, mzMax);

    }

    ERR_RETURN
}
