//
// Created by anichols on 5/9/24.
//

#include "XICPeakManager.h"

#include "EigenKernelUtils.h"
#include "ErrorUtils.h"
#include "MsFrame.h"
#include "PeakIntegratomatic.h"


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
    float stopThresholdFraction
    ) {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(filterLength, 3, ErrorUtilsParam::IncludeThreshold); ree
    e = ErrorUtils::isAboveThreshold(polynomialOrder, 1, ErrorUtilsParam::IncludeThreshold); ree
    e = ErrorUtils::isAboveThreshold(smoothCount, 1, ErrorUtilsParam::IncludeThreshold); ree
    e = ErrorUtils::isAboveThreshold(stopThresholdFraction, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree

    m_filterLength = filterLength;
    m_polynomialOrder = polynomialOrder;
    m_smoothCount = smoothCount;
    m_stopThresholdValue = stopThresholdFraction;

    m_mzHashedVsXicPoints.clear();

    ERR_RETURN
}

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

    e = extractXICs(
        mzValsToExtract,
        ppmTolerance,
        &turboXic
        ); ree;

    ERR_RETURN
}

Err XICPeakManager:: getXIC(
    float mzVal,
    XICPoints *xicPoints
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_mzHashedVsXicPoints); ree;
    e = ErrorUtils::isAboveThreshold(mzVal, 0.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    const MzHashed mzHashed = MathUtils::hashDecimal(mzVal, S_GLOBAL_SETTINGS.HASHING_PRECISION);

    *xicPoints = m_mzHashedVsXicPoints.value(mzHashed);
    ERR_RETURN
}

 Err XICPeakManager::cacheXICPeakManager(const QString &outputFilePath) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_mzHashedVsXicPoints); ree;

    QFile file(outputFilePath);

    if(file.exists()) {
        const bool removed = file.remove();
        e = ErrorUtils::isTrue(removed);
    }

    if (!file.open(QIODevice::WriteOnly)) {
        qDebug() << "Failed to open file!";
        rrr(eFileError);
    }

    QDataStream out(&file);
    out.setVersion(QDataStream::Qt_5_15);

    const QList<XICPoints> &xicPointses = m_mzHashedVsXicPoints.values();

    QVector<QVector<XICPoint>> vec;
    std::transform (
        xicPointses.begin(),
        xicPointses.end(),
        std::back_inserter(vec),
        [](const XICPoints &p){return QVector<XICPoint>(p.begin(), p.end());}
        );

    out << m_filterLength;
    out << m_polynomialOrder;
    out << m_smoothCount;
    out << m_stopThresholdValue;
    out << m_mzHashedVsXicPoints.keys().toVector();
    out << vec;

    file.close();

    ERR_RETURN
}

Err XICPeakManager::loadXICPeakManagerCache(const QString& outputFilePath) {

    ERR_INIT
    e = ErrorUtils::fileExists(outputFilePath); ree;

    QFile file(outputFilePath);
    if (!file.open(QIODevice::ReadOnly)) {
        qDebug() << "Failed to open file!";
        rrr(eFileError);
    }

    QDataStream in(&file);
    in.setVersion(QDataStream::Qt_5_15);

    QVector<MzHashed> mzHashedVec;
    QVector<QVector<XICPoint>> vec;

    in >> m_filterLength;
    in >> m_polynomialOrder;
    in >> m_smoothCount;
    in >> m_stopThresholdValue;
    in >> mzHashedVec;
    in >> vec;

    e = ErrorUtils::isEqual(mzHashedVec.size(), vec.size()); ree;
    for (int i = 0; i < mzHashedVec.size(); i++) {
        const QVector<XICPoint> &v = vec.at(i);
        m_mzHashedVsXicPoints.insert(mzHashedVec.at(i), std::vector<XICPoint>(v.begin(), v.end()));
    }

    file.close();

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

    Err simpleIntegrator(
        const Eigen::VectorX<float> &vec,
        float stopThresholdFraction,
        QVector<QPair<PeakIntegrationIndexes, Intensity>> *peakIntegrationIndexesVsIntensity
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(vec.size() > 0); ree;

    Eigen::VectorX<float> eVec = vec;
    EigenUtils::thresholdVector(static_cast<float>(1.01), &eVec);

    const QMap<int, float> vecApexs = EigenUtils::apexes(eVec);
    if (vecApexs.isEmpty()) {
        ERR_RETURN
    }

    Eigen::VectorX<float> apexes =EigenUtils::convertQMapToEigenVector(vecApexs, vecApexs.lastKey() + 1);
    QVector<QPair<int, float>> apexPairs = EigenUtils::returnTopXIndexAndValues(apexes, vecApexs.size());

    for (const QPair<int, float> &pr : apexPairs) {

        const int apexIndex = pr.first;
        const float apexValue = pr.second;

        if (MathUtils::tZero(apexValue)) {
            continue;
        }

        const float stopThreshold = apexValue * stopThresholdFraction;

        float rightStopVal = apexValue;
        int rightStopIndex = apexIndex;

        int rightCurrentIndex = apexIndex;
        while (rightCurrentIndex < eVec.size()) {

            const float currentValue = eVec(rightCurrentIndex);
            if (currentValue < stopThreshold) {
                rightStopIndex = rightCurrentIndex;
                break;
            }

            if (currentValue <= rightStopVal) {
                rightStopVal = currentValue;
                rightStopIndex = rightCurrentIndex;
                rightCurrentIndex++;
                continue;
            }

            break;
        }

        float leftStopVal = apexValue;
        int leftStopIndex = apexIndex;

        int leftCurrentIndex = apexIndex;
        while (leftCurrentIndex < eVec.size()) {

            const float currentValue = eVec(leftCurrentIndex);
            if (currentValue < stopThreshold) {
                leftStopIndex = leftCurrentIndex;
                break;
            }

            if (currentValue <= leftStopVal) {
                leftStopVal = currentValue;
                leftStopIndex = leftCurrentIndex;
                leftCurrentIndex--;
                continue;
            }

            break;
        }

        peakIntegrationIndexesVsIntensity->push_back({
             {std::max(leftStopIndex, 0), std::min(rightStopIndex, static_cast<int>(vec.size() - 1))},
             apexValue}
        );

        for (int i = leftStopIndex; i <= rightStopIndex; i++) {
            eVec.coeffRef(i) = 0.0;
        }
    }

    ERR_RETURN
}

}//namespace
Err XICPeakManager::extractXICs(
    const QVector<float> &mzValsToExtract,
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

        if (xicPoints.empty()) {
            continue;
        }

        // Eigen::VectorX<float> xicPointsVec = loadXicPointsToVectorX(xicPoints);
        // for (int smoothIter = 0; smoothIter < m_smoothCount; smoothIter++) {
        //     xicPointsVec = EigenKernelUtils::convolveVectorWithKernel(xicPointsVec, savitskyGolayKernel);
        // }

        // QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationIndexesVsIntensity;
        // e = simpleIntegrator(
        //     xicPointsVec,
        //     m_stopThresholdValue,
        //     &peakIntegrationIndexesVsIntensity
        //     ); ree;

        const MzHashed mzHashed = MathUtils::hashDecimal(mzVal, S_GLOBAL_SETTINGS.HASHING_PRECISION);
        m_mzHashedVsXicPoints.insert(mzHashed , xicPoints);

    }

    ERR_RETURN
}
