//
// Created by anichols on 12/13/22.
//

#include "PeakIntegratomatic.h"

#include "ErrorUtils.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"

#include <QDebug>

class Q_DECL_HIDDEN PeakIntegratomatic::Private
{

public:

    Private();
    ~Private();

    Err init(const PeakIntegratomaticParameters &params);

    Err findAllPeaksLimitsInXIC(
            const QVector<float> &intensityVec,
            QVector<PeakIntegrationIndexes> *peakLimits,
            QVector<float> *intensityVecSmoothed
    );

private:

    PeakIntegratomaticParameters m_params;
    Eigen::VectorX<float> m_gaussFilter;
    Eigen::VectorX<float> m_mexicanHatFilter;

};


PeakIntegratomatic::Private::Private() {}
PeakIntegratomatic::Private::~Private() {}


Err PeakIntegratomatic::Private::init(const PeakIntegratomaticParameters &params) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid(), Error::eValueError); ree;
    m_params = params;

    m_gaussFilter = EigenKernelUtils::buildGaussianFilter1D<float>(
            m_params.filterLength,
            m_params.sigma
    );

    m_mexicanHatFilter = EigenKernelUtils::buildMexicanHatFilter1D<float>(
            m_params.filterLength,
            m_params.sigma * 1.5
    );

    ERR_RETURN
}


namespace {

    double maxOfVectorSegment(
            const QVector<float> &vec,
            int startIndex,
            int endIndex
    ) {

        const QVector<float> vecTrunc = vec.mid(startIndex, 1 + endIndex - startIndex);
        return *std::max_element(vecTrunc.begin(), vecTrunc.end());
    }

    //TODO add this to eigen kernel utils
    Eigen::VectorX<float> addPaddingToVector(
            const Eigen::VectorX<float> &vec,
            int paddingCount
            ) {

        const int finalPaddingSize = vec.size() > paddingCount ? paddingCount : 2 * paddingCount;

        Eigen::VectorX<float> vecPadded(vec.size() + finalPaddingSize);
        vecPadded.setZero();

        const int vecInsertPoint = vec.size() > paddingCount ? 0 : paddingCount;
        vecPadded.segment(vecInsertPoint, vec.size()) = vec;

        return vecPadded;
    }

}//namespace
Err PeakIntegratomatic::Private::findAllPeaksLimitsInXIC(
        const QVector<float> &intensityVec,
        QVector<PeakIntegrationIndexes> *peakLimits,
        QVector<float> *intensityVecSmoothed
) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_gaussFilter.size() > 0); ree;
    e = ErrorUtils::isTrue(m_mexicanHatFilter.size() > 0); ree;

    peakLimits->clear();
    intensityVecSmoothed->clear();

    const Eigen::VectorX<float> smoothedVec = EigenUtils::convertQVectorToEigenVector(intensityVec);

    const int paddingMultiplier = 1;
    const int paddingSize = static_cast<int>(m_gaussFilter.size()) * paddingMultiplier;

    Eigen::VectorX<float> smoothedVecPadded = addPaddingToVector(
            smoothedVec,
            paddingSize
            );

    for (int i = 0; i < m_params.smoothCount; i++) {
        smoothedVecPadded = EigenKernelUtils::convolveVectorWithKernel(smoothedVecPadded, m_gaussFilter);
    }

    const double maxVal = smoothedVec.maxCoeff();
    if (MathUtils::tZero(maxVal)) {
        *peakLimits = {{0, intensityVec.size() -1}};
        ERR_RETURN
    }

    smoothedVecPadded /= maxVal;

    *intensityVecSmoothed = EigenUtils::convertEigenVectorToQVector<float>(smoothedVecPadded);

    Eigen::VectorX<float> mexicanHatsmoothedVec
        = EigenKernelUtils::convolveVectorWithKernel(smoothedVecPadded, m_mexicanHatFilter);

//#define PRINT_VECS
#ifdef PRINT_VECS
    qDebug() << intensityVec;
    qDebug() << *intensityVecSmoothed;
    qDebug() << EigenUtils::convertEigenVectorToQVector<double>(mexicanHatsmoothedVec);
#endif

    Eigen::VectorX<float> mexicanHatsmoothedVecApexes = mexicanHatsmoothedVec;
    EigenUtils::thresholdVector(static_cast<float>(0.0), &mexicanHatsmoothedVecApexes);
    const QMap<int, float> maxima = EigenUtils::apexes(mexicanHatsmoothedVecApexes);

    QMap<int, float> minima = EigenUtils::troughs(mexicanHatsmoothedVec);
    minima.insert(-1, 0);
    minima.insert(static_cast<int>(mexicanHatsmoothedVec.size()) - 1, 0);

    const QVector<int> maximaKeys = maxima.keys().toVector();
    const QVector<int> minimaKeys = minima.keys().toVector();

    int startIndex = 0;
    const double signalToNoiseThreshold
        = m_params.signalToNoiseRatio * MathUtils::median(smoothedVecPadded);

    for (int mx : maximaKeys) {

        if (mx <= 1) {
            continue;
        }

        for (int endIndex = startIndex; endIndex < minimaKeys.size(); endIndex++) {

            if (endIndex > minimaKeys.size()) {
                break;
            }

            const int startVal = minimaKeys.at(startIndex);
            const int endVal = minimaKeys.at(endIndex);

            if (!(startVal < mx && mx < endVal)) {
                startIndex = endIndex;
                continue;
            }

            const float maxOfVecSegment = maxOfVectorSegment(*intensityVecSmoothed, startVal, endVal);

            if (maxOfVecSegment > signalToNoiseThreshold) {

                const bool sizeCheck = intensityVec.size() > paddingSize;

                const int startValCorrected = sizeCheck ? startVal : startVal - paddingSize;
                const int endValCorrected = sizeCheck ? endVal : endVal - paddingSize;

                peakLimits->push_back({startValCorrected, endValCorrected});
            }

            startIndex = endIndex;
            break;
        }
    }

    smoothedVecPadded /= smoothedVecPadded.maxCoeff();
    smoothedVecPadded *= maxVal;
    *intensityVecSmoothed = EigenUtils::convertEigenVectorToQVector<float>(smoothedVecPadded);

    ERR_RETURN
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


PeakIntegratomatic::PeakIntegratomatic() : d_ptr(new Private()) {}

PeakIntegratomatic::~PeakIntegratomatic(){}


Err PeakIntegratomatic::init(const PeakIntegratomaticParameters &params) {

    ERR_INIT
    e = d_ptr->init(params); ree;
    ERR_RETURN
}


Err PeakIntegratomatic::findAllPeaksLimitsInXIC(
        const QVector<float> &intensityVec,
        QVector<PeakIntegrationIndexes> *peakLimits,
        QVector<float> *intensityVecSmoothed
){

    ERR_INIT

    e = d_ptr->findAllPeaksLimitsInXIC(
            intensityVec,
            peakLimits,
            intensityVecSmoothed
    ); ree;

    ERR_RETURN
}

namespace {

    Err smoothVector(
            int filterLength,
            int smoothCount,
            Eigen::VectorX<float> *vec
            ) {

        ERR_INIT

        const int filterLengthMin = 5;
        filterLength = std::max(filterLength, filterLengthMin);

        const int order = 2;
        const int derivative = 0;
        const int rate = 1;

        for (int i = 0; i < smoothCount; i++) {
            e = EigenKernelUtils::savitskyGolaySmooth(
                    filterLength,
                    order,
                    derivative,
                    rate,
                    vec
                    ); ree;
        }

        ERR_RETURN
    }


}//namespace
Err PeakIntegratomatic::simpleIntegrator(
        const QVector<float> &vec,
        float stopThresholdFraction,
        int filterLength,
        int smoothCount,
        QVector<PeakIntegrationIndexes> *peakIntegrationIndexes,
        QVector<float> *smoothedVec
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(vec); ree;

    const int minFilterLength = 3;
    e = ErrorUtils::isAboveThreshold(
            filterLength,
            minFilterLength,
            ErrorUtilsParam::IncludeThreshold
    ); ree;

    const float minStopThresholdFraction = 0.005;
    e = ErrorUtils::isAboveThreshold(
            stopThresholdFraction,
            minStopThresholdFraction,
            ErrorUtilsParam::IncludeThreshold
    ); ree;

    Eigen::VectorX<float> eVec = EigenUtils::convertQVectorToEigenVector(vec);
    e = smoothVector(
            filterLength,
            smoothCount,
            &eVec
    ); ree;

    EigenUtils::thresholdVector(static_cast<float>(1.01), &eVec);

    const int topNApexes = 10;
    const QMap<int, float> vecApexs = EigenUtils::apexes(eVec);

    if (vecApexs.isEmpty()) {
        ERR_RETURN
    }

    Eigen::VectorX<float> apexes =EigenUtils::convertQMapToEigenVector(vecApexs, vecApexs.lastKey() + 1);
    QVector<QPair<int, float>> apexPaairs = EigenUtils::returnTopXIndexAndValues(apexes, topNApexes);

    if (eVec.maxCoeff() < 5) {
        ERR_RETURN
    }

    for (auto it : apexPaairs) {
        qDebug() << it.first << it.second;
    }
    qDebug() << vec;
    qDebug() << EigenUtils::convertEigenVectorToQVector(eVec);
    qDebug() << "**********" << vecApexs.size();

    for (auto it = vecApexs.begin(); it != vecApexs.end(); it++) {
        const int apexIndex = it.key();
        const float apexValue = it.value();
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

        peakIntegrationIndexes->push_back({leftStopIndex, rightStopIndex});
    }

    *smoothedVec = EigenUtils::convertEigenVectorToQVector(eVec);

    ERR_RETURN
}