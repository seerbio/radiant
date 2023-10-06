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
            const QVector<double> &intensityVec,
            QVector<PeakIntegrationIndexes> *peakLimits,
            QVector<double> *intensityVecSmoothed
    );

private:

    PeakIntegratomaticParameters m_params;
    Eigen::VectorX<double> m_gaussFilter;
    Eigen::VectorX<double> m_mexicanHatFilter;

};


PeakIntegratomatic::Private::Private() {}
PeakIntegratomatic::Private::~Private() {}


Err PeakIntegratomatic::Private::init(const PeakIntegratomaticParameters &params) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid(), Error::eValueError); ree;
    m_params = params;

    m_gaussFilter = EigenKernelUtils::buildGaussianFilter1D(
            m_params.filterLength,
            m_params.sigma
    );

    m_mexicanHatFilter = EigenKernelUtils::buildMexicanHatFilter1D(
            m_params.filterLength,
            m_params.sigma * 1.5
    );

    ERR_RETURN
}


namespace {

    double maxOfVectorSegment(
            const QVector<double> &vec,
            int startIndex,
            int endIndex
    ) {

        const QVector<double> vecTrunc = vec.mid(startIndex, 1 + endIndex - startIndex);
        return *std::max_element(vecTrunc.begin(), vecTrunc.end());
    }

}//namespace
Err PeakIntegratomatic::Private::findAllPeaksLimitsInXIC(
        const QVector<double> &intensityVec,
        QVector<PeakIntegrationIndexes> *peakLimits,
        QVector<double> *intensityVecSmoothed
) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_gaussFilter.size() > 0); ree;
    e = ErrorUtils::isTrue(m_mexicanHatFilter.size() > 0); ree;

    peakLimits->clear();
    intensityVecSmoothed->clear();

    const double filterToVecRejectionRatio = 1.5;
    if (intensityVec.size() < static_cast<int>(std::round(m_gaussFilter.size() * filterToVecRejectionRatio))) {
        *peakLimits = {{0, intensityVec.size() -1}};
        ERR_RETURN
    }

    Eigen::VectorX<double> smoothedVec = EigenUtils::convertQVectorToEigenVector(intensityVec);

    for (int i = 0; i < m_params.smoothCount; i++) {
        smoothedVec = EigenKernelUtils::convolveVectorWithKernel(smoothedVec, m_gaussFilter);
    }

    const double maxVal = smoothedVec.maxCoeff();
    if (MathUtils::tZero(maxVal)) {
        *peakLimits = {{0, intensityVec.size() -1}};
        ERR_RETURN
    }

    smoothedVec /= maxVal;
    *intensityVecSmoothed = EigenUtils::convertEigenVectorToQVector<double>(smoothedVec);

    Eigen::VectorX<double> mexicanHatsmoothedVec
        = EigenKernelUtils::convolveVectorWithKernel(smoothedVec, m_mexicanHatFilter);

//#define PRINT_VECS
#ifdef PRINT_VECS
    qDebug() << intensityVec;
    qDebug() << *intensityVecSmoothed;
    qDebug() << EigenUtils::convertEigenVectorToQVector<double>(mexicanHatsmoothedVec);
#endif

    Eigen::VectorX<double> mexicanHatsmoothedVecApexes = mexicanHatsmoothedVec;
    EigenUtils::thresholdVector(0.0, &mexicanHatsmoothedVecApexes);
    const QMap<int, double> maxima = EigenUtils::apexes(mexicanHatsmoothedVecApexes);

    QMap<int, double> minima = EigenUtils::troughs(mexicanHatsmoothedVec);
    minima.insert(-1, 0);
    minima.insert(static_cast<int>(mexicanHatsmoothedVec.size()) - 1, 0);

    const QVector<int> maximaKeys = maxima.keys().toVector();
    const QVector<int> minimaKeys = minima.keys().toVector();

    int startIndex = 0;
    const double signalToNoiseThreshold
        = m_params.signalToNoiseRatio * MathUtils::median(smoothedVec);

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

            const double maxOfVecSegment = maxOfVectorSegment(*intensityVecSmoothed, startVal, endVal);

            if (maxOfVecSegment > signalToNoiseThreshold) {
                peakLimits->push_back({startVal, endVal});
            }

            startIndex = endIndex;
            break;
        }
    }

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
        const QVector<double> &intensityVec,
        QVector<PeakIntegrationIndexes> *peakLimits,
        QVector<double> *intensityVecSmoothed
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
            Eigen::VectorX<double> *vec
            ) {

        ERR_INIT

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
        const QVector<double> &vec,
        double stopThresholdFraction,
        int filterLength,
        int smoothCount,
        PeakIntegrationIndexes *peakIntegrationIndexes
        ) {

    ERR_INIT

    QVector<double> smoothedVecNoReturn;
    e = simpleIntegrator(
            vec,
            stopThresholdFraction,
            filterLength,
            smoothCount,
            peakIntegrationIndexes,
            &smoothedVecNoReturn
            ); ree;

    ERR_RETURN
}


Err PeakIntegratomatic::simpleIntegrator(
        const QVector<double> &vec,
        double stopThresholdFraction,
        int filterLength,
        int smoothCount,
        PeakIntegrationIndexes *peakIntegrationIndexes,
        QVector<double> *smoothedVec
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(vec); ree;

    const int minFilterLength = 3;
    e = ErrorUtils::isAboveThreshold(
            filterLength,
            minFilterLength,
            ErrorUtilsParam::IncludeThreshold
    ); ree;

    const double minStopThresholdFraction = 0.005;
    e = ErrorUtils::isAboveThreshold(
            stopThresholdFraction,
            minStopThresholdFraction,
            ErrorUtilsParam::IncludeThreshold
    ); ree;

    Eigen::VectorX<double> eVec = EigenUtils::convertQVectorToEigenVector(vec);
    e = smoothVector(
            filterLength,
            smoothCount,
            &eVec
    ); ree;

    const QVector<QPair<int, double>> vecApex = EigenUtils::returnTopXIndexAndValues(eVec, 1);
    const int apexIndex = vecApex.front().first;
    const double apexValue = vecApex.front().second;
    const double stopThreshold = apexValue * stopThresholdFraction;

    double rightStopVal = apexValue;
    int rightStopIndex = apexIndex;

    int rightCurrentIndex = apexIndex;
    while (rightCurrentIndex < eVec.size()) {

        const double currentValue = eVec(rightCurrentIndex);
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

    double leftStopVal = apexValue;
    int leftStopIndex = apexIndex;

    int leftCurrentIndex = apexIndex;
    while (leftCurrentIndex < eVec.size()) {

        const double currentValue = eVec(leftCurrentIndex);
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

    *smoothedVec = EigenUtils::convertEigenVectorToQVector(eVec);
    *peakIntegrationIndexes = {leftStopIndex, rightStopIndex};

    ERR_RETURN
}