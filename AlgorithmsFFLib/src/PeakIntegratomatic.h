//
// Created by anichols on 12/13/22.
//

#ifndef PYTHIACPP_PEAKINTEGRATOMATIC_H
#define PYTHIACPP_PEAKINTEGRATOMATIC_H


#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "PythiaParameterReader.h"


using namespace Error;


struct ALGORITHMSFFLIB_EXPORTS PeakIntegratomaticParameters {

    int filterLength = 2;
    int smoothCount = 3;
    double sigma = 1.0;
    double signalToNoiseRatio = 3.0;
    double stopThresholdFraction = 0.1;

    [[nodiscard]] bool isValid() const {
        return filterLength >= 2
            && smoothCount >= 0
            && sigma > 0
            && signalToNoiseRatio > 0;
    }

    void printParams() {
        qDebug() << "filterLength" << filterLength;
        qDebug() << "smoothCount" << smoothCount;
        qDebug() << "sigma" << sigma;
        qDebug() << "signalToNoiseRatio" << signalToNoiseRatio;
        qDebug() << "stopThresholdFraction" << stopThresholdFraction;
    }

    static PeakIntegratomaticParameters buildPeakIntegratomaticParams(const PythiaParameters  &pythiaParameters) {

        PeakIntegratomaticParameters params;
        params.smoothCount = pythiaParameters.smoothCount;
        params.filterLength = pythiaParameters.filterLength;
        params.sigma = pythiaParameters.sigma;
        params.signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;
        params.stopThresholdFraction = pythiaParameters.stopThresholdFraction;

        return params;
    }

};


class ALGORITHMSFFLIB_EXPORTS PeakIntegratomatic {

public:

    PeakIntegratomatic();
    ~PeakIntegratomatic();

    Err init(const PeakIntegratomaticParameters &params);

    Err simpleIntegrator(
            const QVector<float> &vec,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
    );

private:

    Q_DISABLE_COPY(PeakIntegratomatic) class Private;
    const QScopedPointer<Private> d_ptr;


};


#endif //PYTHIACPP_PEAKINTEGRATOMATIC_H
