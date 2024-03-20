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
    int smoothCount = 1;
    double sigma = 1.0;
    double signalToNoiseRatio = 3.0;
    double stopThresholdFraction = 0.1;

    /**
    * @brief Checks the validity of filter parameters.
    *
    * This member function checks the validity of filter parameters. It ensures that filterLength is greater than or equal to 2,
    * smoothCount is greater than or equal to 0, sigma is greater than 0, and signalToNoiseRatio is greater than 0.
    *
    * @return A boolean indicating the validity of filter parameters.
    */
    [[nodiscard]] bool isValid() const {
        return filterLength >= 2
            && smoothCount >= 0
            && sigma > 0
            && signalToNoiseRatio > 0;
    }

    void printParams() const {
        qDebug() << "filterLength" << filterLength;
        qDebug() << "smoothCount" << smoothCount;
        qDebug() << "sigma" << sigma;
        qDebug() << "signalToNoiseRatio" << signalToNoiseRatio;
        qDebug() << "stopThresholdFraction" << stopThresholdFraction;
    }

    /**
    * @brief Builds PeakIntegratomaticParameters from PythiaParameters.
    *
    * This static function constructs PeakIntegratomaticParameters by copying relevant parameters from the provided PythiaParameters.
    *
    * @param pythiaParameters The PythiaParameters instance from which to extract parameters.
    * @return The constructed PeakIntegratomaticParameters.
    */
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

    /**
    * @brief Initializes the Private instance with PeakIntegratomaticParameters.
    *
    * This function initializes the Private instance of PeakIntegratomatic with the provided PeakIntegratomaticParameters.
    * It also checks the validity of the parameters and builds the necessary filters.
    *
    * @param params The PeakIntegratomaticParameters to use for initialization.
    * @return An error code indicating the success or failure of the initialization.
    */
    Err init(const PeakIntegratomaticParameters &params);

    /**
    * @brief Performs simple integration on the given vector and identifies peak integration regions.
    *
    * This function performs simple integration on the input vector and identifies peak integration regions.
    * It returns the peak integration indexes along with their corresponding intensity values.
    *
    * @param vec The input vector for integration.
    * @param peakIntegrationIndexesVsIntensity A pointer to the QVector to store the peak integration indexes along with intensity values.
    * @return An error code indicating the success or failure of the integration.
    */
    Err simpleIntegrator(
            const QVector<float> &vec,
            int topNApexes,
            int maxPeakIntegrations,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
    );

private:

    Q_DISABLE_COPY(PeakIntegratomatic) class Private;
    const QScopedPointer<Private> d_ptr;


};


#endif //PYTHIACPP_PEAKINTEGRATOMATIC_H
