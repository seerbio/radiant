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
    float stopThresholdFraction = 0.5;
	float hysteresis = 0.1;

	[[nodiscard]] bool isValid() const {
		return (stopThresholdFraction >= 0.0 && stopThresholdFraction <= 1.0)
				&& (hysteresis >= 0.0 && hysteresis <= 1.0);
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
            const float* vec,
            int vecSize,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
            ) const;

private:

	PeakIntegratomaticParameters m_params;

};


#endif //PYTHIACPP_PEAKINTEGRATOMATIC_H
