//
// Created by anichols on 12/13/22.
//

#include "PeakIntegratomatic.h"

#include "ErrorUtils.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"

#include <QDebug>




PeakIntegratomatic::PeakIntegratomatic() = default;

PeakIntegratomatic::~PeakIntegratomatic() = default;


Err PeakIntegratomatic::init(const PeakIntegratomaticParameters &params) {

    ERR_INIT

	e = ErrorUtils::isTrue(params.isValid(), eValueError);
	m_params = params;

    ERR_RETURN
}

Err PeakIntegratomatic::simpleIntegrator(
	const QVector<int> &apexes,
	const float* vec,
	int vecSize,
	QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_params.isValid()); ree;

	int lastGlobalStopLeft = std::numeric_limits<int>::max();
	int lastGlobalStopRight = 0;

	for (int apex : apexes) {
		const float stopThresold = vec[apex] * m_params.stopThresholdFraction;

		int lastLeftStop = apex;
		float lastLeftValue = vec[apex];

		while (lastLeftStop > 0) {

			const int currentLeftStop = lastLeftStop - 1;
			const float currentLeftValue = vec[currentLeftStop];

			if (currentLeftValue <= stopThresold
				|| currentLeftValue > lastLeftValue * (1.0 + m_params.hysteresis)
				|| lastLeftStop - 1 == lastGlobalStopRight
				) {
				lastLeftStop = currentLeftStop;
				lastGlobalStopLeft = lastLeftStop;
				break;
			}

			lastLeftStop = currentLeftStop;
			lastLeftValue = currentLeftValue;
			lastGlobalStopLeft = lastLeftStop;
		}

		int lastRightStop = apex;
		float lastRightValue = vec[apex];
		while (lastRightStop < vecSize) {

			const int currentRightStop = lastRightStop + 1;
			const float currentRightValue = vec[currentRightStop];

			if (currentRightValue <= stopThresold
				|| currentRightValue > lastRightValue * (1.0 + m_params.hysteresis)
				|| lastRightStop + 1 == lastGlobalStopLeft
				) {
				lastRightStop = currentRightStop;
				lastGlobalStopRight = lastRightStop;
				break;
				}

			lastRightStop = currentRightStop;
			lastRightValue = currentRightValue;
			lastGlobalStopRight = lastRightStop;
		}

		peakIntegrationIndexesVsIntensity->push_back({{lastLeftStop, lastRightStop}, vec[apex]});
	}

    ERR_RETURN
}
