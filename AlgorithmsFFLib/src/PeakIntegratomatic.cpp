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
	QVector<std::tuple<PeakIntegrationIndexes, Intensity, FrameIndex>> *peakIntegrationIndexesVsIntensity
) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_params.isValid()); ree;

	int lastGlobalStopLeft = std::numeric_limits<int>::max();
	int lastGlobalStopRight = 0;

	int previousApex = 0;
	for (int i = 0; i < apexes.size(); i++) {

		const FrameIndex apexFrameIndex = apexes[i];
		const int nextApex = i < apexes.size() - 1 ? apexes[i + 1] : 2 * apexFrameIndex;

		const float apexIntensity = vec[apexFrameIndex];
		const float stopThresold = apexIntensity * m_params.stopThresholdFraction;

		int lastLeftStop = apexFrameIndex;
		float lastLeftValueLowest = apexIntensity;

		while (lastLeftStop > 0) {

			const int currentLeftStop = lastLeftStop - 1;
			const float currentLeftValue = vec[currentLeftStop];

			if (lastLeftStop <= previousApex) {
				lastLeftStop = static_cast<int>((previousApex + apexFrameIndex) / 2);
				lastGlobalStopLeft = lastLeftStop;
				break;
			}

			if (currentLeftValue <= stopThresold
				|| currentLeftValue > lastLeftValueLowest * (1.0 + m_params.hysteresis)
				|| lastLeftStop - 1 == lastGlobalStopRight
				|| currentLeftValue > apexIntensity
				) {
				lastLeftStop = currentLeftStop;
				lastGlobalStopLeft = lastLeftStop;
				break;
			}

			lastLeftStop = currentLeftStop;
			lastLeftValueLowest = std::min(currentLeftValue, lastLeftValueLowest);
			lastGlobalStopLeft = lastLeftStop;
		}

		int lastRightStop = apexFrameIndex;
		float lastRightValueLowest = vec[apexFrameIndex];
		while (lastRightStop < vecSize) {

			const int currentRightStop = lastRightStop + 1;
			const float currentRightValue = vec[currentRightStop];

			if (lastRightStop >= nextApex) {
				lastRightStop = static_cast<int>((nextApex + apexFrameIndex) / 2);
				lastGlobalStopRight = lastRightStop;
				break;
			}

			if (currentRightValue <= stopThresold
				|| currentRightValue > lastRightValueLowest * (1.0 + m_params.hysteresis)
				|| lastRightStop + 1 == lastGlobalStopLeft
				|| currentRightValue > apexIntensity
				) {
				// lastRightStop = currentRightStop;
				lastGlobalStopRight = lastRightStop;
				break;
				}

			lastRightStop = currentRightStop;
			lastRightValueLowest = std::min(currentRightValue, lastRightValueLowest);
			lastGlobalStopRight = lastRightStop;
		}

		previousApex = apexFrameIndex;
		peakIntegrationIndexesVsIntensity->push_back({{lastLeftStop, lastRightStop}, vec[apexFrameIndex], apexFrameIndex});
	}

    ERR_RETURN
}
