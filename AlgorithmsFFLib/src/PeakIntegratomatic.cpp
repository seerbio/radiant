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
	const float* vec,
	int vecSize,
	QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_params.isValid()); ree;

	float lastIntensity = 0.0f;
	float apexIntensity = 0.0f;

	bool inPeak = false;
	bool isBackSide = false;

	int peakStart = 0;
	int peakEnd = 0;

	for (int i = 0; i < vecSize; i++) {

		const float currentIntensity = vec[i];

		if (currentIntensity > lastIntensity && !MathUtils::tZero(currentIntensity) && !inPeak) {
			qDebug() << "PeakStart";
			peakStart = i - 1;
			inPeak = true;
		}

		qDebug() << i << lastIntensity << currentIntensity << apexIntensity << inPeak;

		if (!inPeak) {
			lastIntensity = currentIntensity;
			continue;
		}

		if (currentIntensity > lastIntensity) {
			apexIntensity = currentIntensity;
		}

		if (currentIntensity < m_params.stopThresholdFraction * apexIntensity) {
			peakEnd = i - 1;
			peakIntegrationIndexesVsIntensity->push_back({{peakStart, peakEnd}, apexIntensity});

			inPeak = false;
			isBackSide = false;
			apexIntensity = 0.0f;
			lastIntensity = currentIntensity;
			peakStart = peakEnd;
			qDebug() << "PeakEnd";
		}

		lastIntensity = currentIntensity;
	}





    ERR_RETURN
}
