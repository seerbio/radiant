//
// Created by anichols on 12/11/22.
//

#include "FeatureFinderHill.h"

#include "ErrorUtils.h"
#include "MathUtils.h"


void FeatureFinderHill::addPoint(
        ScanNumberIndex scanNumberIndex,
        ScanNumber scanNumber,
        double mzVal,
        double intensityVal
) {
    m_scanNumberIndexes.push_back(scanNumberIndex);
    m_mzVals.push_back(mzVal);
    m_scanNumbers.push_back(scanNumber);
    m_intensities.push_back(intensityVal);
}

double FeatureFinderHill::mzMean() const {
    return MathUtils::mean(m_mzVals);
}

double FeatureFinderHill::mzStDev() const {
    return MathUtils::stDev(m_mzVals);
}

int FeatureFinderHill::maxIntensityScanNumber() const {
    return m_scanNumbers.at(MathUtils::findMaxIndexInVector(m_intensities));
}

int FeatureFinderHill::scanCount() const {
    return m_scanNumbers.size();
}

QVector<int> FeatureFinderHill::scanNumbers() const {
    return m_scanNumbers;
}

QVector<double> FeatureFinderHill::intensities() const {
    return m_intensities;
}

Err FeatureFinderHill::refineHill(
        int startIndex,
        int endIndex
) {

    ERR_INIT
    e = ErrorUtils::isIndexInRange(m_scanNumbers, startIndex, Error::eValueError); ree;
    e = ErrorUtils::isIndexInRange(m_scanNumbers, endIndex, Error::eValueError); ree;

    const int truncDistance = endIndex - startIndex + 1;

    m_scanNumbers = m_scanNumbers.mid(startIndex, truncDistance);
    m_intensities = m_intensities.mid(startIndex, truncDistance);
    m_mzVals = m_mzVals.mid(startIndex, truncDistance);

    ERR_RETURN
}

void FeatureFinderHill::updateIntensities(const QVector<double> &newIntensities) {
    m_intensities = newIntensities;
}

QPair<double, double> FeatureFinderHill::mzMinMax() const {
    const auto minMaxMz = std::minmax_element(m_mzVals.begin(), m_mzVals.end());
    return {*minMaxMz.first, *minMaxMz.second};
}

QPair<ScanNumber , ScanNumber> FeatureFinderHill::minMaxScanNumber() const {
    return {m_scanNumbers.front(), m_scanNumbers.back()};
}

double FeatureFinderHill::maxIntensityValue() const {
    return *std::max_element(m_intensities.begin(), m_intensities.end());
}

QVector<double> FeatureFinderHill::mzVals() const {
    return m_mzVals;
}
