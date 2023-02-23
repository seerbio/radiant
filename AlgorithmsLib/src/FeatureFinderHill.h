//
// Created by anichols on 12/11/22.
//

#ifndef PYTHIACPP_FEATUREFINDERHILL_H
#define PYTHIACPP_FEATUREFINDERHILL_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"

#include <QVector>


using namespace Error;


class ALGORITHMSLIB_EXPORTS FeatureFinderHill {

public:

    FeatureFinderHill() = default;
    ~FeatureFinderHill() = default;

    void addPoint(
        ScanNumber scanNumber,
        double mzVal,
        double intensityVal
    );

    Err refineHill(
            int startIndex,
            int endIndex
    );

    void updateIntensities(const QVector<double> &newIntensities);

    [[nodiscard]] double mzMean() const;

    [[nodiscard]] double mzStDev() const;

    [[nodiscard]] QPair<double, double> mzMinMax() const;

    [[nodiscard]] double maxIntensityScanNumber() const;

    [[nodiscard]] double maxIntensityValue() const;

    [[nodiscard]] int scanCount() const;

    [[nodiscard]] QVector<int> scanNumbers() const;

    [[nodiscard]] QPair<ScanNumber , ScanNumber> minMaxScanNumber() const;

    [[nodiscard]] QVector<double> intensities() const;


private:

    QVector<double> m_mzVals;
    QVector<int> m_scanNumbers;
    QVector<double> m_intensities;

};


#endif //PYTHIACPP_FEATUREFINDERHILL_H
