//
// Created by anichols on 12/11/22.
//

#ifndef PYTHIACPP_FEATUREFINDERHILL_H
#define PYTHIACPP_FEATUREFINDERHILL_H

#include "AlgorithmsFFLib_Exports.h"
#include "BiophysicalCalcs.h"
#include "CSVReader.h"
#include "Error.h"
#include "GlobalSettings.h"

#include <QVector>


using namespace Error;

class ALGORITHMSFFLIB_EXPORTS FeatureFinderHill {

public:

    FeatureFinderHill() = default;
    ~FeatureFinderHill() = default;

    bool operator == (const FeatureFinderHill &c) const {
        if (mzMean() == c.mzMean() && maxIntensityScanNumberIndex() == c.maxIntensityScanNumberIndex())
            return true;
        return false;
    }

    void addPoint(
        ScanNumberIndex scanNumberIndex,
        ScanNumber scanNumber,
        double mzVal,
        float intensityVal
    );

    [[nodiscard]] double mzMean() const;

    [[nodiscard]] double mzStDev() const;

    [[nodiscard]] QPair<double, double> mzMinMax() const;

    [[nodiscard]] int maxIntensityScanNumber() const;

    [[nodiscard]] int maxIntensityScanNumberIndex() const;

    [[nodiscard]] float intensityValueMax() const;

    [[nodiscard]] int scanCount() const;

    [[nodiscard]] QPair<ScanNumber , ScanNumber> scanNumberMinMax() const;

    [[nodiscard]] QPair<ScanNumber , ScanNumber> scanNumberIndexMinMax() const;

    [[nodiscard]] QVector<int> scanNumbers() const;

    [[nodiscard]] QVector<int> scanNumberIndexes() const;

    [[nodiscard]] QVector<double> mzVals() const;

    [[nodiscard]] QVector<float> intensities() const;


private:

    QVector<double> m_mzVals;
    QVector<ScanNumberIndex> m_scanNumberIndexes;
    QVector<ScanNumber> m_scanNumbers;
    QVector<float> m_intensities;

};


#endif //PYTHIACPP_FEATUREFINDERHILL_H
