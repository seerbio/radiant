//
// Created by anichols on 12/11/22.
//

#ifndef PYTHIACPP_FEATUREFINDERHILL_H
#define PYTHIACPP_FEATUREFINDERHILL_H

#include "AlgorithmsLib_Exports.h"
#include "BiophysicalCalcs.h"
#include "CSVReader.h"
#include "Error.h"
#include "GlobalSettings.h"

#include <QVector>


using namespace Error;


struct ALGORITHMSLIB_EXPORTS FeatureFinderMS1Feature {

    double mzFound = -1.0;
    double cosineSim = -1.0;
    Charge charge = -1;
    int monoIsotopeOffset = 0;
    QPair<ScanNumber, ScanNumber> scanNumberMinMax;
    ScanPoints foundIsotopes;

    [[nodiscard]] double featureMass() const {
        return BiophysicalCalcs::calculateMassFromThomson(
                mzFound,
                charge,
                monoIsotopeOffset
                );
    }

};


class ALGORITHMSLIB_EXPORTS FeatureFinderHill {

public:

    FeatureFinderHill() = default;
    ~FeatureFinderHill() = default;

    bool operator == (const FeatureFinderHill &c) const
    {
        if (mzMean() == c.mzMean() && maxIntensityScanNumberIndex() == c.maxIntensityScanNumber())
            return true;
        return false;
    }

    void addPoint(
        ScanNumberIndex scanNumberIndex,
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

    [[nodiscard]] int maxIntensityScanNumber() const;

    [[nodiscard]] int maxIntensityScanNumberIndex() const;

    [[nodiscard]] double intensityValueMax() const;

    [[nodiscard]] int scanCount() const;

    [[nodiscard]] QPair<ScanNumber , ScanNumber> scanNumberMinMax() const;

    [[nodiscard]] QPair<ScanNumber , ScanNumber> scanNumberIndexMinMax() const;

    [[nodiscard]] QVector<int> scanNumbers() const;

    [[nodiscard]] QVector<int> scanNumberIndexes() const;

    [[nodiscard]] QVector<double> mzVals() const;

    [[nodiscard]] QVector<double> intensities() const;


private:

    QVector<double> m_mzVals;
    QVector<ScanNumberIndex> m_scanNumberIndexes;
    QVector<ScanNumber> m_scanNumbers;
    QVector<double> m_intensities;

};

struct ALGORITHMSLIB_EXPORTS FeatureFinderHillPlus {
    FeatureFinderHill featureFinderHill;
    IonType ionType;
    IonIndex ionIndex = -1;
    double cosineSimToAnchor = -1.0;

    double bestIsotopologueCosineSim = -1.0;
    Charge isotopologueCharge = -1;
    double isotopologueIntensity = -1.0;
};

namespace FeatureFinderHillUtils {

    void ALGORITHMSLIB_EXPORTS sortFeatureFinderHillsPlussesIntensityDesc(
            QVector<FeatureFinderHillPlus> *featureFinderHills
            );

    void ALGORITHMSLIB_EXPORTS sortFeatureFinderHillsPlussesMzAsc(
            QVector<FeatureFinderHillPlus> *featureFinderHills
    );

}


#endif //PYTHIACPP_FEATUREFINDERHILL_H
