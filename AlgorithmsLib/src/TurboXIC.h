//
// Created by anichols on 12/18/22.
//

#ifndef PYTHIACPP_TURBOXIC_H
#define PYTHIACPP_TURBOXIC_H

#include "AlgorithmsLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"


using namespace Error;

class XICPoints {

public:

    XICPoints() = default;
    ~XICPoints() = default;

    QMap<ScanNumber, double> scanNumbersVsIntensityVals;
    QMap<ScanNumber, QVector<double>> scanNumberVsMzVals;
};


class ALGORITHMSLIB_EXPORTS TurboXIC {


public:

    TurboXIC();
    ~TurboXIC();

    Err init(const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber);

    XICPoints extractPointsXIC(
            double mzMin,
            double mzMax,
            ScanNumber scanNumberMin,
            ScanNumber scanNumberMax
    );

    [[nodiscard]] ScanPoints extractSpectrum(
            double mzMin,
            double mzMax,
            ScanNumber scanNumberMin,
            ScanNumber scanNumberMax
            ) const;

    Err getRTreeLimits(
            double *scanNumberMin,
            double *scanNumberMax,
            double *mzMin,
            double *mzMax
            ) const;

    bool isInit();


private:

    Q_DISABLE_COPY(TurboXIC) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_TURBOXIC_H
