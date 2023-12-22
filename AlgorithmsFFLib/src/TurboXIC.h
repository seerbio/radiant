//
// Created by anichols on 12/18/22.
//

#ifndef PYTHIACPP_TURBOXIC_H
#define PYTHIACPP_TURBOXIC_H

#include "AlgorithmsFFLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "PointFF.h"


using namespace Error;

struct ALGORITHMSFFLIB_EXPORTS XICPoint {
    float mz = -1.0;
    float intensity = -1.0;
    ScanNumber scanNumber = -1;
};

class ALGORITHMSFFLIB_EXPORTS XICPoints {

public:

    XICPoints() = default;
    ~XICPoints() = default;

    std::vector<XICPoint> xicPoints;

};


class ALGORITHMSFFLIB_EXPORTS TurboXIC {


public:

    TurboXIC();
    ~TurboXIC();

    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);
    Err init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints);

    XICPoints extractPointsXIC(
            float mzMin,
            float mzMax
    );

    Err getRTreeLimits(
            float *mzMin,
            float *mzMax
            ) const;

    bool isInit();


private:

    Q_DISABLE_COPY(TurboXIC) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_TURBOXIC_H
