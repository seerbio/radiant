//
// Created by anichols on 12/18/22.
//

#ifndef PYTHIACPP_TURBOXIC_H
#define PYTHIACPP_TURBOXIC_H

#include "AlgorithmsLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"


using namespace Error;


class ALGORITHMSLIB_EXPORTS TurboXIC {


public:

    TurboXIC();
    ~TurboXIC();

    Err init(const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber);

    XICPoints extractPoints(
            double mzMin,
            double mzMax,
            ScanNumber scanNumberMin,
            ScanNumber scanNumberMax
    );


private:

    Q_DISABLE_COPY(TurboXIC) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_TURBOXIC_H
