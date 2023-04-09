//
// Created by Drucifer on 4/27/2022.
//

#ifndef DEISOTOPER_TANDEM_H
#define DEISOTOPER_TANDEM_H

#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;


class DeisotoperTandem {


public:

    static Err deisotopeTandemScansParallel(
            const QMap<ScanNumber, ScanPoints> &tandemScans,
            double ppmTolerance,
            bool runParallel,
            QMap<ScanNumber, ScanPoints> *tandemScansDeisotoped
            );

    static ScanPoints deisotopeTandemScan(
            const ScanPoints &tandemScan,
            double ppmTolerance
            );


};


#endif //DEISOTOPER_TANDEM_H
