//
// Created by Drucifer on 4/27/2022.
//

#ifndef DEISOTOPER_TANDEM_H
#define DEISOTOPER_TANDEM_H

#include "GlobalSettings.h"


class DeisotoperTandem {


public:

    static ScanPoints deisotopeTandemScan(
            const ScanPoints &tandemScan,
            double ppmTolerance
            );


};


#endif //DEISOTOPER_TANDEM_H
