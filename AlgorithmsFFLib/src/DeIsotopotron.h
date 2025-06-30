//
// Created by andrewnichols on 6/29/25.
//

#ifndef DEISOTOPOTRON_H
#define DEISOTOPOTRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;

class ALGORITHMSFFLIB_EXPORTS DeIsotopotron {

public:

	static Err deisotopeTandemScan(
		float extractionTolerancePPM,
		ScanPoints *scanPoints
		);

};



#endif //DEISOTOPOTRON_H
