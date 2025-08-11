//
// Created by andrewnichols on 8/11/25.
//

#ifndef APEXCONNECTOR_H
#define APEXCONNECTOR_H

#include "AlgorithmsFFLib_Exports.h"

#include "GlobalSettings.h"
#include "Error.h"

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS ApexConnector {

public:

	static Err connectApexes(
		const QVector<QVector<int>> &apexes,
		QVector<QVector<int>> *connectedApexes
		);

};



#endif //APEXCONNECTOR_H
