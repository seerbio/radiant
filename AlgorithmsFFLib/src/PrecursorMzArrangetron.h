//
// Created by andrewnichols on 9/5/25.
//

#ifndef PRECURSORMZARRANGETRON_H
#define PRECURSORMZARRANGETRON_H

#include "AlgorithmsFFLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "MsReaderPointerAcc.h"

using namespace Error;
using PrecursorMzKey = int;

struct ScanNumberMzIntensity {
	MsScanInfo *scanInfoPntr = nullptr;
	float mzVal = -1;
	float intensityVal = -1;
};

class ALGORITHMSFFLIB_EXPORTS PrecursorMzArrangetron {

	friend class PrecursorMzArrangetronTests;

public:

	PrecursorMzArrangetron();
	~PrecursorMzArrangetron() = default;

	Err init(MsReaderPointerAcc *msReaderPtr);

	[[nodiscard]]bool isInit() const;

	Err trancheMsScanPoints(
		int threadCount,
		QVector<QVector<ScanNumberMzIntensity*>> *scanNumberMzIntensityTranched
		);

private:

	Err extractScansParallel();



private:
	MsReaderPointerAcc *m_msReaderPtr;

	QMap<PrecursorMzKey, QVector<ScanNumberMzIntensity>> m_precursorMzKeyVsScanNumberMzIntensity;

};



#endif //PRECURSORMZARRANGETRON_H
