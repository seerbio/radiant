//
// Created by andrewnichols on 7/29/25.
//

#ifndef MSFRAMEV2_H
#define MSFRAMEV2_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "TurboXIC.h"

//TODO BUILD QTESTS

class ALGORITHMSFFLIB_EXPORTS MsFrameV2 {

public:

	MsFrameV2() = default;
	~MsFrameV2() = default;

	Err init(const QVector<MsScan> &msScans);

	[[nodiscard]] bool isInit() const;

	[[nodiscard]] ScanTime getScanTimeFromFrameIndex(FrameIndex frameIndex) const;
	[[nodiscard]] ScanTime getScanTimeFromScanNumber(ScanNumber scanNumber) const;

	[[nodiscard]] ScanNumber getScanNumber(FrameIndex frameIndex) const;
	[[nodiscard]] ScanNumber getScanNumber(ScanTime scanTime) const;

	[[nodiscard]] FrameIndex getFrameIndex(ScanNumber scanNumber) const;
	[[nodiscard]] FrameIndex getFrameIndex(ScanTime scanTime) const;

	TurboXIC* getTurboXICPntr();

private:

	TurboXIC m_turboXIC;

	QVector<ScanNumber> m_scanNumbers;
	QVector<ScanTime> m_scanTimes;

	QVector<MsScanInfo*> m_msScanInfosTargetKey;
};



#endif //MSFRAMEV2_H
