//
// Created by andrewnichols on 7/29/25.
//

#include "MsFrameV2.h"

Err MsFrameV2::init(const QVector<MsScan> &msScans) {

	ERR_INIT

	e = ErrorUtils::isNotEmpty(msScans);

	std::transform(
		msScans.begin(),
		msScans.end(),
		std::back_inserter(m_msScanInfosTargetKey),
		[](const MsScan &ms){return ms.msScanInfoPntr;}
		);

	std::transform(
		msScans.begin(),
		msScans.end(),
		std::back_inserter(m_scanNumbers),
		[](const MsScan &ms){return ms.msScanInfoPntr->scanNumber;}
		);

	std::transform(
		msScans.begin(),
		msScans.end(),
		std::back_inserter(m_scanTimes),
		[](const MsScan &ms){return ms.msScanInfoPntr->scanTime;}
		);

	e = m_turboXIC.init(msScans); ree;

	ERR_RETURN
}

bool MsFrameV2::isInit() const {
	return m_turboXIC.isInit()
		   && !m_scanNumbers.isEmpty()
		   && !m_scanTimes.isEmpty()
		   && !m_msScanInfosTargetKey.isEmpty();
}

ScanTime MsFrameV2::getScanTimeFromFrameIndex(FrameIndex frameIndex) const {
	return m_scanTimes[frameIndex];
}

ScanTime MsFrameV2::getScanTimeFromScanNumber(ScanNumber scanNumber) const {
	const long frameIndex = getFrameIndex(scanNumber);
	return m_scanTimes[static_cast<FrameIndex>(frameIndex)];
}

ScanNumber MsFrameV2::getScanNumber(ScanTime scanTime) const {
	const long frameIndex = getFrameIndex(scanTime);
	return m_scanNumbers[static_cast<FrameIndex>(frameIndex)];
}

FrameIndex MsFrameV2::getFrameIndex(ScanNumber scanNumber) const {

	auto it = std::lower_bound(m_scanNumbers.begin(), m_scanNumbers.end(), scanNumber);

	if (it == m_scanNumbers.begin() && *it > scanNumber) {
		return std::numeric_limits<FrameIndex>::quiet_NaN();;
	}

	if (it == m_scanNumbers.end() || *it > scanNumber) {
		--it;
	}

	const long frameIndex = it - m_scanNumbers.begin();
	return static_cast<FrameIndex>(frameIndex);
}

FrameIndex MsFrameV2::getFrameIndex(ScanTime scanTime) const {

	auto it = std::lower_bound(m_scanTimes.begin(), m_scanTimes.end(), scanTime);

	if (it == m_scanTimes.begin() && *it > scanTime) {
		return std::numeric_limits<FrameIndex>::quiet_NaN();;
	}

	if (it == m_scanTimes.end() || *it > scanTime) {
		--it;
	}

	const long frameIndex = it - m_scanTimes.begin();
	return static_cast<FrameIndex>(frameIndex);
}

int MsFrameV2::frameIndexSize() const {
	return m_scanNumbers.size();
}

float MsFrameV2::mzMin() const {
	const auto it = std::min_element(
		m_msScanInfosTargetKey.begin(),
		m_msScanInfosTargetKey.end(),
		[](const MsScanInfo *l, const MsScanInfo *r){return l->mzMin < r->mzMin;}
		);
	return (*it)->mzMin;
}

float MsFrameV2::mzMax() const {
	const auto it = std::max_element(
		m_msScanInfosTargetKey.begin(),
		m_msScanInfosTargetKey.end(),
		[](const MsScanInfo *l, const MsScanInfo *r){return l->mzMax < r->mzMax;}
		);
	return (*it)->mzMax;
}

ScanNumber MsFrameV2::getScanNumber(FrameIndex frameIndex) const {
	return m_scanNumbers[frameIndex];
}

TurboXIC * MsFrameV2::getTurboXICPntr() {
	return &m_turboXIC;
}
