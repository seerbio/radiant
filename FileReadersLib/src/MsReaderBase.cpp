//
// Created by Drucifer on 3/20/2022.
//

#include "MsReaderBase.h"

#include "ErrorUtils.h"

#include <QElapsedTimer>
#include <QFile>
#include <QtConcurrent/QtConcurrent>

#include <iostream>

void MsReaderBase::setMsScanInfo(const QMap<ScanNumber, MsScanInfo> &msScanInfos) {
    m_msScanInfo = msScanInfos;
}

void MsReaderBase::setScanPoints(const QMap<ScanNumber, ScanPoints> &scanPoints) {
    m_scanPoints = scanPoints;
}

Err MsReaderBase::getScanInfo(
        ScanNumber scanNumber,
        MsScanInfo *msScanInfo
        ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msScanInfo.contains(scanNumber)); ree;
    *msScanInfo = m_msScanInfo.value(scanNumber);

    ERR_RETURN
}

QMap<ScanNumber, MsScanInfo> MsReaderBase::getMsScanInfos(int msLevel) {

    QMap<ScanNumber, MsScanInfo> msScanInfos;

    for (const MsScanInfo &mi : m_msScanInfo) {

        if (msLevel != mi.msLevel) {
            continue;
        }

        msScanInfos.insert(mi.scanNumber, mi);
    }

    return msScanInfos;
}

namespace {

    constexpr auto  IonsSortMzAsc = [](const ScanPoint &l, const ScanPoint &r){return l.x() < r.x();};

    constexpr auto IonsSortIntensityAsc = [](const ScanPoint &l, const ScanPoint &r){return l.y() < r.y();};

}//namespace
ScanPoints MsReaderBase::sortScanPoints(
        const ScanPoints &scanPoints,
        const ScanPointsSort &sort
        ) {

    ScanPoints vd = scanPoints;

    if (sort == ScanPointsSort::AscMz) {
        std::sort(vd.begin(),  vd.end(), IonsSortMzAsc);
    }
    else if (sort == ScanPointsSort::DescMz) {
        std::sort(vd.rbegin(),  vd.rend(), IonsSortMzAsc);
    }
    else if (sort == ScanPointsSort::AscIntensity) {
        std::sort(vd.begin(),  vd.end(), IonsSortIntensityAsc);
    }
    else if (sort == ScanPointsSort::DescIntensity) {
        std::sort(vd.rbegin(),  vd.rend(), IonsSortIntensityAsc);
    }

    return vd;
}

Err MsReaderBase::openFile(const QString &filePath) {
    return Error::eNoError;
}

Err MsReaderBase::closeFile() {
    return Error::eNoError;
}

int MsReaderBase::getNearestScanNumberFromScanTime(double scanTime) {

    QVector<double> scanTimes;
    for (const MsScanInfo &si : m_msScanInfo) {
        scanTimes.push_back(si.scanTime);
    }

    const int nearestIndex = MathUtils::closest(scanTimes, scanTime);

    const QList<int> &keys = m_msScanInfo.keys();

    return keys.at(nearestIndex);
}

QMap<ScanNumber, ScanTime> MsReaderBase::getScanNumberVsScanTime() const {

    QMap<ScanNumber, double> scanNumberVsScanTime;

    for (const MsScanInfo &msi : m_msScanInfo) {
        scanNumberVsScanTime.insert(msi.scanNumber, msi.scanTime);
    }

    return scanNumberVsScanTime;
}

Err MsReaderBase::splitScanPoints(
        const ScanPoints &scanPoints,
        QVector<double> *mzVals,
        QVector<double> *intensityVals
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;

    for (const ScanPoint &sp : scanPoints) {
        mzVals->push_back(sp.x());
        intensityVals->push_back(sp.y());
    }

    ERR_RETURN
}

Err MsReaderBase::getScanPoints(
        int msLevel,
        QMap<ScanNumber, ScanPoints> *scanPoints
        ) {

    ERR_INIT

    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const ScanPoints &scanPointsVec = it.value();

        MsScanInfo msScanInfo;
        e = getScanInfo(scanNumber, &msScanInfo); ree

        if (msScanInfo.msLevel != msLevel) {
            continue;
        }

        scanPoints->insert(scanNumber, scanPointsVec);
    }

    ERR_RETURN;
}

QMap<ScanNumber, MsScanInfo> MsReaderBase::getMsScanInfos() {
    return m_msScanInfo;
}

QMap<ScanNumber, ScanPoints> MsReaderBase::getScanPoints() {
    return m_scanPoints;
}

Err MsReaderBase::zipScanPoints(
        const QVector<double> &mzVals,
        const QVector<double> &intensityVals,
        ScanPoints *scanPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isEqual(mzVals.size(), intensityVals.size()); ree;

    for (int i = 0; i < mzVals.size(); i++) {
        scanPoints->push_back({mzVals.at(i), intensityVals.at(i)});
    }

    ERR_RETURN
}

Err MsReaderBase::getMsScanInfo(
        ScanNumber scanNumber,
        MsScanInfo *msScanInfo
        ){

    ERR_INIT

    e = ErrorUtils::isTrue(m_msScanInfo.contains(scanNumber)); ree;

    *msScanInfo = m_msScanInfo.value(scanNumber);

    ERR_RETURN
}

Err MsReaderBase::getScanPoints(
        int scanNumber,
        ScanPoints *scanPoints
        ) {

    ERR_INIT

    scanPoints->clear();

    e = ErrorUtils::isTrue(m_msScanInfo.contains(scanNumber)); ree;
    *scanPoints = m_scanPoints.value(scanNumber);

    ERR_RETURN
}

void MsReaderBase::printSize() {
    qDebug() << "MsScanInfo size" << m_msScanInfo.size();
    qDebug() << "ScanPoints size" << m_scanPoints.size();

}

Err MsReaderBase::updateScanPoints(const QMap<ScanNumber, ScanPoints> &scansToUpdate) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scansToUpdate); ree;

    for (auto it = scansToUpdate.begin(); it != scansToUpdate.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const ScanPoints &scanPoints = it.value();

        e = ErrorUtils::isTrue(m_scanPoints.contains(scanNumber)); ree;

        m_scanPoints[scanNumber] = scanPoints;
    }

    ERR_RETURN
}

bool doUniqueKeysCycle(
        const QVector<UniqueMsInfoScanKey> &uniqueInfoKeys,
        const QList<MsScanInfo> &msScanInfos
        ) {

    int cycleCounter = 0;

    int uniqueInfoKeysIndex = 0;
    for (const MsScanInfo &msScanInfo : msScanInfos) {

        if (uniqueInfoKeysIndex == uniqueInfoKeys.size()) {
            cycleCounter++;
            uniqueInfoKeysIndex = 0;
        }

        if (uniqueInfoKeys.at(uniqueInfoKeysIndex++) != msScanInfo.targetScanKey()) {
            return false;
        }

    }

    if (cycleCounter < 5) {
        return false;
    }

    const int expectedCycleSize = msScanInfos.size() / uniqueInfoKeys.size();

    return abs(cycleCounter - expectedCycleSize) <= 1;
}

bool MsReaderBase::isDIA() {

    const int msLevel = 2;

    QMap<UniqueMsInfoScanKey, int> uniqueKeyCounter;

    const QMap<ScanNumber, MsScanInfo> tandemScanInfos = getMsScanInfos(msLevel);

    for (const MsScanInfo &msScanInfo : tandemScanInfos) {
        uniqueKeyCounter[msScanInfo.targetScanKey()]++;
    }

    const bool uniqueKeysCycle = doUniqueKeysCycle(
            uniqueKeyCounter.keys().toVector(),
            tandemScanInfos.values()
            );

    return uniqueKeysCycle;
}

Err MsReaderBase::collateTandemPrecursorTargetsDIA(
        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrame
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;
    e = ErrorUtils::isNotEmpty(m_scanPoints); ree;
    e = ErrorUtils::isTrue(isDIA()); ree;

    const int msLevel = 2;
    const QMap<ScanNumber, MsScanInfo> tandemScanInfos = getMsScanInfos(msLevel);

    for (auto it = tandemScanInfos.begin(); it != tandemScanInfos.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const MsScanInfo msScanInfo = it.value();

        ScanPoints scanPoints;
        e = getScanPoints(
                scanNumber,
                &scanPoints
                ); ree;

        (*diaTargetFrame)[msScanInfo.targetScanKey()].insert(scanNumber, scanPoints);
    }

    ERR_RETURN
}
