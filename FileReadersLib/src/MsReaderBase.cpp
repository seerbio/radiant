//
// Created by Drucifer on 3/20/2022.
//

#include "MsReaderBase.h"

#include "ErrorUtils.h"

#include <QElapsedTimer>
#include <QFile>
#include <QtConcurrent/QtConcurrent>

#include <iostream>

MsReaderBase::MsReaderBase() : m_fileIsCalibrated(false) {}

void MsReaderBase::setMsScanInfo(const QMap<ScanNumber, MsScanInfo> &msScanInfos) {
    m_msScanInfo = msScanInfos;
}

Err MsReaderBase::setScanPoints(const QMap<ScanNumber, ScanPoints> &scanPoints) {

    ERR_INIT

    const int originalScanPointsSize = getMsScanInfos().size();
    e = ErrorUtils::isEqual(scanPoints.size(), originalScanPointsSize); ree;

    m_scanPoints = scanPoints;

    ERR_RETURN
}

Err MsReaderBase::getMsScanInfo(
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
void MsReaderBase::sortScanPoints(
        const ScanPointsSort &sort,
        ScanPoints *scanPoints
        ) {

    if (sort == ScanPointsSort::AscMz) {
        std::sort(scanPoints->begin(),  scanPoints->end(), IonsSortMzAsc);
    }
    else if (sort == ScanPointsSort::DescMz) {
        std::sort(scanPoints->rbegin(),  scanPoints->rend(), IonsSortMzAsc);
    }
    else if (sort == ScanPointsSort::AscIntensity) {
        std::sort(scanPoints->begin(),  scanPoints->end(), IonsSortIntensityAsc);
    }
    else if (sort == ScanPointsSort::DescIntensity) {
        std::sort(scanPoints->rbegin(),  scanPoints->rend(), IonsSortIntensityAsc);
    }

}

Err MsReaderBase::openFile(const QString &filePath) {
    return Error::eNoError;
}

Err MsReaderBase::openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
) {
    return Error::eNoError;
}

Err MsReaderBase::openFile(
        const QString &filePath,
        const QString &columnToFilterBy
) {
    return Error::eNoError;
}

Err MsReaderBase::closeFile() {
    return Error::eNoError;
}

int MsReaderBase::getNearestScanNumberFromScanTime(double scanTime) {

    if (m_scanNumberVsScanTime.isEmpty()) {
        for (const MsScanInfo &msScanInfo : m_msScanInfo) {
            m_scanNumberVsScanTime.insert(msScanInfo.scanNumber, msScanInfo.scanTime);
        }
    }

    const QVector<ScanTime> scanTimes = m_scanNumberVsScanTime.values().toVector();
    const QVector<ScanNumber> scanNumbers = m_scanNumberVsScanTime.keys().toVector();

    const int nearestIndex = MathUtils::closest(scanTimes, scanTime);

    return scanNumbers.at(nearestIndex);
}

Err MsReaderBase::getNearestScanNumberFromScanTime(
        ScanTime scanTime,
        const QVector<ScanNumber> &scanNumbers,
        const QVector<ScanTime> &scanTimes,
        ScanNumber *scanNumber
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanNumbers); ree;
    e = ErrorUtils::isEqual(scanNumbers.size(), scanTimes.size()); ree;

    const int nearestIndex = MathUtils::closest(scanTimes, scanTime);
    *scanNumber = scanNumbers.at(nearestIndex);

    ERR_RETURN
}

Err MsReaderBase::getNearestScanNumberFromScanNumber(
        ScanNumber scanNumber,
        const QVector<ScanNumber> &scanNumbers,
        ScanNumber *closestScanNumber
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanNumbers); ree;

    const int nearestIndex = MathUtils::closest(scanNumbers, scanNumber);

    *closestScanNumber = scanNumbers.at(nearestIndex);

    ERR_RETURN
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
        QVector<float> *intensityVals
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;

    for (const ScanPoint &sp : scanPoints) {
        mzVals->push_back(sp.x());
        intensityVals->push_back(sp.y());
    }

    ERR_RETURN
}

Err MsReaderBase::splitScanPoints(
        ScanPoints* scanPoints,
        QVector<double> *mzVals,
        QVector<float> *intensityVals
) {

    ERR_INIT

    e = ErrorUtils::isFalse(scanPoints->isEmpty()); ree;

    for (const ScanPoint &sp : *scanPoints) {
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
        e = getMsScanInfo(scanNumber, &msScanInfo); ree

        if (msScanInfo.msLevel != msLevel) {
            continue;
        }

        scanPoints->insert(scanNumber, scanPointsVec);
    }

    ERR_RETURN;
}

Err MsReaderBase::getScanPoints(
        int msLevel,
        QMap<ScanNumber, ScanPoints*> *scanPoints
        ) {

    ERR_INIT

    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {

        const ScanNumber scanNumber = it.key();
        ScanPoints *scanPointsVec = &it.value();

        MsScanInfo msScanInfo;
        e = getMsScanInfo(scanNumber, &msScanInfo); ree

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
        const QVector<float> &intensityVals,
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

QPair<Err, ScanPoints*> MsReaderBase::getScanPoints(int scanNumber) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msScanInfo.contains(scanNumber)); rree;
    e = ErrorUtils::isTrue(m_scanPoints.contains(scanNumber)); rree;
    return {e, &m_scanPoints[scanNumber]};
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
        const QVector<MzTargetKey> &uniqueInfoKeys,
        const QList<MsScanInfo> &msScanInfos
        ) {

    int cycleCounter = 0;

    int uniqueInfoKeysIndex = 0;
    for (const MsScanInfo &msScanInfo : msScanInfos) {

        if (uniqueInfoKeysIndex == uniqueInfoKeys.size()) {
            cycleCounter++;
            uniqueInfoKeysIndex = 0;
        }

        if (uniqueInfoKeys.at(uniqueInfoKeysIndex++) != msScanInfo.targetKey()) {
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

    QMap<MzTargetKey, int> uniqueKeyCounter;

    const QMap<ScanNumber, MsScanInfo> tandemScanInfos = getMsScanInfos(msLevel);

    for (const MsScanInfo &msScanInfo : tandemScanInfos) {
        uniqueKeyCounter[msScanInfo.targetKey()]++;
    }

    const bool uniqueKeysCycle = doUniqueKeysCycle(
            uniqueKeyCounter.keys().toVector(),
            tandemScanInfos.values()
            );

    return uniqueKeysCycle;
}

Err MsReaderBase::collateMS2MzTargetFrames(
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints *>> *diaTargetFrame) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;
    e = ErrorUtils::isNotEmpty(m_scanPoints); ree;
//    e = ErrorUtils::isTrue(isDIA()); ree;

    const int msLevel = 2;
    const QMap<ScanNumber, MsScanInfo> tandemScanInfos = getMsScanInfos(msLevel);

    for (auto it = tandemScanInfos.begin(); it != tandemScanInfos.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const MsScanInfo msScanInfo = it.value();

        QPair<Err, ScanPoints*> scanPointsResult = getScanPoints(scanNumber);
        e = scanPointsResult.first; ree;
        (*diaTargetFrame)[msScanInfo.targetKey()].insert(scanNumber, scanPointsResult.second);
    }

    qDebug() << "DIA Target Frames Count:" << diaTargetFrame->size();
    ERR_RETURN
}

QVector<MsScanInfo> MsReaderBase::getUniqueTandemMsScanInfos() {

    const int msLevel = 2;
    const QMap<ScanNumber, MsScanInfo> tandemScanInfos = getMsScanInfos(msLevel);

    QMap<MzTargetKey, MsScanInfo> uniqueMsScanInfos;

    for (const MsScanInfo &info : tandemScanInfos) {
        uniqueMsScanInfos[info.targetKey()] = info;
    }

    return uniqueMsScanInfos.values().toVector();
}

Err MsReaderBase::printFileInfo() {

    ERR_INIT

    int ms1ScanSize = -1;
    int ms2ScanSize = -1;
    int msLevel = 1;
    {
        const QMap<ScanNumber, MsScanInfo> ms1Scans = getMsScanInfos(msLevel);
        ms1ScanSize = ms1Scans.size();

        const QMap<ScanNumber, MsScanInfo> ms2Scans = getMsScanInfos(++msLevel);
        ms2ScanSize = ms2Scans.size();
    }

    const QVector<MsScanInfo> uniqueTandemScanInfos = getUniqueTandemMsScanInfos();

    qDebug() << "MsData FilePath" << m_filePath;
    qDebug() << "MS1 Scan Count" << ms1ScanSize;
    qDebug() << "MS2 Scan Count" << ms2ScanSize;
    qDebug() << "File is DIA" << isDIA();

//#define WRITE_TARGET_CE_FILE
#ifdef WRITE_TARGET_CE_FILE
    qDebug() << "ACHTUNG, ACHTUNG, ACHTUNG!!!! WRITING CF FILE IN:"; einfo;
    writeTargetCollisionEnergyFile();
#endif

    ERR_RETURN
}

QString MsReaderBase::filePath() {
    return m_filePath;
}

int MsReaderBase::getFrameCount() {
    return getUniqueTandemMsScanInfos().size();
}

QPair<double, double> MsReaderBase::scanTimeMinMax() {
    return {m_msScanInfo.first().scanTime, m_msScanInfo.last().scanTime};
}

bool MsReaderBase::isInit() {
    return !m_msScanInfo.isEmpty();
}

void MsReaderBase::reset() {
    QMap<ScanNumber, ScanPoints>().swap(m_scanPoints);
}
