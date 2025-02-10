//
// Created by Drucifer on 3/20/2022.
//

#include "MsReaderBase.h"

#include "ErrorUtils.h"

#include <QElapsedTimer>
#include <QFile>
#include <QtConcurrent/QtConcurrent>

#include <iostream>

MsReaderBase::MsReaderBase()
: m_fileIsCalibrated(false)
, m_isTIMS(false)
, m_mzMs1Min(std::numeric_limits<float>::max())
, m_mzMs1Max(-1.0)
, m_mzMs2Min(std::numeric_limits<float>::max())
, m_mzMs2Max(-1.0)
{}

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

void MsReaderBase::reset() {
    QMap<ScanNumber, ScanPoints>().swap(m_scanPoints);
    QMap<ScanNumber, MsScanInfo>().swap(m_msScanInfo);
}

Err MsReaderBase::openFile(const QString &filePath) {
    return Error::eFunctionNotImplemented;
}

Err MsReaderBase::openFile(
        const QString &filePath,
        const QString &columnToFilterBy,
        const QPair<double, double> &filterRange
) {
    return Error::eFunctionNotImplemented;
}

Err MsReaderBase::openFile(
        const QString &filePath,
        const QString &columnToFilterBy
        ) {
    return Error::eFunctionNotImplemented;
}

Err MsReaderBase::closeFile() {
    QMap<ScanNumber, MsScanInfo>().swap(m_msScanInfo);
    QMap<ScanNumber, ScanPoints>().swap(m_scanPoints);
    QMap<ScanNumber, ScanTime>().swap(m_scanNumberVsScanTime);
    return Error::eNoError;
}

Err MsReaderBase::extractScanPoints(
    const QVector<MsScanInfo*> &msScanInfos,
    QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
    ){
    return eFunctionNotImplemented;
}

Err MsReaderBase::getMzTargetScanPoints(
    const MzTargetKey &targetKey,
    QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
    ){
    return eFunctionNotImplemented;
}

QString MsReaderBase::filePath() {
    return m_filePath;
}


bool MsReaderBase::isDIA() {

    const int msLevel = 2;

    QMap<MzTargetKey, int> uniqueKeyCounter;

    const QMap<ScanNumber, MsScanInfo> tandemScanInfos = getMsScanInfos(msLevel);

    for (const MsScanInfo &msScanInfo : tandemScanInfos) {
        uniqueKeyCounter[msScanInfo.targetKey()]++;
    }

    const double targetKeyCountMean = MathUtils::mean(uniqueKeyCounter.values());
    const double allowableDistance = 2.0;

    QVector<bool> withinCountDistance;
    std::transform(
            uniqueKeyCounter.begin(),
            uniqueKeyCounter.end(),
            std::back_inserter(withinCountDistance),
            [targetKeyCountMean, allowableDistance](int targetKeyCount){
                return std::abs(targetKeyCount - targetKeyCountMean) <= allowableDistance
                            && targetKeyCount > allowableDistance;
            }
            );

    return std::all_of(withinCountDistance.begin(), withinCountDistance.end(), [](bool b){return b;});
}

bool MsReaderBase::isInit() {
    return !m_msScanInfo.isEmpty();
}

QPair<ScanTime, ScanTime > MsReaderBase::scanTimeMinMax() {
    return {m_msScanInfo.first().scanTime, m_msScanInfo.last().scanTime};
}

QMap<ScanNumber, ScanPoints> MsReaderBase::getScanPoints() {
    return m_scanPoints;
}

QMap<ScanNumber, ScanPoints *> MsReaderBase::getScanPointsPntrs() {

    QMap<ScanNumber, ScanPoints *> scanNumberVsScanPointsPtrs;

    for (auto it = m_scanPoints.begin(); it != m_scanPoints.end(); it++) {
        scanNumberVsScanPointsPtrs.insert(it.key(), &it.value());
    }

    return scanNumberVsScanPointsPtrs;
}

Err MsReaderBase::getScanPoints(
        int msLevel,
        QMap<ScanNumber, ScanPoints> *scanPoints
) {

    ERR_INIT

    scanPoints->clear();

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

QPair<Err, ScanPoints*> MsReaderBase::getScanPoints(int scanNumber) {

    ERR_INIT

    e = ErrorUtils::contains(scanNumber, m_msScanInfo); rree;
    e = ErrorUtils::contains(scanNumber, m_scanPoints); rree;
    return {e, &m_scanPoints[scanNumber]};
}

Err MsReaderBase::collateMS2MzTargetFrames(
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints *>> *diaTargetFrame
        ) {

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

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "DIA Target Frames Count:" << diaTargetFrame->size();
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

int MsReaderBase::getFrameCount() {
    return getUniqueTandemMsScanInfos().size();
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

Err MsReaderBase::getMsScanInfos(
    const MzTargetKey &mzTargetKey,
    QVector<MsScanInfo*> *msScanInfos
    ) const {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_mzTargetVsScanInfosPntrs); ree;
    e = ErrorUtils::contains(mzTargetKey, m_mzTargetVsScanInfosPntrs); ree;

    *msScanInfos = m_mzTargetVsScanInfosPntrs.value(mzTargetKey);

    ERR_RETURN
}

QMap<ScanNumber, MsScanInfo> MsReaderBase::getMsScanInfos() {
    return m_msScanInfo;
}

Err MsReaderBase::getMsScanInfo(
        ScanNumber scanNumber,
        MsScanInfo *msScanInfo
) const {

    ERR_INIT

    e = ErrorUtils::contains(scanNumber, m_msScanInfo); ree;

    *msScanInfo = m_msScanInfo.value(scanNumber);

    ERR_RETURN
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

ScanNumber MsReaderBase::getNearestScanNumberFromScanTime(ScanTime scanTime) {

    if (m_scanNumberVsScanTime.isEmpty()) {
        for (const MsScanInfo &msScanInfo : m_msScanInfo) {
            m_scanNumberVsScanTime.insert(msScanInfo.scanNumber, msScanInfo.scanTime);
        }
    }

    const QVector<ScanTime> scanTimes = m_scanNumberVsScanTime.values().toVector();
    const QVector<ScanNumber> scanNumbers = m_scanNumberVsScanTime.keys().toVector();

    const ScanNumber nearestIndex = MathUtils::closest(scanTimes, scanTime);

    return scanNumbers.at(nearestIndex);
}

QMap<ScanNumber, ScanTime> MsReaderBase::getScanNumberVsScanTime() const {

    QMap<ScanNumber, ScanTime > scanNumberVsScanTime;

    for (const MsScanInfo &msi : m_msScanInfo) {
        scanNumberVsScanTime.insert(msi.scanNumber, msi.scanTime);
    }

    return scanNumberVsScanTime;
}

Err MsReaderBase::getHiLoMzPrecursors(QPair<MzMin, MzMax> *precursorMzLoVsMzHi) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;

    const QVector<MsScanInfo> uniqueMsScanInfos = getUniqueTandemMsScanInfos();
    const auto minMaxMsScanInfo = std::minmax_element(
        uniqueMsScanInfos.begin(),
        uniqueMsScanInfos.end(),
        [](const MsScanInfo &l, const MsScanInfo &r) {
            return l.precursorTargetMz < r.precursorTargetMz;
        }
        );

    *precursorMzLoVsMzHi = {
        minMaxMsScanInfo.first->precursorTargetMz - minMaxMsScanInfo.first->isoWindowLower,
        minMaxMsScanInfo.second->precursorTargetMz + minMaxMsScanInfo.second->isoWindowLower
    };

    ERR_RETURN
}

Err MsReaderBase::driftTimeFromIonMobilityIndex(
    const IonMobilityIndex &ionMobilityIndex,
    double *driftTime
    ) const {

    ERR_INIT

    e = ErrorUtils::contains(ionMobilityIndex, m_frameIndexVsDriftTime); ree;

    *driftTime = m_frameIndexVsDriftTime.value(ionMobilityIndex);

    ERR_RETURN
}

Err MsReaderBase::splitScanPoints(
        const ScanPoints &scanPoints,
        QVector<float> *mzVals,
        QVector<float> *intensityVals
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;

    mzVals->clear();
    intensityVals->clear();

    for (const ScanPoint &sp : scanPoints) {
        mzVals->push_back(sp.x());
        intensityVals->push_back(sp.y());
    }

    ERR_RETURN
}

Err MsReaderBase::splitScanPoints(
        ScanPoints* scanPoints,
        QVector<float> *mzVals,
        QVector<float> *intensityVals
) {

    ERR_INIT

    e = ErrorUtils::isFalse(scanPoints->isEmpty()); ree;

    mzVals->clear();
    intensityVals->clear();

    for (const ScanPoint &sp : *scanPoints) {
        mzVals->push_back(sp.x());
        intensityVals->push_back(sp.y());
    }

    ERR_RETURN
}

Err MsReaderBase::zipScanPoints(
        const QVector<float> &mzVals,
        const QVector<float> &intensityVals,
        ScanPoints *scanPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isEqual(mzVals.size(), intensityVals.size()); ree;
    scanPoints->clear();

    for (int i = 0; i < mzVals.size(); i++) {
        scanPoints->push_back({mzVals.at(i), intensityVals.at(i)});
    }

    ERR_RETURN
}

void MsReaderBase::printSize() const {
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsScanInfo size" << m_msScanInfo.size();
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "ScanPoints size" << m_scanPoints.size();
}

Err MsReaderBase::printFileInfo() {

    ERR_INIT

    int ms1ScanSize = -1;
    int ms2ScanSize = -1;
    int msLevel = 1;

    {
        const QMap<ScanNumber, MsScanInfo> ms1Scans = getMsScanInfos(msLevel);
        ms1ScanSize = ms1Scans.size();

        for (const MsScanInfo &msScanInfo : ms1Scans) {
            const ScanPoints &scanPoints = m_scanPoints.value(msScanInfo.scanNumber);
            m_mzMs1Min = std::min(m_mzMs1Min, scanPoints.front().x());
            m_mzMs1Max = std::max(m_mzMs1Max, scanPoints.back().x());
        }

        const QMap<ScanNumber, MsScanInfo> ms2Scans = getMsScanInfos(++msLevel);
        ms2ScanSize = ms2Scans.size();

        for (const MsScanInfo &msScanInfo : ms2Scans) {
            const ScanPoints &scanPoints = m_scanPoints.value(msScanInfo.scanNumber);
            m_mzMs2Min = std::min(m_mzMs1Min, scanPoints.front().x());
            m_mzMs2Max = std::max(m_mzMs1Max, scanPoints.back().x());
        }
    }

    const QVector<MsScanInfo> uniqueTandemScanInfos = getUniqueTandemMsScanInfos();
    const int meanFrameScanCount = static_cast<int>(std::round(static_cast<float>(ms2ScanSize) / uniqueTandemScanInfos.size()));

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MsData FilePath" << m_filePath;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MS1 Scan Range" << m_mzMs1Min << "-" << m_mzMs1Max;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MS1 Scan Count" << ms1ScanSize;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MS2 Scan Range" << m_mzMs2Min << "-" << m_mzMs2Max;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MS2 Scan Count" << ms2ScanSize;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MS2 Target Window Count" << uniqueTandemScanInfos.size();
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "MS2 Target Window Mean Scan Count" << meanFrameScanCount;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "File is DIA" << isDIA();

    ERR_RETURN
}

bool MsReaderBase::isTIMS() const {
    return m_isTIMS;
}

float MsReaderBase::mzMs2Min() const {
    return m_mzMs2Min;
}

float MsReaderBase::mzMs2Max() const {
    return m_mzMs2Max;
}

QMap<FrameNumberTIMS, Ms1FrameTIMS>* MsReaderBase::frameNumberVsMS1FrameTIMSPntr() {
    return &m_frameNumberVsMS1FrameTIMS;
}
