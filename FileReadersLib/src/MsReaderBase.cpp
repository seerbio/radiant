//
// Created by Drucifer on 3/20/2022.
//

#include "MsReaderBase.h"

#include "ErrorUtils.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>

#include <iostream>


Err MsReaderBase::buildScanIonWithScanInfoInputMs2(const QPair<double, double> &mzMinMax,
                                                   int skipEveryNScans,
                                                   QVector<ScanIonWithScanInfo> *scanIonWithScanInfos) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;
    e = ErrorUtils::isNotEmpty(m_scanIons); ree;

    const int minAllowedSkipEveryNScansValue = 1;
    skipEveryNScans = std::max(skipEveryNScans, minAllowedSkipEveryNScansValue);

    const QMap<ScanNumber, bool> allowedScanNumber = buildSkipNAllowableIndexesMs2(skipEveryNScans);

    for (const ScanIon &si : m_scanIons) {

        if (si.mz < mzMinMax.first || si.mz > mzMinMax.second || !allowedScanNumber.value(si.scanNumber)) {
            continue;
        }

        ScanIonWithScanInfo scanIonWithScanInfo;
        scanIonWithScanInfo.scanIon = si;
        scanIonWithScanInfo.msScanInfo = m_msScanInfo.value(si.scanNumber);

        scanIonWithScanInfos->push_back(scanIonWithScanInfo);
    }

    ERR_RETURN
}

QMap<ScanNumber, bool>  MsReaderBase::buildSkipNAllowableIndexesMs2(int skipEveryNScans) {

    const QList<ScanNumber> &scanNumbers = m_msScanInfo.keys();

    QMap<ScanNumber, bool> ms2SkipScanNumbers;
    for (int i = 0; i < scanNumbers.size(); i++) {

        const MsScanInfo &msScanInfo = m_msScanInfo.value(i);

        if (i % skipEveryNScans != 0 || msScanInfo.msLevel == 1) {
            ms2SkipScanNumbers.insert(scanNumbers.value(i), false);
            continue;
        }

        ms2SkipScanNumbers.insert(scanNumbers.value(i), true);
    }

    return ms2SkipScanNumbers;
}

Err MsReaderBase::getScanInfo(ScanNumber scanNumber,
                              MsScanInfo *msScanInfo) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msScanInfo.contains(scanNumber)); ree;
    *msScanInfo = m_msScanInfo.value(scanNumber);

    ERR_RETURN
}

Err MsReaderBase::getScanPoints(ScanNumber scanNumber, QVector<QPointF> *scanPoints) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msScanInfo.contains(scanNumber)); ree;
    *scanPoints = m_scanIonsMapped.value(scanNumber);

    ERR_RETURN
}

Err MsReaderBase::buildScanPointsByScanNumber(QMap<ScanNumber, ScanPoints> *scanPointsByScanNumber) {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_scanIons); ree;

    for (const ScanIon &si : m_scanIons) {
        (*scanPointsByScanNumber)[si.scanNumber].push_back({si.mz, si.intensity});
    }

    ERR_RETURN
}

Err MsReaderBase::buildScanPointsByScanNumber(
        MsLevel msLevel,
        QMap<ScanNumber, ScanPoints> *scanPointsByScanNumber
) {
    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_scanIons); ree;

    for (const ScanIon &si : m_scanIons) {

        MsScanInfo msScanInfo;
        e = getScanInfo(si.scanNumber, &msScanInfo);

        if (msScanInfo.msLevel != msLevel) {
            continue;
        }

        (*scanPointsByScanNumber)[si.scanNumber].push_back({si.mz, si.intensity});
    }

    ERR_RETURN
}

QVector<MsScanInfo> MsReaderBase::getMsScanInfos(int msLevel) {

    QVector<MsScanInfo> msScanInfos;

    for (const MsScanInfo &mi : m_msScanInfo) {

        if (msLevel != mi.msLevel) {
            continue;
        }

        msScanInfos.push_back(mi);
    }

    return msScanInfos;
}

namespace {

    constexpr auto  IonsSortMzAsc = [](const ScanPoint &l, const ScanPoint &r){return l.x() < r.x();};

    constexpr auto IonsSortIntensityAsc = [](const ScanPoint &l, const ScanPoint &r){return l.y() < r.y();};

}//namespace
ScanPoints MsReaderBase::sortScanPoints(const ScanPoints &scanPoints,
                                        const ScanPointsSort &sort) {

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

Err MsReaderBase::openFile(const QString &filePath,
                           double defaultIsolationWindow,
                           int defaultCollisionEnergy) {
    return Error::eNoError;
}

Err MsReaderBase::closeFile() {
    return Error::eNoError;
}


void MsReaderBase::loadScanIonsMapped() {

    for (const ScanIon &si : m_scanIons) {
        m_scanIonsMapped[si.scanNumber].push_back({si.mz, si.intensity});
    }

//    for (const QVector<QPointF> &sis : m_scanIonsMapped) {
//        qDebug() << "SDFSDFDS" << sis.size();
//    }

}


Err MsReaderBase::getBasePeakChromatogram(
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        QMap<ScanNumber, QPointF> *basePeakChromatogram
) {

    ERR_INIT

    basePeakChromatogram->clear();

    e = ErrorUtils::isNotEmpty(m_scanIonsMapped); ree;

    const auto maxPointLogic = [](const QPointF &l, const QPointF &r){
        return l.y() < r.y();
    };

    for (auto it = m_msScanInfo.begin(); it != m_msScanInfo.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const MsScanInfo msScanInfo = it.value();

        const UniqueMsInfoScanKey currentUniqueMsInfoScanKey = msScanInfo.msScanInfoKey();

        if (currentUniqueMsInfoScanKey != uniqueMsInfoScanKey
            && uniqueMsInfoScanKey != S_GLOBAL_SETTINGS.NONE) {
            continue;
        }

        qDebug() << "Loading scan number" << scanNumber;

        const QVector<QPointF> &scanPoints = m_scanIonsMapped.value(scanNumber);
        const QPointF maxScanPoint = *std::max_element(scanPoints.begin(), scanPoints.end(), maxPointLogic);

        basePeakChromatogram->insert(scanNumber, {msScanInfo.scanTime, maxScanPoint.y()});
    }

    ERR_RETURN
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

QVector<QPointF> MsReaderBase::getNearestTandemScan(double scanTime) {

    const ScanNumber scanNumber = getNearestScanNumberFromScanTime(scanTime);
    return m_scanIonsMapped.value(scanNumber);
}

QMap<ScanNumber, QVector<QPointF>> MsReaderBase::scanIonsMapped() {

    if(m_scanIonsMapped.isEmpty()) {
        loadScanIonsMapped();
    }

    return m_scanIonsMapped;
}


QString MsReaderBase::buildUniqueTandemWindowKey(const MsScanInfo &si) {

    QString key = QString::number(si.precursorTargetMz - si.precursorWindowOffsetLower)
                + MsReaderBase::separator()
                + QString::number(si.precursorTargetMz + si.precursorWindowOffsetUpper);

    return key;
}


QMap<IsolationWindowKey, QVector<ScanNumber>> MsReaderBase::getAllUniqueTandemIsolationWindows() {

    QMap<IsolationWindowKey, QVector<ScanNumber>> uniqueTandemIsolationWindows;

    for (const MsScanInfo &si : m_msScanInfo) {

        if (si.msLevel < 2) {
            continue;
        }

        const IsolationWindowKey uniqueKey = buildUniqueTandemWindowKey(si);
        uniqueTandemIsolationWindows[uniqueKey].push_back(si.scanNumber);
    }

    return uniqueTandemIsolationWindows;
}

Err MsReaderBase::getRangesFromIsolationWindowKey(
        const IsolationWindowKey &key,
        QPair<double, double> *mzIsoRange
) {

    ERR_INIT

    const QStringList splitIsolationWindowKey = key.split(
            MsReaderBase::separator(),
            QString::SkipEmptyParts
    );

    const int expectedSize = 2;
    e = ErrorUtils::isEqual(splitIsolationWindowKey.size(), expectedSize); ree;

    double mzIsoRangeStart;
    e = ErrorUtils::toDouble(splitIsolationWindowKey.front(), &mzIsoRangeStart); ree;

    double mzIsoRangeEnd;
    e = ErrorUtils::toDouble(splitIsolationWindowKey.back(), &mzIsoRangeEnd); ree;

    *mzIsoRange = {mzIsoRangeStart, mzIsoRangeEnd};

    ERR_RETURN
}


QStringList MsReaderBase::getUniqueMsInfoScanKeys() {

    QSet<QString> uniqueScanInfoKeys;
    for (const MsScanInfo &msScanInfo : m_msScanInfo) {
        uniqueScanInfoKeys.insert(msScanInfo.msScanInfoKey());
    }

    return uniqueScanInfoKeys.toList();
}


Err MsReaderBase::buildScanNumbersByUniqueMsInfoScanKey() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;

    for (const MsScanInfo msScanInfo : m_msScanInfo) {
        m_scanNumbersByUniqueMsInfoScanKey[msScanInfo.msScanInfoKey()].push_back(msScanInfo.scanNumber);
    }

    ERR_RETURN
}


Err MsReaderBase::getScansByUniqueMsInfoScanKey(
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        QMap<ScanNumber, ScanPoints> *scanPointsMap
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_scanIonsMapped); ree;

    if (m_scanNumbersByUniqueMsInfoScanKey.isEmpty()) {
        e = buildScanNumbersByUniqueMsInfoScanKey(); ree;
    }

    e = ErrorUtils::contains(
            uniqueMsInfoScanKey,
            m_scanNumbersByUniqueMsInfoScanKey.keys().toVector()
    ); ree;

    const QVector<ScanNumber> scanNumbers = m_scanNumbersByUniqueMsInfoScanKey.value(uniqueMsInfoScanKey);

    for (ScanNumber scanNumber : scanNumbers) {
        scanPointsMap->insert(scanNumber, m_scanIonsMapped.value(scanNumber));
    }

    ERR_RETURN
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


QMap<ScanNumber, double> MsReaderBase::retentionTimeByScanNumber() const {

    QMap<ScanNumber, double> retentionTimeByScanNumber;

    for (const MsScanInfo &msi : m_msScanInfo) {
        retentionTimeByScanNumber.insert(msi.scanNumber, msi.scanTime);
    }

    return retentionTimeByScanNumber;
}

namespace {

    using NominalMzMass = int;

    void sortTandemScanIons(QVector<TandemScanIon> &tandemScanIons) {

        const double fudgeFactor = 0.0001;
        const auto sortLogic = [fudgeFactor](const TandemScanIon &l, const TandemScanIon &r){

            const int lmz = static_cast<int>(std::round(l.mz / fudgeFactor));
            const int rmz = static_cast<int>(std::round(r.mz / fudgeFactor));

            if (lmz == rmz) {
                return l.precursorTargetMz < r.precursorTargetMz;
            }

            return lmz < rmz;
        };

        std::sort(tandemScanIons.begin(), tandemScanIons.end(), sortLogic);
    }

    QVector<TandemScanIon> sortTandemScanIonsInChunks(QMap<NominalMzMass, QVector<TandemScanIon>> *scanIonsByNominalMass) {

        QVector<TandemScanIon> sortedTandemScanIons;
        for (auto it = scanIonsByNominalMass->begin(); it != scanIonsByNominalMass->end(); it++) {
            QVector<TandemScanIon> &tandemMsIonsAtNominalMzFrag = it.value();
            sortTandemScanIons(tandemMsIonsAtNominalMzFrag);
        }

        return sortedTandemScanIons;
    }

}//namespace
Err MsReaderBase::createTandemScanIonsCache(const QString &cacheFilePath) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isNotEmpty(cacheFilePath); ree;
    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;
    e = ErrorUtils::isNotEmpty(m_scanIons); ree;

    m_scanIonsMapped = scanIonsMapped();

    QMap<NominalMzMass, QVector<TandemScanIon>> scanIonsByNominalMass;

    for (auto it = m_msScanInfo.begin(); it != m_msScanInfo.end(); it++) {

        const ScanNumber scanNumber = it.key();
        const MsScanInfo &msScanInfo = it.value();
        const ScanPoints &scanPoints = m_scanIonsMapped.value(scanNumber);

        for (const ScanPoint &sp : scanPoints) {

            TandemScanIon tsi;
            tsi.scanNumber = scanNumber;
            tsi.collisionEnergy = msScanInfo.collisionEnergy;
            tsi.mz = sp.x();
            tsi.intensity = sp.y();
            tsi.precursorTargetMz = msScanInfo.precursorTargetMz;
            tsi.precursorTargetLowerWindow = msScanInfo.precursorWindowOffsetLower;
            tsi.precursorTargetUpperWindow = msScanInfo.precursorWindowOffsetUpper;
            tsi.scanTime = msScanInfo.scanTime;

            const auto nominalMzMass = static_cast<NominalMzMass>(std::round(tsi.mz));
            scanIonsByNominalMass[nominalMzMass].push_back(tsi);
        }
    }

    sortTandemScanIonsInChunks(&scanIonsByNominalMass);

    for (const QVector<TandemScanIon> &tandemScanIons : scanIonsByNominalMass) {
        m_tandemScanIons.append(tandemScanIons);
    }

    QFile writeFile(cacheFilePath);
    if (!writeFile.open(QIODevice::WriteOnly)) {
        qDebug() << "Could not write to file:" << cacheFilePath
                 << "Error string:" << writeFile.errorString();
        return Error::eFileError;
    }

    QDataStream out(&writeFile);
    out.setVersion(QDataStream::Qt_5_12);
    out << m_tandemScanIons;

    qDebug() << "Ms File Ions written in" << et.elapsed() << "mSec";
    qDebug() << "Ms FileIons library written to:" << cacheFilePath;

    ERR_RETURN
}

Err MsReaderBase::readFromCache(const QString cacheFileURI) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    m_tandemScanIons.clear();

    QFile readFile(cacheFileURI);
    QDataStream in(&readFile);

    in.setVersion(QDataStream::Qt_5_12);

    if (!readFile.open(QIODevice::ReadOnly)) {
        qDebug() << "Could not read the file:" << cacheFileURI
                 << "Error string:" << readFile.errorString();
        return Error::eFileError;
    }

    in >> m_tandemScanIons;

    e = ErrorUtils::isNotEmpty(m_tandemScanIons); ree;

    qDebug() << "Tandem Scan Ions loaded from" << cacheFileURI
             << "in" <<  et.elapsed() << "mSec";

    ERR_RETURN
}

Err MsReaderBase::buildUniqueTandemScanIons() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_tandemScanIons); ree;

    for (int i = 0; i < m_tandemScanIons.size(); i++) {

        const TandemScanIon &tsi = m_tandemScanIons.at(i);
        m_uniqueTandemScanIons[tsi.hashedKey()].push_back(i);
    }

    ERR_RETURN
}

bool MsReaderBase::cacheExists(const QString cacheFileURI) {

    QFileInfo fi(cacheFileURI);
    return fi.exists();
}
