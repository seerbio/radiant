//
// Created by Drucifer on 3/20/2022.
//

#include "MsReaderBase.h"

#include "ErrorUtils.h"

#include <QElapsedTimer>
#include <QFile>
#include <QtConcurrent/QtConcurrent>

#include <iostream>


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

Err MsReaderBase::getScanInfo(
        ScanNumber scanNumber,
        MsScanInfo *msScanInfo
        ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msScanInfo.contains(scanNumber)); ree;
    *msScanInfo = m_msScanInfo.value(scanNumber);

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

Err MsReaderBase::openFile(
        const QString &filePath,
        bool useCache
        ) {
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

QString MsReaderBase::buildUniqueTandemWindowKey(const MsScanInfo &si) {

    QString key = QString::number(si.precursorTargetMz - si.precursorWindowOffsetLower)
                + MsReaderBase::separator()
                + QString::number(si.precursorTargetMz + si.precursorWindowOffsetUpper);

    return key;
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
Err MsReaderBase::createMsReaderCache(const QString &cacheFilePath) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    e = ErrorUtils::isNotEmpty(cacheFilePath); ree;
    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;
    e = ErrorUtils::isNotEmpty(m_scanIons); ree;

    QFile writeFile(cacheFilePath);
    if (!writeFile.open(QIODevice::WriteOnly)) {
        qDebug() << "Could not write to file:" << cacheFilePath
                 << "Error string:" << writeFile.errorString();
        return Error::eFileError;
    }

    QDataStream out(&writeFile);
    out.setVersion(QDataStream::Qt_5_12);
    out << m_scanIons;
    out << m_msScanInfo;

    qDebug() << "Ms File Ions written in" << et.elapsed() << "mSec";
    qDebug() << "Ms FileIons library written to:" << cacheFilePath;

    ERR_RETURN
}

Err MsReaderBase::readFromCache(const QString &cacheFileURI) {

    ERR_INIT

    QElapsedTimer et;
    et.start();

    m_scanIons.clear();
    m_msScanInfo.clear();

    QFile readFile(cacheFileURI);
    QDataStream in(&readFile);

    in.setVersion(QDataStream::Qt_5_12);

    if (!readFile.open(QIODevice::ReadOnly)) {
        qDebug() << "Could not read the file:" << cacheFileURI
                 << "Error string:" << readFile.errorString();
        return Error::eFileError;
    }

    in >> m_scanIons;
    in >> m_msScanInfo;

    e = ErrorUtils::isNotEmpty(m_scanIons); ree;
    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;

    qDebug() << "Tandem Scan Ions loaded from" << cacheFileURI
             << "in" <<  et.elapsed() << "mSec";

    ERR_RETURN
}

bool MsReaderBase::cacheExists(const QString &cacheFileURI) {

    QFileInfo fi(cacheFileURI);
    return fi.exists();
}

Err MsReaderBase::tandemScanIons(QVector<TandemScanIon> *tandemScanIons) {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_scanIons); ree;
    e = ErrorUtils::isNotEmpty(m_msScanInfo); ree;

    QMap<NominalMzMass, QVector<TandemScanIon>> scanIonsByNominalMass;

    MsScanInfo msScanInfo;
    for (const ScanIon &si : m_scanIons) {

        if (msScanInfo.scanNumber != si.scanNumber) {
            msScanInfo = m_msScanInfo.value(si.scanNumber);
            e = ErrorUtils::isNotEqual(msScanInfo.scanNumber, -1); ree;
        }

        if (msScanInfo.msLevel == 1) {
            continue;
        }

        TandemScanIon tsi;
        tsi.scanNumber = si.scanNumber;
        tsi.mz = si.mz;
        tsi.intensity = si.intensity;
        tsi.precursorTargetMz = msScanInfo.precursorTargetMz;
        tsi.precursorTargetLowerWindow = msScanInfo.precursorWindowOffsetLower;
        tsi.precursorTargetUpperWindow = msScanInfo.precursorWindowOffsetUpper;

        const auto nominalMzMass = static_cast<NominalMzMass>(std::round(tsi.mz));
        scanIonsByNominalMass[nominalMzMass].push_back(tsi);

    }

    sortTandemScanIonsInChunks(&scanIonsByNominalMass);

    for (const QVector<TandemScanIon> &tsi : scanIonsByNominalMass) {
        tandemScanIons->append(tsi);
    }

    ERR_RETURN
}


Err MsReaderBase::sortDIATandemScansByMzTarget(
        const QVector<TandemScanIon> &tandemScanIons,
        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaFrames
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(tandemScanIons); ree;

    for (const TandemScanIon &tsi : tandemScanIons) {
        (*diaFrames)[tsi.targetScanKey()][tsi.scanNumber].push_back({tsi.mz, tsi.intensity});
    }

    ERR_RETURN
}
