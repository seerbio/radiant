//
// Created by Drucifer on 3/20/2022.
//

#include "MsReaderBase.h"

#include "ErrorUtils.h"

#include <QElapsedTimer>
#include <QFile>
#include <QtConcurrent/QtConcurrent>

#include <iostream>


Err MsReaderBase::getScanInfo(
        ScanNumber scanNumber,
        MsScanInfo *msScanInfo
        ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msScanInfo.contains(scanNumber)); ree;
    *msScanInfo = m_msScanInfo.value(scanNumber);

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

QMap<ScanNumber, ScanPoints> MsReaderBase::scanNumberVsScanPoints(int msLevel) {

    QMap<ScanNumber, ScanPoints> pointsOfReturn;


    return {};
}
