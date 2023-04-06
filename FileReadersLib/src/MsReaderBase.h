//
// Created by Drucifer on 3/20/2022.
//

#ifndef MSREADERBASE_H
#define MSREADERBASE_H

#include "Error.h"
#include "FileReadersLib_Exports.h"

#include "GlobalSettings.h"
#include "MathUtils.h"

#include <QDataStream>


using namespace Error;

enum TIME_UNITS {
    MINUTES = 0,
    SECONDS,
    MILLISECONDS
};

enum class ScanPointsSort {
    AscMz,
    AscIntensity,
    DescMz,
    DescIntensity
};

struct FILEREADERSLIB_EXPORTS TandemScanIon {

    ScanNumber scanNumber = -1;
    double mz = -1.0;
    double intensity = -1.0;
    double precursorTargetMz = -1.0;
    double precursorTargetLowerWindow = -1.0;
    double precursorTargetUpperWindow = -1.0;

    [[nodiscard]] QString targetScanKey() const {
        return QString::number(std::round(precursorTargetMz * 1000));
    }

};

struct FILEREADERSLIB_EXPORTS MsScanInfo {

    ScanNumber scanNumber = -1;
    MsLevel msLevel = -1;
    double scanTime = -1.0;
    double TIC = -1.0;
    double basePeakIntensity = -1.0;

    double precursorTargetMz = -1.0;
    int precursorTargetCharge = -1;
    double precursorWindowOffsetLower = -1.0;
    double precursorWindowOffsetUpper = -1.0;
    int collisionEnergy = -1;
    int faimsVoltage = -1;

    [[nodiscard]] QString msScanInfoKey() const {

        if (msLevel < 2) {
            return "MS1";
        }
        const double mzStart = precursorTargetMz - precursorWindowOffsetLower;
        const double mzEnd = precursorTargetMz + precursorWindowOffsetUpper;

        QString key = QString::number(mzStart)
                        + S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP
                        + QString::number(mzEnd);

        return key;
    }

    friend QDataStream &operator <<(QDataStream &stream, const MsScanInfo &msi) {
        stream << msi.scanNumber;
        stream << msi.msLevel;
        stream << msi.scanTime;
        stream << msi.TIC;
        stream << msi.basePeakIntensity;
        stream << msi.precursorTargetMz;
        stream << msi.precursorTargetCharge;
        stream << msi.precursorWindowOffsetLower;
        stream << msi.precursorWindowOffsetUpper;
        stream << msi.collisionEnergy;
        stream << msi.faimsVoltage;

        return stream;
    }

    friend QDataStream &operator >>(QDataStream &stream, MsScanInfo &msi) {

        stream >> msi.scanNumber;
        stream >> msi.msLevel;
        stream >> msi.scanTime;
        stream >> msi.TIC;
        stream >> msi.basePeakIntensity;
        stream >> msi.precursorTargetMz;
        stream >> msi.precursorTargetCharge;
        stream >> msi.precursorWindowOffsetLower;
        stream >> msi.precursorWindowOffsetUpper;
        stream >> msi.collisionEnergy;
        stream >> msi.faimsVoltage;

        return stream;
    }

};

struct FILEREADERSLIB_EXPORTS ScanIon {

    ScanNumber scanNumber = -1;
    double mz = -1.0;
    double intensity = -1.0;

    ScanIon() = default;

    ScanIon(int scanNumber,
            double mz,
            double intensity)
            : scanNumber(scanNumber)
            , mz(mz)
            , intensity(intensity) {}

    friend QDataStream &operator <<(QDataStream &stream, const ScanIon &si) {
        stream << si.scanNumber;
        stream << si.mz;
        stream << si.intensity;

        return stream;
    }

    friend QDataStream &operator >>(QDataStream &stream, ScanIon &si) {

        stream >> si.scanNumber;
        stream >> si.mz;
        stream >> si.intensity;

        return stream;
    }

};

class FILEREADERSLIB_EXPORTS MsReaderBase {

public:

    friend class MsReaderBaseTests;
    friend class MsReaderMZMLTests;

    MsReaderBase() = default;

    ~MsReaderBase() = default;

    virtual Err openFile(const QString &filePath);

    virtual Err closeFile();

    static QChar separator() {return S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP;}

    static Err sortDIATandemScansByMzTarget(
            const QVector<TandemScanIon> &tandemScanIons,
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaFrames
            );

    QMap<ScanNumber, bool> buildSkipNAllowableIndexesMs2(int skipEveryNScans);

    Err getScanInfo(
            ScanNumber scanNumber,
            MsScanInfo *msScanInfo
            ) const;

    QVector<MsScanInfo> getMsScanInfos(int msLevel);

    Err buildScanPointsByScanNumber(QMap<ScanNumber, ScanPoints> *scanPointsByScanNumber);

    Err buildScanPointsByScanNumber(
            MsLevel msLevel,
            QMap<ScanNumber, ScanPoints> *scanPointsByScanNumber
    );

    static ScanPoints sortScanPoints(
            const ScanPoints &scanPoints,
            const ScanPointsSort &sort = ScanPointsSort::AscMz
                    );

    int getNearestScanNumberFromScanTime(double scanTime);

    [[nodiscard]] QMap<ScanNumber, double> retentionTimeByScanNumber() const;

    static QString buildUniqueTandemWindowKey(const MsScanInfo &si);

    static Err splitScanPoints(
            const ScanPoints &scanPoints,
            QVector<double> *mzVals,
            QVector<double> *intensityVals
    );

    Err tandemScanIons(QVector<TandemScanIon> *tandemScanIons);

    static QVector<TandemScanIon> sortTandemScanIonsInChunks(
            QMap<NominalMzMass, QVector<TandemScanIon>> *scanIonsByNominalMass
            );

protected:

    QMap<ScanNumber, MsScanInfo> m_msScanInfo;
    QVector<ScanIon> m_scanIons;

};


#endif //MSREADERBASE_H
