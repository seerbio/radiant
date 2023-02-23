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

    friend QDataStream &operator <<(QDataStream &stream, const TandemScanIon &tsi) {
        stream << tsi.scanNumber;
        stream << tsi.mz;
        stream << tsi.intensity;
        stream << tsi.precursorTargetMz;
        stream << tsi.precursorTargetLowerWindow;
        stream << tsi.precursorTargetUpperWindow;

        return stream;
    }

    friend QDataStream &operator >>(QDataStream &stream, TandemScanIon &tsi) {

        stream >> tsi.scanNumber;
        stream >> tsi.mz;
        stream >> tsi.intensity;
        stream >> tsi.precursorTargetMz;
        stream >> tsi.precursorTargetLowerWindow;
        stream >> tsi.precursorTargetUpperWindow;

        return stream;
    }

    [[nodiscard]] QString hashedKey() const {

        const int hashingPrecision = 3;


        const int mzHashed = MathUtils::hashDecimal(mz, hashingPrecision);
        const int targetHashed = MathUtils::hashDecimal(precursorTargetMz, hashingPrecision);

        return S_GLOBAL_SETTINGS.MODIFICATION_STRING_FORMAT
                .arg(mzHashed)
                .arg(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP)
                .arg(targetHashed);
    }

    static QPair<MZION , TARGETMZ> unhashKey(const QString &hashedKey) {
        const QStringList splitStr = hashedKey.split(
                S_GLOBAL_SETTINGS.MODIFICATION_STRING_FORMAT,
                QString::SkipEmptyParts
                );

        const int hashingPrecision = 3;

        return {MathUtils::unHashDecimal<double>(splitStr.first().toInt(), hashingPrecision),
                MathUtils::unHashDecimal<double>(splitStr.back().toInt(), hashingPrecision)};
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

    QString msScanInfoKey() const {

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

};

struct ScanIonWithScanInfo {
    ScanIon scanIon;
    MsScanInfo msScanInfo;
};


class FILEREADERSLIB_EXPORTS MsReaderBase {

public:

    friend class MsReaderBaseTests;
    friend class MsReaderMZMLTests;

    MsReaderBase() = default;

    ~MsReaderBase() = default;

    virtual Err openFile(
            const QString &filePath,
            bool useCache = true
            );

    virtual Err closeFile();

    static QChar separator() {return S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP;}

    Err  buildScanIonWithScanInfoInputMs2(const QPair<double, double> &mzMinMax,
                                          int skipEveryNScans,
                                          QVector<ScanIonWithScanInfo> *scanIonWithScanInfos);

    QMap<ScanNumber, bool> buildSkipNAllowableIndexesMs2(int skipEveryNScans);

    Err getScanInfo(ScanNumber scanNumber,  MsScanInfo *msScanInfo) const;

    Err getScanPoints(ScanNumber scanNumber,  QVector<QPointF> *scanPoints) const;

    QVector<MsScanInfo> getMsScanInfos(int msLevel);

    Err buildScanPointsByScanNumber(QMap<ScanNumber, ScanPoints> *scanPointsByScanNumber);

    Err buildScanPointsByScanNumber(
            MsLevel msLevel,
            QMap<ScanNumber, ScanPoints> *scanPointsByScanNumber
    );

    static ScanPoints sortScanPoints(const ScanPoints &scanPoints,
                                     const ScanPointsSort &sort = ScanPointsSort::AscMz);

    Err createTandemScanIonsCache(const QString &cacheFilePath);

    Err readFromCache(const QString &cacheFileURI);

    static bool cacheExists(const QString mzMLFileURI);


    static Err getRangesFromIsolationWindowKey(
            const IsolationWindowKey &key,
            QPair<double, double> *mzIsoRange
    );

    void loadScanIonsMapped();

    QMap<ScanNumber, QVector<QPointF>> scanIonsMapped();

    Err getBasePeakChromatogram(
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            QMap<ScanNumber, QPointF> *basePeakChromatogram
    );

    int getNearestScanNumberFromScanTime(double scanTime);

    QVector<QPointF> getNearestTandemScan(double scanTime);

    QMap<IsolationWindowKey, QVector<ScanNumber>> getAllUniqueTandemIsolationWindows();

    QMap<ScanNumber, double> retentionTimeByScanNumber() const;

    QStringList getUniqueMsInfoScanKeys();

    static QString buildUniqueTandemWindowKey(const MsScanInfo &si);

    Err getScansByUniqueMsInfoScanKey(
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            QMap<ScanNumber, ScanPoints> *scanPointsMap
    );

    static Err splitScanPoints(
            const ScanPoints &scanPoints,
            QVector<double> *mzVals,
            QVector<double> *intensityVals
    );

    Err buildUniqueTandemScanIons();

    QVector<TandemScanIon> tandemScanIons();

protected:

    Err buildScanNumbersByUniqueMsInfoScanKey();

protected:

    QMap<ScanNumber, MsScanInfo> m_msScanInfo;
    QMap<ScanNumber, ScanPoints> m_scanIonsMapped;
    QMap<UniqueMsInfoScanKey, QVector<ScanNumber>> m_scanNumbersByUniqueMsInfoScanKey;

    QVector<ScanIon> m_scanIons;
    QVector<TandemScanIon> m_tandemScanIons;
    QMap<UniqueHashedMzAndTarget, QVector<TandemScansIndex>> m_uniqueTandemScanIons;

};


#endif //MSREADERBASE_H
