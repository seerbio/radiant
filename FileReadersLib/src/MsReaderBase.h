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


enum class ScanPointsSort {
    AscMz,
    AscIntensity,
    DescMz,
    DescIntensity
};

struct FILEREADERSLIB_EXPORTS MsScanInfo {

    int msLevel = -1;
    ScanNumber scanNumber = 1;
    double scanTime = -1.0;
    double collisionEnergy = -1.0;
    double precursorTargetMz = -1.0;
    double isoWindowLower = -1.0;
    double isoWindowUpper = -1.0;
    double ionMobilityDriftTime = -1.0;
    IonMobilityIndex ionMobilityIndex = -1;

    [[nodiscard]] QString targetScanKey() const {
        return targetScanKey(
                precursorTargetMz - isoWindowLower,
                precursorTargetMz + isoWindowUpper
                );
    }

    static QString targetScanKey(double mzStart, double mzEnd) {
        return QString::number(std::round(1000 * ((mzStart + mzEnd) / 2)));
    }
};


class FILEREADERSLIB_EXPORTS MsReaderBase {

public:

    friend class MsReaderBaseTests;
    friend class MsReaderMZMLTests;
    friend class MsReaderParquetTests;

    MsReaderBase();

    ~MsReaderBase() = default;

    void setMsScanInfo(const QMap<ScanNumber, MsScanInfo> &msScanInfos);

    Err setScanPoints(const QMap<ScanNumber, ScanPoints> &scanPoints);

    virtual Err openFile(const QString &filePath);

    virtual Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy,
            const QPair<double, double> &filterRange
            );

    virtual Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy
    );

    virtual Err closeFile();

    QString filePath();

    bool isDIA();

    bool isInit();

    QPair<double, double> scanTimeMinMax();

    QMap<ScanNumber, ScanPoints> getScanPoints();

    Err getScanPoints(
            int msLevel,
            QMap<ScanNumber, ScanPoints> *scanPoints
            );

    Err getScanPoints(
            int scanNumber,
            ScanPoints *scanPoints
            );

    Err getMsScanInfo(
            ScanNumber scanNumber,
            MsScanInfo *msScanInfo
            ) const;

    Err collateTandemPrecursorTargetsDIA(
            QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrame
            );

    QVector<MsScanInfo> getUniqueTandemMsScanInfos();

    int getFrameCount();

    Err updateScanPoints(const QMap<ScanNumber, ScanPoints> &scansToUpdate);

    QMap<ScanNumber, MsScanInfo> getMsScanInfos(int msLevel);
    QMap<ScanNumber, MsScanInfo> getMsScanInfos();
    Err getMsScanInfo(
            ScanNumber scanNumber,
            MsScanInfo *msScanInfo
            );

    static ScanPoints sortScanPoints(
            const ScanPoints &scanPoints,
            const ScanPointsSort &sort = ScanPointsSort::AscMz
                    );

    int getNearestScanNumberFromScanTime(double scanTime);

    int getNearestScanNumberFromScanNumber(int scanNumber);

    [[nodiscard]] QMap<ScanNumber, ScanTime> getScanNumberVsScanTime() const;

    static Err splitScanPoints(
            const ScanPoints &scanPoints,
            QVector<double> *mzVals,
            QVector<double> *intensityVals
            );

    static Err zipScanPoints(
            const QVector<double> &mzVals,
            const QVector<double> &intensityVals,
            ScanPoints *scanPoints
            );

    void printSize();

    Err printFileInfo();

protected:

    bool m_fileIsCalibrated;

    QMap<ScanNumber, MsScanInfo> m_msScanInfo;
    QMap<ScanNumber, ScanPoints>  m_scanPoints;

    QVector<ScanNumber> m_scanNumbers;
    QVector<ScanTime> m_scanTimes;

    QString m_filePath;

};


#endif //MSREADERBASE_H
