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

    [[nodiscard]] MzTargetKey mzTargetKey() const {
        return mzTargetKey(
                precursorTargetMz - isoWindowLower,
                precursorTargetMz + isoWindowUpper
                );
    }

    static MzTargetKey mzTargetKey(double mzStart, double mzEnd) {
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
            int msLevel,
            QMap<ScanNumber, ScanPoints*> *scanPoints
    );

    QPair<Err, ScanPoints*> getScanPoints(int scanNumber);

    Err getMsScanInfo(
            ScanNumber scanNumber,
            MsScanInfo *msScanInfo
            ) const;

    Err collateTandemPrecursorTargetsDIA(
            QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrame
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

    static Err getNearestScanNumberFromScanTime(
            ScanTime scanTime,
            const QVector<ScanNumber> &scanNumbers,
            const QVector<ScanTime> &scanTimes,
            ScanNumber *scanNumber
            );

    Err getNearestScanNumberFromScanNumber(
            ScanNumber scanNumber,
            const QVector<ScanNumber> &scanNumbers,
            ScanNumber *closestScanNumber
            );

    [[nodiscard]] QMap<ScanNumber, ScanTime> getScanNumberVsScanTime() const;

    static Err splitScanPoints(
            const ScanPoints &scanPoints,
            QVector<double> *mzVals,
            QVector<float> *intensityVals
            );

    static Err zipScanPoints(
            const QVector<double> &mzVals,
            const QVector<float> &intensityVals,
            ScanPoints *scanPoints
            );

    void printSize();

    Err printFileInfo();

    Err writeTargetCollisionEnergyFile();

protected:

    bool m_fileIsCalibrated;

    QMap<ScanNumber, MsScanInfo> m_msScanInfo;
    QMap<ScanNumber, ScanPoints>  m_scanPoints;

    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;

    QString m_filePath;

};


#endif //MSREADERBASE_H
