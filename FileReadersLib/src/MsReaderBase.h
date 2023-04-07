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
        return QString::number(std::round(precursorTargetMz * 1000));
    }
};


class FILEREADERSLIB_EXPORTS MsReaderBase {

public:

    friend class MsReaderBaseTests;
    friend class MsReaderMZMLTests;
    friend class MsReaderParquetTests;

    MsReaderBase() = default;

    ~MsReaderBase() = default;

    virtual Err openFile(const QString &filePath);

    virtual Err closeFile();

    QMap<ScanNumber, ScanPoints> getScanPoints();
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints(int msLevel);
    Err getScanPoints(
            int scanNumber,
            ScanPoints *scanPoints
            );

    Err getScanInfo(
            ScanNumber scanNumber,
            MsScanInfo *msScanInfo
            ) const;

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

protected:

    QMap<ScanNumber, MsScanInfo> m_msScanInfo;
    QMap<ScanNumber, ScanPoints>  m_scanPoints;

};


#endif //MSREADERBASE_H
