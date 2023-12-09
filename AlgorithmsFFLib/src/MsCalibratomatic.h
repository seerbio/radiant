//
// Created by anichols on 4/16/23.
//

#ifndef PYTHIADIACPP_MSCALIBRATOMATIC_H
#define PYTHIADIACPP_MSCALIBRATOMATIC_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsCalibrationReader.h"
#include "PythiaParameterReader.h"
#include "XYMappermatic.h"

using namespace Error;

class PeptideSequence;

class ALGORITHMSFFLIB_EXPORTS MsCalibratomatic {

public:

    MsCalibratomatic();
    ~MsCalibratomatic() = default;

    Err init(const QString &msCalibrationFilePath);

    Err init(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows);

    // either FrameIndex, or ScanNumber can be key as they are both ints.
    Err recalibrateScanPoints(
            const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints
            );

    [[nodiscard]] double mzStDev();
    [[nodiscard]] double scanTimeStDev();

    Err predictScanTime(
            double iRT,
            double *predictedScanTime
            ) const;

    bool isInit() const;

private:

    Err buildMzCalibrator();
    Err buildIRTCalibrator();

private:

    //Never cleared
    PythiaParameters m_params;
    double m_mzStDev;
    double m_scanTimeStd;
    QString m_msCalibrationFilePath;

    QVector<MsCalibarationReaderRow> m_msCalibarationReaderRows;

    XYMappermatic m_iRTtoScanTimeMapper;
    XYMappermatic m_mzToRecalMz;

    bool m_isInit;

};


#endif //PYTHIADIACPP_MSCALIBRATOMATIC_H
