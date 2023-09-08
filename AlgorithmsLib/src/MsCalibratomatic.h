//
// Created by anichols on 4/16/23.
//

#ifndef PYTHIADIACPP_MSCALIBRATOMATIC_H
#define PYTHIADIACPP_MSCALIBRATOMATIC_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsUtils.h"
#include "PythiaParameterReader.h"
#include "XYMappermatic.h"

using namespace Error;

class PeptideSequence;

class ALGORITHMSLIB_EXPORTS MsCalibratomatic {

public:

    MsCalibratomatic();
    ~MsCalibratomatic() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &msCalibrationFilePath
            );

    // either FrameIndex, or ScanNumber can be key as they are both ints.
    Err recalibrateScanPoints(
            const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints,
            QMap<ScanNumber, ScanPoints> *recalScanNumberVsScanPoints
            );

    [[nodiscard]] double newStDev();

    Err setScanTimeWindowNew(double scanTimeWindow);

private:

    Err buildMzCalibrator();
    Err buildIRTCalibrator();

private:

    //Never cleared
    PythiaParameters m_params;
    double m_stDevNew;
    double m_scanTimeWindowNew;
    QString m_msCalibrationFilePath;

    XYMappermatic m_iRTtoScanTimeMapper;
    XYMappermatic m_mzToRecalMz;

};


#endif //PYTHIADIACPP_MSCALIBRATOMATIC_H
