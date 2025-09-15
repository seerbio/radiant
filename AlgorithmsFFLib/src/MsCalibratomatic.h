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

enum class MSLevelEnum {
    MS1,
    MS2
};

class PeptideSequence;

class ALGORITHMSFFLIB_EXPORTS MsCalibratomatic {

public:

    MsCalibratomatic();
    ~MsCalibratomatic() = default;

    Err buildRTMapper(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows);
    Err buildIMMapper(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows);

    Err setMassCalibrationCoeffs(
        const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows,
        const MSLevelEnum &msLevelEnum
        );

    [[nodiscard]] Err recalibrateScanPoints(
        const MSLevelEnum &msLevel,
        const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints
        ) const;

    Err recalibrateScanPoints(
        const MSLevelEnum &msLevel,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
        ) const;

    Err recalibrateScanPoint(
        const MSLevelEnum &msLevel,
        float mzVal,
        float *mzValRecal
        ) const;

    [[nodiscard]] float mzStDevMS1() const;

    [[nodiscard]] float mzStDevMS2() const;

    [[nodiscard]] float scanTimeMean(float nStdDevs = 1.0f) const;
    [[nodiscard]] float scanTimeStDev(float nStdDevs = 1.0f) const;
    [[nodiscard]] float ionMobilityStDev(float nStdDevs = 1.0f) const;

    void setScanTimeStDev(double val);
    void setIonMobilityStDev(double val);
    void setMzStDevMS2(double val);

    Err predictScanTime(
        float iRT,
        float *predictedScanTime
        ) const;

    [[nodiscard]] bool isInitRT() const;

    [[nodiscard]] bool isInitCalMS1() const;

    [[nodiscard]] bool isInitCalMS2() const;

    Err setCalibrationCoeffsUsingAllMeans();

private:

    PythiaParameters m_params;
    double m_mzStDevMS1;
    double m_mzStDevMS2;
    double m_scanTimeMean;
    double m_scanTimeStd;
    double m_ionMobilityMean;
    double m_ionMobilityStd;

    XYMappermatic m_iRTtoScanTimeMapper;
    XYMappermatic m_iIMtoScanTimeMapper;
    QVector<double> m_calibrationCurveCoeffsMS1;
    QVector<double> m_calibrationCurveCoeffsMS2;

    QVector<QVector<double>> m_calibrationCurveCoEffsAll;

    bool m_isInitRT;
    bool m_isInitIM;
    bool m_isInitMS1;
    bool m_isInitMS2;

    int m_polynomialOrderMassCal;

};

#endif //PYTHIADIACPP_MSCALIBRATOMATIC_H
