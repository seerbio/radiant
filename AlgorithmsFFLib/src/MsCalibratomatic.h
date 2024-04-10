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

    /**
    * @brief Initializes MsCalibratomatic with the provided MS calibration file path.
    *
    * This function initializes MsCalibratomatic by checking the existence of the specified MS calibration file,
    * reading its content using ParquetReader, and performing further initialization based on the read data.
    *
    * @param msCalibrationFilePath A QString representing the file path of the MS calibration file.
    * @return Err indicating the success or failure of the initialization process.
    */
    Err init(const QString &msCalibrationFilePath);

    /**
    * @brief Initializes MsCalibratomatic with the provided MS calibration data.
    *
    * This function initializes MsCalibratomatic using a QVector of MsCalibarationReaderRow, which represents the
    * calibration data read from a source. It ensures that the input data is not empty, sets the internal data
    * member m_msCalibarationReaderRows, builds the internal IRT and mz calibrators, and checks that the computed
    * standard deviations for mz and scan time are greater than 0.0. Finally, it sets the isInit flag to true.
    *
    * @param msCalibarationReaderRows A QVector of MsCalibarationReaderRow containing MS calibration data.
    * @return Err indicating the success or failure of the initialization process.
    */
    Err init(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows);

    /**
    * @brief Initialize the calibration using the provided calibration reader rows.
    *
    * @param msCalibarationReaderRows The calibration reader rows to initialize the calibration.
    * @return An error code indicating the success or failure of the initialization.
    */
    Err initRtOnly(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows);

    /**
    * @brief Initialize the calibration using only the provided Mz calibration reader rows.
    *
    * @param msCalibarationReaderRows The Mz calibration reader rows to initialize the calibration.
    * @return An error code indicating the success or failure of the initialization.
    */
    Err initMzOnly(const QVector<MsCalibarationReaderRow> &msCalibarationReaderRows);

    /**
    * @brief Recalibrate the given scan points using the calibration model.
    *
    * This function recalibrates the provided scan points using the calibration model for Mz values.
    *
    * @param scanNumberVsScanPoints A QMap containing scan numbers as keys and corresponding scan points as values.
    * @return An error code indicating the success or failure of the recalibration process.
    */
    Err recalibrateScanPoints(
            const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints
            );

    /**
    * @brief Recalibrate the given scan points using the calibration model.
    *
    * This function recalibrates the provided scan points using the calibration model for Mz values.
    *
    * @param scanNumberVsScanPoints A QMap containing scan numbers as keys and corresponding scan points as values.
    * @return An error code indicating the success or failure of the recalibration process.
    */
    Err recalibrateScanPoints(QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints);

    /**
    * @brief Get the standard deviation of Mz values used in the calibration.
    *
    * This function returns the standard deviation of Mz values used during the calibration process.
    *
    * @return The standard deviation of Mz values.
    */
    [[nodiscard]] float mzStDev();

    /**
    * @brief Get the standard deviation of scan time values used in the calibration.
    *
    * This function calculates and returns the standard deviation of scan time values used during the calibration process.
    *
    * @param nStdDevs Number of standard deviations to multiply with the calculated standard deviation.
    * @return The calculated standard deviation of scan time values.
    */
    [[nodiscard]] float scanTimeStDev(int nStdDevs = 1);

    /**
    * @brief Set the standard deviation of scan time values used in the calibration.
    *
    * This function sets the standard deviation of scan time values used during the calibration process and prints the updated values to the console.
    *
    * @param val The new value for the standard deviation of scan time.
    */
    void setScanTimeStDev(double val);

    /**
    * @brief Predict scan time for a given iRT value.
    *
    * This function predicts the scan time for a given iRT (indexed retention time) value using the calibrated mapping and returns the result.
    *
    * @param iRT The input iRT value for which to predict the scan time.
    * @param predictedScanTime A pointer to store the predicted scan time result.
    * @return Error code indicating the success or failure of the prediction process.
    */
    Err predictScanTime(
        float iRT,
        float *predictedScanTime
        ) const;

    /**
    * @brief Check if the calibration has been initialized.
    *
    * This function checks whether the calibration has been successfully initialized.
    *
    * @return True if the calibration is initialized, false otherwise.
    */
    [[nodiscard]] bool isInit() const;

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

    QVector<double> m_mzToRecalCalCurveCoeffs;

    bool m_isInit;

};


#endif //PYTHIADIACPP_MSCALIBRATOMATIC_H
