 //
// Created by anichols on 4/8/23.
//

#ifndef PYTHIADIACPP_MSFRAME_H
#define PYTHIADIACPP_MSFRAME_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"
#include "ParquetReader.h"
#include "PythiaParameterReader.h"

using namespace Error;

namespace MsFrameScanPointRowsNamespace {
    const QString FRAME_INDEX = QStringLiteral("frameIndex");
    const QString MZ_VALS = QStringLiteral("mzVals");
    const QString INTENSITY_VALS = QStringLiteral("intensityVals");

    const QStringList keysToCheck = {
         FRAME_INDEX,
         MZ_VALS,
         INTENSITY_VALS
    };
}//namespace

/**
 * See ParquetReaderInputBase for documentation
 */
struct ALGORITHMSFFLIB_EXPORTS MsFrameScanPointRows : public ParquetReaderInputBase {

    FrameIndex frameIndex = -1;
    QVector<float> mzVals;
    QVector<float> intensityVals;

    QMap<QString, QVariant> map() override {

        using namespace MsFrameScanPointRowsNamespace;

        return {
            {FRAME_INDEX, QVariant(frameIndex)},
            {MZ_VALS, QVariant(qVectorToQByteArray(mzVals))},
            {INTENSITY_VALS, QVariant(qVectorToQByteArray(intensityVals))}
        };
    }

    Err initFromRead(const ParquetReaderInputBase &row) override {

        using namespace MsFrameScanPointRowsNamespace;

        ERR_INIT

        const QMap<QString, QVariant> &dataMap = row.dataMap();
        const bool allKeysPresent = ParquetReaderInputBase::checkIfExpectedKeysArePresent(
                dataMap,
                keysToCheck
        );

        e = ErrorUtils::isTrue(allKeysPresent); ree;

        frameIndex = dataMap.value(FRAME_INDEX).toInt();
        mzVals = bytesArrayToQVector<float>(dataMap.value(MZ_VALS).toByteArray());
        intensityVals = bytesArrayToQVector<float>(dataMap.value(INTENSITY_VALS).toByteArray());

        ERR_RETURN
    }

};


class ALGORITHMSFFLIB_EXPORTS MsFrame {

public:

    friend class MsFrameTests;

    MsFrame();
    ~MsFrame();

    /**
    * @brief Initialize the MsFrame with scan points and scan times.
    *
    * This function initializes the MsFrame with the provided scan points and scan times.
    *
    * @param scanPoints A map containing scan numbers and corresponding scan points.
    * @param scanNumberVsScanTime A map containing scan numbers and corresponding scan times.
    * @return An error code indicating the success or failure of the initialization.
    */
    Err init(
            const QMap<ScanNumber, ScanPoints*> &scanPoints,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
            );

    /**
    * @brief Check if the MsFrame is valid.
    *
    * This function checks if the MsFrame is valid by ensuring it contains an adequate number of scans.
    *
    * @return True if the MsFrame is valid, false otherwise.
    */
    [[nodiscard]] bool isValid() const ;

    /**
    * @brief Write MsFrame scans to a Parquet file.
    *
    * This function writes the scan points associated with each frame to a Parquet file.
    *
    * @param framesVsScanPoints A QMap containing FrameIndex as keys and ScanPoints as values.
    * @param outputFilePath The path to the output Parquet file.
    * @return An Err code indicating the success or failure of the operation.
    */
    static Err writeFrameScans(
            const QMap<FrameIndex, ScanPoints> &framesVsScanPoints,
            const QString &outputFilePath
            );

    /**
    * @brief Write MsFrame scans to a Parquet file.
    *
    * This function writes the scan points associated with each frame to a Parquet file.
    *
    * @param framesVsScanPoints A QMap containing FrameIndex as keys and pointers to ScanPoints as values.
    * @param outputFilePath The path to the output Parquet file.
    * @return An Err code indicating the success or failure of the operation.
    */
    static Err writeFrameScans(
            const QMap<FrameIndex, ScanPoints*> &framesVsScanPoints,
            const QString &outputFilePath
    );

    /**
    * @brief Get the total number of scans in the MsFrame.
    *
    * This function returns the count of scans in the MsFrame.
    *
    * @return An integer representing the number of scans in the MsFrame.
    */
    [[nodiscard]] int scanCount() const;

    /**
    * @brief Get a map of FrameIndex to ScanPoints pointers.
    *
    * This function returns a QMap containing the mapping of FrameIndex to ScanPoints pointers
    * for all the scans in the MsFrame.
    *
    * @return A QMap<FrameIndex, ScanPoints*> representing the mapping of FrameIndex to ScanPoints pointers.
    */
    [[nodiscard]] QMap<FrameIndex, ScanPoints*> frameIndexVsScanPoints() const;

    /**
    * @brief Get a map of ScanNumber to ScanPoints pointers.
    *
    * This function returns a QMap containing the mapping of ScanNumber to ScanPoints pointers
    * for all the scans in the MsFrame.
    *
    * @return A QMap<ScanNumber, ScanPoints*> representing the mapping of ScanNumber to ScanPoints pointers.
    */
    [[nodiscard]] QMap<ScanNumber, ScanPoints*> scanNumberVsScanPoints() const;

    /**
    * @brief Get the ScanNumber corresponding to a given FrameIndex.
    *
    * This function retrieves the ScanNumber associated with the provided FrameIndex.
    *
    * @param frameIndex The FrameIndex for which to get the corresponding ScanNumber.
    * @return The ScanNumber corresponding to the given FrameIndex.
    */
    [[nodiscard]] ScanNumber scanNumberFromFrameIndex(FrameIndex frameIndex) const;

    /**
    * @brief Get the ScanTime corresponding to a given ScanNumber.
    *
    * This function retrieves the ScanTime associated with the provided ScanNumber.
    *
    * @param scanNumber The ScanNumber for which to get the corresponding ScanTime.
    * @return The ScanTime corresponding to the given ScanNumber.
    */
    [[nodiscard]] ScanTime scanTimeFromScanNumber(ScanNumber scanNumber) const;

    [[nodiscard]] ScanTime scanTimeFromFrameIndex(FrameIndex frameIndex) const;

    /**
    * @brief Get the ScanNumber corresponding to the closest ScanTime.
    *
    * This function retrieves the ScanNumber associated with the closest ScanTime
    * to the provided ScanTime value.
    *
    * @param scanTime The target ScanTime for which to find the closest ScanNumber.
    * @return The ScanNumber corresponding to the closest ScanTime.
    */
    [[nodiscard]] ScanNumber scanNumberFromScanTime(ScanTime scanTime) const;

    /**
    * @brief Get the FrameIndex corresponding to the given ScanTime.
    *
    * This function retrieves the FrameIndex associated with the provided ScanTime
    * using a KD-Tree search.
    *
    * @param scanTime The target ScanTime for which to find the associated FrameIndex.
    * @param frameIndex Output parameter to store the retrieved FrameIndex.
    * @return An error code indicating the success or failure of the operation.
    */
    [[nodiscard]] Err frameIndexFromScanTime(ScanTime scanTime, FrameIndex *frameIndex) const;

    /**
    * @brief Get the ScanNumber corresponding to the given FrameIndex.
    *
    * This function retrieves the ScanNumber associated with the provided FrameIndex.
    *
    * @param frameIndex The target FrameIndex for which to find the associated ScanNumber.
    * @return The ScanNumber corresponding to the given FrameIndex.
    */
    [[nodiscard]] ScanNumber frameIndexFromScanNumber(ScanNumber scanNumber) const;

    /**
    * @brief Get the ScanPoints pointer associated with the given ScanNumber.
    *
    * This function retrieves the ScanPoints pointer associated with the provided ScanNumber.
    *
    * @param scanNumber The target ScanNumber for which to get the associated ScanPoints pointer.
    * @return The ScanPoints pointer corresponding to the given ScanNumber. Returns nullptr if not found.
    */
    [[nodiscard]] ScanPoints* getScanPointsByScanNumber(ScanNumber scanNumber) const;

	[[nodiscard]] float getScanPointsMedianIntensity(ScanNumber scanNumber) const;


private:

    Err buildFrameIndexVsScanNumber();

	Err buildScanNumberVsMedianIntensity();

    Err initFrameIndexVsScanTimeKDTree();

private:

    QMap<ScanNumber, ScanPoints*> m_framePntrs;
	QMap<ScanNumber, float> m_scanNumberVsMedianIntensity;
    MzTargetKey m_uniqueMsInfoScanKey;
    QMap<FrameIndex, ScanNumber> m_frameIndexVsScanNumber;
    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;

    Q_DISABLE_COPY(MsFrame) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_MSFRAME_H
