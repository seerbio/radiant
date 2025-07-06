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

using FrameNumberTIMS = int;
using Ms1FrameTIMS = QMap<IonMobilityIndex, ScanPoints>;
using Ms1FrameTIMSPntrs = QMap<IonMobilityIndex, ScanPoints*>;

enum class ScanPointsSort {
    AscMz,
    AscIntensity,
    DescMz,
    DescIntensity
};


class FILEREADERSLIB_EXPORTS MsScanInfo {

public:

    int msLevel = -1;
    ScanNumber scanNumber = -1;
    float scanTime = -1.0;
    float collisionEnergy = -1.0;
    float precursorTargetMz = -1.0;
    float isoWindowLower = -1.0;
    float isoWindowUpper = -1.0;
    float ionMobilityDriftTime = -1.0;
	float mzMin = -1.0;
	float mzMax = -1.0;
	float intensityMin = -1.0;
	float intensityMax = -1.0;
	int pointCount = -1;
    IonMobilityIndex ionMobilityIndex = -1;

    long scanOffsetStart = -1;
    long scanOffsetEnd = -1;

    [[nodiscard]] MzTargetKey targetKey() const {
        return targetKey(
                precursorTargetMz - isoWindowLower,
                precursorTargetMz + isoWindowUpper
                );
    }

    static MzTargetKey targetKey(float mzStart, float mzEnd) {
        return QString::number(std::round(1000 * ((mzStart + mzEnd) / 2)));
    }
};

struct MsScan {
	MsScanInfo* msScanInfoPntr = nullptr;
	MzVals mzVals;
	IntensityVals intensityVals;
};

class FILEREADERSLIB_EXPORTS MsReaderBase {

public:

    friend class MsReaderBaseTests;
    friend class MsReaderMZMLTests;
    friend class MsReaderParquetTests;

    MsReaderBase();

    virtual ~MsReaderBase() = default;

    /**
    * @brief Sets the MS scan information.
    *
    * This method is a setter for the m_msScanInfo member variable.
    * It allows replacing the current MS scan information with provided one.
    *
    * @param msScanInfos A const reference to a QMap containing ScanNumber keys and MsScanInfo values.
    */
    void setMsScanInfo(const QMap<ScanNumber, MsScanInfo> &msScanInfos);

    /**
    * @brief Sets the ScanPoints for the MS Reader.
    *
    * This method sets the ScanPoints for the MS Reader. Before setting the ScanPoints,
    * it verifies if the number of provided ScanPoints matches the number of existing ScanPoints.
    *
    * @param scanPoints A constant reference to a QMap containing ScanNumber keys and ScanPoints values.
    *
    * @return Returns an Err object which indicates the success or failure of the operation. If the operation is
    * successful, and the number of provided ScanPoints matches the number of existing ones, an Err object initialized
    * with a success state is returned. If the number of ScanPoints does not match the original number,
    * the function will return an Err object initialized with a failure state.
    */
    // Err setScanPoints(const QMap<ScanNumber, ScanPoints> &scanPoints);

    /**
    * @brief Resets the MS Reader.
    *
    * This method clears all the ScanPoints in the MS Reader by swapping the
    * m_scanPoints member variable with an empty QMap.
    */
    // void reset();

    /**
    * @brief Opens a file in the MS Reader (declaration).
    *
    * This is a pure virtual function that must be overridden in derived MS Reader classes.
    * The function is responsible for opening a file in the MS Reader.
    *
    * @param filePath A QString representing the path to the file to be opened.
    *
    * @return Expected to return an Err object indicating the success or failure of the operation in the actual implementation.
    */
    virtual Err openFile(const QString &filePath);

    /**
    * @brief Opens a file in the MS Reader and applies a filter (declaration).
    *
    * This is a pure virtual function that must be overridden in derived MS Reader classes. The function is
    * responsible for opening a file in the MS Reader and applying a filter based on a specified column
    * and a range of values.
    *
    * @param filePath A QString representing the path to the file to be opened.
    * @param columnToFilterBy A QString deciding the column by which to filter the data.
    * @param filterRange A QPair of doubles representing the lower and upper limit of the range to filter by.
    *
    * @return Expected to return an Err object indicating the success or failure of the operation in the actual implementation.
    */
    virtual Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy,
            const QPair<double, double> &filterRange
            );

    /**
    * @brief Opens a file in the MS Reader and only retreives the values of the column specified.
    *
    * This is a pure virtual function that must be overridden in derived MS Reader classes. The function is responsible
    * for opening a file in the MS Reader and applying a filter based on a specified column.
    *
    * @param filePath A QString representing the path to the file to be opened.
    * @param columnToFilterBy A QString deciding the column by which to filter the data.
    *
    * @return Expected to return an Err object indicating the success or failure of the operation in the actual implementation.
    */
    virtual Err openFile(
            const QString &filePath,
            const QString &columnToFilterBy
    );

    /**
    * @brief Closes the currently opened file in the MS Reader.
    *
    * This method is responsible for closing the currently opened file in the MS Reader. It does this by
    * clearing all the stored data related to the file from the member variables.
    *
    * @return Returns an Err object with a success state indicating that the operation has completed successfully.
    */
    virtual Err closeFile();

    virtual Err extractScanPoints(
            const QVector<MsScanInfo*> &msScanInfos,
            QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
            );

	virtual Err extractScanPoints(
		const QVector<MsScanInfo*> &msScanInfos,
		QVector<MsScan> *msScans
		);

    virtual Err getMzTargetScanPoints(
            const MzTargetKey& targetKey,
            QMap<ScanNumber, ScanPoints>* scanNumberVsScanPoints
            );

	virtual Err getMzTargetScanPoints(
		const MzTargetKey& targetKey,
		QVector<MsScan> *msScans
		);

    /**
    * @brief Gets the file path of the currently opened file in the MS Reader.
    *
    * This method is a getter for the m_filePath member variable and allows to retrieve the
    * file path of the currently opened file in the MS Reader.
    *
    * @return Returns a QString representing the file path of the currently opened file.
    */
    QString filePath();

    /**
    * @brief Checks whether the MS Reader operates in Data Independent Analysis (DIA) mode.
    *
    * This method checks whether the MS Reader operates in Data Independent Analysis (DIA) mode by inspecting
    * whether scanning targets cycle repetitively. It does this by iterating over the ScanInfos of the tandem (MS/MS) scans,
    * counting the occurrences of each MzTargetKey, and then checks whether these keys cycle.
    *
    * @return Returns true if the MS Reader operates in DIA mode (the unique keys cycle), and false otherwise.
    */
    bool isDIA();

    /**
    * @brief Checks if the MS Reader has been initialized.
    *
    * This method determines if the MS Reader has been properly initialized based on
    * whether the m_msScanInfo member variable is empty or not. If m_msScanInfo is not
    * empty, it indicates that the file has been opened and the MS Reader is initialized.
    *
    * @return Returns true if m_msScanInfo is not empty, indicating that the MS Reader has been initialized.
    * If m_msScanInfo is empty, returns false.
    */
    bool isInit();

    /**
    * @brief Retrieves the minimum and maximum scan times from the currently opened file in the MS Reader.
    *
    * This method retrieves the minimum and maximum scan times from the m_msScanInfo variable, which represents
    * the currently opened file in the MS Reader. The minimum and maximum are determined based on the order
    * of the scan times in m_msScanInfo.
    *
    * @return Returns a QPair where the first element is the minimum scan time and the second element is the maximum scan time.
    */
    QPair<ScanTime , ScanTime > scanTimeMinMax();

    /**
    * @brief Retrieves the scan points from the MS Reader.
    *
    * This method is a getter for the m_scanPoints member variable. It allows to retrieve the scan points
    * that were read from the file in the MS Reader.
    *
    * @return Returns a QMap where each key-value pair represents the scan number and the corresponding scan points respectively.
    */
    // QMap<ScanNumber, ScanPoints> getScanPoints();

    /**
    * @brief Retrieves pointers to the scan points in the MS Reader.
    *
    * This method creates a QMap with pointers to the ScanPoints objects stored in the member variable m_scanPoints.
    * It can be used when direct access to these objects is required without copying them.
    *
    * @return Returns a QMap where each key-value pair represents the scan number and the corresponding scan points pointer respectively.
    */
    // QMap<ScanNumber, ScanPoints*> getScanPointsPntrs();

    /**
    * @brief Retrieves scan points of a specific MS level from the MS Reader.
    *
    * This method iterates over the stored scan points in the MS Reader and retrieves only those that
    * correspond to the specified MS level. It does this by checking the MS level of each scan point and
    * including only those of the desired level.
    *
    * @param msLevel The MS level for which to fetch the scan points.
    * @param scanPoints A pointer to a QMap where the fetched scan points should be stored.
    *
    * @return Returns an Err object. If successful, the function returns an Err object initialized with a
    * success state. If any error occurs during the fetching of the scan points, it returns an Err object initialized with a failure state.
    */
    // Err getScanPoints(
    //         int msLevel,
    //         QMap<ScanNumber, ScanPoints> *scanPoints
    //         );

    /**
    * @brief Retrieves pointers to scan points of a specific MS level from the MS Reader.
    *
    * This method iterates over the stored scan points in the MS Reader and retrieves pointers to those scan points that
    * correspond to the specified MS level. It does this by checking the MS level of each scan point and
    * including only those of the desired level.
    *
    * @param msLevel The MS level for which to fetch the scan point pointers.
    * @param scanPoints A pointer to a QMap where the fetched scan points pointers should be stored. The stored QMap will have scan numbers as keys and pointers to corresponding scan points as values.
    *
    * @return Returns an Err object. If successful, the function returns an Err object initialized with a
    * success state. If any error occurs during the retrieval of the scan points, it returns an Err object initialized with a failure state.
    */
    // Err getScanPoints(
    //         int msLevel,
    //         QMap<ScanNumber, ScanPoints*> *scanPoints
    // );

    /**
    * @brief Retrieves scan points of a specific ScanNumber.
    *
    * This method first checks to see that the scanPoints and msScanInfo are present for the given scanNumber.
    * If they are both present, the method will retrieve the scanPoints from the m_scanPoints.
    *
    * @param scanNumber The ScanNumber for which to fetch the scan points.
    *
    * @return Returns a QPair<Err, ScanPoints> object. If successful, the function returns an Err object initialized with a
    * success state as well as the ScanPoints. If any error occurs during the fetching of the scan points, it returns an Err object initialized with a failure state.
    */
    // QPair<Err, ScanPoints*> getScanPoints(int scanNumber);

    /**
    * @brief Collates MS2 MzTarget frames for Data-Independent Acquisition (DIA).
    *
    * This function takes a QMap of MzTargetKey to QMap of ScanNumber to ScanPoints,
    * and populates it with MS2 MzTarget frames by iterating over tandemScanInfos.
    * It checks for the existence and non-empty nature of necessary data structures,
    * retrieves MS scan information for MS level 2, and associates the corresponding
    * ScanPoints with the MzTargetKey and ScanNumber. Finally, it logs the count of
    * DIA Target Frames and returns any encountered errors.
    *
    * @param diaTargetFrame A pointer to a QMap<MzTargetKey, QMap<ScanNumber, ScanPoints *>>,
    *                      representing the target frame structure to be populated.
    * @return Err The error code indicating success or failure of the operation.
    */
    // Err collateMS2MzTargetFrames(
    //         QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrame
    // );

    /**
    * @brief Retrieves unique tandem MS scan information for MS level 2.
    *
    * This function obtains tandem MS scan information for MS level 2 using the
    * `getMsScanInfos` method. It then iterates over the obtained scan infos,
    * keeping track of unique entries based on MzTargetKey. The unique scan infos
    * are stored in a QMap and the corresponding values are converted to a QVector,
    * which is then returned.
    *
    * @return QVector<MsScanInfo> A vector containing unique tandem MS scan information
    *                             for MS level 2.
    */
    QVector<MsScanInfo> getUniqueTandemMsScanInfos();

    /**
    * @brief Retrieves the frame count for tandem MS scans at MS level 2.
    *
    * This function calculates and returns the frame count by obtaining
    * unique tandem MS scan information for MS level 2 using the
    * `getUniqueTandemMsScanInfos` method and retrieving the size of the resulting vector.
    *
    * @return int The count of frames for tandem MS scans at MS level 2.
    */
    int getFrameCount();

    /**
    * @brief Retrieves MS scan information for a specific MS level.
    *
    * This function takes an MS level as a parameter and returns a QMap
    * containing MS scan information (keyed by ScanNumber) corresponding to
    * the specified MS level. It iterates over the existing MsScanInfo entries
    * in the class, filtering out entries that do not match the given MS level,
    * and populates the resulting QMap with the valid entries.
    *
    * @param msLevel The MS level for which scan information is to be retrieved.
    * @return QMap<ScanNumber, MsScanInfo> A map of MS scan information for the specified MS level.
    *
    */
    QMap<ScanNumber, MsScanInfo> getMsScanInfos(int msLevel);
    Err getMsScanInfos(
            const MzTargetKey &mzTargetKey,
            QVector<MsScanInfo*> *msScanInfos
            ) const;

    /**
    * @brief Retrieves all available MS scan information for a particular msLevel.
    *
    * This function returns a QMap containing all available MS scan information
    * keyed by ScanNumber. It provides access to the complete set of stored
    * MsScanInfo entries within the class for a particular msLevel.
    *
    * @return QMap<ScanNumber, MsScanInfo> A map of all available MS scan information.
    *
    */
    QMap<ScanNumber, MsScanInfo> getMsScanInfos();

    /**
    * @brief Retrieves MS scan information for a specific ScanNumber.
    *
    * This function takes a ScanNumber as a parameter and retrieves the
    * corresponding MsScanInfo from the internal storage. It checks whether
    * the provided ScanNumber is present in the stored MS scan information,
    * and if so, it sets the output parameter `msScanInfo` with the retrieved
    * MsScanInfo. Any errors encountered during the process are indicated by
    * the returned Err code.
    *
    * @param scanNumber The ScanNumber for which MS scan information is to be retrieved.
    * @param msScanInfo A pointer to a MsScanInfo object that will be populated with
    *                   the retrieved MS scan information.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    Err getMsScanInfo(
            ScanNumber scanNumber,
            MsScanInfo *msScanInfo
            ) const;

    /**
    * @brief Sorts ScanPoints based on the specified sorting criteria.
    *
    * This function takes a `ScanPointsSort` enum representing the desired sorting order,
    * and a pointer to a `ScanPoints` object that needs to be sorted. It uses standard
    * C++ algorithms to perform the sorting based on the chosen criteria. The supported
    * sorting options include ascending and descending order for both m/z and intensity.
    *
    * @param sort The sorting order specified by the `ScanPointsSort` enum.
    * @param scanPoints A pointer to the `ScanPoints` object that will be sorted in-place.
    *
    */
    static void sortScanPoints(
            const ScanPointsSort &sort,
            ScanPoints *scanPoints
                    );

    /**
    * @brief Retrieves the nearest ScanNumber based on a given ScanTime.
    *
    * This function calculates the nearest ScanNumber to the provided ScanTime
    * by utilizing a mapping between ScanNumber and ScanTime. If the internal
    * mapping is not initialized, it populates it using the available MsScanInfo
    * entries. It then retrieves the ScanTime values and ScanNumber keys from
    * the mapping, finds the closest ScanTime index, and returns the corresponding
    * ScanNumber. The MathUtils::closest function is used for index calculation.
    *
    * @param scanTime The ScanTime for which the nearest ScanNumber is to be retrieved.
    * @return ScanNumber The nearest ScanNumber to the provided ScanTime.
    *
    */
    ScanNumber getNearestScanNumberFromScanTime(ScanTime scanTime);

    /**
    * @brief Retrieves a mapping of ScanNumber to ScanTime for available MS scan information.
    *
    * This function iterates over the stored MsScanInfo entries and creates a QMap
    * containing a mapping between ScanNumber and ScanTime. Each entry in the map
    * corresponds to a ScanNumber with its associated ScanTime value. The resulting
    * map is then returned.
    *
    * @return QMap<ScanNumber, ScanTime> A mapping of ScanNumber to ScanTime for available
    *                                  MS scan information
    */
    [[nodiscard]] QMap<ScanNumber, ScanTime> getScanNumberVsScanTime() const;

    Err getHiLoMzPrecursors(QPair<MzMin, MzMax> *precursorMzLoVsMzHi);

    // Err driftTimeFromIonMobilityIndex(
    //         const IonMobilityIndex &ionMobilityIndex,
    //         double *driftTime
    //         ) const;

    /**
    * @brief Splits ScanPoints into separate vectors for m/z and intensity values.
    *
    * This function takes a ScanPoints object and two QVector<float> pointers
    * for m/z values (`mzVals`) and intensity values (`intensityVals`). It checks
    * whether the input ScanPoints is not empty, clears the output vectors, and then
    * iterates over the ScanPoints, populating the m/z and intensity vectors.
    * Any encountered errors during the process are indicated by the returned Err code.
    *
    * @param scanPoints The input ScanPoints to be split.
    * @param mzVals A pointer to a QVector<float> that will be populated with m/z values.
    * @param intensityVals A pointer to a QVector<float> that will be populated with intensity values.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    static Err splitScanPoints(
            const ScanPoints &scanPoints,
            QVector<float> *mzVals,
            QVector<float> *intensityVals
            );

    /**
    * @brief Splits ScanPoints into separate vectors for m/z and intensity values.
    *
    * This function takes a pointer to a ScanPoints object (`scanPoints`) and two
    * QVector<float> pointers for m/z values (`mzVals`) and intensity values (`intensityVals`).
    * It checks whether the input ScanPoints is not empty, clears the output vectors, and then
    * iterates over the ScanPoints, populating the m/z and intensity vectors.
    * Any encountered errors during the process are indicated by the returned Err code.
    *
    * @param scanPoints A pointer to the input ScanPoints to be split.
    * @param mzVals A pointer to a QVector<float> that will be populated with m/z values.
    * @param intensityVals A pointer to a QVector<float> that will be populated with intensity values.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    static Err splitScanPoints(
            ScanPoints *scanPoints,
            QVector<float> *mzVals,
            QVector<float> *intensityVals
    );

    /**
    * @brief Zips two vectors of m/z and intensity values into ScanPoints.
    *
    * This function takes two QVector<float> objects (`mzVals` and `intensityVals`)
    * representing m/z and intensity values, and a pointer to a ScanPoints object (`scanPoints`).
    * It checks whether the sizes of `mzVals` and `intensityVals` are equal, clears the output
    * ScanPoints, and then iterates over the vectors to create ScanPoints with corresponding
    * m/z and intensity values. Any encountered errors during the process are indicated
    * by the returned Err code.
    *
    * @param mzVals The vector of m/z values to be zipped.
    * @param intensityVals The vector of intensity values to be zipped.
    * @param scanPoints A pointer to the output ScanPoints to store the zipped data.
    * @return Err The error code indicating success or failure of the operation.
    *
    */
    static Err zipScanPoints(
            const QVector<float> &mzVals,
            const QVector<float> &intensityVals,
            ScanPoints *scanPoints
            );

    void printSize() const;

    Err printFileInfo();
    [[nodiscard]] bool isTIMS() const;

    [[nodiscard]] float mzMs2Min() const;
    [[nodiscard]] float mzMs2Max() const;

    // QMap<FrameNumberTIMS, Ms1FrameTIMS>* frameNumberVsMS1FrameTIMSPntr();


protected:

    bool m_fileIsCalibrated;
    bool m_isTIMS;

    QMap<ScanNumber, MsScanInfo> m_msScanInfo;
    QMap<MzTargetKey, QVector<MsScanInfo*>> m_mzTargetVsScanInfosPntrs;
    // QMap<ScanNumber, ScanPoints>  m_scanPoints;
    QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;

    // QMap<FrameNumberTIMS, Ms1FrameTIMS> m_frameNumberVsMS1FrameTIMS;
    // QMap<FrameIndex, double> m_frameIndexVsDriftTime;


    QString m_filePath;

    float m_mzMs1Min;
    float m_mzMs1Max;
    float m_mzMs2Min;
    float m_mzMs2Max;
	float m_intensityMs1Min;
	float m_intensityMs1Max;
	float m_intensityMs2Min;
	float m_intensityMs2Max;

};


#endif //MSREADERBASE_H
