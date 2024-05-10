//
// Created by anichols on 12/18/22.
//

#ifndef PYTHIACPP_TURBOXIC_H
#define PYTHIACPP_TURBOXIC_H

#include "AlgorithmsFFLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "PointFF.h"


using namespace Error;

struct ALGORITHMSFFLIB_EXPORTS XICPoint {
    float mz = -1.0;
    float intensity = -1.0;
    ScanNumber scanNumber = -1;

    friend QDebug operator<<(QDebug dbg, const XICPoint& obj) {
        dbg.nospace() << "XICPoint(" << obj.scanNumber << " " << obj.mz << ", " << obj.intensity << ") ";
        return dbg;
    }

    friend QDataStream &operator <<(QDataStream &stream, const XICPoint &obj) {
        stream << obj.scanNumber;
        stream << obj.mz;
        stream << obj.intensity;

        return stream;
    }

    friend QDataStream &operator >>(QDataStream &stream, XICPoint &obj) {
        stream >> obj.scanNumber;
        stream >> obj.mz;
        stream >> obj.intensity;

        return stream;
    }
};

using XICPoints = std::vector<XICPoint>;


class ALGORITHMSFFLIB_EXPORTS TurboXIC {


public:

    TurboXIC();
    ~TurboXIC();

    /**
    * @brief Initializes the TurboXIC private implementation.
    *
    * This method initializes the TurboXIC private implementation by populating the scan number vs. intensity pointers,
    * creating an R-tree for spatial indexing, and performing necessary checks.
    *
    * @param scanNumberVsScanPoints A map of scan numbers to scan points.
    * @return An error code indicating success or failure.
    */
    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);

    /**
    * @brief Initializes the TurboXIC private implementation.
    *
    * This method initializes the TurboXIC private implementation by populating the scan number vs. intensity pointers,
    * creating an R-tree for spatial indexing, and performing necessary checks.
    *
    * @param scanNumberVsScanPoints A pointer to a map of scan numbers to scan points.
    * @return An error code indicating success or failure.
    */
    Err init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints);

    /**
    * @brief Extracts XIC points within the specified m/z range.
    *
    * This method extracts XIC (Extracted Ion Chromatogram) points from the TurboXIC private implementation
    * within the specified m/z range.
    *
    * @param mzMin The minimum m/z value for the extraction.
    * @param mzMax The maximum m/z value for the extraction.
    * @return XICPoints containing the extracted XIC points.
    */
    XICPoints extractPointsXIC(
            float mzMin,
            float mzMax
    );

    /**
    * @brief Retrieves the limits of the RTree in the TurboXIC private implementation.
    *
    * This method retrieves the minimum and maximum m/z values of the RTree used in the TurboXIC private implementation.
    *
    * @param mzMin Pointer to store the minimum m/z value.
    * @param mzMax Pointer to store the maximum m/z value.
    * @return Error code indicating the success or failure of the operation.
    */
    Err getRTreeLimits(
            float *mzMin,
            float *mzMax
            ) const;

    /**
    * @brief Checks if the TurboXIC private implementation is initialized.
    *
    * This method checks whether the TurboXIC private implementation is initialized by verifying the presence
    * of the RTree and ensuring it is not empty.
    *
    * @return True if the TurboXIC private implementation is initialized, otherwise false.
    */
    bool isInit();


private:

    Q_DISABLE_COPY(TurboXIC) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_TURBOXIC_H
