//
// Created by anichols on 12/18/22.
//

#ifndef PYTHIACPP_TURBOXIC_H
#define PYTHIACPP_TURBOXIC_H

#include "AlgorithmsFFLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"
#include "MsReaderBase.h"
#include "PointFF.h"


using namespace Error;


struct ALGORITHMSFFLIB_EXPORTS XICPoint {
    float mz = -1.0;
    float intensity = -1.0;
    ScanNumber scanNumber = -1;
	FrameIndex frameIndex = -1;
    IonMobilityIndex ionMobilityIndex = -1;

    friend QDebug operator<<(QDebug dbg, const XICPoint& obj) {
        dbg.nospace() << "XICPoint(" << obj.frameIndex << " " << obj.scanNumber << " " << obj.mz << ", " << obj.intensity << ") ";
        return dbg;
    }

    friend QDataStream &operator <<(QDataStream &stream, const XICPoint &obj) {
        stream << obj.scanNumber;
        stream << obj.frameIndex;
        stream << obj.mz;
        stream << obj.intensity;
    	stream << obj.ionMobilityIndex;

        return stream;
    }

    friend QDataStream &operator >>(QDataStream &stream, XICPoint &obj) {
        stream >> obj.scanNumber;
        stream >> obj.frameIndex;
        stream >> obj.mz;
        stream >> obj.intensity;
        stream >> obj.ionMobilityIndex;

        return stream;
    }
};

using XICPoints = QVector<XICPoint>;
using XICPointsPntrs = QVector<XICPoint*>;


class ALGORITHMSFFLIB_EXPORTS TurboXIC {


public:

	friend class TurboXICTests;

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
    [[nodiscard]] Err init(const QVector<MsScan> &msScans);
	Err init(const QString& filePath);

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
    [[nodiscard]] Err extractPointsXIC(
            float mzMin,
            float mzMax,
            XICPointsPntrs *xicPointsPntrs
            );

    /**
    * @brief Checks if the TurboXIC private implementation is initialized.
    *
    * This method checks whether the TurboXIC private implementation is initialized by verifying the presence
    * of the RTree and ensuring it is not empty.
    *
    * @return True if the TurboXIC private implementation is initialized, otherwise false.
    */
    [[nodiscard]] bool isInit() const;

	[[nodiscard]] Err write(const QString& filePath) const;


    static void filterXICPointsByScanNumber(
        ScanNumber scanNumberMin,
        ScanNumber scanNumberMzx,
        XICPointsPntrs *xicPoints
        );


private:

	float* m_mzVals;
	float* m_indexesMz;
	uint32_t* m_indexesI;
	XICPoints m_xicPoints;

	quint64 m_scanPointsCountAlignas;
	quint64 m_alignasPiecesOfEight;

};


#endif //PYTHIACPP_TURBOXIC_H
