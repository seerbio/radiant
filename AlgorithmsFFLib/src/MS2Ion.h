//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_MS2ION_H
#define PYTHIADIACPP_MS2ION_H


#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS MS2Ion {

public:

    float mz = -1.0;
    float intensity = -1.0;
    QString ionLabel;
    int charge = -1;
    float targetDecoyFrequencyRatio = -1.0;

    MS2Ion() = default;
    ~MS2Ion() = default;

    friend QDebug operator<<(QDebug dbg, const MS2Ion& obj) {
        dbg.nospace() << "MS2Ion(" << obj.mz << ", " << obj.intensity << ", " << obj.ionLabel << ") ";
        return dbg;
    }

    /**
    * @brief Gets ion label information for an MS2 ion.
    *
    * This function extracts and returns information about the ion label, including the ion index and type.
    *
    * @param ionInfo Output parameter, a QPair<IonIndex, IonType>* containing the ion index and type.
    * @return An Err enum indicating the success or failure of the operation.
    *
    * NOTE: the index (1 based, not 0) and the type of ion from ionLabel
    * NOTE: An index of 0 is used to denote precursor.
    */
    Err getIonLabelInfo(QPair<IonIndex, IonType> *ionInfo) const;

    /**
    * @brief Sorts MS2 ions in ascending order based on their mz values.
    *
    * This function sorts the provided QVector of MS2 ions in ascending order based on their mz values.
    *
    * @param ms2Ions A QVector<MS2Ion>* representing the MS2 ions to be sorted.
    */
    static void sortMS2IonsMzAsc(QVector<MS2Ion> *ms2Ions);

    /**
    * @brief Sorts MS2 ions in descending order based on their intensity values.
    *
    * This function sorts the provided QVector of MS2 ions in descending order based on their intensity values.
    *
    * @param ms2Ions A QVector<MS2Ion>* representing the MS2 ions to be sorted.
    */
    static void sortMS2IonsIntensityDesc(QVector<MS2Ion> *ms2Ions);

    /**
    * @brief Filters MS2 ions based on mz values within a specified range.
    *
    * This function filters the provided QVector of MS2 ions based on mz values within the specified range.
    *
    * @param mzMin A double representing the minimum mz value for filtering.
    * @param mzMax A double representing the maximum mz value for filtering.
    * @param ms2Ions A QVector<MS2Ion>* representing the MS2 ions to be filtered.
    */
    static void filterMS2IonsByMz(
        double mzMin,
        double mzMax,
        QVector<MS2Ion> *ms2Ions
    );

    /**
    * @brief Converts MS2 ions to ScanPoints.
    *
    * This function converts a QVector of MS2 ions to ScanPoints, where each MS2Ion is converted to a ScanPoint.
    *
    * @param ms2Ions A QVector<MS2Ion> representing the MS2 ions to be converted.
    * @return ScanPoints representing the converted ScanPoints.
    */
    static ScanPoints ms2IonsToScanPoints(const QVector<MS2Ion> &ms2Ions);
};


#endif //PYTHIADIACPP_MS2ION_H
