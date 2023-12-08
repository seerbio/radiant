//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_MS2ION_H
#define PYTHIADIACPP_MS2ION_H


#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "PointDF.h"


using namespace Error;

class ALGORITHMSFFLIB_EXPORTS MS2Ion {

public:

    double mz = -1.0;
    float intensity = -1.0;
    IRT iRT = -1.0;
    QString ionLabel;
    int charge = -1;

    MS2Ion() = default;
    ~MS2Ion() = default;

    friend QDebug operator<<(QDebug dbg, const MS2Ion& obj) {
        dbg.nospace() << "MS2Ion(" << obj.mz << ", " << obj.intensity << ") ";
        return dbg;
    }

    /*!
    * @brief Returns the index (1 based, not 0) and the type of ion from ionLabel
     * NOTE: An index of 0 is used to denote precursor.
    */
    Err getIonLabelInfo(QPair<IonIndex, IonType> *ionInfo) const;

    static void sortMS2IonsMzAsc(QVector<MS2Ion> *ms2Ions);

    static void sortMS2IonsIntensityDesc(QVector<MS2Ion> *ms2Ions);

    static void filterMS2IonsByMz(
        double mzMin,
        double mzMax,
        QVector<MS2Ion> *ms2Ions
    );

    static ScanPoints ms2IonsToScanPoints(const QVector<MS2Ion> &ms2Ions);
};


#endif //PYTHIADIACPP_MS2ION_H
