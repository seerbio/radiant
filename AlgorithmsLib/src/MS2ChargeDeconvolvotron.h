//
// Created by anichols on 7/7/23.
//

#ifndef PYTHIADIACPP_MS2CHARGEDECONVOLVOTRON_H
#define PYTHIADIACPP_MS2CHARGEDECONVOLVOTRON_H


#include "AlgorithmsLib_Exports.h"
#include "CSVReader.h"
#include "Error.h"
#include "GlobalSettings.h"


class ALGORITHMSLIB_EXPORTS MS2ChargeDeconvolvotron {

public:

    friend class MS2ChargeDeconvolvotronTests;

    MS2ChargeDeconvolvotron();
    ~MS2ChargeDeconvolvotron();

    Err init(
            const QString &modelFilePath,
            double ppmTolerance
            );

    Err deisotopeScanPoints(
            const ScanPoints &scanPoints,
            ScanPoints *scanPointsDeisotoped
            );

private:

    QVector<double> buildMzValsForChargeMono(double mzMaxIntensity);

    Err testChargeMonoCaller(
            const ScanPoints &scanPoints,
            double mzVal,
            QVector<float> *predVec
    );

private:

    Q_DISABLE_COPY(MS2ChargeDeconvolvotron) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_MS2CHARGEDECONVOLVOTRON_H
