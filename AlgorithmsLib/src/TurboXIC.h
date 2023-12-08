//
// Created by anichols on 12/18/22.
//

#ifndef PYTHIACPP_TURBOXIC_H
#define PYTHIACPP_TURBOXIC_H

#include "AlgorithmsLib_Exports.h"
#include "GlobalSettings.h"
#include "Error.h"


using namespace Error;

class XICPoints {

public:

    XICPoints() = default;
    ~XICPoints() = default;

    QMap<ScanNumber, ScanPoints> scanNumbersVsScanPoints;
    QMap<ScanNumber, double> scanNumbersVsIntensity;

    [[nodiscard]] Err buildScanNumbersVsIntensityVals() {

        ERR_INIT

        if (scanNumbersVsIntensity.isEmpty()) {

            QMap<ScanNumber, double> scanNumberVsIntensity;

            for (auto it = scanNumbersVsScanPoints.begin(); it != scanNumbersVsScanPoints.end(); it++) {

                const ScanNumber scanNumber = it.key();
                const ScanPoints &scanPoints = it.value();
                for (const ScanPoint &sp : scanPoints) {
                    scanNumberVsIntensity[scanNumber] += sp.y();
                }

            }

            scanNumbersVsIntensity = scanNumberVsIntensity;
        }

        ERR_RETURN
    }

    [[nodiscard]] QMap<ScanNumber, QVector<double>> scanNumberVsMzVals() const {

        QMap<ScanNumber, QVector<double>> output;

        for (auto it = scanNumbersVsScanPoints.begin(); it != scanNumbersVsScanPoints.end(); it++) {

            const ScanNumber scanNumber = it.key();
            const ScanPoints &scanPoints = it.value();
            for (const ScanPoint &sp : scanPoints) {
                output[it.key()].push_back(sp.x());
            }

        }

        return output;
    }

};


class ALGORITHMSLIB_EXPORTS TurboXIC {


public:

    TurboXIC();
    ~TurboXIC();

    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);
    Err init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints);

    XICPoints extractPointsXIC(
            double mzMin,
            double mzMax,
            ScanNumber scanNumberMin,
            ScanNumber scanNumberMax
    );

    [[nodiscard]] ScanPoints extractSpectrum(
            double mzMin,
            double mzMax,
            ScanNumber scanNumberMin,
            ScanNumber scanNumberMax
            ) const;

    Err getRTreeLimits(
            double *scanNumberMin,
            double *scanNumberMax,
            double *mzMin,
            double *mzMax
            ) const;

    bool isInit();


private:

    Q_DISABLE_COPY(TurboXIC) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_TURBOXIC_H
