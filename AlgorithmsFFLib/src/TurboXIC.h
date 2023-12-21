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

class ALGORITHMSFFLIB_EXPORTS XICPoints {

public:

    XICPoints() = default;
    ~XICPoints() = default;

    QMap<ScanNumber, ScanPoints> scanNumbersVsScanPoints;
    QMap<ScanNumber, float> scanNumbersVsIntensity;

    void buildScanNumbersVsIntensityVals() {

        if (scanNumbersVsIntensity.isEmpty()) {

            QMap<ScanNumber, float> scanNumberVsIntensity;

            for (auto it = scanNumbersVsScanPoints.begin(); it != scanNumbersVsScanPoints.end(); it++) {

                const ScanNumber scanNumber = it.key();
                const ScanPoints &scanPoints = it.value();
                for (const ScanPoint &sp : scanPoints) {
                    scanNumberVsIntensity[scanNumber] += sp.y();
                }

            }

            scanNumbersVsIntensity = scanNumberVsIntensity;
        }
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


class ALGORITHMSFFLIB_EXPORTS TurboXIC {


public:

    TurboXIC();
    ~TurboXIC();

    Err init(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);
    Err init(QMap<ScanNumber, ScanPoints*> *scanNumberVsScanPoints);

    XICPoints extractPointsXIC(
            float mzMin,
            float mzMax
    );

    Err getRTreeLimits(
            float *mzMin,
            float *mzMax
            ) const;

    bool isInit();


private:

    Q_DISABLE_COPY(TurboXIC) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIACPP_TURBOXIC_H
