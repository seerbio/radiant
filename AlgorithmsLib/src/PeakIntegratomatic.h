//
// Created by anichols on 12/13/22.
//

#ifndef PYTHIACPP_PEAKINTEGRATOMATIC_H
#define PYTHIACPP_PEAKINTEGRATOMATIC_H


#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;


struct PeakIntegratomaticParameters {

    int filterLength = 7;
    int smoothCount = 3;
    double sigma = 1.0;
    double signalToNoiseRatio = 3;

    bool isValid() const {
        return filterLength > 3
            && smoothCount > 0
            && sigma > 0
            && signalToNoiseRatio > 0;
    }

};


class ALGORITHMSLIB_EXPORTS PeakIntegratomatic {

public:

    PeakIntegratomatic();
    ~PeakIntegratomatic();

    Err init(const PeakIntegratomaticParameters &params);

    Err findAllPeaksLimitsInXIC(
            const QVector<double> &intensityVec,
            QVector<PeakIntegrationIndexes> *peakLimits,
            QVector<double> *intensityVecSmoothed
    );

    /* Integrates the highest peak only.*/
    static Err simpleIntegrator(
            const QVector<double> &vec,
            double stopThresholdFraction,
            int filterLength,
            int smoothCount,
            PeakIntegrationIndexes *peakIntegrationIndexes
            );

private:

    Q_DISABLE_COPY(PeakIntegratomatic) class Private;
    const QScopedPointer<Private> d_ptr;


};


#endif //PYTHIACPP_PEAKINTEGRATOMATIC_H
