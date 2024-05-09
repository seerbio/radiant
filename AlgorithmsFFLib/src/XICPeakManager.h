//
// Created by anichols on 5/9/24.
//

#ifndef XICPEAKMANAGER_H
#define XICPEAKMANAGER_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"


using namespace Error;

class MsFrame;
class PeakIntegratomatic;
class TurboXIC;

class ALGORITHMSFFLIB_EXPORTS XICPeakManager {

public:

    XICPeakManager();
    ~XICPeakManager() = default;

    Err init(
        int filterLength,
        int polynomialOrder,
        int smoothCount,
        double stopThresholdFraction
        );

    Err findPeaks(
        const MsFrame &msFrame,
        const QVector<float> &mzValsToExtract,
        float ppmTolerance
        );

private:

    Err extractPeaks(
        const QVector<float> &mzValsToExtract,
        const PeakIntegratomatic &peakIntegratomatic,
        float ppmTolerance,
        TurboXIC *turboXic
        );

private:
    int m_filterLength;
    int m_polynomialOrder;
    int m_smoothCount;
    double m_stopThresholdValue;


};



#endif //XICPEAKMANAGER_H
