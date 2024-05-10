//
// Created by anichols on 5/9/24.
//

#ifndef XICPEAKMANAGER_H
#define XICPEAKMANAGER_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "TurboXIC.h"

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
        float stopThresholdFraction
        );

    Err findPeaks(
        const MsFrame &msFrame,
        const QVector<float> &mzValsToExtract,
        float ppmTolerance
        );

    Err cacheXICPeakManager(const QString &outputFilePath);
    Err loadXICPeakManagerCache(const QString &outputFilePath);

private:

    Err extractXICs(
        const QVector<float> &mzValsToExtract,
        float ppmTolerance,
        TurboXIC *turboXic
        );

private:
    int m_filterLength;
    int m_polynomialOrder;
    int m_smoothCount;
    float m_stopThresholdValue;

    QHash<MzHashed , XICPoints> m_mzHashedVsXicPoints;

};



#endif //XICPEAKMANAGER_H
