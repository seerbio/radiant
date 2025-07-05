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
class TargetDecoyCandidatePair;

class ALGORITHMSFFLIB_EXPORTS XICPeakManager {

public:

    XICPeakManager();
    ~XICPeakManager();

    Err init(
        const MsFrame &msFrame,
        const QVector<TargetDecoyCandidatePair*> &targetDecoyPointers,
        int topNMs2Ions,
        float ppmTolerance
        );

    Err init(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyPointers,
        int topNMs2Ions,
        float ppmTolerance,
        TurboXIC *turboXic
        );

    [[nodiscard]] bool isValid() const;

    Err getXIC(
        float mzVal,
        XICPointsPntrs *xicPoints
        ) const;


private:
    int m_filterLength;
    int m_polynomialOrder;
    int m_smoothCount;
    float m_ppmTolerance;
    float m_stopThresholdValue;

    bool m_isInit;

    TurboXIC *m_turboXic;
    bool m_deleteTurboXicLocal;

};



#endif //XICPEAKMANAGER_H
