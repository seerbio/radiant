//
// Created by anichols on 6/4/24.
//

#ifndef SPECTRUMCENTRICMZTARGETFRAMESEARCH_H
#define SPECTRUMCENTRICMZTARGETFRAMESEARCH_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"
#include "PythiaParameterReader.h"
#include "MsCalibratomatic.h"
#include "TargetDecoyCandidatePair.h"

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS SpectrumCentricMzTargetFrameSearch {

public:

    SpectrumCentricMzTargetFrameSearch() = default;
    ~SpectrumCentricMzTargetFrameSearch() = default;

    Err init(
        const PythiaParameters &pythiaParameters,
        const MsCalibratomatic &msCalibratomatic,
        const QMap<ScanNumber, ScanPoints*> &diaTargetFrame,
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
        );

    Err assignIdsToScans();

private:
        PythiaParameters m_pythiaParameters;
        MsCalibratomatic m_msCalibratomatic;
        QMap<ScanNumber, ScanPoints*> m_diaTargetFrame;
        QVector<TargetDecoyCandidatePair*> m_targetDecoyCandidatePairs;
        QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;

};



#endif //SPECTRUMCENTRICMZTARGETFRAMESEARCH_H
