//
// Created by anichols on 6/4/24.
//

#ifndef SPECTRUMCENTRICMZTARGETFRAMESEARCH_H
#define SPECTRUMCENTRICMZTARGETFRAMESEARCH_H

#include "AlgorithmsFFLib_Exports.h"
#include "CandidateScores.h"
#include "Error.h"
#include "PythiaParameterReader.h"
#include "MsCalibratomatic.h"

using namespace Error;

class DeconvolvotronResult;

class ALGORITHMSFFLIB_EXPORTS SpectrumCentricMzTargetFrameSearch {

public:

    SpectrumCentricMzTargetFrameSearch() = default;
    ~SpectrumCentricMzTargetFrameSearch() = default;

    Err init(
        const PythiaParameters &pythiaParameters,
        const MsCalibratomatic &msCalibratomatic,
        const QMap<ScanNumber, ScanPoints*> &diaTargetFrame,
        const QVector<CandidateScores*> &candidateScoresPntrs,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime
        );

    Err assignIdsToScans(QVector<QPair<IdStr, DeconvolvotronResult>> *idStrVsScore);

private:
        PythiaParameters m_pythiaParameters;
        MsCalibratomatic m_msCalibratomatic;
        QMap<ScanNumber, ScanPoints*> m_diaTargetFrame;
        QVector<CandidateScores*> m_candidateScoresPntrs;
        QMap<ScanNumber, ScanTime> m_scanNumberVsScanTime;

};



#endif //SPECTRUMCENTRICMZTARGETFRAMESEARCH_H
