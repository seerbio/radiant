//
// Created by anichols on 4/16/23.
//

#ifndef PYTHIADIACPP_MSCALIBRATOMATIC_H
#define PYTHIADIACPP_MSCALIBRATOMATIC_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "MsFrameScoreVectorReader.h"
#include "MsUtils.h"
#include "NearestNeighborsSearch.h"
#include "ProteinDigestomatic.h"
#include "PythiaParameterReader.h"

using namespace Error;

class PeptideSequence;

class ALGORITHMSLIB_EXPORTS MsCalibratomatic {

public:

    MsCalibratomatic();
    ~MsCalibratomatic() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &firstPassSearchFilePath,
            int calPointK
            );

    // either FrameIndex, or ScanNumber can be key as they are both ints.
    Err recalibratePoints(
            const QMap<FrameIndex, ScanPoints> &indexVsScanPoints,
            QMap<FrameIndex, ScanPoints> *recalIndexVsScanPoints
            );

    [[nodiscard]] double newStDev();

private:

    Err buildCalibrator();


private:

    //Never cleared
    PythiaParameters m_params;
    NearestNeighborsSearch m_nnSearch;
    int m_calPointK;
    double m_stDevNew;
    QString m_firstPassSearchFilePath;

};


#endif //PYTHIADIACPP_MSCALIBRATOMATIC_H
