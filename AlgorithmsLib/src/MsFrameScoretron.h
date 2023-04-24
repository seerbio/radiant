//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRON_H
#define PYTHIADIACPP_MSFRAMESCORETRON_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FragLibReader.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "PythiaParameterReader.h"

using namespace Error;

class ALGORITHMSLIB_EXPORTS MsFrameScoretron {

public:

    MsFrameScoretron() = default;
    ~MsFrameScoretron() = default;

    static QPair<Err, QPair<UniqueMsInfoScanKey, QString>> scoreCandidates(
            const PythiaParameters &params,
            const QString &msDataFilePath,
            const QString &fragLibFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const QPair<double, double> &mzTargetStartStop
    );

};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
