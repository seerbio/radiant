//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRON_H
#define PYTHIADIACPP_MSFRAMESCORETRON_H


#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FragLibraryTronDIA.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "PythiaParameterReader.h"

using namespace Error;

class ALGORITHMSLIB_EXPORTS MsFrameScoretron {

public:

    MsFrameScoretron() = default;
    ~MsFrameScoretron() = default;

    static QPair<Err, QPair<UniqueMsInfoScanKey, QString>> scoreCandidatesFrameWrite(
            const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk,
            const PythiaParameters &params,
            const QString &msDataFilePath
    );

};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
