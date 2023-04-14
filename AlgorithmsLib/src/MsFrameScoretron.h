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

    Err init(const PythiaParameters &params);

    Err scoreCandidatesPerFrameParallel(
            const QVector<MsFrame> &msFrames,
            const QString &msDataFilePath,
            FragLibraryTronDIA *fragLibraryTronDia
    );

    Err buildTargetCandidatesForFrame(
            const QVector<MsFrame> &msFrames,
            FragLibraryTronDIA *fragLibraryTronDia,
            QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, QVector<MS2Ion>>> *framePredictions
    ) const;

private:

    PythiaParameters m_params;


};


#endif //PYTHIADIACPP_MSFRAMESCORETRON_H
