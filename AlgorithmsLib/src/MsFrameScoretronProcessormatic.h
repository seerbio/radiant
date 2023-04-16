//
// Created by anichols on 4/14/23.
//

#ifndef PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
#define PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H


#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "FragLibraryTronDIA.h"
#include "GlobalSettings.h"
#include "MsFrame.h"
#include "MsFrameScoreVectorReader.h"
#include "MsUtils.h"
#include "PythiaParameterReader.h"


using namespace Error;


struct ScoredPSM : public ParquetReaderInputBase {
    FrameIndex frameIndex = -1;
    PeptideStringWithMods peptideStringWithMods;
    ExtractPoints extractPoints;

};


class FrameIndexCandidate;

class ALGORITHMSLIB_EXPORTS MsFrameScoretronProcessormatic {

public:

    friend class MsFrameScoretronProcessormaticTests;

    MsFrameScoretronProcessormatic() = default;
    ~MsFrameScoretronProcessormatic() = default;

    Err init(const PythiaParameters &pythiaParameters);

    Err rescoreMsFrame(
            const QVector<MsFrameScanPointRows> &msFrameScanPointRows,
            const QVector<MsFrameScoreVectorReaderRow> &scoreVectors,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &framePredictions
            );

private:

    PythiaParameters m_params;

};


#endif //PYTHIADIACPP_MSFRAMESCORETRONPROCESSORMATIC_H
