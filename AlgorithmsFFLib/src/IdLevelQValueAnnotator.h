#ifndef PYTHIADIACPP_IDLEVELQVALUEANNOTATOR_H
#define PYTHIADIACPP_IDLEVELQVALUEANNOTATOR_H

#include "AlgorithmsFFLib_Exports.h"
#include "Error.h"

#include <QVector>

class CandidateScores;

using namespace Error;

class ALGORITHMSFFLIB_EXPORTS IdLevelQValueAnnotator {

public:

    static Err annotate(
        QVector<CandidateScores*> *candidateScores,
        bool useClassifierScore = true,
        double localRtBinSeconds = 0.0
        );

};

#endif //PYTHIADIACPP_IDLEVELQVALUEANNOTATOR_H
