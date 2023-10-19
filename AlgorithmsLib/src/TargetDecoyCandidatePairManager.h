//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H

#include "AlgorithmsLib_Exports.h"

#include "BoostRTreeWrapper.h"
#include "Error.h"
#include "GlobalSettings.h"

#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePair.h"


using namespace Error;

class FragLibReaderRow;


class ALGORITHMSLIB_EXPORTS TargetDecoyCandidatePairManager {

public:

    TargetDecoyCandidatePairManager();
    ~TargetDecoyCandidatePairManager() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibFileUri
            );

    Err getTargetDecoyCandidatePairPointers(
            double mzMin,
            double mzMax,
            QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
            );

    Err getTargetDecoyCandidatePairPointers(
            double mzMin,
            double mzMax,
            double randomSelectionFraction,
            QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
    );

    bool isInit();

    static Err peptideStringWithModsFromPeptideSequenceChargeKey(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            PeptideStringWithMods *peptideStringWithMods,
            Charge *charge
            );

private:

    Err buildTargetDecoyCandidatePairs(const QVector<FragLibReaderRow> &fragLibReaderRows);

    Err buildIndexVsTargetDecoyCandidatePairPtrs();

    Err initBoostRTreeWrapper();

private:

    QVector<TargetDecoyCandidatePair> m_targetDecoyCandidatePairs;
    QMap<TargetDecoyCandidatePairIndex, TargetDecoyCandidatePair*> m_indexVsTargetDecoyCandidatePairPtrs;

    PythiaParameters m_pythiaParameters;
    BoostRTreeWrapper m_boostRTreeWrapper;

    bool m_isInit;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
