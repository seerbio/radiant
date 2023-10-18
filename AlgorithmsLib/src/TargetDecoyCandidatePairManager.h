//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H

#include "AlgorithmsLib_Exports.h"


#include "Error.h"
#include "GlobalSettings.h"

#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePair.h"


using namespace Error;

class FragLibReaderRow;


class ALGORITHMSLIB_EXPORTS TargetDecoyCandidatePairManager {

public:

    TargetDecoyCandidatePairManager() = default;
    ~TargetDecoyCandidatePairManager() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibFileUri
            );

    static Err peptideStringWithModsFromPeptideSequenceChargeKey(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            PeptideStringWithMods *peptideStringWithMods,
            Charge *charge
            );

private:

    Err buildTargetDecoyCandidatePairs(const QVector<FragLibReaderRow> &fragLibReaderRows);

    Err buildMS2Ions(
            const FragLibReaderRow &flrr,
            QVector<MS2Ion> *ms2Ions
    ) const;

    Err buildIndexVsTargetDecoyCandidatePairPtrs();

private:

    QVector<TargetDecoyCandidatePair> m_targetDecoyCandidatePairs;
    QMap<TargetDecoyCandidatePairIndex, TargetDecoyCandidatePair*> m_indexVsTargetDecoyCandidatePairPtrs;

    PythiaParameters m_pythiaParameters;
};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
