//
// Created by anichols on 10/17/23.
//

#ifndef PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
#define PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H

#include "AlgorithmsFFLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"

#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePair.h"


using namespace Error;

class FragLibReaderRow;


class ALGORITHMSFFLIB_EXPORTS TargetDecoyCandidatePairManager {

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

    static Err peptideStringWithModsFromPeptideSequenceChargeKey(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            PeptideStringWithMods *peptideStringWithMods,
            Charge *charge
            );

private:

    Err buildTargetDecoyCandidatePairs(const QVector<FragLibReaderRow> &fragLibReaderRows);


private:

    QVector<TargetDecoyCandidatePair> m_targetDecoyCandidatePairs;
    PythiaParameters m_pythiaParameters;

    Q_DISABLE_COPY(TargetDecoyCandidatePairManager) class Private;
    const QScopedPointer<Private> d_ptr;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
