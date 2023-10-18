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


class ALGORITHMSLIB_EXPORTS TargetDecoyCandidatePairManager {

public:

    TargetDecoyCandidatePairManager() = default;
    ~TargetDecoyCandidatePairManager() = default;

    Err init(
            const PythiaParameters &pythiaParameters,
            const QString &fragLibFileUri
            );

private:

    QVector<TargetDecoyCandidatePair> m_targetDecoyCandidatePairs;
    QMap<TargetDecoyCandidatePairIndex, TargetDecoyCandidatePair*> m_indexVsTargetDecoyCandidatePairPtr;

};


#endif //PYTHIADIACPP_TARGETDECOYCANDIDATEPAIRMANAGER_H
