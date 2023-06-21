//
// Created by anichols on 6/21/23.
//

#ifndef PYTHIADIACPP_DECOYGENERATOMATIC_H
#define PYTHIADIACPP_DECOYGENERATOMATIC_H

#include "AlgorithmsLib_Exports.h"
#include "Error.h"
#include "ErrorUtils.h"
#include "GlobalSettings.h"


#include <QStringList>

//TODO add function to excise N-term Methionine.
using namespace Error;


class ALGORITHMSLIB_EXPORTS DecoyGeneratomatic {

public:

    static QString midShuffleLogicDecoysLogic(const QString &peptideSeq);

    static QString diannDecoyLogic(const QString &peptideSeq);


};


#endif //PYTHIADIACPP_DECOYGENERATOMATIC_H
