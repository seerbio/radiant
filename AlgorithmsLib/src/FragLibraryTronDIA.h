//
// Created by anichols on 4/10/23.
//

#ifndef PYTHIADIACPP_FRAGLIBRARYTRONDIA_H
#define PYTHIADIACPP_FRAGLIBRARYTRONDIA_H

#include "AlgorithmsLib_Exports.h"

#include "Error.h"
#include "GlobalSettings.h"
#include "MathUtils.h"
#include "PythiaParameterReader.h"




class ALGORITHMSLIB_EXPORTS FragLibraryTronDIA {

public:

    FragLibraryTronDIA() = default;
    ~FragLibraryTronDIA() = default;

    Err init(const PythiaParameters &pythiaParameters);

    Err readFragLibFile(const QString &fragLibFilePath);


private:


    PythiaParameters m_params;

};


#endif //PYTHIADIACPP_FRAGLIBRARYTRONDIA_H
