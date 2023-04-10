//
// Created by anichols on 4/10/23.
//

#include "FragLibraryTronDIA.h"

#include "ErrorUtils.h"


Err FragLibraryTronDIA::readFragLibFile(const QString &fragLibFilePath) {

    ERR_INIT




    ERR_RETURN
}

Err FragLibraryTronDIA::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

    m_params = pythiaParameters;

    ERR_RETURN
}
