//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAFFWorkflow.h"

#include "MsReaderParquet.h"




PythiaDIAFFWorkflow::PythiaDIAFFWorkflow(){}

Err PythiaDIAFFWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &fastaUri
        ) {

    ERR_INIT

    pythiaParameters.print();

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::fileExists(fragLibUri); ree;

    if (!fastaUri.isEmpty()) {
        e = ErrorUtils::fileExists(fastaUri); ree;
        m_fastaUri = fastaUri;
    }

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::processFile(const QString &_msDataFilePath) {

    ERR_INIT



    ERR_RETURN
}


