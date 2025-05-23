//
// Created by andrewnichols on 5/22/25.
//

#include "PythiaDDAFFWorkflow.h"


PythiaDDAFFWorkflow::PythiaDDAFFWorkflow() {}
PythiaDDAFFWorkflow::~PythiaDDAFFWorkflow() {}


Err PythiaDDAFFWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &fastaUri,
        const QString &outputFolderPath
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

        ParallelUtils::printSystemDetails();

        m_pythiaParameters = pythiaParameters;

        ERR_RETURN
}

Err PythiaDDAFFWorkflow::processFile(const QString &msDataFilePath) {

        ERR_INIT

        e = ErrorUtils::fileExists(msDataFilePath); ree;
        MsReaderPointerAcc msReaderPointerAcc;

        msReaderPointerAcc.setUseLazyLoading(false);

        e = msReaderPointerAcc.openFile(msDataFilePath); ree;
        msReaderPointerAcc.ptr->printSize();


        ERR_RETURN
}