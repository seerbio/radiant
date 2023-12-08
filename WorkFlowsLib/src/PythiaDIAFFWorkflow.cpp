//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAFFWorkflow.h"

#include "MsReaderParquet.h"
#include "TargetDecoyCandidatePairManager.h"


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

Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT





    ERR_RETURN
}

Err PythiaDIAFFWorkflow::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
        const QList<MsScanInfo> &msScanInfos,
        double selectionFraction,
        QMap<UniqueMsInfoScanKey, TargetDecoyCandidatePair> *uniqueInfoScanKeyVsTargetDecoyCandidatePointers
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msScanInfos); ree;
    uniqueInfoScanKeyVsTargetDecoyCandidatePointers->clear();

    TargetDecoyCandidatePairManager targetDecoyCandidatePairManager;
    e = targetDecoyCandidatePairManager.init(
            m_pythiaParameters,
            m_fragLibUri
            ); ree;

    for (const MsScanInfo &msScanInfo : msScanInfos) {
        qDebug() << msScanInfo.precursorTargetMz << msScanInfo.isoWindowLower << msScanInfo.isoWindowUpper;
    }

    ERR_RETURN
}


