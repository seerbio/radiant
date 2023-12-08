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

    e = m_targetDecoyCandidatePairManager.init(
            m_pythiaParameters,
            m_fragLibUri
            ); ree;

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT






    ERR_RETURN
}

Err PythiaDIAFFWorkflow::buildCalibration() {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    const bool useNeuralNetworkScores = false;
    const int minTrainingCountTranche = 50;

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;

    const double calibrationTrainingFraction = 0.2;

    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
            m_pythiaParameters.mzMinDataStructure,
            m_pythiaParameters.mzMaxDataStructure,
            calibrationTrainingFraction,
            &scoredTargetDecoyPointers
            ); ree;


    ERR_RETURN
}

Err PythiaDIAFFWorkflow::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
        const QVector<MsScanInfo> &msScanInfos,
        double selectionFraction,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *uniqueInfoScanKeyVsTargetDecoyCandidatePointers
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msScanInfos); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    uniqueInfoScanKeyVsTargetDecoyCandidatePointers->clear();

    for (const MsScanInfo &msScanInfo : msScanInfos) {

        if (msScanInfo.msLevel < 2) {
            continue;
        }

        QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
        e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
                msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower,
                msScanInfo.precursorTargetMz + msScanInfo.isoWindowLower,
                selectionFraction,
                &targetDecoyPointers
                ); ree;

        uniqueInfoScanKeyVsTargetDecoyCandidatePointers->insert(msScanInfo.mzTargetKey(), targetDecoyPointers);
        qDebug() << "MzTargetKey" << msScanInfo.mzTargetKey() << targetDecoyPointers.size() << "targetDecoyPairs found";
    }

    ERR_RETURN
}


