//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

class TargetDecoyPairParallelInput {

public:
    UniqueMsInfoScanKey msInfoScanKey;
    QMap<ScanNumber, ScanPoints>* diaTargetFrame;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    double ppmTol = -1.0;
    int topNMs2Ions = -1.0;
};


Err TargetDecoyCandidatePairScoretron::init(
        const PythiaParameters &pythiaParameters,
        const QVector<MsScanInfo> &msScanInfos,
        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(msScanInfos); ree;
    e = ErrorUtils::isNotEmpty(*diaTargetFrames); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairManager->isInit()); ree;

    m_pythiaParameters = pythiaParameters;
    m_msScanInfos = msScanInfos;
    m_diaTargetFrames = diaTargetFrames;
    m_targetDecoyCandidatePairManager = targetDecoyCandidatePairManager;

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron::scoreTargetDecoyPairs() {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(m_msScanInfos); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;


    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron::buildParallelInput(
        double randomSelectionFraction,
        QVector<TargetDecoyPairParallelInput> *input
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(m_msScanInfos); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;

    for (const MsScanInfo &msScanInfo : m_msScanInfos) {

        const double mzMin = msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower;
        const double mzMax = msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper;

        TargetDecoyPairParallelInput tdppi;
        tdppi.topNMs2Ions = m_pythiaParameters.topNMs2Ions;
        tdppi.ppmTol = m_pythiaParameters.ms2ExtractionWidthPPM;
        tdppi.diaTargetFrame = &(*m_diaTargetFrames)[msScanInfo.targetScanKey()];
        tdppi.msInfoScanKey = msScanInfo.targetScanKey();
        e = m_targetDecoyCandidatePairManager->getTargetDecoyCandidatePairPointers(
                mzMin,
                mzMax,
                &tdppi.targetDecoyPointers
        ); ree;

        input->push_back(tdppi);

    }

    ERR_RETURN
}
