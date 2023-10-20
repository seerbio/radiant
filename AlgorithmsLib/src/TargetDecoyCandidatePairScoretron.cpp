//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "CandidateScorertron.h"
#include "MsCalibratomatic.h"
#include "MsFrame.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>

class TargetDecoyPairParallelInput {

public:
    UniqueMsInfoScanKey msInfoScanKey;
    MsCalibratomatic msCalibratomatic;
    QMap<ScanNumber, ScanPoints>* diaTargetFrame = nullptr;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    int topNMs2Ions = -1.0;
    PythiaParameters pythiaParameters;
};


Err TargetDecoyCandidatePairScoretron::init(
        const PythiaParameters &pythiaParameters,
        MsReaderPointerAcc *msReaderPointerAcc,
        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*diaTargetFrames); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairManager->isInit()); ree;

    m_pythiaParameters = pythiaParameters;
    m_msReaderPointerAcc = msReaderPointerAcc;
    m_diaTargetFrames = diaTargetFrames;
    m_targetDecoyCandidatePairManager = targetDecoyCandidatePairManager;

    ERR_RETURN
}

namespace {

    Err extractMS2Ions(
            const QVector<MS2Ion> &ms2Ions,
            int scanNumberMin,
            int scanNumberMax,
            double ppmTol,
            TurboXIC *turboXic,
            QMap<MzHashed, XICPoints> *cachedPoints,
            QMap<MzHashed, XICPoints> *xicPointMap
            ) {

        ERR_INIT

        for (const MS2Ion &ms2Ion : ms2Ions) {

            const double mzTol = MathUtils::calculatePPM(ms2Ion.mz, ppmTol);
            const double mzMin = ms2Ion.mz - mzTol;
            const double mzMax = ms2Ion.mz + mzTol;

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            if (cachedPoints->contains(mzHashed)) {
                xicPointMap->insert(mzHashed, cachedPoints->value(mzHashed));
                continue;
            }

            XICPoints xicPoints = turboXic->extractPointsXIC(
                    mzMin,
                    mzMax,
                    scanNumberMin,
                    scanNumberMax
            );

            cachedPoints->insert(mzHashed, xicPoints);
            xicPointMap->insert(mzHashed, xicPoints);

        }

        ERR_RETURN
    }

    Err extractMS2Ions(
            const QVector<MS2Ion> &ms2Ions,
            int scanNumberMin,
            int scanNumberMax,
            double ppmTol,
            TurboXIC *turboXic,
            QMap<MzHashed, XICPoints> *xicPointMap
    ) {

        ERR_INIT

        QMap<MzHashed, XICPoints> cachedPoints;
        e = extractMS2Ions(
                ms2Ions,
                scanNumberMin,
                scanNumberMax,
                ppmTol,
                turboXic,
                &cachedPoints,
                xicPointMap
        ); ree;

        ERR_RETURN
    }

    Err extractScores(
            const QVector<MS2Ion> &ms2IonsTargetOrDecoy,
            int topNMs2Ions,
            double ppmTol,
            double iRT,
            MsCalibratomatic *msCalibratomatic,
            TurboXIC *turboXic,
            CandidateScorertron *candidateScorertron,
            QMap<MzHashed, XICPoints> *cachedPoints
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2IonsTargetOrDecoy); ree;

        double scanNumberRtreeMin;
        double scanNumberRtreeMax;
        double mzRtreeMin;
        double mzRtreeMax;
        e = turboXic->getRTreeLimits(
                &scanNumberRtreeMin,
                &scanNumberRtreeMax,
                &mzRtreeMin,
                &mzRtreeMax
        ); ree;

        QVector<MS2Ion> ms2IonsTarget = ms2IonsTargetOrDecoy;
        const int topNTarget = std::min(topNMs2Ions, ms2IonsTarget.size());
        ms2IonsTarget.resize(topNTarget);

        QMap<MzHashed, XICPoints> mzHashedVsXICPoints;

        if (msCalibratomatic->isInit()) {

            ScanTime scanTimePredicted = -1;
            e = msCalibratomatic->predictScanTime(
                    iRT,
                    &scanTimePredicted
            ); ree;

            // TODO set these when calibration is working
            ScanNumber scanNumberPredictedMin = -1;
            ScanNumber scanNumberPredictedMax = -1;

            e = extractMS2Ions(
                    ms2IonsTarget,
                    scanNumberPredictedMin,
                    scanNumberPredictedMax,
                    ppmTol,
                    turboXic,
                    &mzHashedVsXICPoints
            ); ree;

        } else {

            e = extractMS2Ions(
                    ms2IonsTarget,
                    static_cast<ScanNumber>(scanNumberRtreeMin),
                    static_cast<ScanNumber>(scanNumberRtreeMax),
                    ppmTol,
                    turboXic,
                    cachedPoints,
                    &mzHashedVsXICPoints
            ); ree;

        }

        e = candidateScorertron->calculateScores(
                mzHashedVsXICPoints,
                ms2IonsTarget
        ); ree;

        ERR_RETURN
    }

    Err parallelScoreLogic(const TargetDecoyPairParallelInput &pi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        CandidateScorertron candidateScorertron;
        e = candidateScorertron.init(pi.pythiaParameters, pi.topNMs2Ions); ree;

        MsFrame msFrame;
        e = msFrame.init(
                *pi.diaTargetFrame,
                pi.scanNumberVsScanTime
                ); ree;

        TurboXIC turboXic;
        e = turboXic.init(msFrame.frameIndexVsScanPoints()); ree;

        MsCalibratomatic msCalibratomatic = pi.msCalibratomatic;

        //NOTE: this needs to stay outside of loop or code becomes slow because cache is reset.
        QMap<MzHashed, XICPoints> cachedPoints;

        for (const TargetDecoyCandidatePair* targetDecoyPtr : pi.targetDecoyPointers) {

            e = extractScores(
                    targetDecoyPtr->ms2IonsTarget(),
                    pi.pythiaParameters.topNMs2Ions,
                    pi.pythiaParameters.ms2ExtractionWidthPPM,
                    targetDecoyPtr->iRt(),
                    &msCalibratomatic,
                    &turboXic,
                    &candidateScorertron,
                    &cachedPoints
                    ); ree;

            e = extractScores(
                    targetDecoyPtr->ms2IonsDecoy(),
                    pi.pythiaParameters.topNMs2Ions,
                    pi.pythiaParameters.ms2ExtractionWidthPPM,
                    targetDecoyPtr->iRt(),
                    &msCalibratomatic,
                    &turboXic,
                    &candidateScorertron,
                    &cachedPoints
            ); ree;

        }

        qDebug() << "Target key processed in" << pi.msInfoScanKey << et.elapsed() << "mSec";

        ERR_RETURN
    }

}//namespace
Err TargetDecoyCandidatePairScoretron::scoreTargetDecoyPairs(
        int topNMS2Ions,
        double randomSelectionFraction,
        const MsCalibratomatic &msCalibratomatic,
        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;

    QVector<TargetDecoyPairParallelInput> parallelInputs;
    e = buildParallelInput(
            topNMS2Ions,
            randomSelectionFraction,
            msCalibratomatic,
            &parallelInputs
            ); ree;

    //TODO use for loop here to process

#define PARALLEL_SCORE
#ifdef PARALLEL_SCORE
    QFuture<Err> futures = QtConcurrent::mapped(
            parallelInputs,
            parallelScoreLogic
    );
    futures.waitForFinished();

    for (Err res : futures) {
        e = res; ree;
    }
#else
    for(const TargetDecoyPairParallelInput &tdppi : parallelInputs) {
        e = parallelScoreLogic(tdppi); ree;
    }
#endif

    for (const TargetDecoyPairParallelInput &tdppi : parallelInputs) {
        scoredTargetDecoyPointers->append(tdppi.targetDecoyPointers);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron::buildParallelInput(
        int topNMS2Ions,
        double randomSelectionFraction,
        const MsCalibratomatic &msCalibratomatic,
        QVector<TargetDecoyPairParallelInput> *input
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;

    for (const MsScanInfo &msScanInfo : m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos()) {

        const double mzMin = msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower;
        const double mzMax = msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper;

        TargetDecoyPairParallelInput tdppi;
        tdppi.topNMs2Ions = topNMS2Ions;
        tdppi.diaTargetFrame = &(*m_diaTargetFrames)[msScanInfo.targetScanKey()];
        tdppi.msInfoScanKey = msScanInfo.targetScanKey();
        tdppi.msCalibratomatic = msCalibratomatic;
        tdppi.scanNumberVsScanTime = m_msReaderPointerAcc->ptr->getScanNumberVsScanTime();
        tdppi.pythiaParameters = m_pythiaParameters;

        if (randomSelectionFraction < 0) {
            e = m_targetDecoyCandidatePairManager->getTargetDecoyCandidatePairPointers(
                    mzMin,
                    mzMax,
                    &tdppi.targetDecoyPointers
            ); ree;

        } else {
            e = m_targetDecoyCandidatePairManager->getTargetDecoyCandidatePairPointers(
                    mzMin,
                    mzMax,
                    randomSelectionFraction,
                    &tdppi.targetDecoyPointers
            ); ree;

        }

        input->push_back(tdppi);

    }

    ERR_RETURN
}
