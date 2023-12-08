//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "CandidateScorertron.h"
#include "MsCalibratomatic.h"
#include "MsFrame.h"
#include "ParallelUtils.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>


TargetDecoyCandidatePairScoretron::TargetDecoyCandidatePairScoretron()
: m_diaTargetFrames(nullptr)
, m_msReaderPointerAcc(nullptr)
, m_targetDecoyCandidatePairManager(nullptr)
{}


class TargetDecoyPairParallelInput {

public:
    MzTargetKey msInfoScanKey;
    MsCalibratomatic msCalibratomatic;
    QMap<ScanNumber, ScanPoints*> *diaTargetFrame = nullptr;
    QMap<ScanNumber, ScanPoints*> ms1Frame;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    int topNMs2Ions = -1.0;
    PythiaParameters pythiaParameters;
};

Err TargetDecoyCandidatePairScoretron::init(
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanTimeMS1,
        MsReaderPointerAcc *msReaderPointerAcc,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*diaTargetFrames); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTimeMS1); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairManager->isInit()); ree;

    m_pythiaParameters = pythiaParameters;
    m_msReaderPointerAcc = msReaderPointerAcc;
    m_diaTargetFrames = diaTargetFrames;
    m_ms1Frame = scanNumberVsScanTimeMS1;
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
            FrameIndex frameIndexMin,
            FrameIndex frameIndexMax,
            double ppmTol,
            TurboXIC *turboXic,
            QMap<MzHashed, XICPoints> *xicPointMap
    ) {

        ERR_INIT

        QMap<MzHashed, XICPoints> cachedPoints;
        e = extractMS2Ions(
                ms2Ions,
                frameIndexMin,
                frameIndexMax,
                ppmTol,
                turboXic,
                &cachedPoints,
                xicPointMap
        ); ree;

        ERR_RETURN
    }

    Err extractScores(
            const TargetDecoyCandidatePair* targetDecoyCandidatePair,
            const QVector<MS2Ion> &ms2IonsTargetOrDecoyTheoretical,
            int topNMs2Ions,
            double ppmTol,
            double iRT,
            double scanTimeWindowMinutes,
            MsFrame *msFrame,
            MsCalibratomatic *msCalibratomatic,
            TurboXIC *turboXic,
            CandidateScorertron *candidateScorertron,
            QMap<MzHashed, XICPoints> *cachedPoints,
            CandidateScores *candidateScores
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2IonsTargetOrDecoyTheoretical); ree;

        QVector<MS2Ion> ms2IonsTheoretical = ms2IonsTargetOrDecoyTheoretical;
        const int topNTarget = std::min(topNMs2Ions, ms2IonsTheoretical.size());
        ms2IonsTheoretical.resize(topNTarget);

        QMap<MzHashed, XICPoints> mzHashedVsXICPoints;

        QVector<MS2Ion> ms2IonsTheoreticalIsotopeShadows;
        QMap<MzHashed, XICPoints> mzHashedVsXICPointsIsotopeShadows;

        ScanTime scanTimePredicted = -1;
        if (msCalibratomatic->isInit()) {

            std::transform(
                    ms2IonsTargetOrDecoyTheoretical.begin(),
                    ms2IonsTargetOrDecoyTheoretical.end(),
                    std::back_inserter(ms2IonsTheoreticalIsotopeShadows),
                    [](const MS2Ion &ms2Ion){
                        const double isoChargeDiff = S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge;
                        MS2Ion ms2IonNew = ms2Ion;
                        ms2IonNew.mz -= isoChargeDiff;
                        return ms2IonNew;
                    }
                    );

            e = msCalibratomatic->predictScanTime(
                    iRT,
                    &scanTimePredicted
            ); ree;

            FrameIndex frameIndexPredictedMin;
            e = msFrame->frameIndexFromScanTime(
                    scanTimePredicted - scanTimeWindowMinutes,
                    &frameIndexPredictedMin
                    ); ree;

            FrameIndex frameIndexPredictedMax;
            e = msFrame->frameIndexFromScanTime(
                    scanTimePredicted + scanTimeWindowMinutes,
                    &frameIndexPredictedMax
                    ); ree;

            e = extractMS2Ions(
                    ms2IonsTheoretical,
                    frameIndexPredictedMin,
                    frameIndexPredictedMax,
                    ppmTol,
                    turboXic,
                    &mzHashedVsXICPoints
                    ); ree;

            e = extractMS2Ions(
                    ms2IonsTheoreticalIsotopeShadows,
                    frameIndexPredictedMin,
                    frameIndexPredictedMax,
                    ppmTol,
                    turboXic,
                    &mzHashedVsXICPointsIsotopeShadows
            ); ree;

        } else {

            double frameIndexRtreeMin;
            double frameIndexRtreeMax;
            double mzRtreeMin;
            double mzRtreeMax;
            e = turboXic->getRTreeLimits(
                    &frameIndexRtreeMin,
                    &frameIndexRtreeMax,
                    &mzRtreeMin,
                    &mzRtreeMax
                    ); ree;

            e = extractMS2Ions(
                    ms2IonsTheoretical,
                    static_cast<FrameIndex>(frameIndexRtreeMin),
                    static_cast<FrameIndex>(frameIndexRtreeMax),
                    ppmTol,
                    turboXic,
                    cachedPoints,
                    &mzHashedVsXICPoints
                    ); ree;
        }

        e = candidateScorertron->calculateScores(
                targetDecoyCandidatePair,
                ms2IonsTheoretical,
                mzHashedVsXICPoints,
                ms2IonsTheoreticalIsotopeShadows,
                mzHashedVsXICPointsIsotopeShadows,
                scanTimePredicted,
                msFrame,
                candidateScores
                ); ree;

        ERR_RETURN
    }

    Err parallelScoreLogic(const TargetDecoyPairParallelInput &pi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        CandidateScorertron candidateScorertron;
        e = candidateScorertron.init(
                pi.ms1Frame,
                pi.scanNumberVsScanTime,
                pi.pythiaParameters,
                pi.topNMs2Ions
                ); ree;

        MsFrame msFrame;
        e = msFrame.init(
                *pi.diaTargetFrame,
                pi.scanNumberVsScanTime
                ); ree;

        TurboXIC turboXic;
        e = turboXic.init(msFrame.frameIndexVsScanPoints()); ree;

        MsCalibratomatic msCalibratomatic = pi.msCalibratomatic;

        const int tranchSize = std::min(pi.pythiaParameters.trancheSizeMax, pi.targetDecoyPointers.size());
        const int nTranches = pi.targetDecoyPointers.size() / tranchSize;

        QVector<QVector<TargetDecoyCandidatePair*>> scoredTargetDecoyPointersTranched;
        e = ParallelUtils::trancheVectorForParallelization(
                pi.targetDecoyPointers,
                nTranches,
                &scoredTargetDecoyPointersTranched
        ); ree;

        for (const QVector<TargetDecoyCandidatePair*> &tranche : scoredTargetDecoyPointersTranched) {

            //NOTE: this needs to stay outside of the subsequent loop or code becomes slow because cache is reset.
            QMap<MzHashed, XICPoints> cachedPoints;

            for (TargetDecoyCandidatePair* targetDecoyPtr : tranche) {

                CandidateScores candidateScoresTarget;
                e = extractScores(
                        targetDecoyPtr,
                        targetDecoyPtr->ms2IonsTarget(),
                        pi.topNMs2Ions,
                        pi.pythiaParameters.ms2ExtractionWidthPPM,
                        targetDecoyPtr->iRt(),
                        pi.pythiaParameters.scanTimeWindowMinutes,
                        &msFrame,
                        &msCalibratomatic,
                        &turboXic,
                        &candidateScorertron,
                        &cachedPoints,
                        &candidateScoresTarget
                ); ree;

                CandidateScores candidateScoresDecoy;
                e = extractScores(
                        targetDecoyPtr,
                        targetDecoyPtr->ms2IonsDecoy(),
                        pi.topNMs2Ions,
                        pi.pythiaParameters.ms2ExtractionWidthPPM,
                        targetDecoyPtr->iRt(),
                        pi.pythiaParameters.scanTimeWindowMinutes,
                        &msFrame,
                        &msCalibratomatic,
                        &turboXic,
                        &candidateScorertron,
                        &cachedPoints,
                        &candidateScoresDecoy
                ); ree;

                candidateScoresDecoy.isDecoy = true;
                targetDecoyPtr->uniqueInfoScanKeyVsScoresDecoy()->insert(pi.msInfoScanKey, candidateScoresDecoy);
                targetDecoyPtr->uniqueInfoScanKeyVsScoresTarget()->insert(pi.msInfoScanKey, candidateScoresTarget);
            }
        }

        if (pi.pythiaParameters.verbosity > 1) {
            qDebug() << "Target key processed in" << pi.msInfoScanKey << et.elapsed() << "mSec";
        }
        
        ERR_RETURN
    }

    Err reorderParallelInputs(QVector<TargetDecoyPairParallelInput> *parallelInputs) {

        ERR_INIT
        e = ErrorUtils::isFalse(parallelInputs->isEmpty()); ree;

        QVector<QVector<TargetDecoyPairParallelInput>> parallelInputsTranched;
        e = ParallelUtils::trancheVectorForParallelization(
                *parallelInputs,
                ParallelUtils::numberOfAvailableSystemProcessors(),
                &parallelInputsTranched
                ); ree;

        parallelInputs->clear();

        for (const QVector<TargetDecoyPairParallelInput> &tranche : parallelInputsTranched) {

            for (const TargetDecoyPairParallelInput &pi : tranche) {
                parallelInputs->append(pi);
            }
        }

        ERR_RETURN
    }

}//namespace
Err TargetDecoyCandidatePairScoretron::scoreTargetDecoyPairs(
        int topNMS2Ions,
        const MsCalibratomatic &msCalibratomatic,
        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;

    e = m_targetDecoyCandidatePairManager->clearScores(); ree;

    QVector<TargetDecoyPairParallelInput> parallelInputs;
    e = buildParallelInput(
            *scoredTargetDecoyPointers,
            topNMS2Ions,
            msCalibratomatic,
            &parallelInputs
    ); ree;

    e = reorderParallelInputs(&parallelInputs); ree;

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

//    qDebug() << "Tranche" << ++tranchCounter << "Size" << tranche.size() << "mSec" << et.restart();



    ERR_RETURN
}

bool TargetDecoyCandidatePairScoretron::isInit() {

    return m_pythiaParameters.isValid()
            && !m_diaTargetFrames->isEmpty()
            && !m_ms1Frame.isEmpty()
            && m_targetDecoyCandidatePairManager->isInit();
}

Err TargetDecoyCandidatePairScoretron::buildParallelInput(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
        int topNMS2Ions,
        const MsCalibratomatic &msCalibratomatic,
        QVector<TargetDecoyPairParallelInput> *input
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairs); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isNotEmpty(m_ms1Frame); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;

    for (const MsScanInfo &msScanInfo : m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos()) {

        const double mzMin = msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower;
        const double mzMax = msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper;

        TargetDecoyPairParallelInput tdppi;
        tdppi.topNMs2Ions = topNMS2Ions;
        tdppi.diaTargetFrame = &(*m_diaTargetFrames)[msScanInfo.mzTargetKey()];
        tdppi.msInfoScanKey = msScanInfo.mzTargetKey();
        tdppi.msCalibratomatic = msCalibratomatic;
        tdppi.scanNumberVsScanTime = m_msReaderPointerAcc->ptr->getScanNumberVsScanTime();
        tdppi.pythiaParameters = m_pythiaParameters;
        tdppi.ms1Frame = m_ms1Frame;

        for (TargetDecoyCandidatePair *p : targetDecoyCandidatePairs) {

            const double mz = p->mz();

            if (mz < mzMin || mz > mzMax) {
                continue;
            }

            tdppi.targetDecoyPointers.push_back(p);
        }

        input->push_back(tdppi);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron::setPythiaParameters(const PythiaParameters &pythiaParameters) {
    ERR_INIT
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_pythiaParameters = pythiaParameters;
    ERR_RETURN
}
