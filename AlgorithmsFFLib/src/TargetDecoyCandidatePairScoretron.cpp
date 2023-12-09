//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "CandidateScorertron.h"
#include "FeatureFinderHillBuilder.h"
#include "MsCalibratomatic.h"
//#include "MsFrame.h"
#include "ParallelUtils.h"

#include <QtConcurrent/QtConcurrent>

#include <iostream>


TargetDecoyCandidatePairScoretron::TargetDecoyCandidatePairScoretron()
: m_diaTargetFrames(nullptr)
, m_msReaderPointerAcc(nullptr)
{}


class TargetDecoyPairParallelInput {

public:
    MzTargetKey mzTargetKey;
    MsCalibratomatic msCalibratomatic;
    QMap<ScanNumber, ScanPoints*> *diaTargetFrame = nullptr;
    QMap<ScanNumber, ScanPoints> ms1Frame;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    int topNMs2Ions = -1.0;
    PythiaParameters pythiaParameters;
};

Err TargetDecoyCandidatePairScoretron::init(
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanTimeMS1,
        MsReaderPointerAcc *msReaderPointerAcc,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*diaTargetFrames); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTimeMS1); ree;

    m_pythiaParameters = pythiaParameters;
    m_msReaderPointerAcc = msReaderPointerAcc;
    m_diaTargetFrames = diaTargetFrames;
    m_ms1Frame = scanNumberVsScanTimeMS1;

    ERR_RETURN
}

namespace {

//    Err extractMS2Ions(
//            const QVector<MS2Ion> &ms2Ions,
//            int scanNumberMin,
//            int scanNumberMax,
//            double ppmTol,
//            TurboXIC *turboXic,
//            QMap<MzHashed, XICPoints> *cachedPoints,
//            QMap<MzHashed, XICPoints> *xicPointMap
//            ) {
//
//        ERR_INIT
//
//        for (const MS2Ion &ms2Ion : ms2Ions) {
//
//            const double mzTol = MathUtils::calculatePPM(ms2Ion.mz, ppmTol);
//            const double mzMin = ms2Ion.mz - mzTol;
//            const double mzMax = ms2Ion.mz + mzTol;
//
//            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
//            if (cachedPoints->contains(mzHashed)) {
//                xicPointMap->insert(mzHashed, cachedPoints->value(mzHashed));
//                continue;
//            }
//
//            XICPoints xicPoints = turboXic->extractPointsXIC(
//                    mzMin,
//                    mzMax,
//                    scanNumberMin,
//                    scanNumberMax
//            );
//
//            cachedPoints->insert(mzHashed, xicPoints);
//            xicPointMap->insert(mzHashed, xicPoints);
//
//        }
//
//        ERR_RETURN
//    }
//
//    Err extractMS2Ions(
//            const QVector<MS2Ion> &ms2Ions,
//            FrameIndex frameIndexMin,
//            FrameIndex frameIndexMax,
//            double ppmTol,
//            TurboXIC *turboXic,
//            QMap<MzHashed, XICPoints> *xicPointMap
//    ) {
//
//        ERR_INIT
//
//        QMap<MzHashed, XICPoints> cachedPoints;
//        e = extractMS2Ions(
//                ms2Ions,
//                frameIndexMin,
//                frameIndexMax,
//                ppmTol,
//                turboXic,
//                &cachedPoints,
//                xicPointMap
//        ); ree;
//
//        ERR_RETURN
//    }
//
//    Err extractScores(
//            const TargetDecoyCandidatePair* targetDecoyCandidatePair,
//            const QVector<MS2Ion> &ms2IonsTargetOrDecoyTheoretical,
//            int topNMs2Ions,
//            double ppmTol,
//            double iRT,
//            double scanTimeWindowMinutes,
//            MsFrame *msFrame,
//            MsCalibratomatic *msCalibratomatic,
//            TurboXIC *turboXic,
//            CandidateScorertron *candidateScorertron,
//            QMap<MzHashed, XICPoints> *cachedPoints,
//            CandidateScores *candidateScores
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(ms2IonsTargetOrDecoyTheoretical); ree;
//
//        QVector<MS2Ion> ms2IonsTheoretical = ms2IonsTargetOrDecoyTheoretical;
//        const int topNTarget = std::min(topNMs2Ions, ms2IonsTheoretical.size());
//        ms2IonsTheoretical.resize(topNTarget);
//
//        QMap<MzHashed, XICPoints> mzHashedVsXICPoints;
//
//        QVector<MS2Ion> ms2IonsTheoreticalIsotopeShadows;
//        QMap<MzHashed, XICPoints> mzHashedVsXICPointsIsotopeShadows;
//
//        ScanTime scanTimePredicted = -1;
//        if (msCalibratomatic->isInit()) {
//
//            std::transform(
//                    ms2IonsTargetOrDecoyTheoretical.begin(),
//                    ms2IonsTargetOrDecoyTheoretical.end(),
//                    std::back_inserter(ms2IonsTheoreticalIsotopeShadows),
//                    [](const MS2Ion &ms2Ion){
//                        const double isoChargeDiff = S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge;
//                        MS2Ion ms2IonNew = ms2Ion;
//                        ms2IonNew.mz -= isoChargeDiff;
//                        return ms2IonNew;
//                    }
//                    );
//
//            e = msCalibratomatic->predictScanTime(
//                    iRT,
//                    &scanTimePredicted
//            ); ree;
//
//            FrameIndex frameIndexPredictedMin;
//            e = msFrame->frameIndexFromScanTime(
//                    scanTimePredicted - scanTimeWindowMinutes,
//                    &frameIndexPredictedMin
//                    ); ree;
//
//            FrameIndex frameIndexPredictedMax;
//            e = msFrame->frameIndexFromScanTime(
//                    scanTimePredicted + scanTimeWindowMinutes,
//                    &frameIndexPredictedMax
//                    ); ree;
//
//            e = extractMS2Ions(
//                    ms2IonsTheoretical,
//                    frameIndexPredictedMin,
//                    frameIndexPredictedMax,
//                    ppmTol,
//                    turboXic,
//                    &mzHashedVsXICPoints
//                    ); ree;
//
//            e = extractMS2Ions(
//                    ms2IonsTheoreticalIsotopeShadows,
//                    frameIndexPredictedMin,
//                    frameIndexPredictedMax,
//                    ppmTol,
//                    turboXic,
//                    &mzHashedVsXICPointsIsotopeShadows
//            ); ree;
//
//        } else {
//
//            double frameIndexRtreeMin;
//            double frameIndexRtreeMax;
//            double mzRtreeMin;
//            double mzRtreeMax;
//            e = turboXic->getRTreeLimits(
//                    &frameIndexRtreeMin,
//                    &frameIndexRtreeMax,
//                    &mzRtreeMin,
//                    &mzRtreeMax
//                    ); ree;
//
//            e = extractMS2Ions(
//                    ms2IonsTheoretical,
//                    static_cast<FrameIndex>(frameIndexRtreeMin),
//                    static_cast<FrameIndex>(frameIndexRtreeMax),
//                    ppmTol,
//                    turboXic,
//                    cachedPoints,
//                    &mzHashedVsXICPoints
//                    ); ree;
//        }
//
//        e = candidateScorertron->calculateScores(
//                targetDecoyCandidatePair,
//                ms2IonsTheoretical,
//                mzHashedVsXICPoints,
//                ms2IonsTheoreticalIsotopeShadows,
//                mzHashedVsXICPointsIsotopeShadows,
//                scanTimePredicted,
//                msFrame,
//                candidateScores
//                ); ree;
//
//        ERR_RETURN
//    }

    Err parallelScoreLogic(const TargetDecoyPairParallelInput &pi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        MsCalibratomatic msCalibratomatic = pi.msCalibratomatic;

//        const int tranchSize = std::min(pi.pythiaParameters.trancheSizeMax, pi.targetDecoyPointers.size());
//        const int nTranches = pi.targetDecoyPointers.size() / tranchSize;
//
//        QVector<QVector<TargetDecoyCandidatePair*>> scoredTargetDecoyPointersTranched;
//        e = ParallelUtils::trancheVectorForParallelization(
//                pi.targetDecoyPointers,
//                nTranches,
//                &scoredTargetDecoyPointersTranched
//        ); ree;

        FeatureFinderParameters featureFinderParameters(pi.pythiaParameters);

        FeatureFinderHillBuilder featureFinderHillBuilderMS2;
        e = featureFinderHillBuilderMS2.init(featureFinderParameters); ree;
        e = featureFinderHillBuilderMS2.buildHills(*pi.diaTargetFrame); ree;

        QMap<FrameIndex, ScanPoints> ms1Frame = pi.ms1Frame;
        QMap<FrameIndex, ScanPoints*> ms1FramePntrs;
        for (auto it = ms1Frame.begin(); it != ms1Frame.end(); it++) {
            ms1FramePntrs.insert(it.key(), &it.value());
        }

        FeatureFinderHillBuilder featureFinderHillBuilderMS1;
        e = featureFinderHillBuilderMS1.init(featureFinderParameters); ree;
        e = featureFinderHillBuilderMS1.buildHills(ms1FramePntrs); ree;

        CandidateScorertron candidateScorertron;
        e = candidateScorertron.init(
                pi.scanNumberVsScanTime,
                pi.pythiaParameters,
                pi.topNMs2Ions,
                &msCalibratomatic,
                &featureFinderHillBuilderMS2,
                &featureFinderHillBuilderMS2
                ); ree;

        for (TargetDecoyCandidatePair* tdcp : pi.targetDecoyPointers) {

            CandidateScores candidateScoresTarget;
            e = candidateScorertron.calculateScores(
                    tdcp,
                    tdcp->ms2IonsTarget(),
                    &candidateScoresTarget
                    ); ree;

            CandidateScores candidateScoresDecoy;
            e = candidateScorertron.calculateScores(
                    tdcp,
                    tdcp->ms2IonsDecoy(),
                    &candidateScoresDecoy
                    ); ree;

        }

        if (pi.pythiaParameters.verbosity >= 1) {
            qDebug() << "Target key processed in" << pi.mzTargetKey << et.elapsed() << "mSec";
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
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
        QVector<CandidateScores> *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;

    QVector<TargetDecoyPairParallelInput> parallelInputs;
    e = buildParallelInput(
            topNMS2Ions,
            msCalibratomatic,
            mzTargetKeyVsTargetDecoyCandidatePointers,
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

    ERR_RETURN
}

bool TargetDecoyCandidatePairScoretron::isInit() {

    return m_pythiaParameters.isValid()
            && !m_diaTargetFrames->isEmpty()
            && !m_ms1Frame.isEmpty();
}

Err TargetDecoyCandidatePairScoretron::buildParallelInput(
        int topNMS2Ions,
        const MsCalibratomatic &msCalibratomatic,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
        QVector<TargetDecoyPairParallelInput> *input
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(mzTargetKeyVsTargetDecoyCandidatePointers->isEmpty()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isNotEmpty(m_ms1Frame); ree;

    for (const MsScanInfo &msScanInfo : m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos()) {

        TargetDecoyPairParallelInput tdppi;
        tdppi.topNMs2Ions = topNMS2Ions;
        tdppi.diaTargetFrame = &(*m_diaTargetFrames)[msScanInfo.mzTargetKey()];
        tdppi.mzTargetKey = msScanInfo.mzTargetKey();
        tdppi.msCalibratomatic = msCalibratomatic;
        tdppi.scanNumberVsScanTime = m_msReaderPointerAcc->ptr->getScanNumberVsScanTime();
        tdppi.pythiaParameters = m_pythiaParameters;
        tdppi.ms1Frame = m_ms1Frame;
        tdppi.targetDecoyPointers = mzTargetKeyVsTargetDecoyCandidatePointers->value(tdppi.mzTargetKey);

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
