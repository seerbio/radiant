//
// Created by andrewnichols on 9/27/24.
//

#include "MsCalibratomaticSettertron.h"

#include "ErrorUtils.h"
#include "DiscriminantScoretron.h"
#include "FDRCLassifierNeuralNet.h"
#include "InterferingCandidatesEliminatomatic.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
#include "PythiaDIAFFWorkflowSharedMethods.h"
#include "QValueSettertron.h"
#include "TurboXIC.h"

MsCalibratomaticSettertron::MsCalibratomaticSettertron()
: m_msReaderPointerAcc(nullptr)
, m_targetDecoyCandidatePairManager(nullptr)
, m_targetDecoyCandidatePairScoretron(nullptr)
, m_scanNumberVsFeatureFinderHillBuildersPntrsTIMS(nullptr)
, m_pythiaParameters(nullptr)
{}

Err MsCalibratomaticSettertron::init(
    PythiaParameters *pythiaParameters,
    MsReaderPointerAcc *msReaderPointerAcc,
    TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager,
    TargetDecoyCandidatePairScoretron2 *targetDecoyCandidatePairScoretron,
    QMap<ScanNumber, FeatureFinderHillBuilder*> *scanNumberVsFeatureFinderHillBuildersPntrsTIMS
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairManager->isInit()); ree;
    e = ErrorUtils::isTrue(pythiaParameters->isValid()); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairScoretron->isInit()); ree;

    if (msReaderPointerAcc->ptr->isTIMS()) {
        e = ErrorUtils::isFalse(scanNumberVsFeatureFinderHillBuildersPntrsTIMS->isEmpty()); ree;
    }

    m_msReaderPointerAcc = msReaderPointerAcc;
    m_targetDecoyCandidatePairManager = targetDecoyCandidatePairManager;
    m_pythiaParameters = pythiaParameters;
    m_targetDecoyCandidatePairScoretron = targetDecoyCandidatePairScoretron;
    m_scanNumberVsFeatureFinderHillBuildersPntrsTIMS = scanNumberVsFeatureFinderHillBuildersPntrsTIMS;

    e = m_targetDecoyCandidatePairManager->getTargetDecoyCandidatePairPointers(
        &m_targetDecoyPairPntrs
        ); ree

    ERR_RETURN
}

namespace {

    void filterMs1CandidateRowsByCorr(QVector<CandidateScores*> *candidateScoresMS1Cal) {

        constexpr double cosineSimSumMS1Min = 0.9;
        const auto terminatorLogic = [cosineSimSumMS1Min](const CandidateScores *cs){
            return cs->featuresArray[Features::CosineSim100MS1] <= cosineSimSumMS1Min;
        };

        const auto terminator = std::remove_if(
                candidateScoresMS1Cal->begin(),
                candidateScoresMS1Cal->end(),
                terminatorLogic
                );

        candidateScoresMS1Cal->erase(terminator, candidateScoresMS1Cal->end());
    }

}//namespace
Err MsCalibratomaticSettertron::buildCalibration(MsCalibratomatic *msCalibratomatic) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters->isValid()); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron->isInit()); ree;

    const int numberOfTranches = calculateNumberOfTranches();

    const QVector<MsScanInfo> uniqueMsScanInfos = m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    constexpr int maxUniqueScanInfosTrainingCount = 16;
    constexpr int offset = 0;
    QVector<MsScanInfo> uniqueMsScanInfosCalibration;
    e = PythiaDIAFFWorkflowSharedMethods::buildUniqueMsScanInfosForProcessing(
        uniqueMsScanInfos,
        maxUniqueScanInfosTrainingCount,
        offset,
        &uniqueMsScanInfosCalibration
        ); ree;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFramesPntrs
                                        = m_targetDecoyCandidatePairScoretron->diaTargetFrames();
    if (m_pythiaParameters->useLazyLoading) {
        e = PythiaDIAFFWorkflowSharedMethods::buildDiaTargetFrames(
            uniqueMsScanInfosCalibration,
            m_msReaderPointerAcc,
            &m_msCalibratomatic,
            &m_diaTargetFrames
            ); ree;
        for (auto it = m_diaTargetFrames.begin(); it != m_diaTargetFrames.end(); it++) {
            const MzTargetKey &mzTargetKey = it.key();
            QMap<ScanNumber, ScanPoints> &scanNumberVsScanPoints = it.value();
            for (auto itt = scanNumberVsScanPoints.begin(); itt != scanNumberVsScanPoints.end(); itt++) {
                m_diaTargetFramesPntrs[mzTargetKey].insert(itt.key(), &itt.value());
            }
        }
        *diaTargetFramesPntrs = m_diaTargetFramesPntrs;
        e = m_targetDecoyCandidatePairScoretron->buildMzTargetKeyVsMsFrames(); ree;
    }

    QMap<MzTargetKey, TurboXIC*> mzTargetKeyVsTurboXicPntrs;
    e = PythiaDIAFFWorkflowSharedMethods::buildMzTargetKeyVsTurboXicPntrs(
        uniqueMsScanInfosCalibration,
        m_msReaderPointerAcc->ptr->getScanNumberVsScanTime(),
        diaTargetFramesPntrs,
        &mzTargetKeyVsTurboXicPntrs
        ); ree;

    QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePointersTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            m_targetDecoyPairPntrs,
            numberOfTranches,
            &targetDecoyCandidatePointersTranched
            ); ree;

    int batchCounter = 0;
    for (const QVector<TargetDecoyCandidatePair*> &tdcp : targetDecoyCandidatePointersTranched) {

        QElapsedTimer etBatch;
        etBatch.start();

        QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsBatch = m_targetDecoyCandidatePairsTopScores;
        targetDecoyCandidatePairsBatch.append(tdcp);

        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
        e = PythiaDIAFFWorkflowSharedMethods::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
                targetDecoyCandidatePairsBatch,
                *m_pythiaParameters,
                uniqueMsScanInfosCalibration,
                &mzTargetKeyVsTargetDecoyCandidatePointers
                ); ree;

        constexpr int topNMS2IonsCalibration = 6;
        constexpr bool useExtendedScores = false;
        constexpr bool useNeuralNetworkScores = false;
        constexpr bool useTopNIntegrationsParameter = false;

        m_weights = m_weights.isEmpty()
                  ? DiscriminantScoretron::defaultWeights(useExtendedScores, useNeuralNetworkScores)
                  : m_weights;

        constexpr int splitter = 2;
        const int threadCount = uniqueMsScanInfos.size() < m_pythiaParameters->threadCount
                      ? std::min(uniqueMsScanInfos.size() * splitter, m_pythiaParameters->threadCount)
                      : m_pythiaParameters->threadCount;

        constexpr float minPeakCountCalibration = 4.9;
        m_candidateScorePairs.clear();
        e = m_targetDecoyCandidatePairScoretron->scoreTargetDecoyPairs(
                topNMS2IonsCalibration,
                m_msCalibratomatic,
                minPeakCountCalibration,
                threadCount,
                useExtendedScores,
                useNeuralNetworkScores,
                useTopNIntegrationsParameter,
                mzTargetKeyVsTurboXicPntrs,
                m_weights,
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                &m_candidateScorePairs
                ); ree

        if (m_msReaderPointerAcc->ptr->isTIMS()) {
            e = PythiaDIAFFWorkflowSharedMethods::assignIonMobilityValues(
                *m_pythiaParameters,
                &m_candidateScorePairs,
                m_scanNumberVsFeatureFinderHillBuildersPntrsTIMS
                ); ree;
        }

        QVector<CandidateScores*> candidateScoresVecBatchPntrs;
        QMap<int, int> fdrVsCounts;
        QVector<float> weights;
        e = PythiaDIAFFWorkflowSharedMethods::processBatch(
            m_candidateScorePairs,
            *m_pythiaParameters,
            useExtendedScores,
            useNeuralNetworkScores,
            &candidateScoresVecBatchPntrs,
            &fdrVsCounts,
            &weights
            ); ree;

        constexpr int fdrKey = 5;
        constexpr int fdrKeyMassCalMS2 = 2;
        constexpr int fdrKeyMassCalMS1 = 5;

        if (constexpr int minFDRCountForWeightsUpdate = 100; fdrVsCounts.value(fdrKey) >= minFDRCountForWeightsUpdate) {
            m_weights = weights;
        }

        e = honeIRTAndMassCalibration(
            &candidateScoresVecBatchPntrs,
            fdrVsCounts.value(fdrKey),
            fdrVsCounts.value(fdrKeyMassCalMS2)
            ); ree;

        m_scanTimeStDevs.push_back(m_msCalibratomatic.scanTimeStDev());
        m_ms2PPMStDevs.push_back(m_msCalibratomatic.mzStDevMS2());

        QString fdrString;
        e = FDRCLassifierNeuralNet::outPutFDRCounts(fdrVsCounts, &fdrString); ree;

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                << "Processed batch" << ++batchCounter
                << "of"
                << targetDecoyCandidatePointersTranched.size()
                << etBatch.elapsed()
                << "mSec"
                << qPrintable(fdrString);

        if (fdrVsCounts.value(fdrKey) > m_pythiaParameters->calibrationTrainingVolume
            || &tdcp == targetDecoyCandidatePointersTranched.cend() - 1
            ) {

            constexpr int fdrKey50PercentFDR = 50;
            candidateScoresVecBatchPntrs.resize(std::min(candidateScoresVecBatchPntrs.size(), fdrVsCounts.value(fdrKey50PercentFDR)));

            QVector<CandidateScores*> candidateScoresVecBatchPntrsRecal = candidateScoresVecBatchPntrs;
            candidateScoresVecBatchPntrsRecal.resize(std::min(candidateScoresVecBatchPntrsRecal.size(), fdrVsCounts.value(fdrKeyMassCalMS2)));

            QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
            e = PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
                    MSLevelEnum::MS2,
                    candidateScoresVecBatchPntrsRecal,
                    m_pythiaParameters->verbosity,
                    &msCalibrationReaderRows
                    ); ree;

            e = m_msCalibratomatic.setMassCalibrationCoeffs(
                msCalibrationReaderRows,
                MSLevelEnum::MS2
                ); ree;
            e = m_msCalibratomatic.setCalibrationCoeffsUsingAllMeans(); ree;

            m_scanTimeStDevs.push_back(m_msCalibratomatic.scanTimeStDev());
            m_ms2PPMStDevs.push_back(m_msCalibratomatic.mzStDevMS2());

            constexpr int minVecSize = 3;
            std::sort(m_scanTimeStDevs.begin(), m_scanTimeStDevs.end());
            if (m_scanTimeStDevs.size() >= minVecSize) {
                m_scanTimeStDevs.pop_front();
                m_scanTimeStDevs.pop_back();
            }

            std::sort(m_ms2PPMStDevs.begin(), m_ms2PPMStDevs.end());
            if (m_ms2PPMStDevs.size() >= minVecSize) {
                m_ms2PPMStDevs.pop_front();
                m_ms2PPMStDevs.pop_back();
            }

            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "ScanTimeWindow Mean|Median|Min"
                     << MathUtils::mean(m_scanTimeStDevs)
                     << MathUtils::median(m_scanTimeStDevs)
                     << *std::min({m_scanTimeStDevs.begin(), m_scanTimeStDevs.end()});

            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Ms2 ppm Mean|Median" << MathUtils::mean(m_ms2PPMStDevs) << MathUtils::median(m_ms2PPMStDevs);

            m_msCalibratomatic.setScanTimeStDev(m_scanTimeStDevs.front());
            m_msCalibratomatic.setMzStDevMS2(MathUtils::mean(m_ms2PPMStDevs));

            if (!m_pythiaParameters->useLazyLoading) {
                e = recalibrateMzVals(
                        MSLevelEnum::MS2,
                        m_targetDecoyCandidatePairScoretron->diaTargetFrames(),
                        m_targetDecoyCandidatePairScoretron->ms1ScanNumberVsScanPoints()
                        ); ree;
            }

            candidateScoresVecBatchPntrsRecal = candidateScoresVecBatchPntrs;
            candidateScoresVecBatchPntrsRecal.resize(std::min(candidateScoresVecBatchPntrsRecal.size(), fdrVsCounts.value(fdrKeyMassCalMS1)));
            filterMs1CandidateRowsByCorr(&candidateScoresVecBatchPntrsRecal);
            constexpr int recalibrationPointCountMin = 400;
            qDebug()<< qPrintable(S_GLOBAL_TIMER.elapsed())
                << candidateScoresVecBatchPntrsRecal.size()
                << "found for MS1 Recalibration";
            if (candidateScoresVecBatchPntrsRecal.size() < recalibrationPointCountMin) {
                qWarning() << qPrintable(S_GLOBAL_TIMER.elapsed())
                        << "Skipping MS1 recalibration.  Not enough points found";
                for (TurboXIC* turboXic : mzTargetKeyVsTurboXicPntrs) {delete turboXic;}
                *msCalibratomatic = m_msCalibratomatic;
                m_targetDecoyCandidatePairsTopScores.clear();
                m_entered.clear();
                m_candidateScorePairs.clear();
                ERR_RETURN
            }

            QVector<MsCalibarationReaderRow> msCalibrationReaderRowsMS1;
            e = PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
                    MSLevelEnum::MS1,
                    candidateScoresVecBatchPntrsRecal,
                    m_pythiaParameters->verbosity,
                    &msCalibrationReaderRowsMS1
                    ); ree;

            if (msCalibrationReaderRowsMS1.size() < recalibrationPointCountMin) {
                for (TurboXIC* turboXic : mzTargetKeyVsTurboXicPntrs) {delete turboXic;}
                *msCalibratomatic = m_msCalibratomatic;
                m_targetDecoyCandidatePairsTopScores.clear();
                m_entered.clear();
                m_candidateScorePairs.clear();
                ERR_RETURN
            }

            e = m_msCalibratomatic.setMassCalibrationCoeffs(
                msCalibrationReaderRowsMS1,
                MSLevelEnum::MS1
                ); ree;

            e = recalibrateMzVals(
                    MSLevelEnum::MS1,
                    m_targetDecoyCandidatePairScoretron->diaTargetFrames(),
                    m_targetDecoyCandidatePairScoretron->ms1ScanNumberVsScanPoints()
                    ); ree;

            e = m_targetDecoyCandidatePairScoretron->reloadTurboXICMS1(); ree;

            break;
        }

    }

    for (TurboXIC* turboXic : mzTargetKeyVsTurboXicPntrs) {delete turboXic;}
    m_targetDecoyCandidatePairsTopScores.clear();
    m_entered.clear();
    m_candidateScorePairs.clear();
    *msCalibratomatic = m_msCalibratomatic;
    ERR_RETURN
}

int MsCalibratomaticSettertron::calculateNumberOfTranches() const {

    const auto sizePerTranche = static_cast<double>(m_pythiaParameters->trancheSizeMax);
    const int numberOfTranches = std::max(static_cast<int>(m_targetDecoyPairPntrs.size() / sizePerTranche), 1);
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Number of tranches for calibration:" << numberOfTranches
             << "target count:" << m_targetDecoyPairPntrs.size()
             << "sizePerTranche:" << static_cast<int>(sizePerTranche);

    return numberOfTranches;

}


Err MsCalibratomaticSettertron::honeIRTAndMassCalibration(
    QVector<CandidateScores*> *candidateScoresVecScoredPntrs,
    int topNCandidates,
    int topCandidatesMass
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScoresVecScoredPntrs->isEmpty()); ree;

    PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(candidateScoresVecScoredPntrs);

    QVector<CandidateScores*> candidateScoresVecBatchPntrsResized = *candidateScoresVecScoredPntrs;

    constexpr int minTrainingCountTranche = 50;
    candidateScoresVecBatchPntrsResized.resize(std::max(minTrainingCountTranche, topNCandidates));

    e = InterferingCandidatesEliminatomatic::removeInterferingCandidates(
            m_pythiaParameters->ionsSharedToReject,
            m_pythiaParameters->mzMinMS2,
            m_pythiaParameters->mzMaxMS2,
            &candidateScoresVecBatchPntrsResized
            ); ree;

    if (m_pythiaParameters->verbosity > 0) {
        qDebug() << "Using" << candidateScoresVecBatchPntrsResized.size() << "for iRT Estimation";
        qDebug() << candidateScoresVecBatchPntrsResized.size() << "after filtering";
    }

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
            MSLevelEnum::MS2,
            candidateScoresVecBatchPntrsResized,
            m_pythiaParameters->verbosity,
            &msCalibrationReaderRows
            ); ree;

    for (const CandidateScores *cs : candidateScoresVecBatchPntrsResized) {

        if (m_entered.value(cs->targetDecoyCandidatePair)) {
            continue;
        }

        m_targetDecoyCandidatePairsTopScores.push_back(cs->targetDecoyCandidatePair);
        m_entered.insert(cs->targetDecoyCandidatePair, true);
    }

    e = m_msCalibratomatic.buildRTMapper(msCalibrationReaderRows); ree;

    if (m_pythiaParameters->verbosity > 0) {
        qDebug() << "----- scanTimeWindowStDev x" << m_pythiaParameters->scanTimeWindowStDevs
                 <<":" << m_msCalibratomatic.scanTimeStDev(m_pythiaParameters->scanTimeWindowStDevs);
    }

    constexpr int ms2MassRecalCountMin = 200;
    if (topCandidatesMass > ms2MassRecalCountMin) {

        msCalibrationReaderRows.resize(topCandidatesMass);

        e = m_msCalibratomatic.setMassCalibrationCoeffs(
            msCalibrationReaderRows,
            MSLevelEnum::MS2
            ); ree;
    }

    ERR_RETURN
}

namespace {

    Err recalibrationLogic(
        const MsCalibratomatic &msCalibratomatic,
        const QMap<ScanNumber, ScanPoints*> &scanPoints
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msCalibratomatic.isInitCalMS2()); ree;
        e = ErrorUtils::isFalse(scanPoints.isEmpty()); ree;

        e = msCalibratomatic.recalibrateScanPoints(
            MSLevelEnum::MS2,
            scanPoints
            ); ree;

        ERR_RETURN
    }

}//namespace
Err MsCalibratomaticSettertron::recalibrateMzVals(
        const MSLevelEnum &msLevel,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
        ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msCalibratomatic.isInitCalMS2() || m_msCalibratomatic.isInitCalMS1()); ree;
    e = ErrorUtils::isFalse(diaTargetFrames->isEmpty()); ree;
    e = ErrorUtils::isFalse(scanNumberVsScanTimeMS1->isEmpty()); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Recalibrating mz vals";

    QElapsedTimer et;
    et.start();

    if (msLevel == MSLevelEnum::MS2) {

#define RECAL_PARALLEL
#ifdef RECAL_PARALLEL
        const auto reCalBinder = std::bind(
                recalibrationLogic,
                m_msCalibratomatic,
                std::placeholders::_1
        );

        const QList<QMap<ScanNumber, ScanPoints*>> &tandemPointsVec = diaTargetFrames->values();

        QFuture<Err> futures = QtConcurrent::mapped(
                tandemPointsVec,
                reCalBinder
        );
        futures.waitForFinished();

        for (const Err &err: futures) {
            e = err; ree;
        }
#else
        const QList<MzTargetKey> &diaTargetFrameKeys = diaTargetFrames->keys();
        for (const MzTargetKey &k: diaTargetFrameKeys) {
            e = recalibrationLogic(m_msCalibratomatic, diaTargetFrames[k]); ree;
        }
#endif
    }

    if (msLevel == MSLevelEnum::MS1) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Recalibrating MS1 mz vals frame";
        e = m_msCalibratomatic.recalibrateScanPoints(
            MSLevelEnum::MS1,
            scanNumberVsScanTimeMS1
            ); ree;
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Points recalibrated in" << et.elapsed() << "mSec";

    ERR_RETURN
}
