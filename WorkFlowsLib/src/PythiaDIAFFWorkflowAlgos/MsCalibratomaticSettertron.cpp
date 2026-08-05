//
// Created by andrewnichols on 9/27/24.
//

#include "MsCalibratomaticSettertron.h"

#include "ErrorUtils.h"
#include "DiscriminantScoretron.h"
#include "FDRCLassifierNeuralNet.h"
#include "InterferingCandidatesEliminatomatic.h"
#include "MsReaderPointerAcc.h"
#include "ObjectCSVWriters.h"
#include "ParallelUtils.h"
#include "PythiaDIAFFWorkflowSharedMethods.h"
#include "QValueSettertron.h"
#include "TimsbukIndexTypes.h"
#include "TurboXIC.h"

#include <cmath>

MsCalibratomaticSettertron::MsCalibratomaticSettertron()
: m_msReaderPointerAcc(nullptr)
, m_targetDecoyCandidatePairManager(nullptr)
, m_targetDecoyCandidatePairScoretron(nullptr)
, m_pythiaParameters(nullptr)
, m_batchCounter(0)
, m_excludeDecoys(false)
, m_fdrWeightedMean(0)
{}

Err MsCalibratomaticSettertron::init(
    const QVector<Features> &featuresCalibration,
    PythiaParameters *pythiaParameters,
    MsReaderPointerAcc *msReaderPointerAcc,
    TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager,
    TargetDecoyCandidatePairScoretron2 *targetDecoyCandidatePairScoretron,
    bool excludeDecoys
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairManager->isInit()); ree;
    e = ErrorUtils::isTrue(pythiaParameters->isValid()); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairScoretron->isInit()); ree;
    e = ErrorUtils::isNotEmpty(featuresCalibration); ree;

    m_msReaderPointerAcc = msReaderPointerAcc;
    m_targetDecoyCandidatePairManager = targetDecoyCandidatePairManager;
    m_pythiaParameters = pythiaParameters;
    m_targetDecoyCandidatePairScoretron = targetDecoyCandidatePairScoretron;
    m_featuresCalibration = featuresCalibration;

    e = m_targetDecoyCandidatePairManager->getTargetDecoyCandidatePairPointers(
        &m_targetDecoyPairPntrs
        ); ree

    m_excludeDecoys = excludeDecoys;

    ERR_RETURN
}

namespace {

    constexpr double DEFAULT_ION_MOBILITY_TOLERANCE_ONE_OVER_K0 = 0.1;

    void appendIonMobilityStDevIfValid(
        const MsCalibratomatic &msCalibratomatic,
        QVector<float> *ionMobilityStDevs
        ) {
        if (ionMobilityStDevs == nullptr || !msCalibratomatic.isInitIM()) {
            return;
        }

        const float ionMobilityStDev = msCalibratomatic.ionMobilityStDev();
        if (ionMobilityStDev > 0.0f) {
            ionMobilityStDevs->push_back(ionMobilityStDev);
        }
    }

    void filterMs1CandidateRowsByCorr(QVector<CandidateScores*> *candidateScoresMS1Cal) {

        constexpr double cosineSimSumMS1Min = 0.9;
        const auto terminatorLogic = [cosineSimSumMS1Min](const CandidateScores *cs){
            return cs->featuresArray[CosineSim100MS1] <= cosineSimSumMS1Min;
        };

        const auto terminator = std::remove_if(
                candidateScoresMS1Cal->begin(),
                candidateScoresMS1Cal->end(),
                terminatorLogic
                );

        candidateScoresMS1Cal->erase(terminator, candidateScoresMS1Cal->end());
    }

    ScanNumber nearestScanNumberAtOrBefore(
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        ScanTime scanTime
        ) {

        if (scanNumberVsScanTime.isEmpty()) {
            return -1;
        }

        ScanNumber nearestScanNumber = -1;
        ScanTime nearestScanTime = -1.0f;
        for (auto it = scanNumberVsScanTime.constBegin(); it != scanNumberVsScanTime.constEnd(); ++it) {
            const ScanTime currentScanTime = it.value();
            if (currentScanTime > scanTime) {
                break;
            }

            nearestScanNumber = it.key();
            nearestScanTime = currentScanTime;
        }

        if (nearestScanNumber >= 0) {
            return nearestScanNumber;
        }

        return scanNumberVsScanTime.constBegin().key();
    }

}//namespace
Err MsCalibratomaticSettertron::buildCalibration(MsCalibratomatic *msCalibratomatic) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters->isValid()); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron->isInit()); ree;

    const int numberOfTranches = calculateNumberOfTranches(m_pythiaParameters->verbosity);

    const QVector<MsScanInfo> uniqueMsScanInfos = m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    constexpr double calibrationScanInfoFraction = 0.10;
    constexpr int minUniqueScanInfosTrainingCount = 4;
    constexpr int maxUniqueScanInfosTrainingCount = 16;
    const int uniqueScanInfosTrainingCount = std::clamp(
        static_cast<int>(std::ceil(uniqueMsScanInfos.size() * calibrationScanInfoFraction)),
        minUniqueScanInfosTrainingCount,
        maxUniqueScanInfosTrainingCount
        );
    constexpr int offset = 0;
    QVector<MsScanInfo> uniqueMsScanInfosCalibration;
    e = PythiaDIAFFWorkflowSharedMethods::buildUniqueMsScanInfosForProcessing(
        uniqueMsScanInfos,
        uniqueScanInfosTrainingCount,
        offset,
        &uniqueMsScanInfosCalibration
        ); ree;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFramesPntrs
                                        = m_targetDecoyCandidatePairScoretron->diaTargetFrames();
    if (m_msReaderPointerAcc->useLazyLoading()) {
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
                &m_msCalibratomatic,
                &mzTargetKeyVsTargetDecoyCandidatePointers
                ); ree;

        constexpr int topNMS2IonsCalibration = 6;

        constexpr bool useTopNIntegrationsParameter = false;

        m_weights = m_weights.isEmpty()
                  ? DiscriminantScoretron::defaultWeights(m_featuresCalibration)
                  : m_weights;

        constexpr int splitter = 2;
        const int threadCount = uniqueMsScanInfos.size() < m_pythiaParameters->threadCount
                      ? std::min(uniqueMsScanInfos.size() * splitter, m_pythiaParameters->threadCount)
                      : m_pythiaParameters->threadCount;

        constexpr float minPeakCountCalibration = 4.9;
        m_candidateScorePairs.clear();
        e = m_targetDecoyCandidatePairScoretron->scoreTargetDecoyPairs(
                m_featuresCalibration,
                topNMS2IonsCalibration,
                m_msCalibratomatic,
                minPeakCountCalibration,
                threadCount,
                useTopNIntegrationsParameter,
                mzTargetKeyVsTurboXicPntrs,
                m_weights,
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                &m_candidateScorePairs
                ); ree

        QVector<CandidateScores*> candidateScoresVecBatchPntrs;
        QMap<int, int> fdrVsCounts;
        QVector<float> weights;
        e = PythiaDIAFFWorkflowSharedMethods::processBatch(
            m_featuresCalibration,
            m_candidateScorePairs,
            *m_pythiaParameters,
            &candidateScoresVecBatchPntrs,
            &fdrVsCounts,
            &weights,
            false
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
        appendIonMobilityStDevIfValid(m_msCalibratomatic, &m_ionMobilityStDevs);
        m_ms2PPMStDevs.push_back(m_msCalibratomatic.mzStDevMS2());

        QString fdrString;
        e = FDRCLassifierNeuralNet::outPutFDRCounts(fdrVsCounts, &fdrString); ree;
        m_fdrWeightedMean = PythiaDIAFFWorkflowSharedMethods::weightedFDRMean(fdrVsCounts);

        ++m_batchCounter;
        if (m_pythiaParameters->verbosity >= 0) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                    << "Processed batch" << m_batchCounter
                    << "of"
                    << targetDecoyCandidatePointersTranched.size()
                    << etBatch.elapsed()
                    << "mSec"
                    << qPrintable(fdrString)
                    << "FDR Weighted Mean"
                    << m_fdrWeightedMean;
        }

        if (m_pythiaParameters->optimizeMode) {
            continue;
        }

        if (fdrVsCounts.value(fdrKey) > m_pythiaParameters->calibrationTrainingVolume
            || &tdcp == targetDecoyCandidatePointersTranched.cend() - 1
            ) {

            constexpr int fdrKey50PercentFDR = 50;
            candidateScoresVecBatchPntrs.resize(std::min(candidateScoresVecBatchPntrs.size(), fdrVsCounts.value(fdrKey50PercentFDR)));

            QVector<CandidateScores*> candidateScoresVecBatchPntrsRecal = candidateScoresVecBatchPntrs;
            candidateScoresVecBatchPntrsRecal.resize(std::min(candidateScoresVecBatchPntrsRecal.size(), fdrVsCounts.value(fdrKeyMassCalMS2)));
            if (candidateScoresVecBatchPntrsRecal.isEmpty()) {
                candidateScoresVecBatchPntrsRecal = candidateScoresVecBatchPntrs;
                candidateScoresVecBatchPntrsRecal.resize(std::min(candidateScoresVecBatchPntrsRecal.size(), fdrVsCounts.value(fdrKey)));
            }

            if (candidateScoresVecBatchPntrsRecal.isEmpty()) {
                for (TurboXIC* turboXic : mzTargetKeyVsTurboXicPntrs) {delete turboXic;}
                rrr(eValueError);
            }

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
            appendIonMobilityStDevIfValid(m_msCalibratomatic, &m_ionMobilityStDevs);
            m_ms2PPMStDevs.push_back(m_msCalibratomatic.mzStDevMS2());

            e = setMsCalibratomaticMetrics(); ree;

            if (!m_msReaderPointerAcc->useLazyLoading()) {
                e = recalibrateMzVals(
                        MSLevelEnum::MS2,
                        m_targetDecoyCandidatePairScoretron->diaTargetFrames(),
                        m_targetDecoyCandidatePairScoretron->ms1ScanNumberVsScanPoints()
                        ); ree;
            }

            candidateScoresVecBatchPntrsRecal = candidateScoresVecBatchPntrs;
            candidateScoresVecBatchPntrsRecal.resize(std::min(candidateScoresVecBatchPntrsRecal.size(), fdrVsCounts.value(fdrKeyMassCalMS1)));
            filterMs1CandidateRowsByCorr(&candidateScoresVecBatchPntrsRecal);
            constexpr int recalibrationPointCountMinMS1 = 400;

            if (m_pythiaParameters->verbosity > 0) {
                qDebug()<< qPrintable(S_GLOBAL_TIMER.elapsed())
                    << candidateScoresVecBatchPntrsRecal.size()
                    << "found for MS1 Recalibration";
            }

            if (candidateScoresVecBatchPntrsRecal.size() < recalibrationPointCountMinMS1) {

                if (m_pythiaParameters->verbosity > 0) {
                    qWarning() << qPrintable(S_GLOBAL_TIMER.elapsed())
                        << "Skipping MS1 recalibration.  Not enough points found";
                }

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

            if (msCalibrationReaderRowsMS1.size() < recalibrationPointCountMinMS1) {
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

            // e = recalibrateMzVals(
            //         MSLevelEnum::MS1,
            //         m_targetDecoyCandidatePairScoretron->diaTargetFrames(),
            //         m_targetDecoyCandidatePairScoretron->ms1ScanNumberVsScanPoints()
            //         ); ree;

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

int MsCalibratomaticSettertron::batchCounter() const {
    return m_batchCounter;
}

double MsCalibratomaticSettertron::fdrWeightedMean() {
    return m_fdrWeightedMean;
}

int MsCalibratomaticSettertron::calculateNumberOfTranches(int verbosity) const {

    const auto sizePerTranche = static_cast<double>(m_pythiaParameters->trancheSizeMax);
    const int numberOfTranches = std::max(static_cast<int>(m_targetDecoyPairPntrs.size() / sizePerTranche), 1);

    if (verbosity > 0) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "Number of tranches for calibration:" << numberOfTranches
                 << "target count:" << m_targetDecoyPairPntrs.size()
                 << "sizePerTranche:" << static_cast<int>(sizePerTranche);
    }

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

    if (candidateScoresVecBatchPntrsResized.isEmpty()) {
        ERR_RETURN
    }

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
            MSLevelEnum::MS2,
            candidateScoresVecBatchPntrsResized,
            m_pythiaParameters->verbosity,
            &msCalibrationReaderRows
            ); ree;

    for (const CandidateScores *cs : candidateScoresVecBatchPntrsResized) {

        if (m_entered.value(cs->targetDecoyCandidatePair) ) {
            continue;
        }

        if (m_excludeDecoys) {
            if (cs->isDecoy) {
                continue;
            }
        }

        m_targetDecoyCandidatePairsTopScores.push_back(cs->targetDecoyCandidatePair);
        m_entered.insert(cs->targetDecoyCandidatePair, true);
    }

    e = m_msCalibratomatic.buildRTMapper(msCalibrationReaderRows); ree;
    if (m_msReaderPointerAcc->ptr->hasIonMobility()) {
        e = buildIonMobilityCalibrationFromCentroidApex(candidateScoresVecBatchPntrsResized); ree;
    }

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

Err MsCalibratomaticSettertron::buildIonMobilityCalibrationFromCentroidApex(
    const QVector<CandidateScores*> &candidateScores
    ) {

    ERR_INIT

    QVector<CandidateScores*> annotatedCandidateScores = candidateScores;
    e = annotateCentroidIonMobilityForCalibrationCandidates(&annotatedCandidateScores); ree;

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
        MSLevelEnum::MS2,
        annotatedCandidateScores,
        m_pythiaParameters->verbosity,
        &msCalibrationReaderRows
        ); ree;

    e = m_msCalibratomatic.buildIMMapper(msCalibrationReaderRows); ree;

    ERR_RETURN
}

Err MsCalibratomaticSettertron::annotateCentroidIonMobilityForCalibrationCandidates(
    QVector<CandidateScores*> *candidateScores
    ) const {

    ERR_INIT

    if (candidateScores == nullptr || candidateScores->isEmpty()) {
        ERR_RETURN
    }

    if (m_msReaderPointerAcc == nullptr
        || m_msReaderPointerAcc->ptr.isNull()
        || !m_msReaderPointerAcc->ptr->hasIonMobility()) {
        ERR_RETURN
    }

    QMap<ScanNumber, ScanPoints> *ms1ScanNumberVsScanPoints = m_targetDecoyCandidatePairScoretron->ms1ScanNumberVsScanPoints();
    e = ErrorUtils::isTrue(ms1ScanNumberVsScanPoints != nullptr); ree;
    e = ErrorUtils::isFalse(ms1ScanNumberVsScanPoints->isEmpty()); ree;

    const QMap<ScanNumber, ScanTime> scanNumberVsScanTime = m_msReaderPointerAcc->ptr->getScanNumberVsScanTime();
    e = ErrorUtils::isFalse(scanNumberVsScanTime.isEmpty()); ree;

    for (CandidateScores *candidateScoresPntr : *candidateScores) {
        if (candidateScoresPntr == nullptr || candidateScoresPntr->targetDecoyCandidatePair == nullptr) {
            continue;
        }

        const float libraryIonMobility = candidateScoresPntr->targetDecoyCandidatePair->iIM();
        if (libraryIonMobility <= 0.0f) {
            continue;
        }

        const ScanNumber scanNumber = nearestScanNumberAtOrBefore(
            scanNumberVsScanTime,
            candidateScoresPntr->scanTime
            );
        if (scanNumber < 0) {
            continue;
        }

        const auto scanPointsIt = ms1ScanNumberVsScanPoints->constFind(scanNumber);
        if (scanPointsIt == ms1ScanNumberVsScanPoints->constEnd()) {
            continue;
        }

        const TimsbukAlignedPointData *alignedPointData = m_msReaderPointerAcc->ptr->alignedPointDataPntr(scanNumber);
        if (alignedPointData == nullptr || !alignedPointData->isAlignedWith(scanPointsIt.value())) {
            continue;
        }

        const float monoIsotopeMz = candidateScoresPntr->targetDecoyCandidatePair->mz(false);
        const float massTol = MathUtils::calculatePPM(
            monoIsotopeMz,
            static_cast<float>(m_pythiaParameters->ms1ExtractionWidthPPM)
            );
        const float mzMin = monoIsotopeMz - massTol;
        const float mzMax = monoIsotopeMz + massTol;

        float bestIntensity = -1.0f;
        float bestDriftTime = -1.0f;
        IonMobilityIndex bestIonMobilityIndex = -1;

        const ScanPoints &scanPoints = scanPointsIt.value();
        for (int pointIndex = 0; pointIndex < scanPoints.size(); ++pointIndex) {
            const ScanPoint &scanPoint = scanPoints.at(pointIndex);
            const float mz = timsbukMzOf(scanPoint);
            if (mz < mzMin || mz > mzMax) {
                continue;
            }

            const float driftTime = timsbukIonMobilityOf(*alignedPointData, pointIndex);
            if (driftTime <= 0.0f) {
                continue;
            }

            const float intensity = timsbukIntensityOf(scanPoint);
            if (intensity <= bestIntensity) {
                continue;
            }

            bestIntensity = intensity;
            bestDriftTime = driftTime;
            bestIonMobilityIndex = static_cast<IonMobilityIndex>(std::lround(driftTime * 10000.0f));
        }

        if (bestIonMobilityIndex < 0) {
            continue;
        }

        candidateScoresPntr->imDriftTime = bestDriftTime;
        candidateScoresPntr->ionMobilityIndex = bestIonMobilityIndex;
        candidateScoresPntr->ionMobilityIndexStart = bestIonMobilityIndex;
        candidateScoresPntr->ionMobilityIndexEnd = bestIonMobilityIndex;
        candidateScoresPntr->featuresArray[Ms1IntensityFoundApex100IM] = bestIntensity;
        candidateScoresPntr->featuresArray[IonMobilityDelta] = bestDriftTime - libraryIonMobility;
        candidateScoresPntr->featuresArray[IonMobilityDeltaAbs] = std::abs(bestDriftTime - libraryIonMobility);
        candidateScoresPntr->featuresArray[IonMobilityPdAbs]
            = std::sqrt(
                std::min(
                    static_cast<double>(std::abs(bestDriftTime - libraryIonMobility)),
                    DEFAULT_ION_MOBILITY_TOLERANCE_ONE_OVER_K0
                    )
                / DEFAULT_ION_MOBILITY_TOLERANCE_ONE_OVER_K0
                );
    }

    ERR_RETURN
}

Err MsCalibratomaticSettertron::setMsCalibratomaticMetrics() {

    ERR_INIT

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

    if (m_pythiaParameters->verbosity >= 0) {
        qDebug()
        << qPrintable(S_GLOBAL_TIMER.elapsed())
        << "ScanTimeWindow Mean|Median|Min"
        << MathUtils::mean(m_scanTimeStDevs)
        << MathUtils::median(m_scanTimeStDevs)
        << *std::min({m_scanTimeStDevs.begin(), m_scanTimeStDevs.end()});

        qDebug()
        << qPrintable(S_GLOBAL_TIMER.elapsed())
        << "Ms2 ppm Mean|Median"
        << MathUtils::mean(m_ms2PPMStDevs)
        << MathUtils::median(m_ms2PPMStDevs);
    }

    if (m_msReaderPointerAcc->ptr->hasIonMobility()) {
        QVector<float> validIonMobilityStDevs = m_ionMobilityStDevs;
        std::sort(validIonMobilityStDevs.begin(), validIonMobilityStDevs.end());
        if (validIonMobilityStDevs.size() >= minVecSize) {
            validIonMobilityStDevs.pop_front();
            validIonMobilityStDevs.pop_back();
        }

        if (!validIonMobilityStDevs.isEmpty()) {
            if (m_pythiaParameters->verbosity >= 0) {
                qDebug()
                << qPrintable(S_GLOBAL_TIMER.elapsed())
                << "IonMobilityWindow Mean|Median|Min"
                << MathUtils::mean(validIonMobilityStDevs)
                << MathUtils::median(validIonMobilityStDevs)
                << *std::min({validIonMobilityStDevs.begin(), validIonMobilityStDevs.end()});
            }

            m_msCalibratomatic.setIonMobilityStDev(validIonMobilityStDevs.front());
        } else if (m_pythiaParameters->verbosity >= 0) {
            qDebug()
            << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "IonMobilityWindow calibration unavailable; using fallback width"
            << DEFAULT_ION_MOBILITY_TOLERANCE_ONE_OVER_K0
            << "1/K0";
        }
    }

    m_msCalibratomatic.setScanTimeStDev(m_scanTimeStDevs.front());
    m_msCalibratomatic.setMzStDevMS2(MathUtils::mean(m_ms2PPMStDevs));

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

	// e = ErrorUtils::isFalse(scanNumberVsScanTimeMS1->isEmpty()); ree;

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

    if (m_pythiaParameters->verbosity > 0) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Points recalibrated in" << et.elapsed() << "mSec";
    }

    ERR_RETURN
}
