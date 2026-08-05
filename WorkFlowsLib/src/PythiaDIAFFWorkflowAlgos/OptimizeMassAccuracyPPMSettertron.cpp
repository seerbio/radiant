//
// Created by andrewnichols on 9/28/24.
//

#include "OptimizeMassAccuracyPPMSettertron.h"

#include <boost/math/distributions/normal.hpp>
#include <cmath>

#include "EigenUtils.h"
#include "FDRCLassifierNeuralNet.h"
#include "MsCalibratomaticSettertron.h"
#include "ParallelUtils.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>

#include "DiscriminantScoretron.h"
#include "PythiaDIAFFWorkflowSharedMethods.h"

OptimizeMassAccuracyPPMSettertron::OptimizeMassAccuracyPPMSettertron()
: m_msReaderPointerAcc(nullptr)
, m_msCalibratomatic(nullptr)
, m_pythiaParameters(nullptr)
, m_targetDecoyCandidatePairScoretron(nullptr)
, m_targetDecoyPairPntrs(nullptr)
{}

Err OptimizeMassAccuracyPPMSettertron::initExec(
    const QVector<Features> &optimzeFeatures,
    MsReaderPointerAcc *msReaderPointerAcc,
    MsCalibratomatic *msCalibratomatic,
    PythiaParameters *pythiaParameters,
    TargetDecoyCandidatePairScoretron2 *targetDecoyCandidatePairScoretron,
    QVector<TargetDecoyCandidatePair*> *targetDecoyPairPntrs
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(msReaderPointerAcc->isInit()); ree;
    e = ErrorUtils::isTrue(pythiaParameters->isValid()); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairScoretron->isInit()); ree;
    e = ErrorUtils::isFalse(targetDecoyPairPntrs->isEmpty()); ree;
    e = ErrorUtils::isNotEmpty(optimzeFeatures); ree;
    e = ErrorUtils::isTrue(msCalibratomatic->isInitRT()); eee_absorb;

    m_msReaderPointerAcc = msReaderPointerAcc;
    m_msCalibratomatic = msCalibratomatic;
    m_pythiaParameters = pythiaParameters;
    m_targetDecoyCandidatePairScoretron = targetDecoyCandidatePairScoretron;
    m_targetDecoyPairPntrs = targetDecoyPairPntrs;
    m_optimizeFeatures = optimzeFeatures;

    optimizePPM();

    ERR_RETURN
}

QVector<float> OptimizeMassAccuracyPPMSettertron::weights() const {
    return m_weights;
}

namespace {

    constexpr int OPTIMIZATION_SUPPORT_FDR_KEY = 5;
    constexpr int OPTIMIZATION_SUPPORT_FDR_COUNT_MIN = 100;
    constexpr int OPTIMIZATION_SUPPORT_TRANCHE_CAP = 5;
    constexpr int PPM_FIT_POLYNOMIAL_ORDER = 2;
    constexpr double PPM_FIT_SAMPLE_INCREMENT = 0.25;
    constexpr int PPM_FIT_EARLY_STOP_PAST_APEX_POINTS = 3;
    constexpr int PPM_FIT_MIN_RESULTS_FOR_EARLY_STOP = 5;

    struct DOEResult {
        double ppm = -1.0;
        double fdrCount = -1;
    };

    struct PpmFitSummary {
        QVector<double> coeffs;
        QVector<double> xPoints;
        QVector<double> yPoints;
        double apexPpm = -1.0;
        double sampledBestPpm = -1.0;
        double evaluatedPpmMin = -1.0;
        double evaluatedPpmMax = -1.0;
        bool hasUsableApex = false;
    };

    QString ppmCacheKey(double ppm) {
        return QString::number(ppm, 'f', 4);
    }

    int nearestEvaluatedPpmIndex(const QVector<DOEResult> &results, double ppm) {
        if (results.isEmpty()) {
            return -1;
        }

        int bestIndex = 0;
        double bestDistance = std::abs(results.first().ppm - ppm);
        for (int i = 1; i < results.size(); ++i) {
            const double distance = std::abs(results.at(i).ppm - ppm);
            if (distance < bestDistance) {
                bestDistance = distance;
                bestIndex = i;
            }
        }

        return bestIndex;
    }

    int fallbackOptimizationSupportFdrKey(const QMap<int, int> &fdrVsCounts) {
        const QVector<int> fallbackKeys = {10, 20, 50};
        for (int fdrKey : fallbackKeys) {
            if (fdrVsCounts.value(fdrKey) > 0) {
                return fdrKey;
            }
        }

        return -1;
    }

    Err buildDOE(
            const PythiaParameters &pythiaParameters,
            double scanTimeStDev,
            int verbosity,
            QVector<PythiaParameters> *pythiaParametersExperiments
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(scanTimeStDev > 0.0); ree;
        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

        const QVector<double> ppmList = {
            // 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 26, 30, 35, 40, 50
            3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 26, 30, 35, 40, 50
        };

        for (double ppm : ppmList ) {
            PythiaParameters params = pythiaParameters;
            params.ms2ExtractionWidthPPM = ppm;
            pythiaParametersExperiments->push_back(params);
        }

        if (verbosity > 0) {
            qDebug() << "Testing PPM Values";
            for (const PythiaParameters &pp : *pythiaParametersExperiments) {
                qDebug() << "ppmTol" << pp.ms2ExtractionWidthPPM;
            }
        }


        ERR_RETURN
    }

    Err buildPpmFitSummary(
            const QVector<DOEResult> &results,
            int verbosity,
            PpmFitSummary *fitSummary
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(results); ree;
        e = ErrorUtils::isTrue(fitSummary != nullptr, eValueError); ree;

        fitSummary->coeffs.clear();
        fitSummary->xPoints.clear();
        fitSummary->yPoints.clear();
        fitSummary->apexPpm = -1.0;
        fitSummary->sampledBestPpm = -1.0;
        fitSummary->evaluatedPpmMin = results.first().ppm;
        fitSummary->evaluatedPpmMax = results.back().ppm;
        fitSummary->hasUsableApex = false;

        Eigen::MatrixX<double> xyMat(results.size() + 1, 2);
        xyMat.setZero();
        for (int row = 0; row < results.size(); row++) {
            const DOEResult &doeResult = results.at(row);
            xyMat.coeffRef(row + 1, 0) = doeResult.ppm;
            xyMat.coeffRef(row + 1, 1) = static_cast<double>(doeResult.fdrCount);
            if (verbosity > 0) {
                qDebug() << doeResult.ppm << doeResult.fdrCount << xyMat.coeff(row, 1);

            }
        }

        EigenUtils::fitPolynomialQRDecomposition(xyMat, PPM_FIT_POLYNOMIAL_ORDER, &fitSummary->coeffs);

        fitSummary->xPoints = {results.first().ppm};
        while (fitSummary->xPoints.back() < results.back().ppm) {
            fitSummary->xPoints.push_back(fitSummary->xPoints.back() + PPM_FIT_SAMPLE_INCREMENT);
        }

        for (double x : fitSummary->xPoints) {
            double y = 0.0;
            for (int i = 0; i < fitSummary->coeffs.size(); i++) {
                y += fitSummary->coeffs.at(i) * std::pow(x, i);
            }
            if (verbosity > 0) {
                qDebug() << x << y;

            }
            fitSummary->yPoints.push_back(y);
        }

        const auto maxIt = std::max_element(fitSummary->yPoints.begin(), fitSummary->yPoints.end());
        fitSummary->sampledBestPpm
            = fitSummary->xPoints.at(maxIt - fitSummary->yPoints.begin());

        if (fitSummary->coeffs.size() > 2 && !MathUtils::tZero(fitSummary->coeffs.at(2))) {
            fitSummary->apexPpm = -fitSummary->coeffs.at(1) / (2.0 * fitSummary->coeffs.at(2));
            fitSummary->hasUsableApex = std::isfinite(fitSummary->apexPpm);
        }

        ERR_RETURN
    }

    Err getTopFrequencyParameters(
            const PpmFitSummary &fitSummary,
            double *ppmSetting
            ) {

        ERR_INIT
        e = ErrorUtils::isTrue(ppmSetting != nullptr, eValueError); ree;
        e = ErrorUtils::isFalse(fitSummary.coeffs.isEmpty(), eValueError); ree;

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PPM Coeffs" << fitSummary.coeffs;
        *ppmSetting = fitSummary.sampledBestPpm;

        ERR_RETURN
    }

    Err getTopFrequencyParameters(
            const QVector<DOEResult> &results,
            int verbosity,
            double *ppmSetting
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(results); ree;
        e = ErrorUtils::isTrue(ppmSetting != nullptr, eValueError); ree;

        PpmFitSummary fitSummary;
        e = buildPpmFitSummary(results, verbosity, &fitSummary); ree;
        e = getTopFrequencyParameters(fitSummary, ppmSetting); ree;

        ERR_RETURN
    }

}//namespace
Err OptimizeMassAccuracyPPMSettertron::optimizePPM() {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msReaderPointerAcc->isInit()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters->isValid()); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron->isInit()); ree;
    e = ErrorUtils::isTrue(m_msCalibratomatic->isInitRT()); ree;

    constexpr int topNMS2IonsOptimization = 12;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Using top:" << topNMS2IonsOptimization << "fragments for optimization";

    const int numberOfTranches = calculateNumberOfTranches();

    const QVector<MsScanInfo> uniqueMsScanInfos = m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    constexpr int maxUniqueScanInfosTrainingCount = 16;
    constexpr int offset = 1;
    QVector<MsScanInfo> uniqueMsScanInfosOptimization;
    e = PythiaDIAFFWorkflowSharedMethods::buildUniqueMsScanInfosForProcessing(
        uniqueMsScanInfos,
        maxUniqueScanInfosTrainingCount,
        offset,
        &uniqueMsScanInfosOptimization
        ); ree;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFramesPntrs
                                    = m_targetDecoyCandidatePairScoretron->diaTargetFrames();
    if (m_msReaderPointerAcc->useLazyLoading()) {
        e = PythiaDIAFFWorkflowSharedMethods::buildDiaTargetFrames(
            uniqueMsScanInfosOptimization,
            m_msReaderPointerAcc,
            m_msCalibratomatic,
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
        uniqueMsScanInfosOptimization,
        m_msReaderPointerAcc->ptr->getScanNumberVsScanTime(),
        diaTargetFramesPntrs,
        &mzTargetKeyVsTurboXicPntrs
        ); ree;

    QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePointersTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            *m_targetDecoyPairPntrs,
            numberOfTranches,
            &targetDecoyCandidatePointersTranched
            ); ree;

    QVector<PythiaParameters> pythiaParametersExperiments;
    e = buildDOE(
            *m_pythiaParameters,
            m_msCalibratomatic->scanTimeStDev(),
            m_pythiaParameters->verbosity,
            &pythiaParametersExperiments
            ); ree;

    constexpr int splitter = 2;
    const int threadCount = uniqueMsScanInfos.size() < m_pythiaParameters->threadCount
                  ? std::min(uniqueMsScanInfos.size() * splitter, m_pythiaParameters->threadCount)
                  : m_pythiaParameters->threadCount;

    constexpr bool useTopNIntegrationsParameter = true;
    constexpr float minPeakCountOptimization = 3.9;
    constexpr int topNMS2Ions = 12;
    const QVector<float> defaultOptimizationWeights = DiscriminantScoretron::defaultWeights(m_optimizeFeatures);

    auto scoreOptimizationTranche = [&](
        const PythiaParameters &pythiaParams,
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePointers,
        QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> *candidateScorePairs
        ) -> Err {

        ERR_INIT

        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
        e = PythiaDIAFFWorkflowSharedMethods::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            targetDecoyCandidatePointers,
            pythiaParams,
            uniqueMsScanInfosOptimization,
            nullptr,
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

        candidateScorePairs->clear();
        e = m_targetDecoyCandidatePairScoretron->scoreTargetDecoyPairs(
                m_optimizeFeatures,
                topNMS2Ions,
                *m_msCalibratomatic,
                minPeakCountOptimization,
                threadCount,
                useTopNIntegrationsParameter,
                mzTargetKeyVsTurboXicPntrs,
                defaultOptimizationWeights,
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                candidateScorePairs
                ); ree

        ERR_RETURN
    };

    auto buildCandidateScorePairsForSweep = [&](
        const PythiaParameters &pythiaParams,
        int trancheCountToUse,
        QMap<QString, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>> *firstTrancheScorePairsByPpm,
        QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> *candidateScorePairsForBatch
        ) -> Err {

        ERR_INIT

        const int trancheCountBounded = std::max(1, std::min(trancheCountToUse, targetDecoyCandidatePointersTranched.size()));
        const QString cacheKey = ppmCacheKey(pythiaParams.ms2ExtractionWidthPPM);

        candidateScorePairsForBatch->clear();
        if (firstTrancheScorePairsByPpm->contains(cacheKey)) {
            candidateScorePairsForBatch->append(firstTrancheScorePairsByPpm->value(cacheKey));
        }
        else {
            QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> firstTrancheCandidateScorePairs;
            e = scoreOptimizationTranche(
                pythiaParams,
                targetDecoyCandidatePointersTranched.first(),
                &firstTrancheCandidateScorePairs
                ); ree;

            firstTrancheScorePairsByPpm->insert(cacheKey, firstTrancheCandidateScorePairs);
            candidateScorePairsForBatch->append(firstTrancheCandidateScorePairs);
        }

        for (int trancheIndex = 1; trancheIndex < trancheCountBounded; ++trancheIndex) {
            QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> trancheCandidateScorePairs;
            e = scoreOptimizationTranche(
                pythiaParams,
                targetDecoyCandidatePointersTranched.at(trancheIndex),
                &trancheCandidateScorePairs
                ); ree;
            candidateScorePairsForBatch->append(trancheCandidateScorePairs);
        }

        ERR_RETURN
    };

    auto runPpmSweep = [&](
        int trancheCountToUse,
        QMap<QString, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>> *firstTrancheScorePairsByPpm,
        QVector<DOEResult> *resultsOut,
        PpmFitSummary *lastFitSummaryOut,
        QVector<float> *bestWeightsOut,
        int *bestIdsAtFivePercent,
        QMap<int, int> *bestFdrVsCountsOut
        ) -> Err {

        ERR_INIT

        resultsOut->clear();
        if (lastFitSummaryOut != nullptr) {
            *lastFitSummaryOut = PpmFitSummary();
        }
        bestWeightsOut->clear();
        *bestIdsAtFivePercent = 0;
        bestFdrVsCountsOut->clear();

        double bestResultCount = -1.0;
        const int trancheCountBounded = std::max(1, std::min(trancheCountToUse, targetDecoyCandidatePointersTranched.size()));

        for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {

            e = m_targetDecoyCandidatePairScoretron->setPythiaParameters(pythiaParams); ree;

            QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> candidateScorePairsForBatch;
            e = buildCandidateScorePairsForSweep(
                pythiaParams,
                trancheCountBounded,
                firstTrancheScorePairsByPpm,
                &candidateScorePairsForBatch
                ); ree;

            QVector<CandidateScores*> candidateScoresVecBatchPntrs;
            QMap<int, int> fdrVsCounts;
            QVector<float> weights;
            e = PythiaDIAFFWorkflowSharedMethods::processBatch(
                m_optimizeFeatures,
                candidateScorePairsForBatch,
                pythiaParams,
                &candidateScoresVecBatchPntrs,
                &fdrVsCounts,
                &weights,
                false
                ); ree;

            QString fdrString;
            e = FDRCLassifierNeuralNet::outPutFDRCounts(fdrVsCounts, &fdrString); ree;

            const double fdrMean = PythiaDIAFFWorkflowSharedMethods::weightedFDRMean(fdrVsCounts);

            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "ppmTol"
                     << pythiaParams.ms2ExtractionWidthPPM
                     << "tranches"
                     << trancheCountBounded
                     << "fdrMean"
                     << fdrMean
                     << "Finished"
                     << qPrintable(fdrString);

            DOEResult res;
            res.ppm = pythiaParams.ms2ExtractionWidthPPM;
            res.fdrCount = fdrMean;
            resultsOut->push_back(res);

            if (res.fdrCount >= bestResultCount) {
                bestResultCount = res.fdrCount;
                *bestWeightsOut = weights;
                *bestIdsAtFivePercent = fdrVsCounts.value(OPTIMIZATION_SUPPORT_FDR_KEY);
                *bestFdrVsCountsOut = fdrVsCounts;
            }

            if (resultsOut->size() < PPM_FIT_MIN_RESULTS_FOR_EARLY_STOP) {
                continue;
            }

            PpmFitSummary fitSummary;
            e = buildPpmFitSummary(*resultsOut, 0, &fitSummary); ree;
            if (lastFitSummaryOut != nullptr) {
                *lastFitSummaryOut = fitSummary;
            }

            if (!fitSummary.hasUsableApex
                || fitSummary.coeffs.size() <= 2
                || fitSummary.coeffs.at(2) >= 0.0
                || fitSummary.apexPpm < fitSummary.evaluatedPpmMin
                || fitSummary.apexPpm > fitSummary.evaluatedPpmMax) {
                continue;
            }

            const int apexEvaluatedIndex = nearestEvaluatedPpmIndex(*resultsOut, fitSummary.apexPpm);
            if (apexEvaluatedIndex < 0) {
                continue;
            }

            const int pointsPastApex = resultsOut->size() - 1 - apexEvaluatedIndex;
            if (pointsPastApex >= PPM_FIT_EARLY_STOP_PAST_APEX_POINTS) {
                break;
            }
        }

        ERR_RETURN
    };

    QMap<QString, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>> firstTrancheScorePairsByPpm;
    QVector<DOEResult> results;
    PpmFitSummary finalFitSummary;
    QVector<float> bestWeights;
    int bestIdsAtFivePercent = 0;
    QMap<int, int> bestFdrVsCounts;
    e = runPpmSweep(
        1,
        &firstTrancheScorePairsByPpm,
        &results,
        &finalFitSummary,
        &bestWeights,
        &bestIdsAtFivePercent,
        &bestFdrVsCounts
        ); ree;

    // #define ENABLE_SECOND_STAGE_PPM_OPTIMIZATION
    #ifdef ENABLE_SECOND_STAGE_PPM_OPTIMIZATION
    if (bestIdsAtFivePercent < OPTIMIZATION_SUPPORT_FDR_COUNT_MIN
        && targetDecoyCandidatePointersTranched.size() > 1) {

        int expandedTrancheCount = 5;
        int supportCountForExpansion = bestIdsAtFivePercent;
        int supportFdrKeyForExpansion = OPTIMIZATION_SUPPORT_FDR_KEY;

        if (supportCountForExpansion == 0) {
            supportFdrKeyForExpansion = fallbackOptimizationSupportFdrKey(bestFdrVsCounts);
            if (supportFdrKeyForExpansion > OPTIMIZATION_SUPPORT_FDR_KEY) {
                supportCountForExpansion = bestFdrVsCounts.value(supportFdrKeyForExpansion);
            }
        }

        if (supportCountForExpansion > 0) {
            expandedTrancheCount = std::min(
                OPTIMIZATION_SUPPORT_TRANCHE_CAP,
                std::min(
                    targetDecoyCandidatePointersTranched.size(),
                    std::max(
                        2,
                        static_cast<int>(std::ceil(
                            static_cast<double>(OPTIMIZATION_SUPPORT_FDR_COUNT_MIN)
                            / static_cast<double>(supportCountForExpansion)
                            ))
                        )
                    )
                );
        }
        else {
            expandedTrancheCount = std::min(OPTIMIZATION_SUPPORT_TRANCHE_CAP, targetDecoyCandidatePointersTranched.size());
            qWarning() << qPrintable(S_GLOBAL_TIMER.elapsed())
                       << "PPM optimization support remained zero across all tracked FDR levels;"
                       << "rerunning sweep with fallback tranche count"
                       << expandedTrancheCount;
        }

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "PPM optimization rerunning sweep with expanded support"
                 << "support_fdr_pct"
                 << supportFdrKeyForExpansion
                 << "support_count"
                 << supportCountForExpansion
                 << "tranches"
                 << expandedTrancheCount;

        e = runPpmSweep(
            expandedTrancheCount,
            &firstTrancheScorePairsByPpm,
            &results,
            &finalFitSummary,
            &bestWeights,
            &bestIdsAtFivePercent,
            &bestFdrVsCounts
            ); ree;
    }
    #endif

    m_weights = bestWeights;

    if (finalFitSummary.coeffs.isEmpty()) {
        e = buildPpmFitSummary(
                results,
                m_pythiaParameters->verbosity,
                &finalFitSummary
                ); ree;
    }

    e = getTopFrequencyParameters(finalFitSummary, &m_pythiaParameters->ms2ExtractionWidthPPM); ree;

    m_pythiaParameters->ms1ExtractionWidthPPM = m_pythiaParameters->ms2ExtractionWidthPPM;
    e = m_targetDecoyCandidatePairScoretron->setPythiaParameters(*m_pythiaParameters); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Optimal ppm setting:"
             << m_pythiaParameters->ms2ExtractionWidthPPM;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Optimal scanTimeWindow setting:"
             << m_msCalibratomatic->scanTimeStDev(m_pythiaParameters->scanTimeWindowStDevs)
             << "minutes";

    if (m_msReaderPointerAcc->ptr->hasIonMobility()) {
        constexpr float DEFAULT_ION_MOBILITY_TOLERANCE_ONE_OVER_K0 = 0.1f;
        const float calibratedIonMobilityWindow
            = m_msCalibratomatic->ionMobilityStDev(m_pythiaParameters->scanTimeWindowStDevs);
        const float ionMobilityWindowSetting
            = (m_msCalibratomatic->isInitIM() && calibratedIonMobilityWindow > 0.0f)
            ? calibratedIonMobilityWindow
            : DEFAULT_ION_MOBILITY_TOLERANCE_ONE_OVER_K0;
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Optimal ionMobilityWindow setting:"
             << ionMobilityWindowSetting
             << "mSec";
    }

    for (TurboXIC* turboXic : mzTargetKeyVsTurboXicPntrs) {delete turboXic;}

    ERR_RETURN
}

int OptimizeMassAccuracyPPMSettertron::calculateNumberOfTranches() const {

    constexpr int optimizationMultiplicationFactor = 5;
    const auto sizePerTranche = static_cast<double>(m_pythiaParameters->trancheSizeMax * optimizationMultiplicationFactor);
    const int numberOfTranches = std::max(static_cast<int>(m_targetDecoyPairPntrs->size() / sizePerTranche), 1);
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "target count:" << m_targetDecoyPairPntrs->size()
             << "sizePerTranche:" << static_cast<int>(sizePerTranche)
             << "%:" << MathUtils::pRound(static_cast<int>(sizePerTranche) / static_cast<float>(m_targetDecoyPairPntrs->size()) * 100, 1);

    return numberOfTranches;
}
