//
// Created by andrewnichols on 9/28/24.
//

#include "OptimizeMassAccuracyPPMSettertron.h"

#include <boost/math/distributions/normal.hpp>

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

    struct DOEResult {
        double ppm = -1.0;
        double fdrCount = -1;
    };

    QString ppmCacheKey(double ppm) {
        return QString::number(ppm, 'f', 4);
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

    Err getTopFrequencyParameters(
            const QVector<DOEResult> &results,
            int verbosity,
            double *ppmSetting
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(results); ree;

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

        constexpr int polynomialOrder = 2;

        QVector<double> coeffs;
        EigenUtils::fitPolynomialQRDecomposition(xyMat, polynomialOrder, &coeffs);

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PPM Coeffs" << coeffs;

        QVector<double> xPoints = {results.first().ppm};
        while (xPoints.back() < results.back().ppm) {
            constexpr double incrementVal = 0.25;
            xPoints.push_back(xPoints.back() + incrementVal);
        }

        QVector<double> yPoints;
        for (double x : xPoints) {
            double y = 0.0;
            for (int i = 0; i < coeffs.size(); i++) {
                y += coeffs.at(i) * std::pow(x, i);
            }
            if (verbosity > 0) {
                qDebug() << x << y;

            }
            yPoints.push_back(y);
        }

        *ppmSetting = xPoints.at(std::max_element(yPoints.begin(), yPoints.end()) - yPoints.begin());

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
        QVector<float> *bestWeightsOut,
        int *bestIdsAtFivePercent,
        QMap<int, int> *bestFdrVsCountsOut
        ) -> Err {

        ERR_INIT

        resultsOut->clear();
        bestWeightsOut->clear();
        *bestIdsAtFivePercent = 0;
        bestFdrVsCountsOut->clear();

        double bestResultCount = 0;
        int countSinceLastHigh = 0;
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
                m_msReaderPointerAcc->ptr->isTIMS()
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
                countSinceLastHigh = 0;
                *bestWeightsOut = weights;
                *bestIdsAtFivePercent = fdrVsCounts.value(OPTIMIZATION_SUPPORT_FDR_KEY);
                *bestFdrVsCountsOut = fdrVsCounts;
                continue;
            }

            constexpr int maxCountsSinceLastHigh = 3;
            if (++countSinceLastHigh >= maxCountsSinceLastHigh) {
                break;
            }
        }

        ERR_RETURN
    };

    QMap<QString, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>> firstTrancheScorePairsByPpm;
    QVector<DOEResult> results;
    QVector<float> bestWeights;
    int bestIdsAtFivePercent = 0;
    QMap<int, int> bestFdrVsCounts;
    e = runPpmSweep(
        1,
        &firstTrancheScorePairsByPpm,
        &results,
        &bestWeights,
        &bestIdsAtFivePercent,
        &bestFdrVsCounts
        ); ree;

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
            &bestWeights,
            &bestIdsAtFivePercent,
            &bestFdrVsCounts
            ); ree;
    }

    m_weights = bestWeights;

    e = getTopFrequencyParameters(
            results,
            m_pythiaParameters->verbosity,
            &m_pythiaParameters->ms2ExtractionWidthPPM
            ); ree;

    m_pythiaParameters->ms1ExtractionWidthPPM = m_pythiaParameters->ms2ExtractionWidthPPM;
    e = m_targetDecoyCandidatePairScoretron->setPythiaParameters(*m_pythiaParameters); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Optimal ppm setting:"
             << m_pythiaParameters->ms2ExtractionWidthPPM;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Optimal scanTimeWindow setting:"
             << m_msCalibratomatic->scanTimeStDev(m_pythiaParameters->scanTimeWindowStDevs)
             << "minutes";

    if (m_msReaderPointerAcc->ptr->isTIMS()) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Optimal ionMobilityWindow setting:"
             << m_msCalibratomatic->ionMobilityStDev(m_pythiaParameters->scanTimeWindowStDevs)
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
