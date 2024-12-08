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

    e = ErrorUtils::isTrue(msCalibratomatic->isInitRT()); eee_absorb;

    m_msReaderPointerAcc = msReaderPointerAcc;
    m_msCalibratomatic = msCalibratomatic;
    m_pythiaParameters = pythiaParameters;
    m_targetDecoyCandidatePairScoretron = targetDecoyCandidatePairScoretron;
    m_targetDecoyPairPntrs = targetDecoyPairPntrs;

    optimizePPM();

    ERR_RETURN
}

QVector<float> OptimizeMassAccuracyPPMSettertron::weights() const {
    return m_weights;
}

namespace {

    struct DOEResult {
        double ppm = -1.0;
        double fdrCount = -1;
    };

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
    if (m_pythiaParameters->useLazyLoading) {
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

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = PythiaDIAFFWorkflowSharedMethods::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            targetDecoyCandidatePointersTranched.first(),
            *m_pythiaParameters,
            uniqueMsScanInfosOptimization,
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

    QVector<PythiaParameters> pythiaParametersExperiments;
    e = buildDOE(
            *m_pythiaParameters,
            m_msCalibratomatic->scanTimeStDev(),
            m_pythiaParameters->verbosity,
            &pythiaParametersExperiments
            ); ree;

    double bestResultCount = 0;
    int countSinceLastHigh = 0;
    QVector<DOEResult> results;
    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {

        e = m_targetDecoyCandidatePairScoretron->setPythiaParameters(pythiaParams); ree;

        constexpr int splitter = 2;
        const int threadCount = uniqueMsScanInfos.size() < m_pythiaParameters->threadCount
                      ? std::min(uniqueMsScanInfos.size() * splitter, m_pythiaParameters->threadCount)
                      : m_pythiaParameters->threadCount;

        constexpr bool useExtendedScores = true;
        constexpr bool useNeuralNetworkScores = false;
        constexpr bool useTopNIntegrationsParameter = true;
        constexpr float minPeakCountOptimization = 3.9;
        constexpr int topNMS2Ions = 12;
        m_candidateScorePairs.clear();
        e = m_targetDecoyCandidatePairScoretron->scoreTargetDecoyPairs(
                topNMS2Ions,
                *m_msCalibratomatic,
                minPeakCountOptimization,
                threadCount,
                useExtendedScores,
                useNeuralNetworkScores,
                useTopNIntegrationsParameter,
                mzTargetKeyVsTurboXicPntrs,
                DiscriminantScoretron::defaultWeights(useExtendedScores, useNeuralNetworkScores),
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                &m_candidateScorePairs
                ); ree

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

        QString fdrString;
        e = FDRCLassifierNeuralNet::outPutFDRCounts(fdrVsCounts, &fdrString); ree;
        const double fdrMean = MathUtils::mean(fdrVsCounts.values());

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "ppmTol"
                 << pythiaParams.ms2ExtractionWidthPPM
                 << "fdrMean"
                 << fdrMean
                 << "Finished"
                 << qPrintable(fdrString);

        DOEResult res;
        res.ppm = pythiaParams.ms2ExtractionWidthPPM;
        res.fdrCount = fdrMean;
        results.push_back(res);

        if (res.fdrCount >= bestResultCount) {
            bestResultCount = res.fdrCount;
            countSinceLastHigh = 0;
            m_weights = weights;
            continue;
        }

        constexpr int maxCountsSinceLastHigh = 3;
        if (++countSinceLastHigh >= maxCountsSinceLastHigh) {
            break;
        }
    }

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

    for (TurboXIC* turboXic : mzTargetKeyVsTurboXicPntrs) {delete turboXic;}

    ERR_RETURN
}

int OptimizeMassAccuracyPPMSettertron::calculateNumberOfTranches() const {

    constexpr int optimizationMultiplicationFactor = 25;
    const auto sizePerTranche = static_cast<double>(m_pythiaParameters->trancheSizeMax * optimizationMultiplicationFactor);
    const int numberOfTranches = std::max(static_cast<int>(m_targetDecoyPairPntrs->size() / sizePerTranche), 1);
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "target count:" << m_targetDecoyPairPntrs->size()
             << "sizePerTranche:" << static_cast<int>(sizePerTranche)
             << "%:" << MathUtils::pRound(static_cast<int>(sizePerTranche) / static_cast<float>(m_targetDecoyPairPntrs->size()) * 100, 1);

    return numberOfTranches;
}
