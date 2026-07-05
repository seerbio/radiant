//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron.h"

#include "CandidateScores.h"
#include "DiscriminantScoretron.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "IsotopicDistributionBuilder.h"
#include "MsUtils.h"
#include "ObjectCSVWriters.h"
#include "TargetDecoyCandidatePair.h"
#include "TimsMs2IonMobilityIndex.h"
#include "TurboXIC.h"
#include "XICPeakManager.h"

#include <QDebug>

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

class Q_DECL_HIDDEN CandidateScorertron::Private {
public:

    using NominalMass = int;
    using AveragineVec = QVector<float>;

    Private();
    ~Private();

    Err init(const PythiaParameters &pythiaParameters);
    void resetDiagnostics();
    [[nodiscard]] QString scoringDiagnosticsSummary(const MzTargetKey &mzTargetKey) const;

    Eigen::VectorX<float> m_kernelIntegration;
    Eigen::VectorX<float> m_kernelMs2;

    struct ScoringDiagnostics {
        quint64 scoreCalls = 0;
        quint64 targetCalls = 0;
        quint64 decoyCalls = 0;
        quint64 timsIonMobilityCalls = 0;
        quint64 zeroIntensityMatrix = 0;
        quint64 zeroIonCountVector = 0;
        quint64 zeroProductVector = 0;
        quint64 emptyPeakIntegrations = 0;
        quint64 emptyBestCorrelationResults = 0;
        quint64 lowCorrelationRejected = 0;
        quint64 noDiscriminantCandidate = 0;
        quint64 scoredCandidates = 0;
    };

    ScoringDiagnostics m_scoringDiagnostics;

private:
    PythiaParameters m_pythiaParameters;

};

CandidateScorertron::Private::Private() {}
CandidateScorertron::Private::~Private() {}

Err CandidateScorertron::Private::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_pythiaParameters = pythiaParameters;

    constexpr int order = 1;
    constexpr int derivative = 0;
    constexpr int rate = 1;

    //TODO consider writing algo to optimizing filter length parameters as well as smoothCount later in the code.

    Eigen::MatrixX<float> kernel;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        m_pythiaParameters.filterLengthMS2,
        order,
        derivative,
        rate,
        &kernel
        ); ree;
    const Eigen::VectorX<float> kernelVec(kernel);
    m_kernelMs2 = kernelVec;

    Eigen::MatrixX<float> kernelIntegration;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        m_pythiaParameters.filterLengthIntegration,
        order,
        derivative,
        rate,
        &kernelIntegration
        ); ree;
    const Eigen::VectorX<float> kernelIntegrationVec(kernelIntegration);
    m_kernelIntegration = kernelIntegrationVec;

    ERR_RETURN
}

void CandidateScorertron::Private::resetDiagnostics() {
    m_scoringDiagnostics = ScoringDiagnostics();
}

QString CandidateScorertron::Private::scoringDiagnosticsSummary(const MzTargetKey &mzTargetKey) const {

    return QStringLiteral(
        "Radiant candidate scoring diagnostics target_key=%1 score_calls=%2 target_calls=%3 decoy_calls=%4 tims_im_calls=%5 "
        "zero_intensity_matrix=%6 zero_ion_count_vector=%7 zero_product_vector=%8 empty_peak_integrations=%9 "
        "empty_best_correlations=%10 low_correlation_rejected=%11 no_discriminant_candidate=%12 scored_candidates=%13"
        )
        .arg(mzTargetKey)
        .arg(m_scoringDiagnostics.scoreCalls)
        .arg(m_scoringDiagnostics.targetCalls)
        .arg(m_scoringDiagnostics.decoyCalls)
        .arg(m_scoringDiagnostics.timsIonMobilityCalls)
        .arg(m_scoringDiagnostics.zeroIntensityMatrix)
        .arg(m_scoringDiagnostics.zeroIonCountVector)
        .arg(m_scoringDiagnostics.zeroProductVector)
        .arg(m_scoringDiagnostics.emptyPeakIntegrations)
        .arg(m_scoringDiagnostics.emptyBestCorrelationResults)
        .arg(m_scoringDiagnostics.lowCorrelationRejected)
        .arg(m_scoringDiagnostics.noDiscriminantCandidate)
        .arg(m_scoringDiagnostics.scoredCandidates);
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

CandidateScorertron::CandidateScorertron()
: m_topNMS2Ions(-1)
, m_xicPeakManager(nullptr)
, m_msFrameMzTarget(nullptr)
, m_turboXicMS1(nullptr)
, m_msFrameMS1(nullptr)
, m_msReaderPointerAcc(nullptr)
, m_timsMs2IonMobilityIndex(nullptr)
, d_ptr(QScopedPointer<Private>(new Private))
, m_minPeakCount(3.9)
, m_scanTimeRange(0)
, m_useTopNIntegrationsParam(false)
, m_useAdaptiveTimsMobilityCentering(false)
{}

CandidateScorertron::~CandidateScorertron() {}

void CandidateScorertron::setUseAdaptiveTimsMobilityCentering(
    bool useAdaptiveTimsMobilityCentering
    ) {
    m_useAdaptiveTimsMobilityCentering = useAdaptiveTimsMobilityCentering;
}

QString CandidateScorertron::scoringDiagnosticsSummary() const {
    return d_ptr->scoringDiagnosticsSummary(m_mzTargetKey);
}

void CandidateScorertron::printScoringDiagnosticsIfEnabled() const {
    if (!m_pythiaParameters.writeFullCandidateDebug || d_ptr->m_scoringDiagnostics.scoreCalls == 0) {
        return;
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << qPrintable(scoringDiagnosticsSummary());
}

Err CandidateScorertron::init(
    const PythiaParameters &pythiaParameters,
    const MsCalibratomatic &msCalibratomatic,
    const MzTargetKey &mzTargetKey,
    int topNMS2Ions,
    float minPeakCount,
    float scanTimeRange,
    const QMap<NominalMzMass, QVector<float>> &averagineTable,
    const QVector<Features> &features,
    bool useTopNIntegrationsParameter,
    XICPeakManager *xicPeakManager,
    MsFrame *msFrameMzTarget,
    TurboXIC *turboXicMS1,
    MsFrame *msFrameMS1,
    MsReaderPointerAcc *msReaderPointerAcc,
    TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;
    e = ErrorUtils::isTrue(msFrameMzTarget->isValid()); ree;

	if (turboXicMS1) {
		e = ErrorUtils::isTrue(turboXicMS1->isInit()); ree;
		e = ErrorUtils::isTrue(msFrameMS1->isValid()); ree;
	}

    e = ErrorUtils::isNotEmpty(mzTargetKey); ree;
    e = ErrorUtils::isNotEmpty(features); ree;
    e = ErrorUtils::isNotEmpty(averagineTable); ree;
    e = ErrorUtils::isFalse(MathUtils::tZero(scanTimeRange)); ree;
    e = ErrorUtils::isAboveThreshold(
        topNMS2Ions,
        S_GLOBAL_SETTINGS.MIN_MS2_IONS,
        ErrorUtilsParam::IncludeThreshold)
    ; ree;

    e = ErrorUtils::isAboveThreshold(
        minPeakCount,
        1.0f,
        ErrorUtilsParam::ExcludeThreshold
        ); ree;

    m_pythiaParameters = pythiaParameters;
    m_topNMS2Ions = topNMS2Ions;
    m_mzTargetKey = mzTargetKey;
    m_xicPeakManager = xicPeakManager;
    m_msFrameMzTarget = msFrameMzTarget;
    m_turboXicMS1 = turboXicMS1;
    m_scanTimeRange = scanTimeRange;
    m_averagineTable = averagineTable;
    m_msFrameMS1 = msFrameMS1;
    m_msReaderPointerAcc = msReaderPointerAcc;
    m_timsMs2IonMobilityIndex = timsMs2IonMobilityIndex;
    m_ms1FrameNumbersTIMS.clear();
    m_features = features;
    m_minPeakCount = minPeakCount;
    m_useTopNIntegrationsParam = useTopNIntegrationsParameter;

    if (m_msReaderPointerAcc != nullptr
        && !m_msReaderPointerAcc->ptr.isNull()
        && m_msReaderPointerAcc->ptr->isTIMS()) {
        const QMap<FrameNumberTIMS, Ms1FrameTIMS> *frameNumberVsMs1FrameTIMS
            = m_msReaderPointerAcc->ptr->frameNumberVsMS1FrameTIMSPntr();
        if (frameNumberVsMs1FrameTIMS != nullptr && !frameNumberVsMs1FrameTIMS->isEmpty()) {
            m_ms1FrameNumbersTIMS = frameNumberVsMs1FrameTIMS->keys().toVector();
        }
    }

    if (msCalibratomatic.isInitRT()) {
        m_msCalibratomatic = msCalibratomatic;
    }

    e = d_ptr->init(m_pythiaParameters); ree;
    d_ptr->resetDiagnostics();

    ERR_RETURN
}

class MatriciesAndVecs {

public:

    Eigen::MatrixX<float> intensityMatrix100;
    Eigen::MatrixX<float> intensityMatrix100Shadow;
    Eigen::MatrixX<float> intensityMatrix45;

    Eigen::MatrixX<float> mzMatrix100;

    Eigen::VectorX<float> intensityVec;
    Eigen::VectorX<float> ionCountVec;
    Eigen::VectorX<float> integrationVecCosineSim;
    Eigen::VectorX<float> productVec;

    [[nodiscard]] bool intensityMatriciesAreValid() const {
        return intensityMatrix100.size() > 0 && mzMatrix100.size() > 0;
    }

    [[nodiscard]] bool integrationVecIsValid() const {
        return ionCountVec.size() > 0;
    }
};

class BestCorrelationResult {

public:
    QVector<float> peakCorrelations;
    QVector<float> peakCorrelations45;
    QVector<float> peakCorrelationsWindow1p5X;
    QVector<float> peakCorrelationsWindow2X;
    float peakCorrelationsSum = -1.0;
    Eigen::MatrixX<float> matBlockTrimmedIntensity;
    Eigen::MatrixX<float> matBlockTrimmedIntensityWindow1p5X;
    Eigen::MatrixX<float> matBlockTrimmedIntensityWindow2X;
    Eigen::MatrixX<float> matBlockTrimmedIntensity45;
    Eigen::MatrixX<float> matBlockTrimmedIntensityShadows;
    Eigen::MatrixX<float> matBlockTrimmedMz;
    Eigen::VectorX<float> integrationVector;
    int bestAnchorColumnIndex = -1;
    int bestAnchorRowIndex = -1;
    QVector<int> apexStarts;
    PeakIntegrationIndexes peakIntegrationIndexes;

    [[nodiscard]] QVector<int> columnNonZeroPeakLengths() const {

        const int cols = static_cast<int>(matBlockTrimmedIntensity.cols());
        QVector<int> peakLengths(cols);

        for (int col = 0; col < cols; col++) {
            peakLengths[col] = EigenUtils::nonZeros(matBlockTrimmedIntensity.col(col));
        }

        return peakLengths;
    }
};

namespace {

    constexpr double ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0 = 0.1;
    constexpr double ALPHADIA_TARGET_MOBILITY_TOLERANCE_ONE_OVER_K0 = 0.06;
    constexpr double ALPHADIA_MOBILITY_FWHM_ONE_OVER_K0 = 0.01;
    constexpr double ALPHADIA_RT_FWHM_SECONDS = 5.0;

    struct LocalIonMobilityPeak {
        bool isValid = false;
        float centerDriftTime = -1.0f;
        ScanTime centerScanTime = -1.0f;
        float minDriftTime = -1.0f;
        float maxDriftTime = -1.0f;
        IonMobilityIndex centerIonMobilityIndex = -1;
        FrameIndex centerFrameIndex = -1;
    };

    struct LocalIonMobilityRtEvidencePoint {
        IonMobilityIndex ionMobilityIndex = -1;
        FrameIndex frameIndex = -1;
        float driftTime = -1.0f;
        ScanTime scanTime = -1.0f;
        double intensity = 0.0;
    };

    struct SymmetricProfileLimits {
        int startIndex = -1;
        int stopIndex = -1;
    };

    SymmetricProfileLimits alphaDiaStyleSymmetricLimits1d(
        const QVector<double> &profile,
        int centerIndex,
        double f = 0.95,
        double centerFraction = 0.01,
        int minSize = 3,
        int maxSize = 20
        ) {

        SymmetricProfileLimits limits;

        if (profile.isEmpty() || centerIndex < 0 || centerIndex >= profile.size()) {
            return limits;
        }

        if (profile.size() <= 1) {
            limits.startIndex = centerIndex;
            limits.stopIndex = centerIndex;
            return limits;
        }

        const double centerIntensity = profile.at(centerIndex);
        double trailingIntensity = centerIntensity;
        int limit = std::min(minSize, std::max(centerIndex, profile.size() - centerIndex - 1));

        for (int s = minSize + 1; s < maxSize; ++s) {
            const int lowerIndex = std::max(centerIndex - s, 0);
            const int upperIndex = std::min(centerIndex + s, profile.size() - 1);
            const double intensity = (profile.at(lowerIndex) + profile.at(upperIndex)) / 2.0;

            if (intensity < f * trailingIntensity) {
                if (intensity > centerIntensity * centerFraction) {
                    limit = s;
                    trailingIntensity = intensity;
                }
                else {
                    break;
                }
            }
            else {
                break;
            }
        }

        limits.startIndex = std::max(centerIndex - limit, 0);
        limits.stopIndex = std::min(centerIndex + limit, profile.size() - 1);
        return limits;
    }

    bool containsMs2IonMobilityFeature(const QVector<Features> &features) {
        return features.contains(Ms2IonMobilityWeightedDelta)
            || features.contains(Ms2IonMobilityWeightedDeltaAbs)
            || features.contains(Ms2IonMobilityApexDeltaAbsMean)
            || features.contains(Ms2IonMobilityApexDeltaAbsStDev)
            || features.contains(Ms2IonMobilityMatchedIonFraction)
            || features.contains(Ms2IonMobilityFwhmMean)
            || features.contains(Ms2IonMobilityFwhmStDev)
            || features.contains(Ms2IonMobilityRtCosineMean)
            || features.contains(Ms2IonMobilityRtCosineStDev)
            || features.contains(Ms2IonMobilityRtApexAgreementFraction);
    }

    int closestMs1FrameIndexAtOrBefore(
        const QVector<FrameNumberTIMS> &frameNumbers,
        FrameNumberTIMS scanNumber
        ) {

        if (frameNumbers.isEmpty()) {
            return -1;
        }

        const auto lower = std::lower_bound(frameNumbers.constBegin(), frameNumbers.constEnd(), scanNumber);
        if (lower == frameNumbers.constBegin()) {
            return 0;
        }
        if (lower == frameNumbers.constEnd()) {
            return frameNumbers.size() - 1;
        }

        const int lowerIndex = static_cast<int>(lower - frameNumbers.constBegin());
        const int previousIndex = lowerIndex - 1;
        const FrameNumberTIMS lowerFrameNumber = frameNumbers.at(lowerIndex);
        const FrameNumberTIMS previousFrameNumber = frameNumbers.at(previousIndex);

        if (std::abs(previousFrameNumber - scanNumber) <= std::abs(lowerFrameNumber - scanNumber)) {
            return previousIndex;
        }

        return lowerFrameNumber > scanNumber ? previousIndex : lowerIndex;
    }

    Err sortBestCorrelationResult(QVector<BestCorrelationResult> *bestCorrelationResults) {

        ERR_INIT

        e = ErrorUtils::isFalse(bestCorrelationResults->isEmpty()); ree;

        const float bestScore = std::max_element(
            bestCorrelationResults->begin(),
            bestCorrelationResults->end(),
            [](const BestCorrelationResult &l, const BestCorrelationResult &r){return l.peakCorrelationsSum < r.peakCorrelationsSum;}
            )->peakCorrelationsSum;

        const auto terminatorLogic = [bestScore](const BestCorrelationResult &r) {
            constexpr float tolerance = 1.0;
            return (bestScore - r.peakCorrelationsSum) > tolerance;
        };
        const auto terminator = std::remove_if(
            bestCorrelationResults->begin(),
            bestCorrelationResults->end(),
            terminatorLogic
            );

        bestCorrelationResults->erase(terminator, bestCorrelationResults->end());

        std::sort(
            bestCorrelationResults->begin(),
            bestCorrelationResults->end(),
            [](const BestCorrelationResult &l, const BestCorrelationResult &r) {

// #define USE_PEAK_LENGTH
#ifdef USE_PEAK_LENGTH
                const int lLen = l.peakIntegrationIndexes.second - l.peakIntegrationIndexes.first;
                const int rLen = r.peakIntegrationIndexes.second - r.peakIntegrationIndexes.first;

                if (lLen == rLen) {
                    return l.matBlockTrimmedIntensity.coeff(l.bestAnchorRowIndex, l.bestAnchorColumnIndex)
                         > r.matBlockTrimmedIntensity.coeff(r.bestAnchorRowIndex, r.bestAnchorColumnIndex);
                }

                return lLen > rLen;
#else
                return l.matBlockTrimmedIntensity.sum() > r.matBlockTrimmedIntensity.sum();
#endif

            });

        ERR_RETURN;
    }


    PeptideSequenceWithModsChargeAndTargetKey buildPeptideSequenceWithModsChargeAndTargetKey(const CandidateScores* candidateScores) {
        return candidateScores->targetDecoyCandidatePair->peptideStringWithMods()
             + "|"
             + QString::number(candidateScores->targetDecoyCandidatePair->charge())
             + "|"
             + candidateScores->targetKey;
    }


}//namespace
Err CandidateScorertron::calculateScores(
    const QVector<MS2Ion> &ms2Ions,
    const QVector<float> &weights,
    TargetDecoyCandidatePair* targetDecoyCandidatePair,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2Ions); ree;
    const bool collectScoringDiagnostics = m_pythiaParameters.writeFullCandidateDebug;
    if (collectScoringDiagnostics) {
        d_ptr->m_scoringDiagnostics.scoreCalls++;
        if (candidateScores != nullptr && candidateScores->isDecoy) {
            d_ptr->m_scoringDiagnostics.decoyCalls++;
        }
        else {
            d_ptr->m_scoringDiagnostics.targetCalls++;
        }
        if (m_timsMs2IonMobilityIndex != nullptr
            && m_timsMs2IonMobilityIndex->isInit()
            && targetDecoyCandidatePair != nullptr
            && targetDecoyCandidatePair->iIM() > 0.0f) {
            d_ptr->m_scoringDiagnostics.timsIonMobilityCalls++;
        }
    }

    candidateScores->targetDecoyCandidatePair = targetDecoyCandidatePair;
    candidateScores->initFeaturesArray();
    candidateScores->targetKey = m_mzTargetKey;
    candidateScores->proteinGroup = targetDecoyCandidatePair->proteinGroups();

    //Note, target key must be set before peptideSequenceWithModsChargeAndTargetKey
    candidateScores->peptideSequenceWithModsChargeAndTargetKey = buildPeptideSequenceWithModsChargeAndTargetKey(candidateScores);

    FrameIndex frameIndexPredictedMin;
    FrameIndex frameIndexPredictedMax;
    e = setPredictedFrameIndexes(
        targetDecoyCandidatePair->iRt(candidateScores->isDecoy),
        candidateScores,
        &frameIndexPredictedMin,
        &frameIndexPredictedMax
        );

    MatriciesAndVecs matriciesAndVecs;
    e = initMatricesdAndVecs(
        targetDecoyCandidatePair,
        ms2Ions,
        frameIndexPredictedMin,
        frameIndexPredictedMax,
        &matriciesAndVecs
        ); ree;

    if (matriciesAndVecs.intensityMatrix100.size() == 0
        || MathUtils::tZero(matriciesAndVecs.intensityMatrix100.maxCoeff())) {
        if (collectScoringDiagnostics) {
            d_ptr->m_scoringDiagnostics.zeroIntensityMatrix++;
        }
    }
    if (matriciesAndVecs.ionCountVec.size() == 0
        || MathUtils::tZero(matriciesAndVecs.ionCountVec.maxCoeff())) {
        if (collectScoringDiagnostics) {
            d_ptr->m_scoringDiagnostics.zeroIonCountVector++;
        }
    }
    if (matriciesAndVecs.productVec.size() == 0
        || MathUtils::tZero(matriciesAndVecs.productVec.maxCoeff())) {
        if (collectScoringDiagnostics) {
            d_ptr->m_scoringDiagnostics.zeroProductVector++;
        }
    }

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
    e = EigenUtils::simpleIntegrator(
        matriciesAndVecs.productVec,
        m_pythiaParameters.stopThresholdFractionMS2,
        m_pythiaParameters.peakDifferenceFractionThreshold,
        &peakIntegrationsVsIntensities
        ); ree;

    if (peakIntegrationsVsIntensities.isEmpty()) {
        if (collectScoringDiagnostics) {
            d_ptr->m_scoringDiagnostics.emptyPeakIntegrations++;
        }
        ERR_RETURN
    }

    QVector<BestCorrelationResult> bestCorrelationResults;
    e = processIntegrationVectorPeakIntegrations(
        matriciesAndVecs,
        peakIntegrationsVsIntensities,
        &bestCorrelationResults
        ); ree;

    if (bestCorrelationResults.isEmpty()) {
        if (collectScoringDiagnostics) {
            d_ptr->m_scoringDiagnostics.emptyBestCorrelationResults++;
        }
        ERR_RETURN
    }

    constexpr int multiplierForKeySettingByTen = 10;
    const int nominalMass
        = static_cast<int>((std::round(targetDecoyCandidatePair->mass() / multiplierForKeySettingByTen) * multiplierForKeySettingByTen));
    e = ErrorUtils::isTrue(m_averagineTable.contains(nominalMass)); ;
    const QVector<float> ms1Averagine = m_averagineTable.value(nominalMass);

// #define TURBO_SCORE
#ifdef TURBO_SCORE
    e = sortBestCorrelationResult(&bestCorrelationResults); ree;

    e = setCandidateScores(
        targetDecoyCandidatePair,
        ms2Ions,
        bestCorrelationResults,
        ms1Averagine,
        candidateScores
        ); ree;
#else
    QVector<FeaturesArray> candidateScoresFeatureArrays;
    QVector<CandidateScores> candidateScoresFeatures;
    for (const BestCorrelationResult &bcr : bestCorrelationResults) {

    	if (bcr.peakCorrelationsSum < 0.1) {
            if (collectScoringDiagnostics) {
                d_ptr->m_scoringDiagnostics.lowCorrelationRejected++;
            }
    		continue;
		}

        CandidateScores cs = *candidateScores;
        e = setCandidateScores(
            targetDecoyCandidatePair,
            ms2Ions,
            {bcr},
            ms1Averagine,
            &cs
            ); ree;

        const FeaturesArray fa = CandidateScores::selectFeaturesArrayFeatures(
            cs.featuresArray,
            m_features
            );

        candidateScoresFeatures.push_back(cs);
        candidateScoresFeatureArrays.push_back(fa);
    }

    QVector<FeaturesArray*> featuresArraysPntrs;
    std::transform(
        candidateScoresFeatureArrays.begin(),
        candidateScoresFeatureArrays.end(),
        std::back_inserter(featuresArraysPntrs),
        [](FeaturesArray &f){return &f;}
        );

    constexpr int threadCount = 1;
    QVector<float> discScores;
    e = DiscriminantScoretron::applyWeights(
        weights,
        threadCount,
        featuresArraysPntrs,
        &discScores
        );

    QVector<QPair<float, CandidateScores>> candidateScoresPairs;
    for (int i = 0; i < discScores.size(); ++i) {
        candidateScoresPairs.push_back({discScores[i], candidateScoresFeatures[i]});
    }

    if (candidateScoresPairs.isEmpty()) {
        if (collectScoringDiagnostics) {
            d_ptr->m_scoringDiagnostics.noDiscriminantCandidate++;
        }
        ERR_RETURN
    }

    std::sort(
        candidateScoresPairs.rbegin(),
        candidateScoresPairs.rend(),
        [](const QPair<float, CandidateScores> &l, const QPair<float, CandidateScores> &r) {
            return l.first < r.first;
        });

    const float topDiscScore = candidateScoresPairs.front().first;

    *candidateScores = candidateScoresPairs.front().second;
    if (candidateScoresPairs.size() > 1) {
        candidateScores->featuresArray[DiscScore1stRunnerUpDiff]
            = topDiscScore - candidateScoresPairs.at(1).first;
    }
    if (candidateScoresPairs.size() > 2) {
        candidateScores->featuresArray[DiscScore2ndRunnerUpDiff]
            = topDiscScore - candidateScoresPairs.at(2).first;
    }

    QVector<float> discScoresSubbed;
    std::transform(
        discScores.begin(),
        discScores.end(),
        std::back_inserter(discScoresSubbed),
        [topDiscScore](float df){return topDiscScore - df;}
        );
    candidateScores->featuresArray[DiscScoresCount] = discScoresSubbed.size();
    candidateScores->featuresArray[DiscScoresMean] = MathUtils::mean(discScoresSubbed);
    candidateScores->featuresArray[DiscScoresStDev] = MathUtils::stDev(discScoresSubbed);

    e = setFullTheoMs2IonsScores(candidateScores); ree;
    if (collectScoringDiagnostics) {
        d_ptr->m_scoringDiagnostics.scoredCandidates++;
    }


#endif

// #define TROUBLE_SHOOT_INTEGRATION
#ifdef TROUBLE_SHOOT_INTEGRATION
    if (targetDecoyCandidatePair->peptideStringWithMods() == "TVC(UniMod:4)LPDGSFPSGSEC(UniMod:4)HISGWGVTETGK"
        && targetDecoyCandidatePair->charge() == 3
        && !candidateScores->isDecoy
        ) {

        QMap<FrameIndex, BestCorrelationResult> bestCorrelationResultsInds;
        for (const BestCorrelationResult &bcr : bestCorrelationResults) {
            bestCorrelationResultsInds.insert(bcr.peakIntegrationIndexes.first, bcr);
        }

        std::sort(
            candidateScoresPairs.rbegin(),
            candidateScoresPairs.rend(),
            [](const QPair<float, CandidateScores> &l, const QPair<float, CandidateScores> &r){return l.first < r.first;}
            );

        for (QPair<float, CandidateScores> &b : candidateScoresPairs) {
            qDebug()
            << b.second.frameIndexStart << ","
            << b.second.frameIndexEnd
            << b.first
            << bestCorrelationResultsInds.value(b.second.frameIndexStart).apexStarts
            << "SDLKFJSDL";
            //
            // std::cout << bestCorrelationResultsInds.value(b.second.frameIndexStart).matBlockTrimmedIntensity << std::endl;;
            // std::cout << " **************** " << std::endl;
        }

        qDebug() << weights;

        const QString &intensityVecPath = QStringLiteral("/home/andrewnichols/Repos/Graphing/intensity.csv");
        const QString &prodVecPath = QStringLiteral("/home/andrewnichols/Repos/Graphing/prod.csv");

        e = ObjectCSVWriters::writeVectorToFile(
            EigenUtils::convertEigenVectorToQVector(matriciesAndVecs.intensityVec),
            intensityVecPath
            ); ree;

        e = ObjectCSVWriters::writeVectorToFile(
            EigenUtils::convertEigenVectorToQVector(matriciesAndVecs.productVec),
            prodVecPath
            ); ree;
    }
#endif

    ERR_RETURN
}

namespace {

    void filterXICPointsByAccuracyPPM(
        float mzVal,
        float ppmTol,
        XICPoints *xicPoints
        ) {

        const float massTol = MathUtils::calculatePPM(mzVal, ppmTol);

        const auto terminatorLogic = [mzVal, massTol](const XICPoint &p) {
            return !(mzVal - massTol < p.mz && p.mz < mzVal + massTol);
        };

        const auto terminator = std::remove_if(
            xicPoints->begin(),
            xicPoints->end(),
            terminatorLogic
            );

        xicPoints->erase(terminator, xicPoints->end());
    }

    void filterXICPointsByFrameIndex(
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        XICPoints *xicPoints
        ) {

        const auto terminatorLogic = [frameIndexPredictedMin, frameIndexPredictedMax](const XICPoint &p) {
            return !(frameIndexPredictedMin < p.scanNumber && p.scanNumber < frameIndexPredictedMax);
        };

        const auto terminator = std::remove_if(
            xicPoints->begin(),
            xicPoints->end(),
            terminatorLogic
            );

        xicPoints->erase(terminator, xicPoints->end());
    }

    bool canUseLibraryIonMobilityFilteredMs2(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex,
        float ionMobilityCenter
        ) {

        if (targetDecoyCandidatePair == nullptr
            || timsMs2IonMobilityIndex == nullptr
            || !timsMs2IonMobilityIndex->isInit()
            || ionMobilityCenter <= 0.0f) {
            return false;
        }

        return true;
    }

    Err extractLibraryIonMobilityFilteredTimsMs2Xic(
        const TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex,
        float ionMobilityMin,
        float ionMobilityMax,
        float mzVal,
        float ppmTol,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        XICPoints *xicPoints
        ) {

        ERR_INIT

        xicPoints->clear();

        if (timsMs2IonMobilityIndex == nullptr || !timsMs2IonMobilityIndex->isInit()) {
            ERR_RETURN
        }

        const float massTol = MathUtils::calculatePPM(mzVal, ppmTol);
        const float mzMin = mzVal - massTol;
        const float mzMax = mzVal + massTol;

        *xicPoints = timsMs2IonMobilityIndex->extractPointsXIC(
            mzMin,
            mzMax,
            frameIndexPredictedMin,
            frameIndexPredictedMax,
            ionMobilityMin,
            ionMobilityMax
            );

        ERR_RETURN
    }

    bool driftTimeFromIonMobilityIndex(
        const MsReaderPointerAcc *msReaderPointerAcc,
        const TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex,
        IonMobilityIndex ionMobilityIndex,
        double *driftTime
        ) {

        if (driftTime == nullptr) {
            return false;
        }

        float indexedDriftTime = -1.0f;
        if (timsMs2IonMobilityIndex != nullptr
            && timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(ionMobilityIndex, &indexedDriftTime)) {
            *driftTime = indexedDriftTime;
            return true;
        }

        if (msReaderPointerAcc == nullptr || msReaderPointerAcc->ptr.isNull()) {
            return false;
        }

        return msReaderPointerAcc->ptr->driftTimeFromIonMobilityIndex(ionMobilityIndex, driftTime) == eNoError;
    }

    void updateApexFromSortedScanPoints(
        const ScanPoints &scanPoints,
        float mzMin,
        float mzMax,
        IonMobilityIndex ionMobilityIndex,
        double driftTime,
        float *apexIntensity,
        IonMobilityIndex *apexIonMobilityIndex,
        double *apexDriftTime
        ) {

        const auto lower = std::lower_bound(
            scanPoints.constBegin(),
            scanPoints.constEnd(),
            mzMin,
            [](const ScanPoint &scanPoint, float mz) {
                return scanPoint.x() < mz;
            }
            );

        for (auto it = lower; it != scanPoints.constEnd() && it->x() <= mzMax; ++it) {
            const float intensity = it->y();
            if (intensity <= *apexIntensity) {
                continue;
            }

            *apexIntensity = intensity;
            *apexIonMobilityIndex = ionMobilityIndex;
            *apexDriftTime = driftTime;
        }
    }

    LocalIonMobilityPeak selectLocalIonMobilityPeakForTimsMs2(
        const TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex,
        const MsFrame *msFrameMzTarget,
        const QVector<MS2Ion> &ms2Ions,
        float libraryIonMobility,
        float ppmTol,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax
        ) {

        LocalIonMobilityPeak peak;

        if (timsMs2IonMobilityIndex == nullptr
            || !timsMs2IonMobilityIndex->isInit()
            || msFrameMzTarget == nullptr
            || !msFrameMzTarget->isValid()
            || ms2Ions.isEmpty()) {
            return peak;
        }

        const float initialIonMobilityMin = libraryIonMobility - static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
        const float initialIonMobilityMax = libraryIonMobility + static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);

        QMap<QPair<IonMobilityIndex, FrameIndex>, double> mobilityFrameVsIntensity;
        for (const MS2Ion &ms2Ion : ms2Ions) {

            const float massTol = MathUtils::calculatePPM(ms2Ion.mz, ppmTol);
            const XICPoints xicPoints = timsMs2IonMobilityIndex->extractPointsXIC(
                ms2Ion.mz - massTol,
                ms2Ion.mz + massTol,
                frameIndexPredictedMin,
                frameIndexPredictedMax,
                initialIonMobilityMin,
                initialIonMobilityMax
                );

            for (const XICPoint &xicPoint : xicPoints) {
                if (xicPoint.ionMobilityIndex < 0 || xicPoint.intensity <= 0.0f) {
                    continue;
                }
                mobilityFrameVsIntensity[{xicPoint.ionMobilityIndex, xicPoint.scanNumber}] += xicPoint.intensity;
            }
        }

        if (mobilityFrameVsIntensity.isEmpty()) {
            return peak;
        }

        QVector<LocalIonMobilityRtEvidencePoint> evidencePoints;
        evidencePoints.reserve(mobilityFrameVsIntensity.size());
        for (auto it = mobilityFrameVsIntensity.constBegin(); it != mobilityFrameVsIntensity.constEnd(); ++it) {

            float driftTime = -1.0f;
            if (!timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(it.key().first, &driftTime)) {
                continue;
            }

            LocalIonMobilityRtEvidencePoint evidencePoint;
            evidencePoint.ionMobilityIndex = it.key().first;
            evidencePoint.frameIndex = it.key().second;
            evidencePoint.driftTime = driftTime;
            evidencePoint.scanTime = msFrameMzTarget->scanTimeFromFrameIndex(it.key().second);
            evidencePoint.intensity = it.value();
            evidencePoints.push_back(evidencePoint);
        }

        if (evidencePoints.isEmpty()) {
            return peak;
        }

        const double mobilitySigma = ALPHADIA_MOBILITY_FWHM_ONE_OVER_K0 / 2.3548;
        const double rtSigma = ALPHADIA_RT_FWHM_SECONDS / 2.3548;
        const double twoMobilitySigmaSquared = 2.0 * mobilitySigma * mobilitySigma;
        const double twoRtSigmaSquared = 2.0 * rtSigma * rtSigma;
        double bestSmoothedIntensity = 0.0;

        for (const LocalIonMobilityRtEvidencePoint &centerPoint : evidencePoints) {

            double smoothedIntensity = 0.0;
            for (const LocalIonMobilityRtEvidencePoint &evidencePoint : evidencePoints) {

                const double mobilityDelta = evidencePoint.driftTime - centerPoint.driftTime;
                if (std::abs(mobilityDelta) > (3.0 * mobilitySigma)) {
                    continue;
                }

                const double rtDelta = evidencePoint.scanTime - centerPoint.scanTime;
                if (std::abs(rtDelta) > (3.0 * rtSigma)) {
                    continue;
                }

                const double mobilityWeight = std::exp(-(mobilityDelta * mobilityDelta) / twoMobilitySigmaSquared);
                const double rtWeight = std::exp(-(rtDelta * rtDelta) / twoRtSigmaSquared);
                smoothedIntensity += evidencePoint.intensity * mobilityWeight * rtWeight;
            }

            if (smoothedIntensity <= bestSmoothedIntensity) {
                continue;
            }

            bestSmoothedIntensity = smoothedIntensity;
            peak.isValid = true;
            peak.centerDriftTime = centerPoint.driftTime;
            peak.centerScanTime = centerPoint.scanTime;
            peak.centerIonMobilityIndex = centerPoint.ionMobilityIndex;
            peak.centerFrameIndex = centerPoint.frameIndex;
        }

        if (!peak.isValid) {
            return peak;
        }

        peak.minDriftTime = peak.centerDriftTime - static_cast<float>(ALPHADIA_TARGET_MOBILITY_TOLERANCE_ONE_OVER_K0);
        peak.maxDriftTime = peak.centerDriftTime + static_cast<float>(ALPHADIA_TARGET_MOBILITY_TOLERANCE_ONE_OVER_K0);

        return peak;
    }

    LocalIonMobilityPeak selectMobilityProfilePeakForTimsMs2(
        const TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex,
        const QVector<MS2Ion> &ms2Ions,
        float libraryIonMobility,
        float ppmTol,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax
        ) {

        LocalIonMobilityPeak peak;

        if (timsMs2IonMobilityIndex == nullptr
            || !timsMs2IonMobilityIndex->isInit()
            || ms2Ions.isEmpty()
            || libraryIonMobility <= 0.0f) {
            return peak;
        }

        constexpr int maxProfileFragments = 6;
        constexpr double mobilityPriorSigmaOneOverK0 = 0.04;
        constexpr double minWeightedProfileIntensity = 1.0;
        constexpr double adaptiveExtractionHalfWidthOneOverK0 = ALPHADIA_TARGET_MOBILITY_TOLERANCE_ONE_OVER_K0;

        const int ionCount = std::min(maxProfileFragments, ms2Ions.size());
        const float broadIonMobilityMin = libraryIonMobility - static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
        const float broadIonMobilityMax = libraryIonMobility + static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);

        QMap<IonMobilityIndex, double> summedMobilityProfile;
        for (int i = 0; i < ionCount; ++i) {
            const MS2Ion &ms2Ion = ms2Ions.at(i);
            const float massTol = MathUtils::calculatePPM(ms2Ion.mz, ppmTol);

            QMap<IonMobilityIndex, double> fragmentProfile;
            float apexIntensity = 0.0f;
            float apexDeltaAbs = static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
            if (!timsMs2IonMobilityIndex->extractMobilityProfile(
                    ms2Ion.mz - massTol,
                    ms2Ion.mz + massTol,
                    frameIndexPredictedMin,
                    frameIndexPredictedMax,
                    broadIonMobilityMin,
                    broadIonMobilityMax,
                    libraryIonMobility,
                    &fragmentProfile,
                    &apexIntensity,
                    &apexDeltaAbs
                    )) {
                continue;
            }

            for (auto it = fragmentProfile.constBegin(); it != fragmentProfile.constEnd(); ++it) {
                summedMobilityProfile[it.key()] += it.value();
            }
        }

        if (summedMobilityProfile.isEmpty()) {
            return peak;
        }

        const double mobilitySigma = ALPHADIA_MOBILITY_FWHM_ONE_OVER_K0 / 2.3548;
        const double twoMobilitySigmaSquared = 2.0 * mobilitySigma * mobilitySigma;
        const double twoPriorSigmaSquared = 2.0 * mobilityPriorSigmaOneOverK0 * mobilityPriorSigmaOneOverK0;
        double bestWeightedIntensity = 0.0;

        for (auto centerIt = summedMobilityProfile.constBegin(); centerIt != summedMobilityProfile.constEnd(); ++centerIt) {
            float centerDriftTime = -1.0f;
            if (!timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(centerIt.key(), &centerDriftTime)) {
                continue;
            }

            double smoothedIntensity = 0.0;
            for (auto profileIt = summedMobilityProfile.constBegin(); profileIt != summedMobilityProfile.constEnd(); ++profileIt) {
                float driftTime = -1.0f;
                if (!timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(profileIt.key(), &driftTime)) {
                    continue;
                }

                const double mobilityDelta = driftTime - centerDriftTime;
                if (std::abs(mobilityDelta) > (3.0 * mobilitySigma)) {
                    continue;
                }

                const double mobilityWeight = std::exp(-(mobilityDelta * mobilityDelta) / twoMobilitySigmaSquared);
                smoothedIntensity += profileIt.value() * mobilityWeight;
            }

            const double libraryDelta = centerDriftTime - libraryIonMobility;
            const double priorWeight = std::exp(-(libraryDelta * libraryDelta) / twoPriorSigmaSquared);
            const double weightedIntensity = smoothedIntensity * priorWeight;
            if (weightedIntensity <= bestWeightedIntensity) {
                continue;
            }

            bestWeightedIntensity = weightedIntensity;
            peak.isValid = true;
            peak.centerDriftTime = centerDriftTime;
            peak.centerIonMobilityIndex = centerIt.key();
        }

        if (!peak.isValid || bestWeightedIntensity < minWeightedProfileIntensity) {
            return LocalIonMobilityPeak();
        }

        peak.minDriftTime = peak.centerDriftTime - static_cast<float>(adaptiveExtractionHalfWidthOneOverK0);
        peak.maxDriftTime = peak.centerDriftTime + static_cast<float>(adaptiveExtractionHalfWidthOneOverK0);

        return peak;
    }

    Err getLibraryIonMobilityFilteredTimsMs2XICs(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const QVector<MS2Ion> &ms2Ions,
        const TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex,
        const MsFrame *msFrameMzTarget,
        float ionMobilityCenter,
        float ppmTol,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        float targetedIonMobilityWindowHalfWidth,
        bool useAdaptiveTimsMobilityCentering,
        QVector<XICPoints> *xicPointsVec100,
        QVector<XICPoints> *xicPointsVec100Shadows,
        QVector<XICPoints> *xicPointsVec45
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;

        xicPointsVec100->reserve(ms2Ions.size());
        xicPointsVec100Shadows->reserve(ms2Ions.size());
        xicPointsVec45->reserve(ms2Ions.size());

        LocalIonMobilityPeak localIonMobilityPeak;
        localIonMobilityPeak.isValid = true;
        localIonMobilityPeak.centerDriftTime = ionMobilityCenter;
        localIonMobilityPeak.centerScanTime = -1.0f;
        localIonMobilityPeak.minDriftTime = ionMobilityCenter - targetedIonMobilityWindowHalfWidth;
        localIonMobilityPeak.maxDriftTime = ionMobilityCenter + targetedIonMobilityWindowHalfWidth;

        if (useAdaptiveTimsMobilityCentering) {
            const LocalIonMobilityPeak observedMobilityPeak = selectMobilityProfilePeakForTimsMs2(
                timsMs2IonMobilityIndex,
                ms2Ions,
                ionMobilityCenter,
                ppmTol,
                frameIndexPredictedMin,
                frameIndexPredictedMax
                );

            constexpr float maxAdaptiveCenterShiftOneOverK0 = 0.04f;
            if (observedMobilityPeak.isValid
                && std::abs(observedMobilityPeak.centerDriftTime - ionMobilityCenter)
                       <= maxAdaptiveCenterShiftOneOverK0) {
                localIonMobilityPeak.centerDriftTime = observedMobilityPeak.centerDriftTime;
                localIonMobilityPeak.centerIonMobilityIndex = observedMobilityPeak.centerIonMobilityIndex;
                localIonMobilityPeak.minDriftTime = observedMobilityPeak.centerDriftTime - targetedIonMobilityWindowHalfWidth;
                localIonMobilityPeak.maxDriftTime = observedMobilityPeak.centerDriftTime + targetedIonMobilityWindowHalfWidth;
            }
        }

        for (const MS2Ion &ms2Ion : ms2Ions) {

            XICPoints xicPoints;
            e = extractLibraryIonMobilityFilteredTimsMs2Xic(
                timsMs2IonMobilityIndex,
                localIonMobilityPeak.minDriftTime,
                localIonMobilityPeak.maxDriftTime,
                ms2Ion.mz,
                ppmTol,
                frameIndexPredictedMin,
                frameIndexPredictedMax,
                &xicPoints
                ); ree;

            XICPoints xicPointsShadows;
            const float isotopeDistanceThomsons = S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge;
            e = extractLibraryIonMobilityFilteredTimsMs2Xic(
                timsMs2IonMobilityIndex,
                localIonMobilityPeak.minDriftTime,
                localIonMobilityPeak.maxDriftTime,
                ms2Ion.mz - isotopeDistanceThomsons,
                ppmTol,
                frameIndexPredictedMin,
                frameIndexPredictedMax,
                &xicPointsShadows
                ); ree;

            xicPointsVec100->push_back(xicPoints);
            xicPointsVec100Shadows->push_back(xicPointsShadows);

            filterXICPointsByAccuracyPPM(
                ms2Ion.mz,
                ppmTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
                &xicPoints
                );
            xicPointsVec45->push_back(xicPoints);
        }

        ERR_RETURN
    }

    Err getXICs(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const QVector<MS2Ion> &ms2Ions,
        float ppmTol,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        XICPeakManager *xicPeakManager,
        const MsFrame *msFrameMzTarget,
        const TimsMs2IonMobilityIndex *timsMs2IonMobilityIndex,
        float ionMobilityCenter,
        float targetedIonMobilityWindowHalfWidth,
        bool useAdaptiveTimsMobilityCentering,
        QVector<XICPoints> *xicPointsVec100,
        QVector<XICPoints> *xicPointsVec100Shadows,
        QVector<XICPoints> *xicPointsVec45
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;

        if (canUseLibraryIonMobilityFilteredMs2(
                targetDecoyCandidatePair,
                timsMs2IonMobilityIndex,
                ionMobilityCenter
                )) {

            e = getLibraryIonMobilityFilteredTimsMs2XICs(
                targetDecoyCandidatePair,
                ms2Ions,
                timsMs2IonMobilityIndex,
                msFrameMzTarget,
                ionMobilityCenter,
                ppmTol,
                frameIndexPredictedMin,
                frameIndexPredictedMax,
                targetedIonMobilityWindowHalfWidth,
                useAdaptiveTimsMobilityCentering,
                xicPointsVec100,
                xicPointsVec100Shadows,
                xicPointsVec45
                ); ree;

            ERR_RETURN
        }

        xicPointsVec100->reserve(ms2Ions.size());
        xicPointsVec100Shadows->reserve(ms2Ions.size());
        xicPointsVec45->reserve(ms2Ions.size());

        for (int i = 0; i < ms2Ions.size(); i++) {

            const MS2Ion &ms2Ion = ms2Ions.at(i);

            XICPoints xicPoints;
            e = xicPeakManager->getXIC(ms2Ion.mz, &xicPoints); ree;

            if (frameIndexPredictedMax > 0) {
                filterXICPointsByFrameIndex(
                frameIndexPredictedMin,
                frameIndexPredictedMax,
                &xicPoints
                );
            }

            if (xicPoints.empty()) {
                xicPointsVec100->push_back({});
                xicPointsVec45->push_back({});
                xicPointsVec100Shadows->push_back({});
                continue;
            }

            XICPoints xicPointsShadows;
            const float isotopeDistanceThomsons = S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge;
            e = xicPeakManager->getXIC(ms2Ion.mz - isotopeDistanceThomsons, &xicPointsShadows); ree;
            if (xicPointsShadows.empty()) {
                xicPointsVec100Shadows->push_back({});
            }
            else {
                if (frameIndexPredictedMax > 0) {
                    filterXICPointsByFrameIndex(
                        frameIndexPredictedMin,
                        frameIndexPredictedMax,
                        &xicPointsShadows
                        );
                }
                xicPointsVec100Shadows->push_back(xicPointsShadows);
            }

            xicPointsVec100->push_back(xicPoints);

            filterXICPointsByAccuracyPPM(
                ms2Ion.mz,
                ppmTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
                &xicPoints
                );
            xicPointsVec45->push_back(xicPoints);

        }

        ERR_RETURN
    }

    FrameIndex findFrameIndexMaxXICPointsVec (const QVector<XICPoints> &xicPointsVec100) {

        FrameIndex frameIndexMax = -1;

        for (const XICPoints &xicPoints : xicPointsVec100) {

            if (xicPoints.empty()) {
                continue;
            }

            const FrameIndex frameIndexMaxXICPoints = std::max_element(
                xicPoints.begin(),
                xicPoints.end(),
                [](const XICPoint &l, const XICPoint &r){return l.scanNumber < r.scanNumber;}
                )->scanNumber;
            frameIndexMax = std::max(frameIndexMax, frameIndexMaxXICPoints);
        }

        return frameIndexMax;
    }

    Err buildEigenMatrix(
        const QVector<XICPoints> &xicPointsVec,
        const Eigen::VectorX<float> &kernelMs2,
        FrameIndex frameIndexMax,
        bool buildMzMatrix,
        int smoothCount,
        Eigen::MatrixX<float> *matIntensity,
        Eigen::MatrixX<float> *matMz
        ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(xicPointsVec); ree;

        const FrameIndex frameIndexBuffer = 2;

        const int rows = frameIndexMax + frameIndexBuffer;

        matIntensity->resize(rows, xicPointsVec.size());
        matIntensity->setZero();

        matMz->resize(rows, xicPointsVec.size());
        matMz->setZero();

        for (int col = 0; col < xicPointsVec.size(); col++) {

            const XICPoints &xicPointsCol = xicPointsVec.at(col);
            for (const XICPoint &p : xicPointsCol) {

                if (p.scanNumber >= rows) {
                    continue;
                }

                matIntensity->coeffRef(p.scanNumber, col) += p.intensity;
                if (buildMzMatrix) {

                    if (matMz->coeff(p.scanNumber, col) > 0) {
                        matMz->coeffRef(p.scanNumber, col) += p.mz;
                        matMz->coeffRef(p.scanNumber, col) /= 2.0;
                        continue;
                    }

                    matMz->coeffRef(p.scanNumber, col) = p.mz;
                }
            }
        }

        for (int smooth = 0; smooth < smoothCount; smooth++) {
            // *matIntensity= EigenKernelUtils::applyKernelToEachColumnInMatrix(*matIntensity, kernelMs2);
        }

        ERR_RETURN
    }

    Err buildIntegrationVector(
        const MatriciesAndVecs &matriciesAndVecs,
        const Eigen::VectorX<float> &kernelIntegration,
        float minPeakCount,
        int maxAnchorColumnIndex,
        Eigen::VectorX<float> *ionCountVec
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;

        Eigen::MatrixX<float> matCount = matriciesAndVecs.intensityMatrix100.leftCols(
                std::min(maxAnchorColumnIndex,
                static_cast<int>(matriciesAndVecs.intensityMatrix100.cols())
                ));

        constexpr float intensityThresholdVal = 0.1;
        constexpr float countValue = 1.0;

        matCount = (matCount.array() > intensityThresholdVal).select(countValue, matCount);
        EigenUtils::thresholdMatrix(0.0f, &matCount);

        Eigen::VectorX<float> integrationVecLocal = matCount.rowwise().sum();
        EigenUtils::thresholdVector(minPeakCount, &integrationVecLocal);

        *ionCountVec = EigenKernelUtils::convolveVectorWithKernel(
	        integrationVecLocal,
	        kernelIntegration
	        );

        ERR_RETURN
    }

    Err buildIntegrationVectorCosineSim(
        const MatriciesAndVecs &matriciesAndVecs,
        const QVector<MS2Ion> &ms2IonsTheo,
        const Eigen::VectorX<float> &kernelIntegration,
        int maxAnchorColumnIndex,
        Eigen::VectorX<float> *integrationVecCosineSim
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;
        e = ErrorUtils::isNotEmpty(ms2IonsTheo); ree;
        e = ErrorUtils::isEqual(ms2IonsTheo.size(), static_cast<int>(matriciesAndVecs.intensityMatrix100.cols())); ree;

        const int columnCount = std::min(maxAnchorColumnIndex, static_cast<int>(matriciesAndVecs.intensityMatrix100.cols()));

        Eigen::MatrixX<float> matIntensityVals = matriciesAndVecs.intensityMatrix100.leftCols(columnCount);

        const int matRows = static_cast<int>(matIntensityVals.rows());
        Eigen::MatrixX<Intensity> matMs2IonsIntensityVals(matRows, matIntensityVals.cols());

        for (int i = 0; i < columnCount; i++) {
            const MS2Ion &ms2Ion = ms2IonsTheo.at(i);
            const QVector<Intensity> ms2IonIntensityCol(matRows, ms2Ion.intensity);

            const Eigen::VectorX<Intensity> v = EigenUtils::convertQVectorToEigenVector(ms2IonIntensityCol);
            matMs2IonsIntensityVals.col(i) = v;
        }

        Eigen::VectorX<float> integrationVecCosineSimLocal;
        e = EigenUtils::rowWiseCosineSimilarOfMatrices(
            matIntensityVals,
            matMs2IonsIntensityVals,
            &integrationVecCosineSimLocal
            ); ree;

        integrationVecCosineSimLocal = integrationVecCosineSimLocal.array().pow(1.0f/3.0f);
        *integrationVecCosineSim = EigenKernelUtils::convolveVectorWithKernel(
            integrationVecCosineSimLocal,
            kernelIntegration
            );

        ERR_RETURN

    }

}//namespace
Err CandidateScorertron::initMatricesdAndVecs(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const QVector<MS2Ion> &ms2Ions,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        MatriciesAndVecs *matriciesAndVecs
        ) const {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;
        e = ErrorUtils::isAboveThreshold(m_topNMS2Ions, 0, ErrorUtilsParam::ExcludeThreshold); ree;

        QVector<MS2Ion> ms2IonsResized = ms2Ions;
        ms2IonsResized.resize(std::min(m_topNMS2Ions, ms2Ions.size()));

        const bool ms2IonsAreSorted = std::is_sorted(
            ms2IonsResized.rbegin(),
            ms2IonsResized.rend(),
            [](const MS2Ion &l, const MS2Ion &r){return l.intensity < r.intensity;}
            );
        e = ErrorUtils::isTrue(ms2IonsAreSorted); ree;

        QVector<XICPoints> xicPointsVec100;
        QVector<XICPoints> xicPointsVec100Shadow;
        QVector<XICPoints> xicPointsVec45;
        const float calibratedIonMobilityCenter = ionMobilityCenter(targetDecoyCandidatePair);
        e = getXICs(
            targetDecoyCandidatePair,
            ms2IonsResized,
            static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM),
            frameIndexPredictedMin,
            frameIndexPredictedMax,
            m_xicPeakManager,
            m_msFrameMzTarget,
            m_timsMs2IonMobilityIndex,
            calibratedIonMobilityCenter,
            static_cast<float>(m_pythiaParameters.timsTargetedMs2IonMobilityWindow),
            m_useAdaptiveTimsMobilityCentering,
            &xicPointsVec100,
            &xicPointsVec100Shadow,
            &xicPointsVec45
            ); ree;

        e = ErrorUtils::isEqual(xicPointsVec100.size(), xicPointsVec45.size()); ree;
        e = ErrorUtils::isEqual(xicPointsVec100.size(), xicPointsVec100Shadow.size()); ree;

        const FrameIndex frameIndexMax = findFrameIndexMaxXICPointsVec(xicPointsVec100);

        constexpr int smoothCount = 1;
        e = buildEigenMatrix(
            xicPointsVec100,
            d_ptr->m_kernelMs2,
            frameIndexMax,
            true,
            smoothCount,
            &matriciesAndVecs->intensityMatrix100,
            &matriciesAndVecs->mzMatrix100
            ); ree;
        // matriciesAndVecs->intensityVec = EigenKernelUtils::convolveVectorWithKernel(
        //     matriciesAndVecs->intensityVec,
        //     d_ptr->m_kernelIntegration
        //     );
        // matriciesAndVecs->intensityVec = matriciesAndVecs->intensityVec.array().log10();
        // matriciesAndVecs->intensityVec /= matriciesAndVecs->intensityVec.maxCoeff();
        // EigenUtils::replaceInf(0.0f, &matriciesAndVecs->intensityVec);
        // matriciesAndVecs->intensityVec *= matriciesAndVecs->intensityVec.maxCoeff();

        Eigen::MatrixX<float> unused;
        e = buildEigenMatrix(
            xicPointsVec100Shadow,
            d_ptr->m_kernelMs2,
            frameIndexMax,
            false,
            smoothCount,
            &matriciesAndVecs->intensityMatrix100Shadow,
            &unused
            ); ree;

        matriciesAndVecs->intensityMatrix100 = m_pythiaParameters.subtractShadows
                          ? matriciesAndVecs->intensityMatrix100 - matriciesAndVecs->intensityMatrix100Shadow
                          : matriciesAndVecs->intensityMatrix100;
        EigenUtils::thresholdMatrix(0.0f, &matriciesAndVecs->intensityMatrix100);

        matriciesAndVecs->intensityVec = matriciesAndVecs->intensityMatrix100.rowwise().sum();

        e = buildIntegrationVector(
            *matriciesAndVecs,
            d_ptr->m_kernelMs2,
            m_minPeakCount,
            m_pythiaParameters.maxAnchorColumnIndex,
            &matriciesAndVecs->ionCountVec
            ); ree;

        matriciesAndVecs->productVec = matriciesAndVecs->ionCountVec.array()
                                     * matriciesAndVecs->intensityVec.array();

		constexpr int smoothCountOverride = 3;
		for (int i = 0; i < smoothCountOverride; ++i) {
			matriciesAndVecs->productVec = EigenKernelUtils::convolveVectorWithKernel(
				matriciesAndVecs->productVec,
				d_ptr->m_kernelMs2
				);
		}

        constexpr int noSmooths = 0;
        e = buildEigenMatrix(
            xicPointsVec45,
            d_ptr->m_kernelMs2,
            frameIndexMax,
            false,
            noSmooths,
            &matriciesAndVecs->intensityMatrix45,
            &unused
            ); ree;

        ERR_RETURN
    }

Err CandidateScorertron::setPredictedFrameIndexes(
    float iRT,
    CandidateScores *candidateScores,
    FrameIndex *frameIndexPredictedMin,
    FrameIndex *frameIndexPredictedMax
    ) const {

    ERR_INIT

    if (m_msCalibratomatic.isInitRT()) {

        const float scanTimeWindow
            = m_msCalibratomatic.scanTimeStDev(static_cast<float>(m_pythiaParameters.scanTimeWindowStDevs));

        e = m_msCalibratomatic.predictScanTime(
                iRT,
                &candidateScores->scanTimePredicted
        ); ree;

        e = m_msFrameMzTarget->frameIndexFromScanTime(
                candidateScores->scanTimePredicted - scanTimeWindow,
                frameIndexPredictedMin
        ); ree;

        e = m_msFrameMzTarget->frameIndexFromScanTime(
                candidateScores->scanTimePredicted + scanTimeWindow,
                frameIndexPredictedMax
        ); ree;

        ERR_RETURN
    }

    *frameIndexPredictedMin = 0;
    *frameIndexPredictedMax = 0;

    ERR_RETURN
}

float CandidateScorertron::ionMobilityCenter(
    const TargetDecoyCandidatePair *targetDecoyCandidatePair
    ) const {

    if (targetDecoyCandidatePair == nullptr) {
        return -1.0f;
    }

    const float libraryIonMobility = targetDecoyCandidatePair->iIM();
    if (libraryIonMobility <= 0.0f) {
        return libraryIonMobility;
    }

    if (m_msReaderPointerAcc == nullptr
        || m_msReaderPointerAcc->ptr.isNull()
        || !m_msReaderPointerAcc->ptr->isTIMS()) {
        return libraryIonMobility;
    }

    if (!m_msCalibratomatic.isInitIM()) {
        return libraryIonMobility;
    }

    float calibratedIonMobility = -1.0f;
    const Err e = m_msCalibratomatic.predictIonMobility(
        libraryIonMobility,
        &calibratedIonMobility
        );
    if (e != eNoError || calibratedIonMobility <= 0.0f) {
        return libraryIonMobility;
    }

    return calibratedIonMobility;
}

namespace {

    QPair<PeakIntegrationIndexes, Intensity> correctPeakIntegrationForSingleRow(
        const QPair<PeakIntegrationIndexes, Intensity> &pii,
        const Eigen::VectorX<float> &ionCountVec
        ) {
        QPair<PeakIntegrationIndexes, Intensity> piiWorking = pii;
        if (piiWorking.first.first == piiWorking.first.second) {
            constexpr int bufferDistance = 1;
            piiWorking.first.first = std::max(0, piiWorking.first.first - bufferDistance);
            piiWorking.first.second = std::min(
                static_cast<int>(ionCountVec.size()) - 1,
                piiWorking.first.second + bufferDistance
                );
        }
        return piiWorking;
    }

    QVector<QVector<int>> getMatrxColumnApexes(const Eigen::MatrixX<float> &matBlock) {

        QVector<QVector<int>> apexIndexesByColumn(matBlock.cols());
        for (int col = 0; col < matBlock.cols(); col++) {
            const Eigen::VectorX<float> &column = matBlock.col(col);
            apexIndexesByColumn[col] = EigenUtils::apexesIndexesOnly(column);
        }

        return apexIndexesByColumn;
    }

    QVector<int> findStartApexes(
            const QVector<QVector<int>> &apexIndexes,
            int apexIndex
            ) {

            QVector<int> apexesStartsToUse(apexIndexes.size(), -1);
            for (int col = 0; col < apexIndexes.size(); col++) {

                const QVector<int> &apexesColumn = apexIndexes.at(col);
                if(apexesColumn.isEmpty()) {
                    continue;
                }

                apexesStartsToUse[col] = apexesColumn.at(MathUtils::closest(apexesColumn, apexIndex));
            }

            return apexesStartsToUse;
        }

    QPair<int, int> simpleIntegratorTrimmer(
        const Eigen::VectorX<float> &vec,
        int apexRowIndex,
        float stopThresholdFraction
        ) {

        const float nearZero = 0.001;
        const float apexIndexValue = vec.coeff(apexRowIndex);
        const float stopThresholdValue = apexIndexValue * stopThresholdFraction;

        int rightStopIndex = apexRowIndex;
        int rightCurrentIndex = apexRowIndex;
        while (rightCurrentIndex < vec.size()) {

            const float currentValue = vec.coeff(rightCurrentIndex);
            if (currentValue < nearZero || currentValue <= stopThresholdValue) {
                rightStopIndex = rightCurrentIndex;
                break;
            }

            if (currentValue > nearZero) {
                rightStopIndex = rightCurrentIndex;
                rightCurrentIndex++;
                continue;
            }

            break;
        }

        int leftStopIndex = apexRowIndex;
        int leftCurrentIndex = apexRowIndex;
        while (leftCurrentIndex >= 0) {

            const float currentValue = vec(leftCurrentIndex);
            if (currentValue < nearZero || currentValue <= stopThresholdValue) {
                leftStopIndex = leftCurrentIndex;
                break;
            }

            if (currentValue > nearZero) {
                leftStopIndex = leftCurrentIndex;
                leftCurrentIndex--;
                continue;
            }

            break;
        }

        return {std::max(leftStopIndex, 0), std::min(rightStopIndex, static_cast<int>(vec.size()) - 1)};
    }

    Eigen::MatrixX<float> trimMatrixBlock(
        const Eigen::MatrixX<float> &matBlock,
        const QVector<int> &apexStarts,
        float stopThresholdFraction
        ) {

        Eigen::MatrixX<float> matBlockTrimmedColumns(matBlock.rows(), matBlock.cols());
        matBlockTrimmedColumns.setZero();

        for (int col = 0; col < matBlock.cols(); col++) {

            const Eigen::VectorX<float> &colVec = matBlock.col(col);
            const QPair<int, int> peakLimits = simpleIntegratorTrimmer(
                colVec,
                apexStarts.at(col),
                stopThresholdFraction
                );

            for(int row = peakLimits.first; row <= peakLimits.second; row++) {
                matBlockTrimmedColumns.coeffRef(row, col) = colVec.coeff(row);
            }
        }
        return matBlockTrimmedColumns;
    }

    Err findBestAnchorColumn(
        const Eigen::VectorX<float> &integrationVecSeg,
        const Eigen::MatrixX<float> &matBlockTrimmed,
        int maxAnchorColumnIndex,
        QVector<float> *peakCorrelations,
        int *bestAnchorColumnIndex,
        int *bestAnchorRowIndex
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matBlockTrimmed.size() > 0); ree;
        e = ErrorUtils::isEqual(integrationVecSeg.size(), matBlockTrimmed.rows()); ree;

        const int colCount = matBlockTrimmed.cols();

    	peakCorrelations->resize(colCount);
    	float bestCorrelation = 0;
        for (int col = 0; col < colCount; col++) {

            const Eigen::VectorX<float> &anchor = matBlockTrimmed.col(col);

        	float cosineSim;
        	e = EigenUtils::cosineSimilarity(anchor, integrationVecSeg, &cosineSim); ree;
        	(*peakCorrelations)[col] = cosineSim;

        	if (cosineSim > bestCorrelation) {
				*bestAnchorColumnIndex = col;
        		bestCorrelation = cosineSim;
        	}
        }

    	*bestAnchorRowIndex = MathUtils::closest(
    		EigenUtils::convertEigenVectorToQVector(integrationVecSeg),
    		integrationVecSeg.maxCoeff()
    		);

        ERR_RETURN
    }

    Err calculatePeakCorrelations(
		const Eigen::VectorX<float> &integrationVecSeg,
        const Eigen::MatrixX<float> &matBlockTrimmed,
        QVector<float> *peakCorrelations
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matBlockTrimmed.size() > 0); ree;

        peakCorrelations->clear();

        peakCorrelations->reserve(static_cast<int>(matBlockTrimmed.cols()));
        for (int col = 0; col < matBlockTrimmed.cols(); col++) {

            const Eigen::VectorX<float> &matCol = matBlockTrimmed.col(col);

            float cosineSim;
            e = EigenUtils::cosineSimilarity(matCol, integrationVecSeg, &cosineSim); ree;
            peakCorrelations->push_back(cosineSim);
        }

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::processIntegrationVectorPeakIntegrations(
    const MatriciesAndVecs &matriciesAndVecs,
    const QVector<QPair<PeakIntegrationIndexes, Intensity>> &peakIntegrationsVsIntensity,
    QVector<BestCorrelationResult> *bestCorrelationResults
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;
    e = ErrorUtils::isTrue(matriciesAndVecs.integrationVecIsValid()); ree;

    const int maxRows = static_cast<int>(matriciesAndVecs.intensityMatrix100.rows());
    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensityResized = peakIntegrationsVsIntensity;
    // if (m_useTopNIntegrationsParam) {
    //     peakIntegrationsVsIntensityResized.resize(std::min(
    //         m_pythiaParameters.topNIntegrations,
    //         peakIntegrationsVsIntensityResized.size()
    //         ));
    // }

    for (const QPair<PeakIntegrationIndexes, Intensity> &pii : peakIntegrationsVsIntensityResized) {

        const QPair<PeakIntegrationIndexes, Intensity> piiWorking = correctPeakIntegrationForSingleRow(
            pii,
            matriciesAndVecs.ionCountVec
            );

        const int ogPeakLength = piiWorking.first.second - piiWorking.first.first + 1;

        Eigen::MatrixX<float> matBlock = matriciesAndVecs.intensityMatrix100.block(
              piiWorking.first.first,
              0,
              ogPeakLength,
              matriciesAndVecs.intensityMatrix100.cols()
              ).eval();

        const QVector<QVector<int>> apexIndexesByColumn = getMatrxColumnApexes(matBlock);

        const Eigen::VectorX<float> integrationVecSegment = matriciesAndVecs.productVec.segment(
            piiWorking.first.first,
            ogPeakLength
            ).eval();

        const QPair<int, float> apexIndex = EigenUtils::returnTopIndexAndValue(integrationVecSegment);

        const QVector<int> apexStarts = findStartApexes(apexIndexesByColumn, apexIndex.first);
        e = ErrorUtils::isEqual(apexStarts.size(), apexIndexesByColumn.size()); ree;

        constexpr float stopThresholdFraction = 0.0;
        // const Eigen::MatrixX<float> matBlockTrimmed = trimMatrixBlock(
        //     matBlock,
        //     apexStarts,
        //     stopThresholdFraction
        //     );

        QVector<float> peakCorrelations;
        int bestAnchorColumnIndex = -1;
        int bestAnchorRowIndex = -1;
        e = findBestAnchorColumn(
			integrationVecSegment,
            matBlock,
            m_pythiaParameters.maxAnchorColumnIndex,
            &peakCorrelations,
            &bestAnchorColumnIndex,
            &bestAnchorRowIndex
            ); ree;

        const int topN = std::min(6, peakCorrelations.size());
        const float peakCorrelationsSum
                = std::accumulate(peakCorrelations.begin(), peakCorrelations.begin() + topN, 0.0f);

        BestCorrelationResult bestCorrelationResultPii;

        bestCorrelationResultPii.peakCorrelations = peakCorrelations;
        bestCorrelationResultPii.peakCorrelationsSum = peakCorrelationsSum;
        bestCorrelationResultPii.matBlockTrimmedIntensity = matBlock;
        bestCorrelationResultPii.bestAnchorColumnIndex = bestAnchorColumnIndex;
        bestCorrelationResultPii.bestAnchorRowIndex = bestAnchorRowIndex;
        bestCorrelationResultPii.apexStarts = apexStarts;
        bestCorrelationResultPii.peakIntegrationIndexes = piiWorking.first;
    	bestCorrelationResultPii.integrationVector = integrationVecSegment;

        constexpr float windowMultiplier1p5X = 0.5;
        const auto frameIndex1p5XMin
            = std::max(static_cast<FrameIndex>(std::round(piiWorking.first.first - (windowMultiplier1p5X * ogPeakLength))), 0);
        const auto frameIndex1p5XMax
            = std::min(static_cast<FrameIndex>(std::round(piiWorking.first.first + (windowMultiplier1p5X * ogPeakLength))), maxRows);
        const int peakLength1p5X = frameIndex1p5XMax - frameIndex1p5XMin;
        Eigen::MatrixX<float> matBlock1p5X = matriciesAndVecs.intensityMatrix100.block(
              frameIndex1p5XMin,
              0,
              peakLength1p5X,
              matriciesAndVecs.intensityMatrix100.cols()
              ).eval();

        constexpr float windowMultiplier2X = 1.0;
        const auto frameIndex2XMin
            = std::max(static_cast<FrameIndex>(std::round(piiWorking.first.first - (windowMultiplier2X * ogPeakLength))), 0);
        const auto frameIndex2XMax
            = std::min(static_cast<FrameIndex>(std::round(piiWorking.first.first + (windowMultiplier2X * ogPeakLength))), maxRows);
        const int peakLength2X = frameIndex2XMax - frameIndex2XMin;
        Eigen::MatrixX<float> matBlock2X = matriciesAndVecs.intensityMatrix100.block(
                      frameIndex2XMin,
                      0,
                      peakLength2X,
                      matriciesAndVecs.intensityMatrix100.cols()
                      ).eval();

        // bestCorrelationResultPii.matBlockTrimmedIntensityWindow1p5X = trimMatrixBlock(
        //     matBlock1p5X,
        //     apexStarts,
        //     stopThresholdFraction
        //     );

    	bestCorrelationResultPii.matBlockTrimmedIntensityWindow1p5X = matBlock1p5X;
    	const Eigen::VectorX<float> integrationVecSegment1p5X = matriciesAndVecs.productVec.segment(
			frameIndex1p5XMin,
			peakLength1p5X
			).eval();

        e = calculatePeakCorrelations(
			integrationVecSegment1p5X,
            bestCorrelationResultPii.matBlockTrimmedIntensityWindow1p5X,
            &bestCorrelationResultPii.peakCorrelationsWindow1p5X
            ); ree;

    	bestCorrelationResultPii.matBlockTrimmedIntensityWindow2X = matBlock2X;
		const Eigen::VectorX<float> integrationVecSegment2X = matriciesAndVecs.productVec.segment(
			frameIndex2XMin,
			peakLength2X
			).eval();

        // bestCorrelationResultPii.matBlockTrimmedIntensityWindow2X = trimMatrixBlock(
        //             matBlock2X,
        //             apexStarts,
        //             stopThresholdFraction
        //             );

        e = calculatePeakCorrelations(
			integrationVecSegment2X,
            bestCorrelationResultPii.matBlockTrimmedIntensityWindow2X,
            &bestCorrelationResultPii.peakCorrelationsWindow2X
            ); ree;

        const PeakIntegrationIndexes &p = piiWorking.first;
        const int pSize = p.second - p.first + 1;
        bestCorrelationResultPii.matBlockTrimmedMz = matriciesAndVecs.mzMatrix100.block(
            p.first,
            0,
            pSize,
            matBlock.cols()
            ).eval();

        bestCorrelationResultPii.matBlockTrimmedIntensityShadows = matriciesAndVecs.intensityMatrix100Shadow.block(
            p.first,
            0,
            pSize,
            matBlock.cols()
            ).eval();

        bestCorrelationResultPii.matBlockTrimmedIntensity45 = matriciesAndVecs.intensityMatrix45.block(
            p.first,
            0,
            pSize,
            matBlock.cols()
            ).eval();

        e = calculatePeakCorrelations(
			integrationVecSegment,
            bestCorrelationResultPii.matBlockTrimmedIntensity45,
            &bestCorrelationResultPii.peakCorrelations45
            ); ree;

        bestCorrelationResults->push_back(bestCorrelationResultPii);

// #define OUTPUT_MATS
#ifdef OUTPUT_MATS
        qDebug() << peakCorrelations << std::accumulate(peakCorrelations.begin(), peakCorrelations.end(), 0.0f) << "SLDFJSDLFLJ";
        for (int as : apexStarts) {
            std::cout << as << " " << apexIndex.first << std::endl;
        }
        std::cout << piiWorking.first.first << " " << piiWorking.first.second << std::endl;
        std::cout << matBlock << std::endl;
        std::cout << bestAnchorColumnIndex << " **********" << std::endl;
#endif
    }

    ERR_RETURN
}

namespace {

    float calculatedCosineSimSumGreaterThan80(const QVector<float> &peakCorrelations) {

        const int top6 = std::min(6, peakCorrelations.size());
        constexpr float cosineSimToAnchorThreshold = 0.8;
        const float scoreGreater = std::accumulate(
            peakCorrelations.begin(),
            peakCorrelations.begin() + top6,
            0.0f,
            [cosineSimToAnchorThreshold](float sum, float f){return f > cosineSimToAnchorThreshold ? f + sum : 0.0f + sum;}
            );

        return scoreGreater;
    }

    Err setCosineSimilarityMetrics(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const BestCorrelationResult& bestCorrelationResult,
        int topNMS2Ions,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        //TODO add kldiv and pearsons corr.

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.size() > 0); ree;
        QVector<MS2Ion> ms2IonsTheo = candidateScores->isDecoy
                                    ? targetDecoyCandidatePair->ms2IonsDecoy()
                                    : targetDecoyCandidatePair->ms2IonsTarget();

    	ms2IonsTheo.resize(std::min(ms2IonsTheo.size(), topNMS2Ions));

        QVector<float> intensitiesTheo;
        std::transform(
            ms2IonsTheo.begin(),
            ms2IonsTheo.end(),
            std::back_inserter(intensitiesTheo),
            [](const MS2Ion &ms2Ion){return ms2Ion.intensity;}
            );

        const Eigen::VectorX<float> intensitiesTheoVec = EigenUtils::convertQVectorToEigenVector(intensitiesTheo);

        Eigen::MatrixX<float> intensitiesTheoMat(
            bestCorrelationResult.matBlockTrimmedIntensity.rows(),
            intensitiesTheoVec.size()
        );
        for (int row = 0; row < intensitiesTheoMat.rows(); row++) {
            intensitiesTheoMat.row(row) = intensitiesTheoVec;
        }

    	Eigen::MatrixX<float> matBlockTrimmedIntensityResized(
			bestCorrelationResult.matBlockTrimmedIntensity.rows(),
			intensitiesTheoVec.size()
		);
		for (int col = 0; col < bestCorrelationResult.matBlockTrimmedIntensity.cols(); col++) {
			matBlockTrimmedIntensityResized.col(col) = bestCorrelationResult.matBlockTrimmedIntensity.col(col);
		}

    	e = ErrorUtils::isEqual(
    		intensitiesTheoMat.cols(),
    		bestCorrelationResult.matBlockTrimmedIntensity.cols()
    		); ree;

        Eigen::VectorX<float> cosineSimsByRow;
         e = EigenUtils::rowWiseCosineSimilarOfMatrices(
             matBlockTrimmedIntensityResized,
             intensitiesTheoMat,
             &cosineSimsByRow
             ); ree;

        Eigen::MatrixX<float> matSummed = bestCorrelationResult.matBlockTrimmedIntensity;
        for (int row = 0; row < matSummed.rows(); row++) {
            Eigen::VectorX<float> r = matSummed.row(row).array() / std::max(matSummed.row(row).sum(), 0.001f);
            matSummed.row(row) = r;
        }

        for (int row = 0; row < intensitiesTheoMat.rows(); row++) {
            Eigen::VectorX<float> r = intensitiesTheoMat.row(row).array() / std::max(intensitiesTheoMat.row(row).sum(), 0.001f);
            intensitiesTheoMat.row(row) = r;
        }

        Eigen::VectorX<float> klDivByRow(intensitiesTheoMat.rows());
        klDivByRow.setZero();
        for (int row = 0; row < intensitiesTheoMat.rows(); row++) {
            float klDiv;
            e = EigenUtils::klDivergence(
                Eigen::VectorX<float>(matSummed.row(row)),
                Eigen::VectorX<float>(intensitiesTheoMat.row(row)),
                &klDiv);
            klDivByRow.coeffRef(row) = klDiv;
        }

        for (int i = 0; i < bestCorrelationResult.peakCorrelations.size(); i++) {
            candidateScores->featuresArray[Features::CosineSimToAnchor1 + i] = bestCorrelationResult.peakCorrelations.at(i);
        }

        const float cosineSimMax = cosineSimsByRow.coeff(bestCorrelationResult.bestAnchorRowIndex);
        candidateScores->featuresArray[CosineSimSpectrum] = cosineSimMax;
        candidateScores->featuresArray[CosineSimSpectrumCubed]
                                                    = static_cast<float>(std::pow(cosineSimMax, 3));

        const float klDivMax = klDivByRow.coeff(bestCorrelationResult.bestAnchorRowIndex);
        candidateScores->featuresArray[KlDivSpectrum] = klDivMax;
        candidateScores->featuresArray[KlDivSpectrumCubeRoot]
                                                    = static_cast<float>(std::pow(cosineSimMax, 1.0f/3.0f));

        const float cosineSimRowsSummed = cosineSimsByRow.sum();
        if (MathUtils::tZero(cosineSimRowsSummed)) {
            candidateScores->featuresArray[CosineSimSpectrumOverTime] = 0.0f;
            candidateScores->featuresArray[CosineSimSpectrumOverTimeCubed] = 0.0f;

            ERR_RETURN
        }

        const int nonZeroRows = EigenUtils::nonZeros(cosineSimsByRow);
        const float cosineSimSpectrumOverTime = cosineSimRowsSummed / static_cast<float>(nonZeroRows);

        candidateScores->featuresArray[CosineSimSpectrumOverTime] = cosineSimSpectrumOverTime;
        candidateScores->featuresArray[CosineSimSpectrumOverTimeCubed]
                                                                = static_cast<float>(std::pow(cosineSimSpectrumOverTime, 3));

        QVector<float> cosineSimsByRowSansZeros;
        cosineSimsByRowSansZeros.reserve(nonZeroRows);
        for (int i = 0; i < cosineSimsByRow.size(); i++) {
            const float f = cosineSimsByRow.coeff(i);
            if (MathUtils::tZero(f)) {
                continue;
            }
            cosineSimsByRowSansZeros.push_back(f);
        }
        candidateScores->featuresArray[CosineSimSpectrumStDev]
                                                        = static_cast<float>(MathUtils::stDev(cosineSimsByRowSansZeros));

        const QVector<float> &cosineSimToAnchorVec = bestCorrelationResult.peakCorrelations;

        const int top6 = std::min(6, cosineSimToAnchorVec.size());
        const float cosineSimSumTop = std::accumulate(
            cosineSimToAnchorVec.begin(),
            cosineSimToAnchorVec.begin() + top6,
            0.0f
    );
        candidateScores->featuresArray[CosineSimSumTop] = cosineSimSumTop;

        float cosineSimSumBottom = 0.1;
        if (cosineSimToAnchorVec.size() > top6) {
            cosineSimSumBottom = std::accumulate(
                    cosineSimToAnchorVec.begin() + top6 + 1,
                    cosineSimToAnchorVec.end(),
                    std::numeric_limits<float>::min()
            );
        }

        candidateScores->featuresArray[CosineSimSumBottom] = cosineSimSumBottom;

        candidateScores->featuresArray[TopBottomRatio]
                = std::log(std::max(1.0f, cosineSimSumTop) / (cosineSimSumTop + cosineSimSumBottom + 1.0f));

        candidateScores->featuresArray[TopBottomRatioNorm]
                = cosineSimSumBottom / static_cast<float>(candidateScores->targetDecoyCandidatePair->totalFragmentCount());

        ERR_RETURN
    }

    constexpr int arraySizeMax = 12;

    Err buildScanTimeVectorFromFrameIndexIntegrationLimits(
        const PeakIntegrationIndexes &peakIntegrationIndexes,
        MsFrame *msFrameMzTarget,
        QVector<ScanTime> *scanTimes
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msFrameMzTarget->isValid()); ree;
        e = ErrorUtils::isTrue(peakIntegrationIndexes.second > peakIntegrationIndexes.first); ree;

        scanTimes->reserve(peakIntegrationIndexes.second - peakIntegrationIndexes.first + 1);
        for (int frameIndex = peakIntegrationIndexes.first; frameIndex <= peakIntegrationIndexes.second; ++frameIndex) {
            scanTimes->push_back(msFrameMzTarget->scanTimeFromFrameIndex(frameIndex));
        }

        ERR_RETURN
    }

    Err calculateTrapezoidalArea(
        const Eigen::MatrixX<float> &xicMat,
        const QVector<float> &scanTimes,
        Eigen::VectorX<float> *areas
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scanTimes); ree;
        e = ErrorUtils::isEqual(scanTimes.size(), static_cast<int>(xicMat.rows())); ree;

        const Eigen::VectorX<float> scanTimesVec = EigenUtils::convertQVectorToEigenVector(scanTimes);

        const int vSize = scanTimes.size() - 1;
        const Eigen::VectorX<float> dScanTimes = scanTimesVec.tail(vSize) - scanTimesVec.head(vSize);
        const Eigen::MatrixX<float> xicMatSums = xicMat.bottomRows(vSize) + xicMat.topRows(vSize);

        *areas = (dScanTimes.transpose() * xicMatSums) * 0.5;

        ERR_RETURN
    }

    Err setFoundMs2Ions(
        const QVector<BestCorrelationResult> &bestCorrelationResults,
        int topNMS2Ions,
        MsFrame *msFrameMzTarget,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        const BestCorrelationResult &bestCorrelationResult = bestCorrelationResults.front();

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedMz.size() > 0); ree;

        const Eigen::MatrixX<float> &matMz = bestCorrelationResult.matBlockTrimmedMz;

        QVector<float> mzMeanValsFound;
        mzMeanValsFound.reserve(static_cast<int>(matMz.cols()));
        QVector<float> mzMeanValsFoundPPM;
        mzMeanValsFoundPPM.reserve(static_cast<int>(matMz.cols()));
        QVector<float> stdMeanValsFound;
        stdMeanValsFound.reserve(static_cast<int>(matMz.cols()));

        for (int col = 0; col < matMz.cols(); col++) {

            const Eigen::VectorX<float> &colMz = matMz.col(col);

            if(MathUtils::tZero(colMz.maxCoeff())) {
                mzMeanValsFound.push_back(-1.0f);
                stdMeanValsFound.push_back(-1.0f);
                continue;
            }

            const Eigen::VectorX<float> &colMzNonZero = EigenUtils::removeZeroElements(colMz);

            mzMeanValsFound.push_back(colMzNonZero.mean());
            stdMeanValsFound.push_back(static_cast<float>(EigenUtils::calculateStDevOfVector(colMzNonZero)));
        }

        const QVector<MS2Ion> &ms2Ions = candidateScores->isDecoy
                                       ? candidateScores->targetDecoyCandidatePair->ms2IonsDecoy()
                                       : candidateScores->targetDecoyCandidatePair->ms2IonsTarget();

        int foundB = 0;
        int foundY = 0;
        for (int i = 0; i < std::min(mzMeanValsFound.size(), arraySizeMax); i++) {

        	const MS2Ion &ms2Ion = ms2Ions.at(i);

            candidateScores->featuresArray[MzFoundMean1 + i] = mzMeanValsFound.at(i);
            const float ppm = mzMeanValsFound.at(i) > 1.0f
                            ? std::min(MathUtils::calculateMassAccuracyPPM(ms2Ion.mz, mzMeanValsFound.at(i)), 100.0f)
                            : 100.0f;
            // candidateScores->featuresArray[MzFoundMean1PPM + i] = ppm;

            if (ppm > 99.9) {
                continue;
            }
            mzMeanValsFoundPPM.push_back(ppm);
            if (ms2Ion.ionLabel.contains('y')) {
                foundY++;
            }
            if (ms2Ion.ionLabel.contains('b')) {
                foundB++;
            }

	       	const float ms2IonMass = static_cast<float>(ms2Ion.mz * ms2Ion.charge) - (ms2Ion.charge * 1.0072);
        	if (ms2IonMass > 650) {
        		candidateScores->featuresArray[MzFoundOverCount650]++;
        	}
        	else {
        		candidateScores->featuresArray[MzFoundUnderCount650]++;
        	}
        }

        candidateScores->featuresArray[MzPPMMean] = MathUtils::mean(mzMeanValsFoundPPM);
        candidateScores->featuresArray[MzPPMMeanAbs] = std::abs(MathUtils::mean(mzMeanValsFoundPPM));
        candidateScores->featuresArray[MzPPMStd] = MathUtils::stDev(mzMeanValsFoundPPM);
        candidateScores->featuresArray[MzPPMStd] = std::isinf(candidateScores->featuresArray[MzPPMStd]) || std::isnan(candidateScores->featuresArray[MzPPMStd])
                                                                ? -1.0f
                                                                : candidateScores->featuresArray[MzPPMStd];

        candidateScores->featuresArray[FoundB] = foundB / static_cast<float>(topNMS2Ions);
        candidateScores->featuresArray[FoundY] = foundY / static_cast<float>(topNMS2Ions);
        candidateScores->featuresArray[FoundPercent] = (foundB + foundY) / static_cast<float>(topNMS2Ions);


        for (int i = 0; i < std::min(stdMeanValsFound.size(), arraySizeMax / 2); i++) {
            candidateScores->featuresArray[MzFoundStDev1 + i] = stdMeanValsFound.at(i);
        }

        QVector<ScanTime> scanTimes;
        e = buildScanTimeVectorFromFrameIndexIntegrationLimits(
            bestCorrelationResult.peakIntegrationIndexes,
            msFrameMzTarget,
            &scanTimes
            );

        const Eigen::VectorX<float> matBlockTrimmedIntensityIntegrationSums = bestCorrelationResult.matBlockTrimmedIntensity.colwise().sum();

        const Eigen::VectorX<float> intensitySumsNormalized
                = matBlockTrimmedIntensityIntegrationSums / std::max(static_cast<float>(bestCorrelationResult.matBlockTrimmedIntensity.maxCoeff()), 1.0f);

        for (int i = 0; i < std::min(static_cast<int>(matBlockTrimmedIntensityIntegrationSums.size()), arraySizeMax); i++) {
            candidateScores->featuresArray[Features::IntensityFoundMax1 + i] = matBlockTrimmedIntensityIntegrationSums.coeff(i);
        }

        for (int i = 0; i < std::min(static_cast<int>(intensitySumsNormalized.size()), arraySizeMax); i++) {
            candidateScores->featuresArray[IntensityFoundMaxNorm1 + i] = intensitySumsNormalized.coeff(i);
        }

        ERR_RETURN
    }

    Err setPeakShapeRatios(
        const BestCorrelationResult &bestCorrelationResult,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.size() > 0);

        constexpr double chunkDivision = 3.0;

        const int bestColumnIndex = bestCorrelationResult.bestAnchorColumnIndex;
        const Eigen::VectorX<float> &bestAnchorColumn = bestCorrelationResult.integrationVector;

        QVector<float> bestAnchorColumnVec = EigenUtils::convertEigenVectorToQVector(bestAnchorColumn);

        const auto terminatorLogic = [](double d){return d < 1.0;};
        const auto terminator = std::remove_if(bestAnchorColumnVec.begin(), bestAnchorColumnVec.end(), terminatorLogic);
        bestAnchorColumnVec.erase(terminator, bestAnchorColumnVec.end());

        const int chunkSize = std::max(1, static_cast<int>(std::round(bestAnchorColumnVec.size() / chunkDivision)));
        const double bestAnchorColumnVecSum = std::accumulate(bestAnchorColumnVec.begin(), bestAnchorColumnVec.end(), 0.0001);

        if (bestAnchorColumnVec.size() < chunkDivision) {
            candidateScores->featuresArray[PeakShapeRatio1] = std::numeric_limits<float>::min();
            candidateScores->featuresArray[PeakShapeRatio2] = 1.0;
            candidateScores->featuresArray[PeakShapeRatio3] = std::numeric_limits<float>::min();
        }
        else {

            candidateScores->featuresArray[PeakShapeRatio1] = std::accumulate(
                    bestAnchorColumnVec.begin(),
                    bestAnchorColumnVec.begin() + chunkSize,
                    std::numeric_limits<float>::min()
            ) / bestAnchorColumnVecSum;

            candidateScores->featuresArray[PeakShapeRatio2] = std::accumulate(
                    bestAnchorColumnVec.begin() + chunkSize,
                    bestAnchorColumnVec.begin() + (chunkSize * 2),
                    std::numeric_limits<float>::min()
            ) / bestAnchorColumnVecSum;

            candidateScores->featuresArray[PeakShapeRatio3] = std::accumulate(
                    bestAnchorColumnVec.begin() + (chunkSize * 2),
                    bestAnchorColumnVec.end(),
                    std::numeric_limits<float>::min()
            ) / bestAnchorColumnVecSum;

        }

        ERR_RETURN
    }

    Err setAminoAcidFrequencies(CandidateScores *candidateScores) {

        ERR_INIT

        QMap<QChar, int> aminoAcidCounts = {
            {'A', 0},
            {'C', 0},
            {'D', 0},
            {'E', 0},
            {'F', 0},
            {'G', 0},
            {'H', 0},
            {'I', 0},
            {'K', 0},
            {'L', 0},
            {'M', 0},
            {'N', 0},
            {'P', 0},
            {'Q', 0},
            {'R', 0},
            {'S', 0},
            {'T', 0},
            {'V', 0},
            {'W', 0},
            {'Y', 0},
            {'B', 0},
            {'J', 0},
            {'O', 0},
            {'U', 0},
            {'X', 0},
            {'Z', 0}
        };

        const QString peptideString = candidateScores->targetDecoyCandidatePair->peptideString();

        for (const QChar aminoAcid : peptideString) {

            if (!aminoAcidCounts.contains(aminoAcid)) {
                qDebug() << peptideString << "missing amino acid" << aminoAcid;
                continue;
            }

            e = ErrorUtils::isTrue(aminoAcidCounts.contains(aminoAcid)); ree;
            aminoAcidCounts[aminoAcid]++;
        }

        candidateScores->featuresArray[AminoAcidCountA] = static_cast<float>(aminoAcidCounts['A']);
        candidateScores->featuresArray[AminoAcidCountC] = static_cast<float>(aminoAcidCounts['C']);
        candidateScores->featuresArray[AminoAcidCountD] = static_cast<float>(aminoAcidCounts['D']);
        candidateScores->featuresArray[AminoAcidCountE] = static_cast<float>(aminoAcidCounts['E']);
        candidateScores->featuresArray[AminoAcidCountF] = static_cast<float>(aminoAcidCounts['F']);
        candidateScores->featuresArray[AminoAcidCountG] = static_cast<float>(aminoAcidCounts['G']);
        candidateScores->featuresArray[AminoAcidCountH] = static_cast<float>(aminoAcidCounts['H']);
        candidateScores->featuresArray[AminoAcidCountI] = static_cast<float>(aminoAcidCounts['I']);
        candidateScores->featuresArray[AminoAcidCountK] = static_cast<float>(aminoAcidCounts['K']);
        candidateScores->featuresArray[AminoAcidCountL] = static_cast<float>(aminoAcidCounts['L']);
        candidateScores->featuresArray[AminoAcidCountM] = static_cast<float>(aminoAcidCounts['M']);
        candidateScores->featuresArray[AminoAcidCountN] = static_cast<float>(aminoAcidCounts['N']);
        candidateScores->featuresArray[AminoAcidCountP] = static_cast<float>(aminoAcidCounts['P']);
        candidateScores->featuresArray[AminoAcidCountQ] = static_cast<float>(aminoAcidCounts['Q']);
        candidateScores->featuresArray[AminoAcidCountR] = static_cast<float>(aminoAcidCounts['R']);
        candidateScores->featuresArray[AminoAcidCountS] = static_cast<float>(aminoAcidCounts['S']);
        candidateScores->featuresArray[AminoAcidCountT] = static_cast<float>(aminoAcidCounts['T']);
        candidateScores->featuresArray[AminoAcidCountV] = static_cast<float>(aminoAcidCounts['V']);
        candidateScores->featuresArray[AminoAcidCountW] = static_cast<float>(aminoAcidCounts['W']);
        candidateScores->featuresArray[AminoAcidCountY] = static_cast<float>(aminoAcidCounts['Y']);

        candidateScores->featuresArray[AminoAcidCountB] = static_cast<float>(aminoAcidCounts['B']);
        candidateScores->featuresArray[AminoAcidCountJ] = static_cast<float>(aminoAcidCounts['J']);
        candidateScores->featuresArray[AminoAcidCountO] = static_cast<float>(aminoAcidCounts['O']);
        candidateScores->featuresArray[AminoAcidCountU] = static_cast<float>(aminoAcidCounts['U']);
        candidateScores->featuresArray[AminoAcidCountX] = static_cast<float>(aminoAcidCounts['X']);
        candidateScores->featuresArray[AminoAcidCountZ] = static_cast<float>(aminoAcidCounts['Z']);

        ERR_RETURN
    }

    Err setShadowCorrelations(
        const BestCorrelationResult &bestCorrelationResult,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

    	if (!bestCorrelationResult.matBlockTrimmedIntensity.sum() > 0) {
    		qDebug()
    		<< "SLDFJSDLKJ"
    		<< bestCorrelationResult.matBlockTrimmedIntensity.rows()
    		<< bestCorrelationResult.matBlockTrimmedIntensity.cols()
    		<< bestCorrelationResult.peakCorrelationsSum
    		<< bestCorrelationResult.peakIntegrationIndexes;
    	}

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.sum() > 0); ree;

        for (int col = 0; col < bestCorrelationResult.matBlockTrimmedIntensity.cols(); col++) {

            const Eigen::VectorX<float> &v1 = bestCorrelationResult.matBlockTrimmedIntensity.col(col);
            const Eigen::VectorX<float> &v2 = bestCorrelationResult.matBlockTrimmedIntensityShadows.col(col);
            const float v1Max = v1.maxCoeff();
            const float v2Max = v2.maxCoeff();

            float cosineSim;
            e = EigenUtils::cosineSimilarity(v1, v2, &cosineSim); ree;
            if (col < 6 && cosineSim > 0.80) {
                candidateScores->featuresArray[ShadowsCosineSimSum] += cosineSim;
            }
            // candidateScores->featuresArray[CosineSimShadowsToAnchor1 + col] = std::max(cosineSim, 0.0f);

        }

        ERR_RETURN
    }

    Err setMs1Averagine(
        const QVector<float> &ms1Averagine,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        const Eigen::VectorX<float> averagineVec = EigenUtils::convertQVectorToEigenVector(ms1Averagine);

        const QVector<float> ms1IsoDisActualVec = {
            candidateScores->featuresArray[CosineSim100MS1PreMono],
            candidateScores->featuresArray[CosineSim100MS1],
            candidateScores->featuresArray[CosineSim100MS1Iso1],
            candidateScores->featuresArray[CosineSim100MS1Iso2]
            };
        const Eigen::VectorX<float> ms1IsoDistActual = EigenUtils::convertQVectorToEigenVector(ms1IsoDisActualVec);

        float cosineSimAveragine;
        e = EigenUtils::cosineSimilarity(
            averagineVec,
            ms1IsoDistActual,
            &cosineSimAveragine
            ); ree;

        candidateScores->featuresArray[MS1Averagine] = cosineSimAveragine;

        ERR_RETURN
    }

    Err setIntegrations(
        const BestCorrelationResult &bestCorrelationResult,
        int peakCenter,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        constexpr int top12 = 12;

        if (peakCenter > 0) {
            const int startScan = std::max(bestCorrelationResult.bestAnchorRowIndex - static_cast<int>(std::round(peakCenter / 2.0)), 0);
            const int distance = startScan + peakCenter > bestCorrelationResult.matBlockTrimmedIntensity.rows()
                               ? static_cast<int>(bestCorrelationResult.matBlockTrimmedIntensity.rows()) - startScan
                               : peakCenter;

            const Eigen::MatrixX<float> matBlockTrimmedIntensityIntegration
                                                      = peakCenter < 0 || bestCorrelationResult.matBlockTrimmedIntensity.rows() <= peakCenter
                                                      ? bestCorrelationResult.matBlockTrimmedIntensity
                                                      : bestCorrelationResult.matBlockTrimmedIntensity.block(
                                                          startScan,
                                                          0,
                                                          distance,
                                                          bestCorrelationResult.matBlockTrimmedIntensity.cols()
                                                          ).eval();

            const Eigen::VectorX<float> matBlockTrimmedIntensityIntegrationSum = matBlockTrimmedIntensityIntegration.colwise().sum();
            QVector<float> vec = EigenUtils::convertEigenVectorToQVector(matBlockTrimmedIntensityIntegrationSum);

            vec.resize(top12);

            candidateScores->integrations = vec;

            ERR_RETURN
        }

        candidateScores->integrations = QVector<float>(top12, 0.0f);
        for (int i = 0; i < top12; i++) {
            candidateScores->integrations[i] = candidateScores->featuresArray[Features::IntensityFoundMax1 + i];
        }

        ERR_RETURN
    }

    Err setNormPeakLengthScores(
        const BestCorrelationResult &bestCorrelationResult,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        e = ErrorUtils::isAboveThreshold(
            bestCorrelationResult.bestAnchorColumnIndex,
            0,
            ErrorUtilsParam::IncludeThreshold
            ); ree;

        const int columnCount = bestCorrelationResult.matBlockTrimmedIntensity.cols();
        QVector<int> peakLengths(columnCount, 0);

        for (int col = 0; col < columnCount; col++) {
            peakLengths[col] = EigenUtils::nonZeros(bestCorrelationResult.matBlockTrimmedIntensity.col(col));
        }

        const int anchorPeakLength = peakLengths[bestCorrelationResult.bestAnchorColumnIndex];
        e = ErrorUtils::isFalse(anchorPeakLength <= 0); ree;

        for (int i = 0; i < peakLengths.size(); i++) {
            candidateScores->featuresArray[Features::MzPeakLengthsNorm1 + i] = peakLengths[i] / static_cast<float>(anchorPeakLength);
        }

        candidateScores->featuresArray[MzPeakLengthsMean] = MathUtils::mean(peakLengths);
        candidateScores->featuresArray[MzPeakLengthsStd] = MathUtils::stDev(peakLengths);

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::setCandidateScores(
    const TargetDecoyCandidatePair *targetDecoyCandidatePair,
    const QVector<MS2Ion> &ms2Ions,
    const QVector<BestCorrelationResult> &bestCorrelationResults,
    const QVector<float> &ms1Averagine,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    const BestCorrelationResult &bestCorrelationResult = bestCorrelationResults.front();

    e = ErrorUtils::isNotEmpty(bestCorrelationResult.peakCorrelations); ree;
    e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.size() > 0); ree;

    candidateScores->initFeaturesArray();
    candidateScores->proteinGroup = targetDecoyCandidatePair->proteinGroups();

    candidateScores->frameIndex = bestCorrelationResult.peakIntegrationIndexes.first
                                + bestCorrelationResult.apexStarts.at(bestCorrelationResult.bestAnchorColumnIndex);

    candidateScores->frameIndexStart = bestCorrelationResult.peakIntegrationIndexes.first;
    candidateScores->frameIndexEnd = bestCorrelationResult.peakIntegrationIndexes.second;

    candidateScores->scanNumber = m_msFrameMzTarget->scanNumberFromFrameIndex(candidateScores->frameIndex);
    candidateScores->scanNumberStart = m_msFrameMzTarget->scanNumberFromFrameIndex(candidateScores->frameIndexStart);
    candidateScores->scanNumberEnd = m_msFrameMzTarget->scanNumberFromFrameIndex(candidateScores->frameIndexEnd);

    candidateScores->scanTime = m_msFrameMzTarget->scanTimeFromScanNumber(candidateScores->scanNumber);

	if (m_msCalibratomatic.isInitRT()) {
		e = m_msCalibratomatic.predictIRT(
				candidateScores->scanTime,
				&candidateScores->empiricalIRT
				); ree;
	}

    candidateScores->scanTimeStart = m_msFrameMzTarget->scanTimeFromScanNumber(candidateScores->scanNumberStart);
    candidateScores->scanTimeEnd = m_msFrameMzTarget->scanTimeFromScanNumber(candidateScores->scanNumberEnd);

    const int top6 = std::min(6, bestCorrelationResult.peakCorrelations.size());
    candidateScores->featuresArray[CosineSimSum100] = std::accumulate(
            bestCorrelationResult.peakCorrelations.begin(),
            bestCorrelationResult.peakCorrelations.begin() + top6,
            std::numeric_limits<float>::min()
            );

    candidateScores->featuresArray[CosineSimSum100Top12] = std::accumulate(
            bestCorrelationResult.peakCorrelations.begin(),
            bestCorrelationResult.peakCorrelations.end(),
            std::numeric_limits<float>::min()
            );

    candidateScores->featuresArray[CosineSimSum100Window1p5X] = std::accumulate(
        bestCorrelationResult.peakCorrelationsWindow1p5X.begin(),
        bestCorrelationResult.peakCorrelationsWindow1p5X.begin() + top6,
        std::numeric_limits<float>::min()
        );

    candidateScores->featuresArray[CosineSimSum100Window2X] = std::accumulate(
        bestCorrelationResult.peakCorrelationsWindow2X.begin(),
        bestCorrelationResult.peakCorrelationsWindow2X.begin() + top6,
        std::numeric_limits<float>::min()
        );

    candidateScores->featuresArray[ScanIonCount]
        = static_cast<float>(m_msFrameMzTarget->getScanPointsByScanNumber(candidateScores->scanNumber)->size());

    candidateScores->featuresArray[CosineSimSum100GreaterThan80]
                                            = calculatedCosineSimSumGreaterThan80(bestCorrelationResult.peakCorrelations);

    candidateScores->featuresArray[TheoFragmentCount]
                                                            = static_cast<float>(targetDecoyCandidatePair->totalFragmentCount());

    candidateScores->featuresArray[TotalIntensityLog]
                                        = std::log(std::max(bestCorrelationResult.matBlockTrimmedIntensity.sum(),
                                                    std::numeric_limits<float>::min()));

     candidateScores->featuresArray[TotalIntensityPeakHeights]
                                = bestCorrelationResult.matBlockTrimmedIntensity.colwise().maxCoeff().sum();

    candidateScores->featuresArray[TotalIntensityRaw]
                                            = std::exp(candidateScores->featuresArray[TotalIntensityLog]);

    candidateScores->featuresArray[Charge] = static_cast<float>(candidateScores->targetDecoyCandidatePair->charge());

    const float scanTimeDelta = candidateScores->scanTime - candidateScores->scanTimePredicted;
    candidateScores->featuresArray[ScanTimeDelta] = scanTimeDelta;
    candidateScores->featuresArray[ScanTimeDeltaAbs] = std::abs(scanTimeDelta);

    candidateScores->featuresArray[ScanTimePredicted] = candidateScores->scanTimePredicted;

    const double pdScanTime = std::sqrt(std::min(std::abs(scanTimeDelta), m_scanTimeRange) / m_scanTimeRange);
    candidateScores->featuresArray[ScanTimePdAbs] = static_cast<float>(pdScanTime);
    candidateScores->featuresArray[ScanTimePd] = scanTimeDelta < 0
                                                         ? -static_cast<float>(std::abs(pdScanTime))
                                                         : static_cast<float>(std::abs(pdScanTime));

    const double pepLength = (-10.0 + candidateScores->targetDecoyCandidatePair->peptideString().size()) / 10.0;
    candidateScores->featuresArray[PeptideLengthNorm] = static_cast<float>(pepLength);

    const auto mz = candidateScores->targetDecoyCandidatePair->mz(false);
    candidateScores->featuresArray[MzNorm] = (mz - 600.0f) * 0.002f;
    candidateScores->featuresArray[IRTPredicted] = candidateScores->targetDecoyCandidatePair->iRt(candidateScores->isDecoy);
    candidateScores->featuresArray[Mass] = candidateScores->targetDecoyCandidatePair->mass();

    const float mzTargetKey = MathUtils::unHashDecimal<float>(m_mzTargetKey.toInt(), S_GLOBAL_SETTINGS.HASHING_PRECISION);
    candidateScores->featuresArray[TargetWindowLocation] = mzTargetKey - mz;
    candidateScores->featuresArray[TargetWindowLocationAbs] = std::abs(mzTargetKey - mz);

    candidateScores->featuresArray[CosineSimSum45]
        = std::max(std::accumulate(bestCorrelationResult.peakCorrelations45.begin(), bestCorrelationResult.peakCorrelations45.begin() + top6, 0.0f), std::numeric_limits<float>::min());

    int bestAlignmentMatrixRowIndex = bestCorrelationResult.bestAnchorRowIndex;
    candidateScores->featuresArray[AllignedMaxIndexesCount] = static_cast<float>(std::count_if(
        bestCorrelationResult.apexStarts.begin(),
        bestCorrelationResult.apexStarts.end(),
        [bestAlignmentMatrixRowIndex](int i){return i == bestAlignmentMatrixRowIndex;}
        ));

    const float length = bestCorrelationResult.peakIntegrationIndexes.second - bestCorrelationResult.peakIntegrationIndexes.first;
    candidateScores->featuresArray[AlignmentIndexMean] = MathUtils::mean(bestCorrelationResult.apexStarts) / std::max(length, 1.0f);
    candidateScores->featuresArray[AlignmentIndexStDev] = MathUtils::stDev(bestCorrelationResult.apexStarts);
    candidateScores->featuresArray[AlignmentCombinedScore]
                        = candidateScores->featuresArray[AlignmentIndexMean]
                        * candidateScores->featuresArray[AllignedMaxIndexesCount];

    candidateScores->featuresArray[MatrixZeroPercentage] = static_cast<float>((bestCorrelationResult.matBlockTrimmedIntensity.array() < 1).count())
                                                                                    / bestCorrelationResult.matBlockTrimmedIntensity.size();
    e = setAminoAcidFrequencies(candidateScores); ree;

    e = setCosineSimilarityMetrics(
        targetDecoyCandidatePair,
        bestCorrelationResult,
        m_topNMS2Ions,
        candidateScores
        );ree

	if (m_turboXicMS1) {
		e = setMs1RelatedScores(
			targetDecoyCandidatePair,
			bestCorrelationResult,
			static_cast<float>(m_pythiaParameters.ms1ExtractionWidthPPM),
			candidateScores
			); ree;
	}

    e = setLibraryIonMobilityRelatedScores(
        targetDecoyCandidatePair,
        candidateScores
        ); ree;

    const bool needsMs2IonMobilityScores = containsMs2IonMobilityFeature(m_features);
    if (needsMs2IonMobilityScores) {
        e = setMs2IonMobilityRelatedScores(
            targetDecoyCandidatePair,
            ms2Ions,
            candidateScores
            ); ree;
    }

    e = setFoundMs2Ions(
        bestCorrelationResults,
        m_topNMS2Ions,
        m_msFrameMzTarget,
        candidateScores
        ); ree;

    e = setPeakShapeRatios(bestCorrelationResult, candidateScores); ree;

    e = setShadowCorrelations(bestCorrelationResult, candidateScores); ree;

    e = setMs1Averagine(ms1Averagine, candidateScores); ree;

    e = setIntegrations(
        bestCorrelationResult,
        m_pythiaParameters.peakCenter,
        candidateScores
        ); ree;

    e = setNormPeakLengthScores(bestCorrelationResult, candidateScores); ree;

    ERR_RETURN
}

namespace {

    Eigen::VectorX<float> buildXicPointsVector(
        const XICPoints &xicPoints,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax
        ) {

        Eigen::VectorX<float> vec(frameIndexMax + 1);
        vec.setZero();

        for (const auto &[mz, intensity, scanNumber, ionMobilityIndex] : xicPoints) {

            if (!(frameIndexMin <= scanNumber && scanNumber <= frameIndexMax)) {
                continue;
            }
            vec.coeffRef(scanNumber) += intensity;
        }

        return vec;
    }

    Err calculateMs1Scores(
        const Eigen::VectorX<float> &kernel,
        const Eigen::VectorX<float> &_anchorColumn,
        float mzToExtract,
        float massTol,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        TurboXIC *turboXicMS1,
        float *cosineSimMS1,
        float *mzMS1Mean,
        float *mzMS1StDev,
        float *ppmDiff,
        float *ms1IntensitySum,
        float *ms1ApexIntensity
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(turboXicMS1->isInit()); ree;
        e = ErrorUtils::isTrue(_anchorColumn.size() > 0); ree;

        Eigen::VectorX<float> anchorColumn = _anchorColumn;

        *cosineSimMS1 = 0.01;
        *mzMS1Mean = 0.0;
        *mzMS1StDev = 10.0;
        *ppmDiff = 50.0;
        *ms1IntensitySum = 0.0;

        XICPoints xicPoints = turboXicMS1->extractPointsXIC(
            mzToExtract - massTol,
            mzToExtract + massTol
            );

        TurboXIC::filterXICPointsByScanNumber(
            frameIndexMin,
            frameIndexMax,
            &xicPoints
            );

        if (xicPoints.empty()) {
            ERR_RETURN
        }

        Eigen::VectorX<float> xicVec = buildXicPointsVector(
            xicPoints,
            frameIndexMin,
            frameIndexMax
            ).segment(frameIndexMin, frameIndexMax - frameIndexMin + 1).eval();

        if (xicVec.size() < anchorColumn.size()) {
            Eigen::VectorX<float> xicVecResized(anchorColumn.size());
            xicVecResized.setZero();
            const int stepSize = static_cast<int>(std::round(anchorColumn.size() / xicVec.size()));

            int ogVecIndex = 0;
            for(int i = 0; i >= anchorColumn.size(); i += stepSize) {

                if (ogVecIndex < xicVec.size()) {
                    break;
                }

                xicVecResized.coeffRef(i) = xicVec.coeff(ogVecIndex++);
            }

            xicVec = xicVecResized;
            if (stepSize >= 1) {
                xicVec = EigenKernelUtils::convolveVectorWithKernel(xicVec, kernel);
                anchorColumn = EigenKernelUtils::convolveVectorWithKernel(anchorColumn, kernel);
            }
        }
        else if (xicVec.size() > anchorColumn.size()) {

            Eigen::VectorX<float> anchorColumnResized(xicVec.size());
            anchorColumnResized.setZero();
            const int stepSize = static_cast<int>(std::round(xicVec.size() / anchorColumn.size()));

            int ogVecIndex = 0;
            for(int i = 0; i < xicVec.size(); i += stepSize) {

                if (ogVecIndex >= anchorColumn.size()) {
                    break;
                }

                anchorColumnResized.coeffRef(i) = anchorColumn.coeff(ogVecIndex++);
            }

            anchorColumn = anchorColumnResized;
            if (stepSize > 1) {
                xicVec = EigenKernelUtils::convolveVectorWithKernel(xicVec, kernel);
                anchorColumn = EigenKernelUtils::convolveVectorWithKernel(anchorColumn, kernel);
            }
        }

        e = EigenUtils::cosineSimilarity(
            anchorColumn,
            xicVec,
            cosineSimMS1
            ); ree;

        const int xicPointsSize = static_cast<float>(xicPoints.size());

        const auto meanLogic = [](float sum, const XICPoint &p){return sum + p.mz;};
        *mzMS1Mean = std::accumulate(
            xicPoints.begin(),
            xicPoints.end(),
            0.0f,
            meanLogic
            ) / xicPointsSize;

        const auto stDevLogic = [mzMS1Mean](float sum, const XICPoint &p){return sum + std::pow((p.mz - *mzMS1Mean) , 2);};
        *mzMS1StDev = std::sqrt(std::accumulate(xicPoints.begin(), xicPoints.end(), 0.0f, stDevLogic)
                    / xicPointsSize);

        *ppmDiff = MathUtils::calculateMassAccuracyPPM(mzToExtract, *mzMS1Mean);
        *ms1IntensitySum = xicVec.sum();
        *ms1ApexIntensity = xicVec.maxCoeff();

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::setMs1RelatedScores(
    const TargetDecoyCandidatePair *targetDecoyCandidatePair,
    const BestCorrelationResult &bestCorrelationResult,
    float ppmTol,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_turboXicMS1->isInit()); ree;

    FrameIndex frameIndexMinMS1;
    e = m_msFrameMS1->frameIndexFromScanTime(candidateScores->scanTimeStart, &frameIndexMinMS1); ree;

    FrameIndex frameIndexMaxMS1;
    e = m_msFrameMS1->frameIndexFromScanTime(candidateScores->scanTimeEnd, &frameIndexMaxMS1); ree;

    const float isotopeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / targetDecoyCandidatePair->charge();

    const float monoIsotopeMz = targetDecoyCandidatePair->mz(false);
    const float monoIsotopeShadowMz = monoIsotopeMz - isotopeDistance;
    const float c13isotopeMz1 = monoIsotopeMz + isotopeDistance;
    const float c13isotopeMz2 = monoIsotopeMz + (isotopeDistance * 2);

    const float massTol = MathUtils::calculatePPM(monoIsotopeMz, ppmTol);

    const Eigen::VectorX<float> anchorColumn
            = bestCorrelationResult.matBlockTrimmedIntensity.col(bestCorrelationResult.bestAnchorColumnIndex);

    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        monoIsotopeMz,
        massTol,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CosineSim100MS1],
        &candidateScores->featuresArray[Ms1MzMeanFound100],
        &candidateScores->featuresArray[Ms1MzStDevFound100],
        &candidateScores->featuresArray[Ms1MzMeanFound100PPM],
        &candidateScores->featuresArray[Ms1IntensityFound100],
        &candidateScores->featuresArray[Ms1IntensityFoundApex100]
        ); ree;

    float unusedApexIntensity;
    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        monoIsotopeMz,
        massTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CosineSim45MS1],
        &candidateScores->featuresArray[Ms1MzMeanFound45],
        &candidateScores->featuresArray[Ms1MzStDevFound45],
        &candidateScores->featuresArray[Ms1MzMeanFound45PPM],
        &candidateScores->featuresArray[Ms1IntensityFound45],
        &unusedApexIntensity
        ); ree;

    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        monoIsotopeShadowMz,
        massTol,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CosineSim100MS1PreMono],
        &candidateScores->featuresArray[Ms1MzMeanFoundPreMono],
        &candidateScores->featuresArray[Ms1MzStDevFoundPreMono],
        &candidateScores->featuresArray[Ms1MzMeanFoundPreMonoPPM],
        &candidateScores->featuresArray[Ms1IntensityFoundPreMono],
        &unusedApexIntensity
        ); ree;

    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        c13isotopeMz1,
        massTol,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CosineSim100MS1Iso1],
        &candidateScores->featuresArray[Ms1MzMeanFoundIso1],
        &candidateScores->featuresArray[Ms1MzStDevFoundIso1],
        &candidateScores->featuresArray[Ms1MzMeanFoundIso1PPM],
        &candidateScores->featuresArray[Ms1IntensityFoundIso1],
        &unusedApexIntensity
        ); ree;
    
    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        c13isotopeMz2,
        massTol,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CosineSim100MS1Iso2],
        &candidateScores->featuresArray[Ms1MzMeanFoundIso2],
        &candidateScores->featuresArray[Ms1MzStDevFoundIso2],
        &candidateScores->featuresArray[Ms1MzMeanFoundIso2PPM],
        &candidateScores->featuresArray[Ms1IntensityFoundIso2],
        &unusedApexIntensity
        ); ree;

    candidateScores->featuresArray[CosineSimSum100MS1]
        = candidateScores->featuresArray[CosineSim100MS1]
        + candidateScores->featuresArray[CosineSim100MS1Iso1]
        + candidateScores->featuresArray[CosineSim100MS1Iso2];
        // + candidateScores->featuresArray[CosineSim100MS1Iso3];
        // - candidateScores->featuresArray[CosineSim100MS1PreMono];

    candidateScores->featuresArray[MonoPreMonoRatio] = std::max(candidateScores->featuresArray[CosineSim100MS1PreMono]
                                                     / std::max(candidateScores->featuresArray[CosineSim100MS1], 1.0f), 1.0f);

    ERR_RETURN
}

Err CandidateScorertron::setLibraryIonMobilityRelatedScores(
    const TargetDecoyCandidatePair *targetDecoyCandidatePair,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    if (m_msReaderPointerAcc == nullptr || m_msReaderPointerAcc->ptr.isNull()) {
        ERR_RETURN
    }

    if (!m_msReaderPointerAcc->ptr->isTIMS()) {
        ERR_RETURN
    }

    const float mobilityCenter = ionMobilityCenter(targetDecoyCandidatePair);
    if (mobilityCenter <= 0.0f) {
        ERR_RETURN
    }

    const QMap<FrameNumberTIMS, Ms1FrameTIMS> *frameNumberVsMs1FrameTIMS
        = m_msReaderPointerAcc->ptr->frameNumberVsMS1FrameTIMSPntr();
    if (frameNumberVsMs1FrameTIMS == nullptr || frameNumberVsMs1FrameTIMS->isEmpty()) {
        ERR_RETURN
    }

    const QVector<FrameNumberTIMS> &frameNumbers = m_ms1FrameNumbersTIMS;
    const int closestIndex = closestMs1FrameIndexAtOrBefore(frameNumbers, candidateScores->scanNumber);
    if (closestIndex < 0 || closestIndex >= frameNumbers.size()) {
        ERR_RETURN
    }

    FrameNumberTIMS ms1FrameNumber = frameNumbers.at(closestIndex);
    if (ms1FrameNumber > candidateScores->scanNumber && closestIndex > 0) {
        ms1FrameNumber = frameNumbers.at(closestIndex - 1);
    }

    const auto ms1FrameIt = frameNumberVsMs1FrameTIMS->constFind(ms1FrameNumber);
    if (ms1FrameIt == frameNumberVsMs1FrameTIMS->constEnd() || ms1FrameIt.value().isEmpty()) {
        ERR_RETURN
    }

    const float monoIsotopeMz = targetDecoyCandidatePair->mz(false);
    const float massTol = MathUtils::calculatePPM(
        monoIsotopeMz,
        static_cast<float>(m_pythiaParameters.ms1ExtractionWidthPPM)
        );
    const float mzMin = monoIsotopeMz - massTol;
    const float mzMax = monoIsotopeMz + massTol;

    float apexIntensity = 0.0f;
    IonMobilityIndex apexIonMobilityIndex = -1;
    double apexDriftTime = -1.0;
    IonMobilityIndex ionMobilityIndexStart = std::numeric_limits<IonMobilityIndex>::max();
    IonMobilityIndex ionMobilityIndexEnd = -1;

    const Ms1FrameTIMS &ms1FrameTIMS = ms1FrameIt.value();
    for (auto frameIt = ms1FrameTIMS.constBegin(); frameIt != ms1FrameTIMS.constEnd(); ++frameIt) {

        const IonMobilityIndex ionMobilityIndex = frameIt.key();

        double driftTime = -1.0;
        if (!driftTimeFromIonMobilityIndex(
            m_msReaderPointerAcc,
            m_timsMs2IonMobilityIndex,
            ionMobilityIndex,
            &driftTime
            )) {
            continue;
        }

        if (std::abs(driftTime - mobilityCenter) > ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0) {
            continue;
        }

        ionMobilityIndexStart = std::min(ionMobilityIndexStart, ionMobilityIndex);
        ionMobilityIndexEnd = std::max(ionMobilityIndexEnd, ionMobilityIndex);

        const ScanPoints &scanPoints = frameIt.value();
        updateApexFromSortedScanPoints(
            scanPoints,
            mzMin,
            mzMax,
            ionMobilityIndex,
            driftTime,
            &apexIntensity,
            &apexIonMobilityIndex,
            &apexDriftTime
            );
    }

    if (ionMobilityIndexEnd >= 0) {
        candidateScores->ionMobilityIndexStart = ionMobilityIndexStart;
        candidateScores->ionMobilityIndexEnd = ionMobilityIndexEnd;
    }

    if (apexIonMobilityIndex < 0) {
        ERR_RETURN
    }

    candidateScores->featuresArray[Ms1IntensityFoundApex100IM] = apexIntensity;
    candidateScores->ionMobilityIndex = apexIonMobilityIndex;
    candidateScores->imDriftTime = static_cast<float>(apexDriftTime);

    const float ionMobilityDelta = candidateScores->imDriftTime - mobilityCenter;
    candidateScores->featuresArray[IonMobilityDelta] = ionMobilityDelta;
    candidateScores->featuresArray[IonMobilityDeltaAbs] = std::abs(ionMobilityDelta);
    candidateScores->featuresArray[IonMobilityPdAbs] = std::sqrt(
        std::min(
            static_cast<double>(std::abs(ionMobilityDelta)),
            ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0
            ) / ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0
        );

    ERR_RETURN
}

Err CandidateScorertron::setMs2IonMobilityRelatedScores(
    const TargetDecoyCandidatePair *targetDecoyCandidatePair,
    const QVector<MS2Ion> &ms2Ions,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    if (m_timsMs2IonMobilityIndex == nullptr || !m_timsMs2IonMobilityIndex->isInit()) {
        ERR_RETURN
    }

    if (m_msReaderPointerAcc == nullptr || m_msReaderPointerAcc->ptr.isNull()) {
        ERR_RETURN
    }

    const float mobilityCenter = ionMobilityCenter(targetDecoyCandidatePair);
    if (mobilityCenter <= 0.0f || ms2Ions.isEmpty()) {
        ERR_RETURN
    }

    constexpr int maxMs2IonMobilityFragmentCount = 6;
    const int topIonCount = std::min({m_topNMS2Ions, ms2Ions.size(), maxMs2IonMobilityFragmentCount});
    if (topIonCount <= 0 || candidateScores->frameIndexEnd < candidateScores->frameIndexStart) {
        ERR_RETURN
    }

    const FrameIndex frameIndexMin = std::max(0, candidateScores->frameIndexStart - 1);
    const FrameIndex frameIndexMax = candidateScores->frameIndexEnd + 1;
    const float targetedIonMobilityWindowHalfWidth
        = static_cast<float>(m_pythiaParameters.timsTargetedMs2IonMobilityWindow);
    const float ionMobilityMin = mobilityCenter - static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
    const float ionMobilityMax = mobilityCenter + static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
    const IonMobilityIndex broadIonMobilityIndexStart = candidateScores->ionMobilityIndexStart;
    const IonMobilityIndex broadIonMobilityIndexEnd = candidateScores->ionMobilityIndexEnd;

    QMap<IonMobilityIndex, double> summedMobilityProfile;
    int matchedIonCount = 0;
    QVector<float> apexDeltaAbsValues;
    QVector<float> mobilityFwhmValues;
    QVector<float> mobilityFwhmWeights;
    QVector<QMap<IonMobilityIndex, double>> fragmentMobilityProfiles;
    apexDeltaAbsValues.reserve(topIonCount);
    mobilityFwhmValues.reserve(topIonCount);
    mobilityFwhmWeights.reserve(topIonCount);
    fragmentMobilityProfiles.reserve(topIonCount);

    constexpr float timsRtMobilityCoelutionMinSpectrumOverTime = 0.08f;
    constexpr float timsRtMobilityCoelutionMinTotalIntensityLog = 8.0f;
    const bool computeRtMobilityCoelutionFeatures
        = (m_features.contains(Ms2IonMobilityRtCosineMean)
           || m_features.contains(Ms2IonMobilityRtCosineStDev)
           || m_features.contains(Ms2IonMobilityRtApexAgreementFraction))
          && candidateScores->featuresArray[CosineSimSpectrumOverTimeCubed]
             >= timsRtMobilityCoelutionMinSpectrumOverTime
          && candidateScores->featuresArray[TotalIntensityLog]
             >= timsRtMobilityCoelutionMinTotalIntensityLog;

    if (!computeRtMobilityCoelutionFeatures) {
        for (int i = 0; i < topIonCount; ++i) {
            const MS2Ion &ms2Ion = ms2Ions.at(i);
            const float massTol = MathUtils::calculatePPM(
                ms2Ion.mz,
                static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
                );

            float apexIntensity = 0.0f;
            float apexDeltaAbs = static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
            QMap<IonMobilityIndex, double> fragmentMobilityProfile;
            const bool hasMobilityProfile = m_timsMs2IonMobilityIndex->extractMobilityProfile(
                ms2Ion.mz - massTol,
                ms2Ion.mz + massTol,
                frameIndexMin,
                frameIndexMax,
                ionMobilityMin,
                ionMobilityMax,
                mobilityCenter,
                &fragmentMobilityProfile,
                &apexIntensity,
                &apexDeltaAbs
                );

            if (!hasMobilityProfile) {
                continue;
            }

            matchedIonCount++;

            for (auto profileIt = fragmentMobilityProfile.constBegin();
                 profileIt != fragmentMobilityProfile.constEnd();
                 ++profileIt) {
                summedMobilityProfile[profileIt.key()] += profileIt.value();
            }

            apexDeltaAbsValues.push_back(apexDeltaAbs);
            fragmentMobilityProfiles.push_back(fragmentMobilityProfile);
            mobilityFwhmWeights.push_back(std::max(ms2Ion.intensity, 0.0f));
        }
    }
    else {
        using RtMobilityKey = quint64;
        constexpr RtMobilityKey invalidRtMobilityKey = std::numeric_limits<RtMobilityKey>::max();
        const auto makeRtMobilityKey = [](FrameIndex frameIndex, IonMobilityIndex ionMobilityIndex) {
            return (static_cast<RtMobilityKey>(static_cast<quint32>(frameIndex)) << 32)
                   | static_cast<quint32>(ionMobilityIndex);
        };
        const auto keyFrameIndex = [](RtMobilityKey key) {
            return static_cast<FrameIndex>(key >> 32);
        };
        const auto keyIonMobilityIndex = [](RtMobilityKey key) {
            return static_cast<IonMobilityIndex>(key & 0xffffffffu);
        };

        std::unordered_map<RtMobilityKey, double> summedRtMobilityProfile;
        QVector<std::unordered_map<RtMobilityKey, double>> fragmentRtMobilityProfiles;
        QVector<RtMobilityKey> fragmentRtMobilityApexes;
        fragmentRtMobilityProfiles.reserve(topIonCount);
        fragmentRtMobilityApexes.reserve(topIonCount);

        for (int i = 0; i < topIonCount; ++i) {
            const MS2Ion &ms2Ion = ms2Ions.at(i);
            const float massTol = MathUtils::calculatePPM(
                ms2Ion.mz,
                static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
                );

            float apexIntensity = 0.0f;
            float apexDeltaAbs = static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
            QMap<IonMobilityIndex, double> fragmentMobilityProfile;
            std::unordered_map<RtMobilityKey, double> fragmentRtMobilityProfile;
            RtMobilityKey fragmentApexKey = invalidRtMobilityKey;

            const XICPoints xicPoints = m_timsMs2IonMobilityIndex->extractPointsXIC(
                ms2Ion.mz - massTol,
                ms2Ion.mz + massTol,
                frameIndexMin,
                frameIndexMax,
                ionMobilityMin,
                ionMobilityMax
                );
            fragmentRtMobilityProfile.reserve(static_cast<size_t>(xicPoints.size()));

            for (const XICPoint &xicPoint : xicPoints) {
                if (xicPoint.intensity <= 0.0f || xicPoint.ionMobilityIndex < 0) {
                    continue;
                }

                float driftTime = -1.0f;
                if (!m_timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(
                        xicPoint.ionMobilityIndex,
                        &driftTime
                        )) {
                    continue;
                }

                const double intensity = std::max(0.0f, xicPoint.intensity);
                fragmentMobilityProfile[xicPoint.ionMobilityIndex] += intensity;

                const RtMobilityKey rtMobilityKey = makeRtMobilityKey(
                    xicPoint.scanNumber,
                    xicPoint.ionMobilityIndex
                    );
                fragmentRtMobilityProfile[rtMobilityKey] += intensity;

                if (xicPoint.intensity > apexIntensity) {
                    apexIntensity = xicPoint.intensity;
                    apexDeltaAbs = std::abs(driftTime - mobilityCenter);
                    fragmentApexKey = rtMobilityKey;
                }
            }

            if (fragmentMobilityProfile.isEmpty() || fragmentRtMobilityProfile.empty()) {
                continue;
            }

            matchedIonCount++;

            for (auto profileIt = fragmentMobilityProfile.constBegin();
                 profileIt != fragmentMobilityProfile.constEnd();
                 ++profileIt) {
                summedMobilityProfile[profileIt.key()] += profileIt.value();
            }

            apexDeltaAbsValues.push_back(apexDeltaAbs);
            fragmentMobilityProfiles.push_back(fragmentMobilityProfile);
            mobilityFwhmWeights.push_back(std::max(ms2Ion.intensity, 0.0f));
            fragmentRtMobilityProfiles.push_back(fragmentRtMobilityProfile);
            fragmentRtMobilityApexes.push_back(fragmentApexKey);

            for (const auto &profileEntry : fragmentRtMobilityProfile) {
                summedRtMobilityProfile[profileEntry.first] += profileEntry.second;
            }
        }

        if (!summedRtMobilityProfile.empty() && fragmentRtMobilityProfiles.size() > 1) {
            RtMobilityKey consensusApex = invalidRtMobilityKey;
            double consensusApexIntensity = 0.0;
            for (const auto &profileEntry : summedRtMobilityProfile) {
                if (profileEntry.second > consensusApexIntensity) {
                    consensusApexIntensity = profileEntry.second;
                    consensusApex = profileEntry.first;
                }
            }

            if (consensusApex != invalidRtMobilityKey && !fragmentRtMobilityApexes.isEmpty()) {
                if (computeRtMobilityCoelutionFeatures) {
                    QVector<float> rtMobilityCosines;
                    rtMobilityCosines.reserve(fragmentRtMobilityProfiles.size());

                    double totalNormSquared = 0.0;
                    for (const auto &totalEntry : summedRtMobilityProfile) {
                        totalNormSquared += totalEntry.second * totalEntry.second;
                    }

                    for (const std::unordered_map<RtMobilityKey, double> &fragmentProfile : fragmentRtMobilityProfiles) {
                        if (fragmentProfile.empty()) {
                            continue;
                        }

                        double fragmentTotalDotProduct = 0.0;
                        double fragmentNormSquared = 0.0;

                        for (const auto &fragmentEntry : fragmentProfile) {
                            const double fragmentIntensity = fragmentEntry.second;
                            const auto totalIt = summedRtMobilityProfile.find(fragmentEntry.first);
                            const double totalIntensity = totalIt == summedRtMobilityProfile.end()
                                                              ? 0.0
                                                              : totalIt->second;
                            fragmentTotalDotProduct += fragmentIntensity * totalIntensity;
                            fragmentNormSquared += fragmentIntensity * fragmentIntensity;
                        }

                        const double dotProduct = fragmentTotalDotProduct - fragmentNormSquared;
                        const double consensusNormSquared = totalNormSquared
                                                            - (2.0 * fragmentTotalDotProduct)
                                                            + fragmentNormSquared;
                        if (fragmentNormSquared <= 0.0 || consensusNormSquared <= 0.0) {
                            continue;
                        }

                        const double cosine = dotProduct / std::sqrt(fragmentNormSquared * consensusNormSquared);
                        rtMobilityCosines.push_back(static_cast<float>(std::clamp(cosine, 0.0, 1.0)));
                    }

                    if (!rtMobilityCosines.isEmpty()) {
                        candidateScores->featuresArray[Ms2IonMobilityRtCosineMean] = MathUtils::mean(rtMobilityCosines);
                        candidateScores->featuresArray[Ms2IonMobilityRtCosineStDev] = MathUtils::stDev(rtMobilityCosines);
                    }
                }

                int agreeingApexCount = 0;
                const FrameIndex consensusFrameIndex = keyFrameIndex(consensusApex);
                const IonMobilityIndex consensusIonMobilityIndex = keyIonMobilityIndex(consensusApex);
                for (const RtMobilityKey &fragmentApex : fragmentRtMobilityApexes) {
                    if (fragmentApex == invalidRtMobilityKey) {
                        continue;
                    }

                    if (std::abs(keyFrameIndex(fragmentApex) - consensusFrameIndex) <= 1
                        && std::abs(keyIonMobilityIndex(fragmentApex) - consensusIonMobilityIndex) <= 2) {
                        agreeingApexCount++;
                    }
                }

                candidateScores->featuresArray[Ms2IonMobilityRtApexAgreementFraction]
                    = agreeingApexCount / static_cast<float>(fragmentRtMobilityApexes.size());

            }
        }
    }

    candidateScores->featuresArray[Ms2IonMobilityMatchedIonFraction]
        = matchedIonCount / static_cast<float>(topIonCount);

    double intensitySum = 0.0;
    double weightedDelta = 0.0;
    double weightedDeltaAbs = 0.0;
    if (!summedMobilityProfile.isEmpty()) {
        QVector<IonMobilityIndex> mobilityIndices = summedMobilityProfile.keys().toVector();
        std::sort(mobilityIndices.begin(), mobilityIndices.end());

        QVector<double> summedProfile;
        summedProfile.reserve(mobilityIndices.size());
        double bestIntensity = -1.0;
        int bestIndex = -1;

        for (int i = 0; i < mobilityIndices.size(); ++i) {
            const IonMobilityIndex ionMobilityIndex = mobilityIndices.at(i);
            const double intensity = summedMobilityProfile.value(ionMobilityIndex);
            summedProfile.push_back(intensity);

            float driftTime = -1.0f;
            if (!m_timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(ionMobilityIndex, &driftTime)) {
                continue;
            }

            const double delta = driftTime - mobilityCenter;
            intensitySum += intensity;
            weightedDelta += intensity * delta;
            weightedDeltaAbs += intensity * std::abs(delta);

            if (intensity > bestIntensity) {
                bestIntensity = intensity;
                bestIndex = i;
            }
        }

        SymmetricProfileLimits limits;
        float observedMobilityForWindow = -1.0f;
        if (bestIndex >= 0) {
            const IonMobilityIndex observedIonMobilityIndex = mobilityIndices.at(bestIndex);
            float observedDriftTime = -1.0f;
            if (m_timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(
                    observedIonMobilityIndex,
                    &observedDriftTime
                    )) {

                candidateScores->ionMobilityIndex = observedIonMobilityIndex;
                candidateScores->imDriftTime = observedDriftTime;
                observedMobilityForWindow = observedDriftTime;

                const float ionMobilityDelta = observedDriftTime - mobilityCenter;
                candidateScores->featuresArray[IonMobilityDelta] = ionMobilityDelta;
                candidateScores->featuresArray[IonMobilityDeltaAbs] = std::abs(ionMobilityDelta);
                candidateScores->featuresArray[IonMobilityPdAbs] = std::sqrt(
                    std::min(
                        static_cast<double>(std::abs(ionMobilityDelta)),
                        ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0
                        ) / ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0
                    );
            }

            limits = alphaDiaStyleSymmetricLimits1d(summedProfile, bestIndex);
            if (limits.startIndex >= 0
                && limits.stopIndex >= limits.startIndex
                && limits.stopIndex < mobilityIndices.size()) {
                candidateScores->ionMobilityIndexStart = mobilityIndices.at(limits.startIndex);
                candidateScores->ionMobilityIndexEnd = mobilityIndices.at(limits.stopIndex);
            }
        }

        QVector<IonMobilityIndex> fwhmMobilityIndices;
        if (observedMobilityForWindow > 0.0f) {
            const IonMobilityIndex indexStart = broadIonMobilityIndexStart >= 0 && broadIonMobilityIndexEnd >= 0
                                                    ? std::min(broadIonMobilityIndexStart, broadIonMobilityIndexEnd)
                                                    : mobilityIndices.front();
            const IonMobilityIndex indexEnd = broadIonMobilityIndexStart >= 0 && broadIonMobilityIndexEnd >= 0
                                                  ? std::max(broadIonMobilityIndexStart, broadIonMobilityIndexEnd)
                                                  : mobilityIndices.back();

            for (IonMobilityIndex ionMobilityIndex = indexStart; ionMobilityIndex <= indexEnd; ++ionMobilityIndex) {
                float driftTime = -1.0f;
                if (!m_timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(ionMobilityIndex, &driftTime)) {
                    continue;
                }

                if (std::abs(driftTime - observedMobilityForWindow)
                    <= targetedIonMobilityWindowHalfWidth) {
                    fwhmMobilityIndices.push_back(ionMobilityIndex);
                }
            }
        }

        if (fwhmMobilityIndices.isEmpty()
            && limits.startIndex >= 0
            && limits.stopIndex >= limits.startIndex
            && limits.stopIndex < mobilityIndices.size()) {
            for (int i = limits.startIndex; i <= limits.stopIndex; ++i) {
                fwhmMobilityIndices.push_back(mobilityIndices.at(i));
            }
        }

        if (!fwhmMobilityIndices.isEmpty()) {

            float mobilityStart = -1.0f;
            float mobilityStop = -1.0f;
            if (m_timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(
                    fwhmMobilityIndices.front(),
                    &mobilityStart
                    )
                && m_timsMs2IonMobilityIndex->driftTimeFromIonMobilityIndex(
                    fwhmMobilityIndices.back(),
                    &mobilityStop
                    )) {

                candidateScores->ionMobilityIndexStart = fwhmMobilityIndices.front();
                candidateScores->ionMobilityIndexEnd = fwhmMobilityIndices.back();

                const float mobilityWidth = std::abs(mobilityStop - mobilityStart);
                const int mobilityWindowBinCount = std::max(1, fwhmMobilityIndices.size());
                const QVector<float> fragmentMobilityProfileWeights = mobilityFwhmWeights;
                mobilityFwhmWeights.clear();
                for (int fragmentIndex = 0; fragmentIndex < fragmentMobilityProfiles.size(); ++fragmentIndex) {
                    const QMap<IonMobilityIndex, double> &fragmentMobilityProfile
                        = fragmentMobilityProfiles.at(fragmentIndex);
                    if (fragmentMobilityProfile.isEmpty()) {
                        continue;
                    }

                    double apexIntensity = 0.0;
                    for (IonMobilityIndex ionMobilityIndex : fwhmMobilityIndices) {
                        apexIntensity = std::max(
                            apexIntensity,
                            fragmentMobilityProfile.value(ionMobilityIndex, 0.0)
                            );
                    }

                    if (apexIntensity <= 0.0) {
                        continue;
                    }

                    const double halfMaxIntensity = apexIntensity / 2.0;
                    int valuesAboveHalfMax = 0;
                    for (IonMobilityIndex ionMobilityIndex : fwhmMobilityIndices) {
                        if (fragmentMobilityProfile.value(ionMobilityIndex, 0.0) > halfMaxIntensity) {
                            valuesAboveHalfMax++;
                        }
                    }

                    const float fractionAboveHalfMax = valuesAboveHalfMax
                                                       / static_cast<float>(mobilityWindowBinCount);
                    mobilityFwhmValues.push_back(fractionAboveHalfMax * mobilityWidth);
                    mobilityFwhmWeights.push_back(fragmentMobilityProfileWeights.at(fragmentIndex));
                }
            }
        }

    }

    if (intensitySum > 0.0) {
        candidateScores->featuresArray[Ms2IonMobilityWeightedDelta]
            = static_cast<float>(weightedDelta / intensitySum);
        candidateScores->featuresArray[Ms2IonMobilityWeightedDeltaAbs]
            = static_cast<float>(weightedDeltaAbs / intensitySum);
    }
    else {
        candidateScores->featuresArray[Ms2IonMobilityWeightedDelta] = 0.0f;
        candidateScores->featuresArray[Ms2IonMobilityWeightedDeltaAbs]
            = static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
    }

    if (!apexDeltaAbsValues.isEmpty()) {
        candidateScores->featuresArray[Ms2IonMobilityApexDeltaAbsMean]
            = MathUtils::mean(apexDeltaAbsValues);
        candidateScores->featuresArray[Ms2IonMobilityApexDeltaAbsStDev]
            = MathUtils::stDev(apexDeltaAbsValues);
    }
    else {
        candidateScores->featuresArray[Ms2IonMobilityApexDeltaAbsMean]
            = static_cast<float>(ALPHADIA_MOBILITY_TOLERANCE_ONE_OVER_K0);
        candidateScores->featuresArray[Ms2IonMobilityApexDeltaAbsStDev] = 0.0f;
    }

    if (!mobilityFwhmValues.isEmpty()) {
        double fwhmWeightSum = std::accumulate(
            mobilityFwhmWeights.begin(),
            mobilityFwhmWeights.end(),
            0.0
            );

        if (fwhmWeightSum <= 0.0) {
            fwhmWeightSum = mobilityFwhmValues.size();
            mobilityFwhmWeights.fill(1.0f);
        }

        double weightedFwhmSum = 0.0;
        for (int i = 0; i < mobilityFwhmValues.size(); ++i) {
            weightedFwhmSum += mobilityFwhmValues.at(i) * mobilityFwhmWeights.at(i);
        }

        const float weightedFwhmMean = static_cast<float>(weightedFwhmSum / fwhmWeightSum);
        QVector<float> fwhmResiduals;
        fwhmResiduals.reserve(mobilityFwhmValues.size());
        for (const float mobilityFwhmValue : mobilityFwhmValues) {
            fwhmResiduals.push_back(mobilityFwhmValue - weightedFwhmMean);
        }

        candidateScores->featuresArray[Ms2IonMobilityFwhmMean] = weightedFwhmMean;
        candidateScores->featuresArray[Ms2IonMobilityFwhmStDev] = MathUtils::stDev(fwhmResiduals);
    }
    else {
        candidateScores->featuresArray[Ms2IonMobilityFwhmMean] = 0.0f;
        candidateScores->featuresArray[Ms2IonMobilityFwhmStDev] = 0.0f;
    }

    ERR_RETURN
}

namespace {

    Err closest(
        const QVector<MS2Ion> &ms2Ions,
        double mzVal,
        MS2Ion *ms2IonFound
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;

        auto it = std::min_element(ms2Ions.begin(), ms2Ions.end(), [mzVal] (const MS2Ion &a, const MS2Ion &b) {
            return std::abs(a.mz - mzVal) <= std::abs(b.mz - mzVal);
        });

        if(it == ms2Ions.end()) {
            ERR_RETURN
        }

        *ms2IonFound = *it;

        ERR_RETURN
    }

    Err extractFullTheoreticalPointsFromScan(
        const ScanPoints* scanPoints,
        const QVector<MS2Ion> &ms2IonsTheoritical,
        double ms2ExtractionWidthPPM,
        QVector<QPair<QPointF, MS2Ion>> *foundPointVsMS2Ions
        ) {

        ERR_INIT

        QVector<QPointF> scanPointsQF;
        std::transform(
            scanPoints->begin(),
            scanPoints->end(),
            std::back_inserter(scanPointsQF),
            [](const ScanPoint& scanPoint){return QPointF(static_cast<double>(scanPoint.x()), static_cast<double>(scanPoint.y()));}
            );

        QVector<double> mzVals;
        std::transform(
            ms2IonsTheoritical.begin(),
            ms2IonsTheoritical.end(),
            std::back_inserter(mzVals),
            [](const MS2Ion& ms2Ion){return static_cast<double>(ms2Ion.mz);}
            );

        const QVector<QPointF> foundPoints = MsUtils::extractPointsFromPoints(
            scanPointsQF,
            mzVals,
            ms2ExtractionWidthPPM,
            true
            );


        foundPointVsMS2Ions->reserve(foundPoints.size());
        for (const QPointF &fp : foundPoints) {
            MS2Ion ms2Ion;
            e = closest(
                ms2IonsTheoritical,
                fp.x(),
                &ms2Ion
                ); ree;

            foundPointVsMS2Ions->push_back({fp, ms2Ion});
        }

        ERR_RETURN
    }

    QVector<int> findLongestSeriesInSet(const QSet<int> &set) {

        QVector<int> vec = {set.begin(), set.end()};
        std::sort(vec.begin(), vec.end());

        QVector<int> longestSeriesCountVec;
        int longestSeriesCountPassing = 0;

        int currentVal = -2;
        for (int v : vec) {

            if (currentVal == v - 1) {
                longestSeriesCountPassing++;
                currentVal = v;
                continue;
            }

            currentVal = v;

            if (longestSeriesCountPassing > 0) {
                longestSeriesCountVec.push_back(longestSeriesCountPassing + 1);
            }
            longestSeriesCountPassing = 0;
        }

        if (longestSeriesCountPassing > 0) {
            longestSeriesCountVec.push_back(longestSeriesCountPassing + 1);
        }

        return longestSeriesCountVec;
    }

    Err buildFoundIonSeriesSets(
        const QVector<QPair<QPointF, MS2Ion>> &foundPointVsMS2Ions,
        QSet<int> *yIonsSeries,
        QSet<int> *bIonsSeries
        ) {

        ERR_INIT

        for (const QPair<QPointF, MS2Ion> &pr : foundPointVsMS2Ions) {
            QPair<IonIndex, IonType> ionInfo;
            e = pr.second.getIonLabelInfo(&ionInfo); ree;

            if (ionInfo.second.contains("y")) {
                yIonsSeries->insert(ionInfo.first);
            }
            else if (ionInfo.second.contains("b")) {
                bIonsSeries->insert(ionInfo.first);
            }
            else {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Ion other than b or y found";
            }
        }

        ERR_RETURN
    }

    Err buildTheoIonSeriesSets(
	    const QVector<MS2Ion> &ms2IonsTheo,
	    QSet<int> *yIonsSeriesTheo,
	    QSet<int> *bIonsSeriesTheo
	    ) {

        ERR_INIT

        for (const MS2Ion &ms2Ion : ms2IonsTheo) {
            QPair<IonIndex, IonType> ionInfo;
            e = ms2Ion.getIonLabelInfo(&ionInfo); ree;

            if (ionInfo.second.contains("y")) {
                yIonsSeriesTheo->insert(ionInfo.first);
            }
            else if (ionInfo.second.contains("b")) {
                bIonsSeriesTheo->insert(ionInfo.first);
            }
            else {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Ion other than b or y found";
            }
        }

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::setFullTheoMs2IonsScores(CandidateScores *candidateScores) const {

    ERR_INIT

    const ScanPoints* scanPoints = m_msFrameMzTarget->getScanPointsByScanNumber(candidateScores->scanNumber);

    const QVector<MS2Ion> ms2IonsTheoritical = candidateScores->isDecoy
                                     ? candidateScores->targetDecoyCandidatePair->ms2IonsDecoy()
                                     : candidateScores->targetDecoyCandidatePair->ms2IonsTarget();

	constexpr int ionLabelsSize = 12;
	candidateScores->ionLabels.resize(ionLabelsSize);
	for (const MS2Ion &msi : ms2IonsTheoritical) {
		if (msi.rank >= ionLabelsSize) continue;
		e = ErrorUtils::isIndexInRange(candidateScores->ionLabels, msi.rank); ree;
		candidateScores->ionLabels[msi.rank] = msi.ionLabel;
	}

    QVector<QPair<QPointF, MS2Ion>> foundPointVsMS2Ions;
    e = extractFullTheoreticalPointsFromScan(
        scanPoints,
        ms2IonsTheoritical,
        m_pythiaParameters.ms2ExtractionWidthPPM,
        &foundPointVsMS2Ions
        ); ree;

    Eigen::VectorX<float> v1(foundPointVsMS2Ions.size());
    Eigen::VectorX<float> v2(foundPointVsMS2Ions.size());
    for (int i = 0; i < foundPointVsMS2Ions.size(); i++) {
        v1.coeffRef(i) = foundPointVsMS2Ions.at(i).first.y();
        v2.coeffRef(i) = foundPointVsMS2Ions.at(i).second.intensity;
    }

    float cosineSimFullTheo;
    e = EigenUtils::cosineSimilarity(v1, v2, &cosineSimFullTheo); ree;
    candidateScores->featuresArray[CosineSimFullTheo] = cosineSimFullTheo;

    candidateScores->featuresArray[IonsFoundFractionFull] = foundPointVsMS2Ions.size() / static_cast<float>(ms2IonsTheoritical.size());
    candidateScores->featuresArray[CosineSimFullTheoXIonsFoundFractionFull] = cosineSimFullTheo * candidateScores->featuresArray[IonsFoundFractionFull];

    QSet<int> yIonsSeries;
    QSet<int> bIonsSeries;
    e = buildFoundIonSeriesSets(foundPointVsMS2Ions, &yIonsSeries, &bIonsSeries);
    const QVector<int> yIonsSeriesLongest = findLongestSeriesInSet(yIonsSeries);
    const QVector<int> bIonsSeriesLongest = findLongestSeriesInSet(bIonsSeries);

    QSet<int> yIonsSeriesTheo;
    QSet<int> bIonsSeriesTheo;
    e = buildTheoIonSeriesSets(ms2IonsTheoritical, &yIonsSeriesTheo, &bIonsSeriesTheo);
    const QVector<int> yIonsSeriesLongestTheo = findLongestSeriesInSet(yIonsSeriesTheo);
    const QVector<int> bIonsSeriesLongestTheo = findLongestSeriesInSet(bIonsSeriesTheo);

    const int yIonsSeriesLongestMax =*std::max_element(yIonsSeriesLongest.begin(), yIonsSeriesLongest.end());
    const int yIonsSeriesLongestTheoMax =*std::max_element(yIonsSeriesLongestTheo.begin(), yIonsSeriesLongestTheo.end());

    const int bIonsSeriesLongestMax =*std::max_element(bIonsSeriesLongest.begin(), bIonsSeriesLongest.end());
    const int bIonsSeriesLongestTheoMax =*std::max_element(bIonsSeriesLongestTheo.begin(), bIonsSeriesLongestTheo.end());

    candidateScores->featuresArray[YIonSeriesMax] = yIonsSeriesLongestMax;
    candidateScores->featuresArray[YIonSeriesCount] = yIonsSeriesLongest.size();
    candidateScores->featuresArray[YIonSeriesMean] = MathUtils::mean(yIonsSeriesLongest);
    candidateScores->featuresArray[YIonSeriesStd] =  MathUtils::stDev(yIonsSeriesLongest);
    candidateScores->featuresArray[YIonSeriesStd] = std::isinf(candidateScores->featuresArray[YIonSeriesStd]) || std::isnan(candidateScores->featuresArray[YIonSeriesStd])
                                                            ? -1.0f
                                                            : candidateScores->featuresArray[YIonSeriesStd];

    candidateScores->featuresArray[YIonSeriesTheoMax] = yIonsSeriesLongestTheoMax;
    candidateScores->featuresArray[YIonSeriesTheoCount] = yIonsSeriesLongestTheo.size();
    candidateScores->featuresArray[YIonSeriesTheoMean] = MathUtils::mean(yIonsSeriesLongestTheo);
    candidateScores->featuresArray[YIonSeriesTheoStd] = MathUtils::stDev(yIonsSeriesLongestTheo);
    candidateScores->featuresArray[YIonSeriesTheoStd] = std::isinf(candidateScores->featuresArray[YIonSeriesTheoStd]) || std::isnan(candidateScores->featuresArray[YIonSeriesTheoStd])
                                                        ? -1.0f
                                                        : candidateScores->featuresArray[YIonSeriesTheoStd];
    candidateScores->featuresArray[YIonSeriesMaxFoundToTheoFraction] = yIonsSeriesLongestMax / std::max(static_cast<float>(yIonsSeriesLongestTheoMax), 1.0f);
    candidateScores->featuresArray[YIonSeriesCountRatio] = candidateScores->featuresArray[YIonSeriesCount]
                                                                   / std::max(candidateScores->featuresArray[YIonSeriesTheoCount], 1.0f);

    candidateScores->featuresArray[BIonSeriesMax] = bIonsSeriesLongestMax;
    candidateScores->featuresArray[BIonSeriesCount] = bIonsSeriesLongest.size();
    candidateScores->featuresArray[BIonSeriesMean] = MathUtils::mean(bIonsSeriesLongest);
    candidateScores->featuresArray[BIonSeriesStd] = MathUtils::stDev(bIonsSeriesLongest);
    candidateScores->featuresArray[BIonSeriesStd] = std::isinf(candidateScores->featuresArray[BIonSeriesStd]) || std::isnan(candidateScores->featuresArray[BIonSeriesStd])
                                                            ? -1.0f
                                                            : candidateScores->featuresArray[BIonSeriesStd];
    candidateScores->featuresArray[BIonSeriesTheoMax] = bIonsSeriesLongestTheoMax;
    candidateScores->featuresArray[BIonSeriesTheoCount] = bIonsSeriesLongestTheo.size();
    candidateScores->featuresArray[BIonSeriesTheoMean] = MathUtils::mean(bIonsSeriesLongestTheo);
    candidateScores->featuresArray[BIonSeriesTheoStd] = MathUtils::stDev(bIonsSeriesLongestTheo);
    candidateScores->featuresArray[BIonSeriesTheoStd] = std::isinf(candidateScores->featuresArray[BIonSeriesTheoStd]) || std::isnan(candidateScores->featuresArray[BIonSeriesTheoStd])
                                                            ? -1.0f
                                                            : candidateScores->featuresArray[BIonSeriesTheoStd];

    candidateScores->featuresArray[BIonSeriesMaxFoundToTheoFraction] = bIonsSeriesLongestMax / std::max(static_cast<float>(bIonsSeriesLongestTheoMax), 1.0f);
    candidateScores->featuresArray[BIonSeriesCountRatio] = candidateScores->featuresArray[BIonSeriesCount]
                                                                   / std::max(candidateScores->featuresArray[BIonSeriesTheoCount], 1.0f);

	const bool featureVectorHasNanOrInfVal = MathUtils::vectorContainsInfOrNaN(candidateScores->featuresArray);
    if (featureVectorHasNanOrInfVal) {
        qDebug() << "FeaturesArray contains at least one NaN or Inf value" << candidateScores->featuresArray;
    }
    e = ErrorUtils::isFalse(featureVectorHasNanOrInfVal); ree;

    ERR_RETURN
}
