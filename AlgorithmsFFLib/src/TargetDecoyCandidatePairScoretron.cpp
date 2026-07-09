//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "CandidateScores.h"
#include "CandidateScorertron.h"
#include "IsotopicDistributionBuilder.h"
#include "MsCalibratomatic.h"
#include "ParallelUtils.h"
#include "TimsMs2IonMobilityIndex.h"
#include "XICPeakManager.h"

#include <QtConcurrent/QtConcurrent>

#include <QSet>

#include <algorithm>
#include <cmath>


TargetDecoyCandidatePairScoretron2::TargetDecoyCandidatePairScoretron2()
: m_msReaderPointerAcc(nullptr)
, m_turboXICMS1(nullptr)
, m_msFrameMS1(nullptr)
, m_turboXIC2DMS1(nullptr)
{}

TargetDecoyCandidatePairScoretron2::~TargetDecoyCandidatePairScoretron2() {
    delete m_turboXICMS1;
    delete m_msFrameMS1;
    for (MsFrame *msFrame : m_mzTargetKeyVsMsFramePntr) {
        delete msFrame;
    }
}

class TargetDecoyPairParallelInput {

public:

    MzTargetKey targetKey;
    MsScanInfo msScanInfo;
    MsCalibratomatic msCalibratomatic;
    QMap<ScanNumber, ScanPoints*> diaTargetFrame;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    MsFrame *msFrameMzTarget = nullptr;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    int topNMs2Ions = -1.0;
    PythiaParameters pythiaParameters;
    QPair<double, double> scanTimeMinMax;
    TurboXIC *turboXicMS1 = nullptr;
    TurboXIC *turboXicMS2 = nullptr;
    MsFrame *msFrameMS1 = nullptr;
    float minPeakCount = -1.0;
    QMap<int, QVector<float>> averagineTable;
    QVector<float> weights;
    QVector<Features> features;
    bool useTopNIntegrationsParameter = false;
    bool useAdaptiveTimsMobilityCentering = false;
    MsReaderPointerAcc *msReaderPointerAcc = nullptr;
    QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePointersAllPntr = nullptr;
    bool splitMzTargetKey = false;
    bool isBottomSplit = false;
};


Err TargetDecoyCandidatePairScoretron2::init(
        const PythiaParameters &pythiaParameters,
        MsReaderPointerAcc *msReaderPointerAcc
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;

    m_msReaderPointerAcc = msReaderPointerAcc;
    m_pythiaParameters = pythiaParameters;

    m_scanNumberVsScanTime = m_msReaderPointerAcc->ptr->getScanNumberVsScanTime();
    m_scanTimeMinMax = m_msReaderPointerAcc->ptr->scanTimeMinMax();
    m_uniqueTandemMsScanInfos = m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    if (!msReaderPointerAcc->useLazyLoading()) {
        e = m_msReaderPointerAcc->ptr->collateMS2MzTargetFrames(&diaTargetFrames); ree;
    }

    if (msReaderPointerAcc->useLazyLoading()) {
        e = m_msReaderPointerAcc->ptr->getMzTargetScanPoints(
            S_GLOBAL_SETTINGS.MS1Key,
            &m_ms1ScanNumberVsScanPoints
            ); ree;
    }
    else {
        constexpr int msLevel = 1;
        e = m_msReaderPointerAcc->ptr->getScanPoints(msLevel, &m_ms1ScanNumberVsScanPoints); ree;
    }

	if (!m_ms1ScanNumberVsScanPoints.isEmpty()) {

		QMap<ScanNumber, ScanPoints*> ms1FramePtrs;

		for (auto it = m_ms1ScanNumberVsScanPoints.begin(); it != m_ms1ScanNumberVsScanPoints.end(); ++it) {
			ms1FramePtrs.insert(it.key(), &it.value());
		}

		m_msFrameMS1 = new MsFrame;
		e = m_msFrameMS1->init(ms1FramePtrs, m_msReaderPointerAcc->ptr->getScanNumberVsScanTime()); ree;

		m_turboXICMS1 = new TurboXIC();
		e = m_turboXICMS1->init(m_msFrameMS1->frameIndexVsScanPoints()); ree;
	}

	if (!diaTargetFrames.isEmpty()) {
		m_diaTargetFrames = diaTargetFrames;
		e = buildMzTargetKeyVsMsFrames(); ree;
	}

    e = buildAveragineTable(); ree;

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron2::buildMzTargetKeyVsMsFrames() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_diaTargetFrames); ree;
    e = ErrorUtils::isNotEmpty(m_scanNumberVsScanTime); ree;

    for (auto it = m_mzTargetKeyVsMsFramePntr.begin(); it != m_mzTargetKeyVsMsFramePntr.end(); ++it) {
        delete it.value();
        m_mzTargetKeyVsMsFramePntr[it.key()] = nullptr;
    }

    for (auto it = m_diaTargetFrames.begin(); it != m_diaTargetFrames.end(); ++it) {
        auto *msFrame = new MsFrame();
        e = msFrame->init(it.value(), m_scanNumberVsScanTime); ree;
        m_mzTargetKeyVsMsFramePntr.insert(it.key(), msFrame);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron2::buildAveragineTable() {

    ERR_INIT

    for (int nominalMass = 300; nominalMass < 10000; nominalMass += 10) {
        QVector<float> preMonoIncluded = {0.0f};
        const QVector<double> isoDis = IsotopicDistributionBuilder::getIsotopicDistribution(static_cast<double>(nominalMass));
        constexpr int maxAveragineVecLength = 3;
        preMonoIncluded.append({isoDis.begin(), isoDis.begin() + maxAveragineVecLength});
        m_averagineTable.insert(nominalMass, preMonoIncluded);
    }

    e = ErrorUtils::isNotEmpty(m_averagineTable); ree;

    ERR_RETURN
}

QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>>* TargetDecoyCandidatePairScoretron2::diaTargetFrames() {
    return &m_diaTargetFrames;
}

QMap<ScanNumber, ScanPoints>* TargetDecoyCandidatePairScoretron2::ms1ScanNumberVsScanPoints() {
    return &m_ms1ScanNumberVsScanPoints;
}

QMap<MzTargetKey, MsFrame*> TargetDecoyCandidatePairScoretron2::mzTargetKeyVsMsFramePntr() {
    return m_mzTargetKeyVsMsFramePntr;
}

Err TargetDecoyCandidatePairScoretron2::reloadTurboXICMS1() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_ms1ScanNumberVsScanPoints); ree;

    delete m_turboXICMS1;
    delete m_msFrameMS1;

    QMap<ScanNumber, ScanPoints*> ms1FramePtrs;
    for (auto it = m_ms1ScanNumberVsScanPoints.begin(); it != m_ms1ScanNumberVsScanPoints.end(); ++it) {
        ms1FramePtrs.insert(it.key(), &it.value());
    }

    m_msFrameMS1 = new MsFrame;
    e = m_msFrameMS1->init(ms1FramePtrs, m_msReaderPointerAcc->ptr->getScanNumberVsScanTime()); ree;

    m_turboXICMS1 = new TurboXIC();
    e = m_turboXICMS1->init(m_msFrameMS1->frameIndexVsScanPoints()); ree;

    ERR_RETURN
}

namespace {

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

    constexpr float TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0 = 0.03f;
    constexpr float TIMS_LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0 = 0.15f;
    constexpr int TIMS_MAIN_TOP_N_MS2_IONS = 12;
    constexpr int TIMS_MAIN_EVIDENCE_TOP_N_MS2_IONS = 4;
    constexpr float TIMS_MAIN_EVIDENCE_BUDGET_FRACTION = 0.75f;
    constexpr int TIMS_MAIN_CANDIDATE_BUDGET_BIN_COUNT = 24;

    struct TimsEvidenceCandidate {
        TargetDecoyCandidatePair *targetDecoyPair = nullptr;
        double evidenceScore = 0.0;
    };

    bool isLibraryIonMobilityInAcquisitionWindow(
        const TargetDecoyCandidatePair *candidate,
        const MsScanInfo &msScanInfo
        ) {

        if (candidate == nullptr) {
            return false;
        }

        const float libraryIonMobility = candidate->iIM();
        if (libraryIonMobility <= 0.0f || msScanInfo.ionMobilityDriftTime <= 0.0f) {
            return true;
        }

        float windowLower = msScanInfo.ionMobilityWindowLower;
        float windowUpper = msScanInfo.ionMobilityWindowUpper;
        if (windowLower <= 0.0f || windowUpper <= 0.0f || windowUpper < windowLower) {
            windowLower = msScanInfo.ionMobilityDriftTime - TIMS_LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0;
            windowUpper = msScanInfo.ionMobilityDriftTime + TIMS_LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0;
        }

        windowLower -= TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0;
        windowUpper += TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0;
        return windowLower <= libraryIonMobility && libraryIonMobility <= windowUpper;
    }

    void sortByIonMobilityCenterDistance(
        float ionMobilityCenter,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
        ) {

        std::sort(
            targetDecoyPointers->begin(),
            targetDecoyPointers->end(),
            [ionMobilityCenter](const TargetDecoyCandidatePair *left, const TargetDecoyCandidatePair *right) {
                return std::abs(left->iIM() - ionMobilityCenter)
                     < std::abs(right->iIM() - ionMobilityCenter);
            }
            );
    }

    int clampedBinIndex(float value, float lower, float width, int binCount) {

        if (binCount <= 1 || width <= 0.0f) {
            return 0;
        }

        int binIndex = static_cast<int>(
            ((value - lower) / width) * static_cast<float>(binCount)
            );
        return std::max(0, std::min(binCount - 1, binIndex));
    }

    void applyTimsMainIonMobilityCandidateBudget(
        const MsScanInfo &msScanInfo,
        const PythiaParameters &pythiaParameters,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
        ) {

        const int candidateBudget = pythiaParameters.timsMainCandidateBudgetPerTargetKey;
        if (candidateBudget <= 0
            || targetDecoyPointers->size() <= candidateBudget) {
            return;
        }

        const float ionMobilityCenter = msScanInfo.ionMobilityDriftTime;
        if (!pythiaParameters.timsStratifyCandidateBudget
            || ionMobilityCenter <= 0.0f
            || msScanInfo.ionMobilityWindowLower <= 0.0f
            || msScanInfo.ionMobilityWindowUpper <= msScanInfo.ionMobilityWindowLower) {

            sortByIonMobilityCenterDistance(ionMobilityCenter, targetDecoyPointers);
            targetDecoyPointers->resize(candidateBudget);
            return;
        }

        const float windowLower = msScanInfo.ionMobilityWindowLower - TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0;
        const float windowUpper = msScanInfo.ionMobilityWindowUpper + TIMS_LIBRARY_IM_FILTER_PAD_ONE_OVER_K0;
        const float windowWidth = windowUpper - windowLower;
        if (windowWidth <= 0.0f) {
            sortByIonMobilityCenterDistance(ionMobilityCenter, targetDecoyPointers);
            targetDecoyPointers->resize(candidateBudget);
            return;
        }

        QVector<QVector<TargetDecoyCandidatePair*>> bins(TIMS_MAIN_CANDIDATE_BUDGET_BIN_COUNT);
        for (TargetDecoyCandidatePair *tdcp : *targetDecoyPointers) {
            if (tdcp == nullptr) {
                continue;
            }
            const int binIndex = clampedBinIndex(
                tdcp->iIM(),
                windowLower,
                windowWidth,
                TIMS_MAIN_CANDIDATE_BUDGET_BIN_COUNT
                );
            bins[binIndex].push_back(tdcp);
        }

        int nonEmptyBins = 0;
        for (const QVector<TargetDecoyCandidatePair*> &bin : bins) {
            if (!bin.isEmpty()) {
                ++nonEmptyBins;
            }
        }
        if (nonEmptyBins == 0) {
            sortByIonMobilityCenterDistance(ionMobilityCenter, targetDecoyPointers);
            targetDecoyPointers->resize(candidateBudget);
            return;
        }

        const int quotaPerNonEmptyBin = std::max(1, candidateBudget / nonEmptyBins);
        QVector<TargetDecoyCandidatePair*> selected;
        selected.reserve(candidateBudget);
        QSet<TargetDecoyCandidatePair*> selectedSet;
        selectedSet.reserve(candidateBudget);

        const float precursorTargetMz = msScanInfo.precursorTargetMz;
        for (QVector<TargetDecoyCandidatePair*> &bin : bins) {
            if (bin.isEmpty()) {
                continue;
            }
            std::sort(
                bin.begin(),
                bin.end(),
                [precursorTargetMz](const TargetDecoyCandidatePair *left, const TargetDecoyCandidatePair *right) {
                    const float leftMzDelta = std::abs(left->mz(false) - precursorTargetMz);
                    const float rightMzDelta = std::abs(right->mz(false) - precursorTargetMz);
                    if (MathUtils::tSame(leftMzDelta, rightMzDelta, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                        return left->totalFragmentCount() > right->totalFragmentCount();
                    }
                    return leftMzDelta < rightMzDelta;
                }
                );

            const int takeCount = std::min(quotaPerNonEmptyBin, bin.size());
            for (int i = 0; i < takeCount && selected.size() < candidateBudget; ++i) {
                selected.push_back(bin.at(i));
                selectedSet.insert(bin.at(i));
            }
        }

        if (selected.size() < candidateBudget) {
            sortByIonMobilityCenterDistance(ionMobilityCenter, targetDecoyPointers);
            for (TargetDecoyCandidatePair *tdcp : *targetDecoyPointers) {
                if (selected.size() >= candidateBudget) {
                    break;
                }
                if (selectedSet.contains(tdcp)) {
                    continue;
                }
                selected.push_back(tdcp);
                selectedSet.insert(tdcp);
            }
        }

        *targetDecoyPointers = selected;
    }

    void applyTimsMainCandidateBudget(
        const MsScanInfo &msScanInfo,
        const PythiaParameters &pythiaParameters,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPointers
        ) {

        applyTimsMainIonMobilityCandidateBudget(
            msScanInfo,
            pythiaParameters,
            targetDecoyPointers
            );
    }

    double maxIntensityInFrameWindow(
        const TurboXIC &turboXic,
        const MS2Ion &ms2Ion,
        float ppmTolerance,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax
        ) {

        if (ms2Ion.mz <= 0.0f || ppmTolerance <= 0.0f) {
            return 0.0;
        }

        if (frameIndexMax < frameIndexMin) {
            std::swap(frameIndexMin, frameIndexMax);
        }

        const float massTol = MathUtils::calculatePPM(ms2Ion.mz, ppmTolerance);
        const XICPoints xicPoints = turboXic.extractPointsXIC(ms2Ion.mz - massTol, ms2Ion.mz + massTol);

        double maxIntensity = 0.0;
        for (const XICPoint &xicPoint : xicPoints) {
            if (xicPoint.scanNumber < frameIndexMin || xicPoint.scanNumber > frameIndexMax) {
                continue;
            }
            if (xicPoint.intensity <= maxIntensity) {
                continue;
            }
            maxIntensity = static_cast<double>(xicPoint.intensity);
        }

        return maxIntensity;
    }

    bool predictedFrameWindow(
        const MsCalibratomatic &msCalibratomatic,
        const PythiaParameters &pythiaParameters,
        const MsFrame &msFrameMzTarget,
        const TargetDecoyCandidatePair *targetDecoyPair,
        FrameIndex *frameIndexMin,
        FrameIndex *frameIndexMax
        ) {

        if (targetDecoyPair == nullptr
            || !msCalibratomatic.isInitRT()
            || !msFrameMzTarget.isValid()) {
            return false;
        }

        float predictedScanTime = -1.0f;
        Err e = msCalibratomatic.predictScanTime(targetDecoyPair->iRt(false), &predictedScanTime);
        if (e != eNoError) {
            return false;
        }

        const float scanTimeWindow
            = msCalibratomatic.scanTimeStDev(static_cast<float>(pythiaParameters.scanTimeWindowStDevs));

        e = msFrameMzTarget.frameIndexFromScanTime(predictedScanTime - scanTimeWindow, frameIndexMin);
        if (e != eNoError) {
            return false;
        }

        e = msFrameMzTarget.frameIndexFromScanTime(predictedScanTime + scanTimeWindow, frameIndexMax);
        if (e != eNoError) {
            return false;
        }

        if (*frameIndexMax < *frameIndexMin) {
            std::swap(*frameIndexMin, *frameIndexMax);
        }

        return true;
    }

    double calculateCheapTimsEvidenceScore(
        const TurboXIC &turboXic,
        const MsCalibratomatic &msCalibratomatic,
        const PythiaParameters &pythiaParameters,
        const MsFrame &msFrameMzTarget,
        const TargetDecoyCandidatePair *targetDecoyPair
        ) {

        if (targetDecoyPair == nullptr) {
            return 0.0;
        }

        FrameIndex frameIndexMin = 0;
        FrameIndex frameIndexMax = 0;
        if (!predictedFrameWindow(
                msCalibratomatic,
                pythiaParameters,
                msFrameMzTarget,
                targetDecoyPair,
                &frameIndexMin,
                &frameIndexMax
                )) {
            return 0.0;
        }

        QVector<MS2Ion> ms2Ions = targetDecoyPair->ms2IonsTarget();
        if (ms2Ions.isEmpty()) {
            return 0.0;
        }

        MS2Ion::sortMS2IonsIntensityDesc(&ms2Ions);

        const int ionCount = std::min(TIMS_MAIN_EVIDENCE_TOP_N_MS2_IONS, ms2Ions.size());
        double weightedSignal = 0.0;
        int matchedIonCount = 0;
        for (int ionIndex = 0; ionIndex < ionCount; ++ionIndex) {
            const MS2Ion &ms2Ion = ms2Ions.at(ionIndex);
            const double maxIntensity = maxIntensityInFrameWindow(
                turboXic,
                ms2Ion,
                static_cast<float>(pythiaParameters.ms2ExtractionWidthPPM),
                frameIndexMin,
                frameIndexMax
                );
            if (maxIntensity <= 0.0) {
                continue;
            }

            const double libraryWeight = std::sqrt(std::max(1.0, static_cast<double>(ms2Ion.intensity)));
            const double signalWeight = std::log1p(maxIntensity) * libraryWeight;
            weightedSignal += signalWeight;
            ++matchedIonCount;
        }

        if (matchedIonCount == 0) {
            return 0.0;
        }

        return weightedSignal * static_cast<double>(matchedIonCount);
    }

    void applyTimsMainEvidenceCandidateBudget(
        const MsScanInfo &msScanInfo,
        const MsCalibratomatic &msCalibratomatic,
        const PythiaParameters &pythiaParameters,
        const MsFrame &msFrameMzTarget,
        const TurboXIC &turboXicMS2,
        QVector<TargetDecoyCandidatePair*> *targetDecoyPointers,
        int *evidenceSelectedCount,
        int *fallbackSelectedCount
        ) {

        if (evidenceSelectedCount != nullptr) {
            *evidenceSelectedCount = 0;
        }
        if (fallbackSelectedCount != nullptr) {
            *fallbackSelectedCount = 0;
        }

        const int candidateBudget = pythiaParameters.timsMainCandidateBudgetPerTargetKey;
        if (candidateBudget <= 0
            || targetDecoyPointers->size() <= candidateBudget
            || !msCalibratomatic.isInitRT()
            || !msFrameMzTarget.isValid()) {
            applyTimsMainCandidateBudget(msScanInfo, pythiaParameters, targetDecoyPointers);
            return;
        }

        QVector<TimsEvidenceCandidate> rankedCandidates;
        rankedCandidates.reserve(targetDecoyPointers->size());
        for (TargetDecoyCandidatePair *tdcp : *targetDecoyPointers) {
            const double evidenceScore = calculateCheapTimsEvidenceScore(
                turboXicMS2,
                msCalibratomatic,
                pythiaParameters,
                msFrameMzTarget,
                tdcp
                );
            rankedCandidates.push_back({tdcp, evidenceScore});
        }

        const auto sortByEvidence = [](const TimsEvidenceCandidate &left, const TimsEvidenceCandidate &right) {
            if (MathUtils::tSame(left.evidenceScore, right.evidenceScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                return left.targetDecoyPair->totalFragmentCount() > right.targetDecoyPair->totalFragmentCount();
            }
            return left.evidenceScore > right.evidenceScore;
        };

        std::sort(rankedCandidates.begin(), rankedCandidates.end(), sortByEvidence);

        const int evidenceBudget = std::max(
            0,
            std::min(
                candidateBudget,
                static_cast<int>(std::round(static_cast<float>(candidateBudget) * TIMS_MAIN_EVIDENCE_BUDGET_FRACTION))
                )
            );

        QVector<TargetDecoyCandidatePair*> selected;
        selected.reserve(candidateBudget);
        QSet<TargetDecoyCandidatePair*> selectedSet;
        selectedSet.reserve(candidateBudget);

        for (const TimsEvidenceCandidate &rankedCandidate : rankedCandidates) {
            if (selected.size() >= evidenceBudget) {
                break;
            }
            if (rankedCandidate.evidenceScore <= 0.0) {
                break;
            }
            selected.push_back(rankedCandidate.targetDecoyPair);
            selectedSet.insert(rankedCandidate.targetDecoyPair);
            if (evidenceSelectedCount != nullptr) {
                ++(*evidenceSelectedCount);
            }
        }

        QVector<TargetDecoyCandidatePair*> fallbackCandidates = *targetDecoyPointers;
        applyTimsMainCandidateBudget(msScanInfo, pythiaParameters, &fallbackCandidates);
        for (TargetDecoyCandidatePair *tdcp : fallbackCandidates) {
            if (selected.size() >= candidateBudget) {
                break;
            }
            if (selectedSet.contains(tdcp)) {
                continue;
            }
            selected.push_back(tdcp);
            selectedSet.insert(tdcp);
            if (fallbackSelectedCount != nullptr) {
                ++(*fallbackSelectedCount);
            }
        }

        for (const TimsEvidenceCandidate &rankedCandidate : rankedCandidates) {
            if (selected.size() >= candidateBudget) {
                break;
            }
            if (selectedSet.contains(rankedCandidate.targetDecoyPair)) {
                continue;
            }
            selected.push_back(rankedCandidate.targetDecoyPair);
            selectedSet.insert(rankedCandidate.targetDecoyPair);
        }

        *targetDecoyPointers = selected;
    }


    QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> parallelScoreLogic(
            const QVector<TargetDecoyPairParallelInput> &inputs
            ) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        e = ErrorUtils::isNotEmpty(inputs); rree;

        QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> outputs;

        for (const TargetDecoyPairParallelInput &pi : inputs) {

        	if (pi.turboXicMS1) {
        		e = ErrorUtils::isTrue(pi.turboXicMS1->isInit()); rree;
        	}

            QVector<TargetDecoyCandidatePair*> targetDecoyPointers = pi.targetDecoyPointers;
            bool builtTargetDecoyPointersFromAllCandidates = false;
            int allCandidatePointerCount = 0;

            if (targetDecoyPointers.isEmpty()) {
                if (pi.targetDecoyCandidatePointersAllPntr == nullptr) {
                    continue;
                }

                const float mzMin
                    = pi.msScanInfo.precursorTargetMz - (pi.msScanInfo.isoWindowLower + pi.pythiaParameters.precursorExtractionWindowThomsons);

                const float mzMax
                    = pi.msScanInfo.precursorTargetMz + (pi.msScanInfo.isoWindowUpper + pi.pythiaParameters.precursorExtractionWindowThomsons);

                const bool useIonMobilityFilter = pi.msScanInfo.ionMobilityDriftTime > 0.0f;
                const auto terminatorLogic = [
                    mzMin,
                    mzMax,
                    useIonMobilityFilter,
                    &pi
                ](const TargetDecoyCandidatePair *tdcp) {
                    const float mzPrecursorTargetDecoyPair = tdcp->mz(false);
                    if (!(mzMin <= mzPrecursorTargetDecoyPair && mzPrecursorTargetDecoyPair <= mzMax)) {
                        return true;
                    }
                    if (useIonMobilityFilter && !isLibraryIonMobilityInAcquisitionWindow(tdcp, pi.msScanInfo)) {
                        return true;
                    }
                    return false;
                };

                targetDecoyPointers = *pi.targetDecoyCandidatePointersAllPntr;
                builtTargetDecoyPointersFromAllCandidates = true;
                allCandidatePointerCount = targetDecoyPointers.size();
                const auto terminator = std::remove_if(
                    targetDecoyPointers.begin(),
                    targetDecoyPointers.end(),
                    terminatorLogic
                    );

                targetDecoyPointers.erase(terminator, targetDecoyPointers.end());
            }

            if (pi.splitMzTargetKey) {
                const int midPoint = targetDecoyPointers.size() / 2;
                targetDecoyPointers = pi.isBottomSplit
                                    ? targetDecoyPointers.mid(0, midPoint)
                                    : targetDecoyPointers.mid(midPoint, targetDecoyPointers.size() - midPoint);
            }

            const bool isTimsReader = pi.msReaderPointerAcc != nullptr
                && !pi.msReaderPointerAcc->ptr.isNull()
                && pi.msReaderPointerAcc->ptr->isTIMS();

            QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
            MsFrame msFrameMzTarget;
            const bool shouldLoadLocalMs2Frame = pi.msFrameMzTarget == nullptr
                && pi.msReaderPointerAcc != nullptr
                && !pi.msReaderPointerAcc->ptr.isNull();

            if (shouldLoadLocalMs2Frame) {

                e = pi.msReaderPointerAcc->ptr->getMzTargetScanPoints(pi.targetKey, &scanNumberVsScanPoints); rtee;

                if (pi.msCalibratomatic.isInitCalMS2()) {
                    e = pi.msCalibratomatic.recalibrateScanPoints(
                        MSLevelEnum::MS2,
                        &scanNumberVsScanPoints
                        ); rtee;
                }

                QMap<ScanNumber, ScanPoints*> scanNumberVsScanPointsPntrs;
                for(auto it = scanNumberVsScanPoints.begin(); it != scanNumberVsScanPoints.end(); it++) {
                    scanNumberVsScanPointsPntrs.insert(it.key(), &it.value());
                }

                e = msFrameMzTarget.init(
                    scanNumberVsScanPointsPntrs,
                    pi.scanNumberVsScanTime
                    ); rtee;
            }

            MsFrame *msFrameMzTargetPntr = msFrameMzTarget.isValid()
                                         ? &msFrameMzTarget
                                         : pi.msFrameMzTarget;

            TurboXIC turboXicMS2Local;
            TurboXIC *turboXicMS2Pntr = pi.turboXicMS2;
            if (turboXicMS2Pntr == nullptr
                && msFrameMzTargetPntr != nullptr
                && msFrameMzTargetPntr->isValid()) {
                e = turboXicMS2Local.init(msFrameMzTargetPntr->frameIndexVsScanPoints()); rree;
                turboXicMS2Pntr = &turboXicMS2Local;
            }

            int preBudgetCandidateCount = targetDecoyPointers.size();
            int evidenceSelectedCandidateCount = 0;
            int fallbackSelectedCandidateCount = 0;
            if (builtTargetDecoyPointersFromAllCandidates
                && isTimsReader
                && pi.msScanInfo.ionMobilityDriftTime > 0.0f
                && targetDecoyPointers.size() > pi.pythiaParameters.timsMainCandidateBudgetPerTargetKey) {

                if (turboXicMS2Pntr != nullptr && turboXicMS2Pntr->isInit()) {
                    applyTimsMainEvidenceCandidateBudget(
                        pi.msScanInfo,
                        pi.msCalibratomatic,
                        pi.pythiaParameters,
                        *msFrameMzTargetPntr,
                        *turboXicMS2Pntr,
                        &targetDecoyPointers,
                        &evidenceSelectedCandidateCount,
                        &fallbackSelectedCandidateCount
                        );
                }
                else {
                    applyTimsMainCandidateBudget(
                    pi.msScanInfo,
                    pi.pythiaParameters,
                    &targetDecoyPointers
                    );
                }
            }

            if (builtTargetDecoyPointersFromAllCandidates && isTimsReader) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMS main target candidate filter"
                         << "target_key" << pi.targetKey
                         << "all_candidates" << allCandidatePointerCount
                         << "pre_budget_candidates" << preBudgetCandidateCount
                         << "filtered_candidates" << targetDecoyPointers.size()
                         << "budget" << pi.pythiaParameters.timsMainCandidateBudgetPerTargetKey
                         << "stratified_budget" << pi.pythiaParameters.timsStratifyCandidateBudget
                         << "evidence_selected" << evidenceSelectedCandidateCount
                         << "fallback_selected" << fallbackSelectedCandidateCount
                         << "im_center" << pi.msScanInfo.ionMobilityDriftTime
                         << "im_window_lower" << pi.msScanInfo.ionMobilityWindowLower
                         << "im_window_upper" << pi.msScanInfo.ionMobilityWindowUpper;
            }

            if (targetDecoyPointers.isEmpty()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << pi.targetKey << "Target key is empty";
                continue;
            }

            const int scoringTopNMs2Ions = builtTargetDecoyPointersFromAllCandidates && isTimsReader
                ? std::min(pi.topNMs2Ions, TIMS_MAIN_TOP_N_MS2_IONS)
                : pi.topNMs2Ions;

            QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> allCandidateScores;
            allCandidateScores.reserve(targetDecoyPointers.size() * 2);
            quint64 zeroTargetDecoyScorePairs = 0;

            MsCalibratomatic msCalibratomatic = pi.msCalibratomatic;

            XICPeakManager xicPeakManager;
            if (turboXicMS2Pntr != nullptr && turboXicMS2Pntr->isInit()) {
                e = xicPeakManager.init(
                    targetDecoyPointers,
                    scoringTopNMs2Ions,
                    static_cast<float>(pi.pythiaParameters.ms2ExtractionWidthPPM),
                    turboXicMS2Pntr
                    ); rree;
            }
            else {
                e = xicPeakManager.init(
                    msFrameMzTarget.isValid() ? msFrameMzTarget : *pi.msFrameMzTarget,
                    targetDecoyPointers,
                    scoringTopNMs2Ions,
                    static_cast<float>(pi.pythiaParameters.ms2ExtractionWidthPPM)
                    ); rree;
            }

            TimsMs2IonMobilityIndex timsMs2IonMobilityIndex;
            TimsMs2IonMobilityIndex *timsMs2IonMobilityIndexPntr = nullptr;

            const bool hasLibraryIonMobility = isTimsReader
                && std::any_of(
                    targetDecoyPointers.constBegin(),
                    targetDecoyPointers.constEnd(),
                    [](const TargetDecoyCandidatePair *tdcp) {
                        return tdcp != nullptr && tdcp->iIM() > 0.0f;
                    }
                    );
            const bool needsMs2IonMobilityIndex = hasLibraryIonMobility
                && containsMs2IonMobilityFeature(pi.features);

            if (needsMs2IonMobilityIndex
                && msFrameMzTargetPntr != nullptr
                && msFrameMzTargetPntr->isValid()) {

                QElapsedTimer timsMs2IndexTimer;
                timsMs2IndexTimer.start();
                MzTargetKeyVsMs2FrameTIMS *ms2FrameTims
                    = pi.msReaderPointerAcc->ptr->mzTargetKeyVsFrameNumberVsMS2FrameTIMSPntr();
                const QMap<FrameIndex, double> *ionMobilityIndexVsDriftTime
                    = pi.msReaderPointerAcc->ptr->frameIndexVsDriftTimePntr();

                if (ms2FrameTims != nullptr && ionMobilityIndexVsDriftTime != nullptr) {
                    const auto targetIt = ms2FrameTims->constFind(pi.targetKey);
                    if (targetIt != ms2FrameTims->constEnd()) {
                        e = timsMs2IonMobilityIndex.init(
                            targetIt.value(),
                            *msFrameMzTargetPntr,
                            *ionMobilityIndexVsDriftTime
                            ); rree;

                        if (timsMs2IonMobilityIndex.isInit()) {
                            timsMs2IonMobilityIndexPntr = &timsMs2IonMobilityIndex;
                        }
                    }
                }

                if (builtTargetDecoyPointersFromAllCandidates) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMS MS2 mobility index"
                             << "target_key" << pi.targetKey
                             << "needed" << needsMs2IonMobilityIndex
                             << "points" << timsMs2IonMobilityIndex.pointCount()
                             << "msec" << timsMs2IndexTimer.elapsed();
                }
            }

            const float scanTimeRange = pi.scanTimeMinMax.second - pi.scanTimeMinMax.first;
            CandidateScorertron candidateScorertron;
            e = candidateScorertron.init(
                pi.pythiaParameters,
                pi.msCalibratomatic,
                pi.targetKey,
                scoringTopNMs2Ions,
                pi.minPeakCount,
                scanTimeRange,
                pi.averagineTable,
                pi.features,
                pi.useTopNIntegrationsParameter,
                &xicPeakManager,
                msFrameMzTargetPntr,
                pi.turboXicMS1,
                pi.msFrameMS1,
                pi.msReaderPointerAcc,
                timsMs2IonMobilityIndexPntr
                ); rree;
            candidateScorertron.setUseAdaptiveTimsMobilityCentering(
                pi.useAdaptiveTimsMobilityCentering
                );

            QElapsedTimer scoreTimer;
            scoreTimer.start();
            QElapsedTimer scoreProgressTimer;
            scoreProgressTimer.start();
            int processedTargetDecoyPointers = 0;
            for (TargetDecoyCandidatePair* tdcp : targetDecoyPointers) {
                ++processedTargetDecoyPointers;

                QVector<MS2Ion> ms2TargetIons = tdcp->ms2IonsTarget();

                if (ms2TargetIons.isEmpty()) {
                    continue;
                }

                QVector<MS2Ion> ms2DecoyIons;
                if (tdcp->isDecoy()) {
                    // Important: if the library entry was marked as a decoy, we must be
                    // careful to handle it correctly in this scoring step. Specifically,
                    // the "target" ions should be treated as a decoy while the "decoy"
                    // ions should be treated as a target.
                    ms2DecoyIons = ms2TargetIons;
                    ms2TargetIons = tdcp->ms2IonsDecoy();
                } else {
                    ms2DecoyIons = tdcp->ms2IonsDecoy();
                }

                CandidateScores candidateScoresTarget;
                candidateScoresTarget.isDecoy = false;
                e = candidateScorertron.calculateScores(
                        ms2TargetIons,
                        pi.weights,
                        tdcp,
                        &candidateScoresTarget
                        ); rree;

                CandidateScores candidateScoresDecoy;
                candidateScoresDecoy.isDecoy = true;
                e = candidateScorertron.calculateScores(
                        ms2DecoyIons,
                        pi.weights,
                        tdcp,
                        &candidateScoresDecoy
                        ); rree;

                if (
                    MathUtils::tZero(candidateScoresTarget.featuresArray[CosineSimSum100])
                    && MathUtils::tZero(candidateScoresDecoy.featuresArray[CosineSimSum100])
                    ) {
                    zeroTargetDecoyScorePairs++;
                    continue;
                }

                allCandidateScores.push_back({candidateScoresTarget, candidateScoresDecoy});

                if (builtTargetDecoyPointersFromAllCandidates
                    && isTimsReader
                    && scoreProgressTimer.elapsed() >= 30000) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMS main target progress"
                             << "target_key" << pi.targetKey
                             << "processed_candidates" << processedTargetDecoyPointers
                             << "filtered_candidates" << targetDecoyPointers.size()
                             << "emitted_pairs" << allCandidateScores.size()
                             << "zero_pairs" << zeroTargetDecoyScorePairs
                             << "elapsed_msec" << scoreTimer.elapsed();
                    scoreProgressTimer.restart();
                }
            }

            if (builtTargetDecoyPointersFromAllCandidates && isTimsReader) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMS main target scored"
                         << "target_key" << pi.targetKey
                         << "filtered_candidates" << targetDecoyPointers.size()
                         << "emitted_pairs" << allCandidateScores.size()
                         << "zero_pairs" << zeroTargetDecoyScorePairs
                         << "msec" << scoreTimer.elapsed();
            }

            if (pi.pythiaParameters.writeFullCandidateDebug) {
                candidateScorertron.printScoringDiagnosticsIfEnabled();
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "Radiant candidate scoring diagnostics target_key="
                         << pi.targetKey
                         << "zero_target_decoy_score_pairs="
                         << zeroTargetDecoyScorePairs
                         << "emitted_pairs="
                         << allCandidateScores.size();
            }

            if (pi.pythiaParameters.verbosity > 0) {
                qDebug() << "Target key processed in" << pi.targetKey << et.restart() << "mSec";
            }

            outputs.push_back({e, allCandidateScores});
        }

        return outputs;
    }

}//namespace
Err TargetDecoyCandidatePairScoretron2::scoreTargetDecoyPairs(
        const QVector<Features> &features,
        int topNMS2Ions,
        const MsCalibratomatic &msCalibratomatic,
        float minPeakCount,
        int threadCount,
        bool useTopNIntegrationsParameter,
        const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXicPntrs,
        const QVector<float> &weights,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
        QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> *candidateScoresPairsVec
        ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isFalse(m_diaTargetFrames.isEmpty()); ree;
    e = ErrorUtils::isNotEmpty(m_scanNumberVsScanTime); ree;
    e = ErrorUtils::isNotEmpty(features); ree;

    candidateScoresPairsVec->clear();

    QVector<TargetDecoyPairParallelInput> parallelInputs;
    e = buildParallelInput(
            features,
            topNMS2Ions,
            m_scanTimeMinMax,
            msCalibratomatic,
            minPeakCount,
            useTopNIntegrationsParameter,
            mzTargetKeyVsTurboXicPntrs,
            weights,
            mzTargetKeyVsTargetDecoyCandidatePointers,
            &parallelInputs
            ); ree;

    QVector<QVector<TargetDecoyPairParallelInput>> parallelInputsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            parallelInputs,
            threadCount,
            &parallelInputsTranched
            ); ree;

    const auto filterEmptyTranchesTerminatorLogic = [](const QVector<TargetDecoyPairParallelInput> &input) {
        return input.isEmpty();
    };
    const auto terminator = std::remove_if(
        parallelInputsTranched.begin(),
        parallelInputsTranched.end(),
        filterEmptyTranchesTerminatorLogic
        );
    parallelInputsTranched.erase(terminator, parallelInputsTranched.end());

    e = ErrorUtils::isNotEmpty(parallelInputsTranched); ree;

#define PARALLEL_SCORE
#ifdef PARALLEL_SCORE
    QFuture<QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>>> futures = QtConcurrent::mapped(
            parallelInputsTranched,
            parallelScoreLogic
            );
    futures.waitForFinished();

    for (const QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> &results : futures) {

        const QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> &result = results;
        for (const QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>> &r : result) {
            e = r.first; ree;
            const QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> &candidateScoresTargetMz = r.second;
            candidateScoresPairsVec->append(candidateScoresTargetMz);
        }
    }
#else
    for(const QVector<TargetDecoyPairParallelInput> &tdppis : parallelInputsTranched) {

        const QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> results = parallelScoreLogic(
                tdppis
                ); ree;

        for (const QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>> &res : results) {
            e = res.first; ree;
            candidateScoresPairsVec->append(res.second);
        }
    }
#endif

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron2::scoreTargetDecoyPairs(
        const QVector<Features> &features,
        int topNMS2Ions,
        const MsCalibratomatic &msCalibratomatic,
        float minPeakCount,
        int threadCount,
        bool useTopNIntegrationsParameter,
        const QVector<MsScanInfo> &msScanInfos,
        const QVector<float> &weights,
        QVector<TargetDecoyCandidatePair*> *targetDecoyCandidateAllPntrs,
        QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> *candidateScoresPairsVec
        ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isFalse(m_diaTargetFrames.isEmpty()); ree;
    e = ErrorUtils::isNotEmpty(m_scanNumberVsScanTime); ree;
    e = ErrorUtils::isNotEmpty(features); ree;

    candidateScoresPairsVec->clear();

    QVector<TargetDecoyPairParallelInput> parallelInputs;
    e = buildParallelInput(
            features,
            topNMS2Ions,
            m_scanTimeMinMax,
            msCalibratomatic,
            minPeakCount,
            useTopNIntegrationsParameter,
            msScanInfos,
            weights,
            targetDecoyCandidateAllPntrs,
            &parallelInputs
            ); ree;

    QVector<QVector<TargetDecoyPairParallelInput>> parallelInputsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            parallelInputs,
            threadCount,
            &parallelInputsTranched
            ); ree;

    const auto filterEmptyTranchesTerminatorLogic = [](const QVector<TargetDecoyPairParallelInput> &input) {
        return input.isEmpty();
    };
    const auto terminator = std::remove_if(
        parallelInputsTranched.begin(),
        parallelInputsTranched.end(),
        filterEmptyTranchesTerminatorLogic
        );
    parallelInputsTranched.erase(terminator, parallelInputsTranched.end());

    e = ErrorUtils::isNotEmpty(parallelInputsTranched); ree;

#define PARALLEL_SCORE2
#ifdef PARALLEL_SCORE2
    QFuture<QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>>> futures = QtConcurrent::mapped(
            parallelInputsTranched,
            parallelScoreLogic
            );
    futures.waitForFinished();

    for (const QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> &results : futures) {

        const QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> &result = results;
        for (const QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>> &r : result) {
            e = r.first; ree;
            QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> candidateScoresTargetMz = r.second;
            candidateScoresPairsVec->append(candidateScoresTargetMz);
        }
    }
#else
    for(const QVector<TargetDecoyPairParallelInput> &tdppis : parallelInputsTranched) {

        const QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> results = parallelScoreLogic(
                tdppis
                ); ree;

        for (const QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>> &res : results) {
            e = res.first; ree;
            candidateScoresPairsVec->append(res.second);
        }
    }
#endif

    ERR_RETURN
}

bool TargetDecoyCandidatePairScoretron2::isInit() const {
    return m_pythiaParameters.isValid();
}

void TargetDecoyCandidatePairScoretron2::setUseAdaptiveTimsMobilityCentering(
    bool useAdaptiveTimsMobilityCentering
    ) {
    m_useAdaptiveTimsMobilityCentering = useAdaptiveTimsMobilityCentering;
}

Err TargetDecoyCandidatePairScoretron2::buildParallelInput(
        const QVector<Features> &features,
        int topNMS2Ions,
        const QPair<double, double> &scanTimeMinMax,
        const MsCalibratomatic &msCalibratomatic,
        float minPeakCount,
        bool useTopNIntegrationsParameter,
        const QMap<MzTargetKey, TurboXIC*> &mzTargetKeyVsTurboXicPntrs,
        const QVector<float> &weights,
        const QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
        QVector<TargetDecoyPairParallelInput> *input
        ) const {

    ERR_INIT

    e = ErrorUtils::isFalse(mzTargetKeyVsTargetDecoyCandidatePointers->isEmpty()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(m_diaTargetFrames); ree;
    // e = ErrorUtils::isNotEmpty(m_ms1ScanNumberVsScanPoints); ree;
    e = ErrorUtils::isNotEmpty(m_mzTargetKeyVsMsFramePntr); ree;
    e = ErrorUtils::isAboveThreshold(minPeakCount, 1.0f, ErrorUtilsParam::ExcludeThreshold); ree;
    e = ErrorUtils::isNotEmpty(features); ree;

    input->reserve(mzTargetKeyVsTargetDecoyCandidatePointers->size());
    const bool avoidTargetKeySplit = m_msReaderPointerAcc != nullptr
        && !m_msReaderPointerAcc->ptr.isNull()
        && m_msReaderPointerAcc->ptr->isTIMS()
        && containsMs2IonMobilityFeature(features);

    for (const MzTargetKey &mzTargetKey : mzTargetKeyVsTargetDecoyCandidatePointers->keys()) {

        const QVector<TargetDecoyCandidatePair*> &tdcpPntrs
                            = mzTargetKeyVsTargetDecoyCandidatePointers->value(mzTargetKey);

        const int bufferOddEvenSize = tdcpPntrs.size() % 2 == 1 ? 1 : 0;

        const int midSize = tdcpPntrs.size() / 2;

        TargetDecoyPairParallelInput tdppi1;
        tdppi1.topNMs2Ions = topNMS2Ions;
        tdppi1.targetKey = mzTargetKey;
        tdppi1.msCalibratomatic = msCalibratomatic;
        tdppi1.pythiaParameters = m_pythiaParameters;
        tdppi1.targetDecoyPointers = tdcpPntrs.mid(0, midSize);
        tdppi1.scanTimeMinMax = scanTimeMinMax;
        tdppi1.diaTargetFrame = m_diaTargetFrames.value(tdppi1.targetKey);
        tdppi1.turboXicMS1 = m_turboXICMS1;
        tdppi1.minPeakCount = minPeakCount;
        tdppi1.averagineTable = m_averagineTable;
        tdppi1.msFrameMS1 = m_msFrameMS1;
        tdppi1.weights = weights;
        tdppi1.features = features;
        tdppi1.useTopNIntegrationsParameter = useTopNIntegrationsParameter;
        tdppi1.useAdaptiveTimsMobilityCentering = m_useAdaptiveTimsMobilityCentering;
        tdppi1.msReaderPointerAcc = m_msReaderPointerAcc;
        tdppi1.scanNumberVsScanTime = m_scanNumberVsScanTime;

        if (!m_msReaderPointerAcc->useLazyLoading()) {
            e = ErrorUtils::contains(tdppi1.targetKey, m_mzTargetKeyVsMsFramePntr); ree;
            tdppi1.msFrameMzTarget = m_mzTargetKeyVsMsFramePntr.value(tdppi1.targetKey);
        }

        if (!mzTargetKeyVsTurboXicPntrs.isEmpty()) {
            e = ErrorUtils::contains(tdppi1.targetKey, mzTargetKeyVsTurboXicPntrs); ree;
            tdppi1.turboXicMS2 = mzTargetKeyVsTurboXicPntrs.value(tdppi1.targetKey);
            e = ErrorUtils::contains(tdppi1.targetKey, m_mzTargetKeyVsMsFramePntr); ree;
            tdppi1.msFrameMzTarget = m_mzTargetKeyVsMsFramePntr.value(tdppi1.targetKey);
        }

        if (avoidTargetKeySplit) {
            tdppi1.targetDecoyPointers = tdcpPntrs;
            input->push_back(tdppi1);
            continue;
        }

        TargetDecoyPairParallelInput tdppi2 = tdppi1;
        tdppi2.targetDecoyPointers = tdcpPntrs.mid(midSize, midSize + bufferOddEvenSize);

        input->push_back(tdppi1);
        input->push_back(tdppi2);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron2::buildParallelInput(
        const QVector<Features> &features,
        int topNMS2Ions,
        const QPair<double, double> &scanTimeMinMax,
        const MsCalibratomatic &msCalibratomatic,
        float minPeakCount,
        bool useTopNIntegrationsParameter,
        const QVector<MsScanInfo> &msScanInfos,
        const QVector<float> &weights,
        QVector<TargetDecoyCandidatePair*> *targetDecoyCandidateAllPntrs,
        QVector<TargetDecoyPairParallelInput> *input
        ) const {

    ERR_INIT

    e = ErrorUtils::isFalse(targetDecoyCandidateAllPntrs->isEmpty()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(m_diaTargetFrames); ree;

	if (m_msFrameMS1) {
		e = ErrorUtils::isTrue(m_msFrameMS1->isValid()); ree;
		e = ErrorUtils::isNotEmpty(m_ms1ScanNumberVsScanPoints); ree;
	}

    e = ErrorUtils::isNotEmpty(m_mzTargetKeyVsMsFramePntr); ree;
    e = ErrorUtils::isNotEmpty(features); ree;
    e = ErrorUtils::isAboveThreshold(minPeakCount, 1.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    const bool avoidTargetKeySplit = m_msReaderPointerAcc != nullptr
        && !m_msReaderPointerAcc->ptr.isNull()
        && m_msReaderPointerAcc->ptr->isTIMS()
        && containsMs2IonMobilityFeature(features);
    const bool splitMzTargetKey = !avoidTargetKeySplit
        && m_pythiaParameters.threadCount > msScanInfos.size();

    for (const MsScanInfo &msi : msScanInfos) {

        TargetDecoyPairParallelInput tdppi1;
        tdppi1.topNMs2Ions = topNMS2Ions;
        tdppi1.targetKey = msi.targetKey();
        tdppi1.msScanInfo = msi;
        tdppi1.msCalibratomatic = msCalibratomatic;
        tdppi1.pythiaParameters = m_pythiaParameters;
        tdppi1.scanTimeMinMax = scanTimeMinMax;
        tdppi1.diaTargetFrame = m_diaTargetFrames.value(tdppi1.targetKey);
        tdppi1.turboXicMS1 = m_turboXICMS1;
        tdppi1.minPeakCount = minPeakCount;
        tdppi1.averagineTable = m_averagineTable;
        tdppi1.msFrameMS1 = m_msFrameMS1;
        tdppi1.weights = weights;
        tdppi1.features = features;
        tdppi1.useTopNIntegrationsParameter = useTopNIntegrationsParameter;
        tdppi1.useAdaptiveTimsMobilityCentering = m_useAdaptiveTimsMobilityCentering;
        tdppi1.msReaderPointerAcc = m_msReaderPointerAcc;
        tdppi1.scanNumberVsScanTime = m_scanNumberVsScanTime;
        tdppi1.targetDecoyCandidatePointersAllPntr = targetDecoyCandidateAllPntrs;
        tdppi1.splitMzTargetKey = splitMzTargetKey;

        if (!m_msReaderPointerAcc->useLazyLoading()) {
            e = ErrorUtils::contains(tdppi1.targetKey, m_mzTargetKeyVsMsFramePntr); ree;
            tdppi1.msFrameMzTarget = m_mzTargetKeyVsMsFramePntr.value(tdppi1.targetKey);
        }

        input->push_back(tdppi1);

        if (splitMzTargetKey) {
            tdppi1.isBottomSplit = true;
            input->push_back(tdppi1);
        }

    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron2::setPythiaParameters(const PythiaParameters &pythiaParameters) {
    ERR_INIT
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_pythiaParameters = pythiaParameters;
    ERR_RETURN
}
