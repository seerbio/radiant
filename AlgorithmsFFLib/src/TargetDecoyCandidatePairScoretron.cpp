//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "CandidateScores.h"
#include "CandidateScorertron.h"
#include "CentroidMs2IonMobilityIndex.h"
#include "IsotopicDistributionBuilder.h"
#include "MsCalibratomatic.h"
#include "Ms2IonMobilityIndexBase.h"
#include "ParallelUtils.h"
#include "XICPeakManager.h"

#include <QtConcurrent/QtConcurrent>
#include <QSharedPointer>

#include <algorithm>
#include <cmath>

class TargetDecoyPairParallelInput;

namespace {

    struct TargetKeyScoringContext {
        QMap<ScanNumber, ScanPoints> ownedScanPoints;
        QSharedPointer<MsFrame> ownedMsFrameMzTarget;
        QSharedPointer<TurboXIC> ownedTurboXicMS2;
        QSharedPointer<CentroidMs2IonMobilityIndex> ownedMs2IonMobilityIndex;
        MsFrame *msFrameMzTarget = nullptr;
        TurboXIC *turboXicMS2 = nullptr;
        Ms2IonMobilityIndexBase *ms2IonMobilityIndex = nullptr;
    };

    bool readerHasIonMobility(const MsReaderPointerAcc *msReaderPointerAcc) {
        return msReaderPointerAcc != nullptr
            && !msReaderPointerAcc->ptr.isNull()
            && msReaderPointerAcc->ptr->hasIonMobility();
    }

    struct Ms2IonMobilityIndexStorage {
        CentroidMs2IonMobilityIndex centroidIndex;
        Ms2IonMobilityIndexBase *indexPntr = nullptr;
    };

    Err ensureTargetScanPointsLoaded(
        const TargetDecoyPairParallelInput &pi,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
        );

    Err initLocalMs2Frame(
        const TargetDecoyPairParallelInput &pi,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints,
        MsFrame *msFrameMzTarget
        );

    Err buildCentroidMs2IonMobilityIndex(
        const TargetDecoyPairParallelInput &pi,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints,
        const MsFrame &msFrameMzTarget,
        Ms2IonMobilityIndexStorage *indexStorage
        );

    Err buildMs2IonMobilityIndex(
        const TargetDecoyPairParallelInput &pi,
        bool needsMs2IonMobilityIndex,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints,
        const MsFrame *msFrameMzTarget,
        Ms2IonMobilityIndexStorage *indexStorage
        );

    bool containsMs2IonMobilityFeature(const QVector<Features> &features);

    Err buildTargetKeyScoringContext(
        const TargetDecoyPairParallelInput &pi,
        TargetKeyScoringContext *context
        );

}


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
    bool useAdaptiveIonMobilityCentering = false;
    MsReaderPointerAcc *msReaderPointerAcc = nullptr;
    QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePointersAllPntr = nullptr;
    bool splitMzTargetKey = false;
    bool isBottomSplit = false;
    const TargetKeyScoringContext *targetKeyContext = nullptr;
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

    Err ensureTargetScanPointsLoaded(
        const TargetDecoyPairParallelInput &pi,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints
        ) {

        ERR_INIT

        if (!scanNumberVsScanPoints->isEmpty()) {
            ERR_RETURN
        }

        e = ErrorUtils::isTrue(pi.msReaderPointerAcc != nullptr, eValueError); ree;
        e = ErrorUtils::isTrue(!pi.msReaderPointerAcc->ptr.isNull(), eValueError); ree;

        e = pi.msReaderPointerAcc->ptr->getMzTargetScanPoints(pi.targetKey, scanNumberVsScanPoints); ree;

        if (pi.msCalibratomatic.isInitCalMS2()) {
            e = pi.msCalibratomatic.recalibrateScanPoints(
                MSLevelEnum::MS2,
                scanNumberVsScanPoints
                ); ree;
        }

        ERR_RETURN
    }

    Err initLocalMs2Frame(
        const TargetDecoyPairParallelInput &pi,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints,
        MsFrame *msFrameMzTarget
        ) {

        ERR_INIT

        if (pi.msFrameMzTarget != nullptr) {
            ERR_RETURN
        }

        e = ensureTargetScanPointsLoaded(pi, scanNumberVsScanPoints); ree;

        QMap<ScanNumber, ScanPoints*> scanNumberVsScanPointsPntrs;
        for (auto it = scanNumberVsScanPoints->begin(); it != scanNumberVsScanPoints->end(); ++it) {
            scanNumberVsScanPointsPntrs.insert(it.key(), &it.value());
        }

        e = msFrameMzTarget->init(
            scanNumberVsScanPointsPntrs,
            pi.scanNumberVsScanTime
            ); ree;

        ERR_RETURN
    }

    Err buildCentroidMs2IonMobilityIndex(
        const TargetDecoyPairParallelInput &pi,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints,
        const MsFrame &msFrameMzTarget,
        Ms2IonMobilityIndexStorage *indexStorage
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(pi.msReaderPointerAcc != nullptr, eValueError); ree;
        e = ErrorUtils::isTrue(!pi.msReaderPointerAcc->ptr.isNull(), eValueError); ree;

        e = ensureTargetScanPointsLoaded(pi, scanNumberVsScanPoints); ree;

        QMap<ScanNumber, const TimsbukAlignedPointData*> scanNumberVsAlignedPointData;
        e = pi.msReaderPointerAcc->ptr->getMzTargetAlignedPointData(
            pi.targetKey,
            &scanNumberVsAlignedPointData
            ); ree;

        e = indexStorage->centroidIndex.init(
            *scanNumberVsScanPoints,
            scanNumberVsAlignedPointData,
            msFrameMzTarget
            ); ree;

        if (indexStorage->centroidIndex.isInit()) {
            indexStorage->indexPntr = &indexStorage->centroidIndex;
        }

        ERR_RETURN
    }

    Err buildMs2IonMobilityIndex(
        const TargetDecoyPairParallelInput &pi,
        bool needsMs2IonMobilityIndex,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanPoints,
        const MsFrame *msFrameMzTarget,
        Ms2IonMobilityIndexStorage *indexStorage
        ) {

        ERR_INIT

        if (!needsMs2IonMobilityIndex
            || msFrameMzTarget == nullptr
            || !msFrameMzTarget->isValid()) {
            ERR_RETURN
        }

        if (pi.msReaderPointerAcc == nullptr
            || pi.msReaderPointerAcc->ptr.isNull()
            || !pi.msReaderPointerAcc->ptr->hasIonMobility()) {
            ERR_RETURN
        }

        e = buildCentroidMs2IonMobilityIndex(
            pi,
            scanNumberVsScanPoints,
            *msFrameMzTarget,
            indexStorage
            ); ree;

        ERR_RETURN
    }

    Err buildTargetKeyScoringContext(
        const TargetDecoyPairParallelInput &pi,
        TargetKeyScoringContext *context
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(context != nullptr, eValueError); ree;

        context->ownedScanPoints.clear();
        context->ownedMsFrameMzTarget.clear();
        context->ownedTurboXicMS2.clear();
        context->ownedMs2IonMobilityIndex.clear();
        context->msFrameMzTarget = pi.msFrameMzTarget;
        context->turboXicMS2 = pi.turboXicMS2;
        context->ms2IonMobilityIndex = nullptr;

        const bool needsMs2IonMobilityIndex = readerHasIonMobility(pi.msReaderPointerAcc)
            && containsMs2IonMobilityFeature(pi.features);
        const bool needsOwnedScanPoints = context->msFrameMzTarget == nullptr
            || context->turboXicMS2 == nullptr
            || needsMs2IonMobilityIndex;

        if (needsOwnedScanPoints) {
            if (!pi.diaTargetFrame.isEmpty()) {
                for (auto it = pi.diaTargetFrame.constBegin(); it != pi.diaTargetFrame.constEnd(); ++it) {
                    if (it.value() != nullptr) {
                        context->ownedScanPoints.insert(it.key(), *it.value());
                    }
                }
            }
            else if (pi.msReaderPointerAcc != nullptr && !pi.msReaderPointerAcc->ptr.isNull()) {
                e = ensureTargetScanPointsLoaded(
                    pi,
                    &context->ownedScanPoints
                    ); ree;
            }
        }

        if (context->msFrameMzTarget == nullptr) {
            if (context->ownedScanPoints.isEmpty()) {
                ERR_RETURN
            }

            QMap<ScanNumber, ScanPoints*> scanNumberVsScanPointsPntrs;
            for (auto it = context->ownedScanPoints.begin(); it != context->ownedScanPoints.end(); ++it) {
                scanNumberVsScanPointsPntrs.insert(it.key(), &it.value());
            }

            context->ownedMsFrameMzTarget.reset(new MsFrame);
            e = context->ownedMsFrameMzTarget->init(
                scanNumberVsScanPointsPntrs,
                pi.scanNumberVsScanTime
                ); ree;
            context->msFrameMzTarget = context->ownedMsFrameMzTarget.data();
        }

        if (context->turboXicMS2 == nullptr
            && context->msFrameMzTarget != nullptr
            && context->msFrameMzTarget->isValid()) {
            context->ownedTurboXicMS2.reset(new TurboXIC);
            e = context->ownedTurboXicMS2->init(context->msFrameMzTarget->frameIndexVsScanPoints()); ree;
            context->turboXicMS2 = context->ownedTurboXicMS2.data();
        }

        if (needsMs2IonMobilityIndex
            && context->msFrameMzTarget != nullptr
            && context->msFrameMzTarget->isValid()
            && !context->ownedScanPoints.isEmpty()) {
            QMap<ScanNumber, const TimsbukAlignedPointData*> scanNumberVsAlignedPointData;
            e = pi.msReaderPointerAcc->ptr->getMzTargetAlignedPointData(
                pi.targetKey,
                &scanNumberVsAlignedPointData
                ); ree;

            context->ownedMs2IonMobilityIndex.reset(new CentroidMs2IonMobilityIndex);
            e = context->ownedMs2IonMobilityIndex->init(
                context->ownedScanPoints,
                scanNumberVsAlignedPointData,
                *context->msFrameMzTarget
                ); ree;

            if (context->ownedMs2IonMobilityIndex->isInit()) {
                context->ms2IonMobilityIndex = context->ownedMs2IonMobilityIndex.data();
            }
        }

        ERR_RETURN
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

    constexpr float LIBRARY_IM_FILTER_PAD_ONE_OVER_K0 = 0.03f;
    constexpr float LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0 = 0.15f;

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
            windowLower = msScanInfo.ionMobilityDriftTime - LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0;
            windowUpper = msScanInfo.ionMobilityDriftTime + LIBRARY_IM_FILTER_FALLBACK_HALF_WIDTH_ONE_OVER_K0;
        }

        windowLower -= LIBRARY_IM_FILTER_PAD_ONE_OVER_K0;
        windowUpper += LIBRARY_IM_FILTER_PAD_ONE_OVER_K0;
        return windowLower <= libraryIonMobility && libraryIonMobility <= windowUpper;
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

            MsFrame *msFrameMzTargetPntr = pi.msFrameMzTarget;
            TurboXIC *turboXicMS2Pntr = pi.turboXicMS2;
            Ms2IonMobilityIndexBase *ms2IonMobilityIndexPntr = nullptr;
            if (pi.targetKeyContext != nullptr) {
                msFrameMzTargetPntr = pi.targetKeyContext->msFrameMzTarget;
                turboXicMS2Pntr = pi.targetKeyContext->turboXicMS2;
                ms2IonMobilityIndexPntr = pi.targetKeyContext->ms2IonMobilityIndex;
            }

            if (targetDecoyPointers.isEmpty()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << pi.targetKey << "Target key is empty";
                continue;
            }

            const int scoringTopNMs2Ions = pi.topNMs2Ions;

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

            const bool hasReaderIonMobility = readerHasIonMobility(pi.msReaderPointerAcc);
            const bool hasLibraryIonMobility = hasReaderIonMobility
                && std::any_of(
                    targetDecoyPointers.constBegin(),
                    targetDecoyPointers.constEnd(),
                    [](const TargetDecoyCandidatePair *tdcp) {
                        return tdcp != nullptr && tdcp->iIM() > 0.0f;
                    }
                    );
            if (!hasLibraryIonMobility) {
                ms2IonMobilityIndexPntr = nullptr;
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
                ms2IonMobilityIndexPntr
                ); rree;
            candidateScorertron.setUseAdaptiveIonMobilityCentering(
                pi.useAdaptiveIonMobilityCentering
                );

            for (TargetDecoyCandidatePair* tdcp : targetDecoyPointers) {
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

    QMap<MzTargetKey, QSharedPointer<TargetKeyScoringContext>> targetKeyContexts;
    for (TargetDecoyPairParallelInput &parallelInput : parallelInputs) {
        if (!targetKeyContexts.contains(parallelInput.targetKey)) {
            auto context = QSharedPointer<TargetKeyScoringContext>::create();
            e = buildTargetKeyScoringContext(parallelInput, context.data()); ree;
            targetKeyContexts.insert(parallelInput.targetKey, context);
        }

        parallelInput.targetKeyContext = targetKeyContexts.value(parallelInput.targetKey).data();
    }

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

    QMap<MzTargetKey, QSharedPointer<TargetKeyScoringContext>> targetKeyContexts;
    for (TargetDecoyPairParallelInput &parallelInput : parallelInputs) {
        if (!targetKeyContexts.contains(parallelInput.targetKey)) {
            auto context = QSharedPointer<TargetKeyScoringContext>::create();
            e = buildTargetKeyScoringContext(parallelInput, context.data()); ree;
            targetKeyContexts.insert(parallelInput.targetKey, context);
        }

        parallelInput.targetKeyContext = targetKeyContexts.value(parallelInput.targetKey).data();
    }

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

void TargetDecoyCandidatePairScoretron2::setUseAdaptiveIonMobilityCentering(
    bool useAdaptiveIonMobilityCentering
    ) {
    m_useAdaptiveIonMobilityCentering = useAdaptiveIonMobilityCentering;
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
        tdppi1.useAdaptiveIonMobilityCentering = m_useAdaptiveIonMobilityCentering;
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

    const bool splitMzTargetKey = m_pythiaParameters.threadCount > msScanInfos.size();

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
        tdppi1.useAdaptiveIonMobilityCentering = m_useAdaptiveIonMobilityCentering;
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
