//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "CandidateScores.h"
#include "CandidateScorertron.h"
#include "IsotopicDistributionBuilder.h"
#include "MsCalibratomatic.h"
#include "ParallelUtils.h"
#include "XICPeakManager.h"

#include <QtConcurrent/QtConcurrent>


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

    QMap<ScanNumber, ScanPoints*> ms1FramePtrs;
    for (auto it = m_ms1ScanNumberVsScanPoints.begin(); it != m_ms1ScanNumberVsScanPoints.end(); ++it) {
        ms1FramePtrs.insert(it.key(), &it.value());
    }

    m_msFrameMS1 = new MsFrame;
    e = m_msFrameMS1->init(ms1FramePtrs, m_msReaderPointerAcc->ptr->getScanNumberVsScanTime()); ree;

    m_turboXICMS1 = new TurboXIC();
    e = m_turboXICMS1->init(m_msFrameMS1->frameIndexVsScanPoints()); ree;

    if (!diaTargetFrames.isEmpty()) {
        e = ErrorUtils::isNotEmpty(m_ms1ScanNumberVsScanPoints); ree;
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


    QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> parallelScoreLogic(
            const QVector<TargetDecoyPairParallelInput> &inputs
            ) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        e = ErrorUtils::isNotEmpty(inputs); rree;

        QVector<QPair<Err, QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>>>> outputs;

        for (const TargetDecoyPairParallelInput &pi : inputs) {
            
            e = ErrorUtils::isTrue(pi.turboXicMS1->isInit()); rree;

            QVector<TargetDecoyCandidatePair*> targetDecoyPointers = pi.targetDecoyPointers;

            if (targetDecoyPointers.isEmpty()) {

                const float mzMin
                    = pi.msScanInfo.precursorTargetMz - (pi.msScanInfo.isoWindowLower + pi.pythiaParameters.precursorExtractionWindowThomsons);

                const float mzMax
                    = pi.msScanInfo.precursorTargetMz + (pi.msScanInfo.isoWindowUpper + pi.pythiaParameters.precursorExtractionWindowThomsons);

                const auto terminatorLogic = [mzMin, mzMax](const TargetDecoyCandidatePair *tdcp) {
                    const float mzPrecursorTargetDecoyPair = tdcp->mz(false);
                    return !(mzMin <= mzPrecursorTargetDecoyPair && mzPrecursorTargetDecoyPair <= mzMax);
                };

                targetDecoyPointers = *pi.targetDecoyCandidatePointersAllPntr;
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

            if (targetDecoyPointers.isEmpty()) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << pi.targetKey << "Target key is empty";
                continue;
            }

            QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> allCandidateScores;
            allCandidateScores.reserve(targetDecoyPointers.size() * 2);

            MsCalibratomatic msCalibratomatic = pi.msCalibratomatic;

            QMap<ScanNumber, ScanPoints> scanNumberVsScanPoints;
            MsFrame msFrameMzTarget;
            if (pi.msFrameMzTarget == nullptr) {

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

            XICPeakManager xicPeakManager;
            if (pi.turboXicMS2 != nullptr) {
                e = xicPeakManager.init(
                    targetDecoyPointers,
                    pi.topNMs2Ions,
                    static_cast<float>(pi.pythiaParameters.ms2ExtractionWidthPPM),
                    pi.turboXicMS2
                    ); rree;
            }
            else {
                e = xicPeakManager.init(
                    msFrameMzTarget.isValid() ? msFrameMzTarget : *pi.msFrameMzTarget,
                    targetDecoyPointers,
                    pi.topNMs2Ions,
                    static_cast<float>(pi.pythiaParameters.ms2ExtractionWidthPPM)
                    ); rree;
            }

            const float scanTimeRange = pi.scanTimeMinMax.second - pi.scanTimeMinMax.first;
            CandidateScorertron candidateScorertron;
            e = candidateScorertron.init(
                pi.pythiaParameters,
                pi.msCalibratomatic,
                pi.targetKey,
                pi.topNMs2Ions,
                pi.minPeakCount,
                scanTimeRange,
                pi.averagineTable,
                pi.features,
                pi.useTopNIntegrationsParameter,
                &xicPeakManager,
                msFrameMzTarget.isValid() ? &msFrameMzTarget : pi.msFrameMzTarget,
                pi.turboXicMS1,
                pi.msFrameMS1
                ); rree;

            for (TargetDecoyCandidatePair* tdcp : targetDecoyPointers) {
                // Skip library entries without any ions
                if (tdcp->ms2IonsTarget().isEmpty()) {
                    continue;
                }

                // Score the target candidate for the pair
                CandidateScores candidateScoresTarget;
                candidateScoresTarget.isDecoy = false;
                e = candidateScorertron.calculateScores(
                    pi.weights,
                    tdcp,
                    &candidateScoresTarget
                ); rree;

                // Score the decoy candidate for the pair
                CandidateScores candidateScoresDecoy;
                candidateScoresDecoy.isDecoy = true;
                e = candidateScorertron.calculateScores(
                    pi.weights,
                    tdcp,
                    &candidateScoresDecoy
                ); rree;

                if (
                    MathUtils::tZero(candidateScoresTarget.featuresArray[CosineSimSum100])
                    && MathUtils::tZero(candidateScoresDecoy.featuresArray[CosineSimSum100])
                    ) {
                    continue;
                }

                allCandidateScores.push_back({candidateScoresTarget, candidateScoresDecoy});
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
    return m_pythiaParameters.isValid() && !m_ms1ScanNumberVsScanPoints.isEmpty();
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
    e = ErrorUtils::isTrue(m_msFrameMS1->isValid()); ree;
    e = ErrorUtils::isNotEmpty(m_ms1ScanNumberVsScanPoints); ree;
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
    e = ErrorUtils::isTrue(m_msFrameMS1->isValid()); ree;
    e = ErrorUtils::isNotEmpty(m_ms1ScanNumberVsScanPoints); ree;
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
