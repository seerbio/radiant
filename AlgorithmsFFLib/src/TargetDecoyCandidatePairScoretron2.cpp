//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron2.h"

#include "CandidateScores.h"
#include "CandidateScorertron2.h"
#include "IsotopicDistributionBuilder.h"
#include "MsCalibratomatic.h"
#include "ParallelUtils.h"
#include "XICPeakManager.h"

#include <QtConcurrent/QtConcurrent>


TargetDecoyCandidatePairScoretron2::TargetDecoyCandidatePairScoretron2()
: m_msReaderPointerAcc(nullptr)
, m_turboXICMS1(nullptr)
{}

TargetDecoyCandidatePairScoretron2::~TargetDecoyCandidatePairScoretron2() {
    delete m_turboXICMS1;
    for (MsFrame *msFrame : m_mzTargetKeyVsMsFrame) {
        delete msFrame;
    }
}

class TargetDecoyPairParallelInput {

public:
    MzTargetKey targetKey;
    MsCalibratomatic msCalibratomatic;
    QMap<ScanNumber, ScanPoints*> diaTargetFrame;
    MsFrame *msFrameMzTarget = nullptr;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    int topNMs2Ions = -1.0;
    PythiaParameters pythiaParameters;
    QPair<double, double> scanTimeMinMax;
    TurboXIC *turboXicMS1 = nullptr;
    float minPeakCount = -1.0;
    QMap<int, QVector<float>> averagineTable;
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
    e = m_msReaderPointerAcc->ptr->collateMS2MzTargetFrames(&diaTargetFrames); ree;

    constexpr int msLevel = 1;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsMS1;
    e = m_msReaderPointerAcc->ptr->getScanPoints(msLevel, &scanNumberVsScanPointsMS1); ree;

    QMap<ScanNumber, ScanPoints*> ms1FramePtrs;
    for (auto it = scanNumberVsScanPointsMS1.begin(); it != scanNumberVsScanPointsMS1.end(); it++) {
        ms1FramePtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrameMS1;
    e = msFrameMS1.init(ms1FramePtrs, m_msReaderPointerAcc->ptr->getScanNumberVsScanTime()); ree;

    m_turboXICMS1 = new TurboXIC();
    e = m_turboXICMS1->init(msFrameMS1.frameIndexVsScanPoints()); ree;

    if (!diaTargetFrames.isEmpty()) {
        e = ErrorUtils::isNotEmpty(scanNumberVsScanPointsMS1); ree;
        m_diaTargetFrames = diaTargetFrames;
        m_ms1ScanNumberVsScanPoints = scanNumberVsScanPointsMS1;
        e = buildMzTargetKeyVsMsFrames(); ree;
    }

    e = buildAveragineTable(); ree;

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron2::buildMzTargetKeyVsMsFrames() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_diaTargetFrames); ree;
    e = ErrorUtils::isNotEmpty(m_scanNumberVsScanTime); ree;

    for (auto it = m_diaTargetFrames.begin(); it != m_diaTargetFrames.end(); ++it) {
        auto *msFrame = new MsFrame();
        e = msFrame->init(it.value(), m_scanNumberVsScanTime); ree;
        m_mzTargetKeyVsMsFrame.insert(it.key(), msFrame);
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

namespace {

    Err buildMzHashedVsCount(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyPointers,
            int topNFragIons,
            QMap<MzHashed, int> *mzHashedVsCount
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyPointers); eee_absorb;
        for (TargetDecoyCandidatePair* tdcp : targetDecoyPointers) {

            int counter = 0;
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsTarget()) {

                if (counter++ >= topNFragIons) {
                    break;
                }

                const MzHashed mzHashed = MathUtils::hashDecimal(
                    ms2Ion.mz,
                    S_GLOBAL_SETTINGS.HASHING_PRECISION
                    );
                (*mzHashedVsCount)[mzHashed]++;

                const float isotopeDistanceThomsons = S_GLOBAL_SETTINGS.ISO_DIFF / static_cast<float>(ms2Ion.charge);
                const MzHashed mzHashedShadow = MathUtils::hashDecimal(
                    ms2Ion.mz - isotopeDistanceThomsons,
                    S_GLOBAL_SETTINGS.HASHING_PRECISION
                    );
                (*mzHashedVsCount)[mzHashedShadow]++;
            }

            counter = 0;
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsDecoy()) {

                if (counter++ >= topNFragIons) {
                    break;
                }

                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                (*mzHashedVsCount)[mzHashed]++;

                const float isotopeDistanceThomsons = S_GLOBAL_SETTINGS.ISO_DIFF / static_cast<float>(ms2Ion.charge);
                const MzHashed mzHashedShadow = MathUtils::hashDecimal(
                    ms2Ion.mz - isotopeDistanceThomsons,
                    S_GLOBAL_SETTINGS.HASHING_PRECISION
                    );

                (*mzHashedVsCount)[mzHashedShadow]++;

            }
        }

        ERR_RETURN
    }

    QVector<QPair<Err, QVector<CandidateScores>>> parallelScoreLogic(
            const QVector<TargetDecoyPairParallelInput> &inputs,
            TurboXIC *turboXICMS1
            ) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        e = ErrorUtils::isNotEmpty(inputs); rree;
        e = ErrorUtils::isTrue(turboXICMS1->isInit()); rree;

        QVector<QPair<Err, QVector<CandidateScores>>> outputs;

        for (const TargetDecoyPairParallelInput &pi : inputs) {

            if (pi.targetDecoyPointers.isEmpty()) {
                qDebug() << pi.targetKey << "Target key is empty";
                continue;
            }

            QVector<CandidateScores> allCandidateScores;

            MsCalibratomatic msCalibratomatic = pi.msCalibratomatic;

            QMap<MzHashed, int> mzHashedVsCount;
            e = buildMzHashedVsCount(
                pi.targetDecoyPointers,
                pi.topNMs2Ions,
                &mzHashedVsCount
                ); rree;

            const QList<int> &mzHashedVsCountKeys = mzHashedVsCount.keys();
            QVector<float> mzValsToExtract;
            std::transform(
                mzHashedVsCountKeys.begin(),
                mzHashedVsCountKeys.end(),
                std::back_inserter(mzValsToExtract),
                [](int mzHashed){return MathUtils::unHashDecimal<float>(mzHashed, S_GLOBAL_SETTINGS.HASHING_PRECISION);}
                );

            XICPeakManager xicPeakManager;
            e = xicPeakManager.init(
                *pi.msFrameMzTarget,
                mzValsToExtract,
                static_cast<float>(pi.pythiaParameters.ms2ExtractionWidthPPM)
                ); rree;

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
                &xicPeakManager,
                pi.msFrameMzTarget,
                pi.turboXicMS1
                ); rree;

            for (TargetDecoyCandidatePair* tdcp : pi.targetDecoyPointers) {

                CandidateScores candidateScoresTarget;
                candidateScoresTarget.isDecoy = false;
                e = candidateScorertron.calculateScores(
                        tdcp->ms2IonsTarget(),
                        tdcp,
                        &candidateScoresTarget
                        ); rree;
                allCandidateScores.push_back(candidateScoresTarget);

                CandidateScores candidateScoresDecoy;
                candidateScoresDecoy.isDecoy = true;
                e = candidateScorertron.calculateScores(
                        tdcp->ms2IonsDecoy(),
                        tdcp,
                        &candidateScoresDecoy
                ); rree;
                allCandidateScores.push_back(candidateScoresDecoy);
            }

            if (pi.pythiaParameters.verbosity > 1) {
                qDebug() << "Target key processed in" << pi.targetKey << et.restart() << "mSec";
            }

            outputs.push_back({e, allCandidateScores});
        }

        return outputs;
    }

}//namespace
Err TargetDecoyCandidatePairScoretron2::scoreTargetDecoyPairs(
        int topNMS2Ions,
        const MsCalibratomatic &msCalibratomatic,
        float minPeakCount,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
        QVector<CandidateScores> *candidateScoresVec
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isFalse(m_diaTargetFrames.isEmpty()); ree;
    e = ErrorUtils::isNotEmpty(m_scanNumberVsScanTime); ree;

    candidateScoresVec->clear();

    QVector<TargetDecoyPairParallelInput> parallelInputs;
    e = buildParallelInput(
            topNMS2Ions,
            m_scanTimeMinMax,
            msCalibratomatic,
            minPeakCount,
            mzTargetKeyVsTargetDecoyCandidatePointers,
            &parallelInputs
            ); ree;

    QVector<QVector<TargetDecoyPairParallelInput>> parallelInputsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            parallelInputs,
            ParallelUtils::numberOfAvailableSystemProcessors(),
            &parallelInputsTranched
            ); ree;

#define PARALLEL_SCORE
#ifdef PARALLEL_SCORE
    const auto batchScoreProcessLogicBinder = std::bind(
            parallelScoreLogic,
            std::placeholders::_1,
            m_turboXICMS1
            );

    QFuture<QVector<QPair<Err, QVector<CandidateScores>>>> futures = QtConcurrent::mapped(
            parallelInputsTranched,
            batchScoreProcessLogicBinder
            );
    futures.waitForFinished();

    for (const QVector<QPair<Err, QVector<CandidateScores>>> &results : futures) {

        const QVector<QPair<Err, QVector<CandidateScores>>> &result = results;
        for (const QPair<Err, QVector<CandidateScores>> &r : result) {
            e = r.first; ree;
            const QVector<CandidateScores> &candidateScoresTargetMz = r.second;
            candidateScoresVec->append(candidateScoresTargetMz);
        }
    }
#else
    for(const QVector<TargetDecoyPairParallelInput> &tdppis : parallelInputsTranched) {

        const QVector<QPair<Err, QVector<CandidateScores>>> results = parallelScoreLogic(
                tdppis,
                m_turboXICMS1
                ); ree;

        for (const QPair<Err, QVector<CandidateScores>> &res : results) {
            e = res.first; ree;
            candidateScoresVec->append(res.second);
        }
    }
#endif

    ERR_RETURN
}

bool TargetDecoyCandidatePairScoretron2::isInit() {

    return m_pythiaParameters.isValid()
       && !m_diaTargetFrames.isEmpty()
       && !m_ms1ScanNumberVsScanPoints.isEmpty();
}

Err TargetDecoyCandidatePairScoretron2::buildParallelInput(
        int topNMS2Ions,
        const QPair<double, double> &scanTimeMinMax,
        const MsCalibratomatic &msCalibratomatic,
        float minPeakCount,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
        QVector<TargetDecoyPairParallelInput> *input
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(mzTargetKeyVsTargetDecoyCandidatePointers->isEmpty()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(m_diaTargetFrames); ree;
    e = ErrorUtils::isNotEmpty(m_ms1ScanNumberVsScanPoints); ree;
    e = ErrorUtils::isNotEmpty(m_uniqueTandemMsScanInfos); ree;
    e = ErrorUtils::isNotEmpty(m_mzTargetKeyVsMsFrame); ree;
    e = ErrorUtils::isAboveThreshold(minPeakCount, 1.0f, ErrorUtilsParam::ExcludeThreshold); ree;

    input->reserve(m_uniqueTandemMsScanInfos.size());

    for (const MsScanInfo &msScanInfo : m_uniqueTandemMsScanInfos) {

        TargetDecoyPairParallelInput tdppi;
        tdppi.topNMs2Ions = topNMS2Ions;
        tdppi.targetKey = msScanInfo.targetKey();
        tdppi.msCalibratomatic = msCalibratomatic;
        tdppi.pythiaParameters = m_pythiaParameters;
        tdppi.targetDecoyPointers = mzTargetKeyVsTargetDecoyCandidatePointers->value(tdppi.targetKey);
        tdppi.scanTimeMinMax = scanTimeMinMax;
        tdppi.diaTargetFrame = m_diaTargetFrames.value(msScanInfo.targetKey());
        tdppi.msFrameMzTarget = m_mzTargetKeyVsMsFrame.value(msScanInfo.targetKey());
        tdppi.turboXicMS1 = m_turboXICMS1;
        tdppi.minPeakCount = minPeakCount;
        tdppi.averagineTable = m_averagineTable;

        input->push_back(tdppi);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron2::setPythiaParameters(const PythiaParameters &pythiaParameters) {
    ERR_INIT
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_pythiaParameters = pythiaParameters;
    ERR_RETURN
}
