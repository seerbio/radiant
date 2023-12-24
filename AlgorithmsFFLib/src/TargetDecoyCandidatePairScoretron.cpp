//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "CandidateScores.h"
#include "CandidateScorertron.h"
#include "FeatureFinderHillBuilder.h"
#include "MsCalibratomatic.h"
#include "ParallelUtils.h"

#include <QtConcurrent/QtConcurrent>


TargetDecoyCandidatePairScoretron::TargetDecoyCandidatePairScoretron()
: m_diaTargetFrames(nullptr)
, m_msReaderPointerAcc(nullptr)
{}


class TargetDecoyPairParallelInput {

public:
    MzTargetKey targetKey;
    MsCalibratomatic msCalibratomatic;
    QMap<ScanNumber, ScanPoints*> diaTargetFrame;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    int topNMs2Ions = -1.0;
    PythiaParameters pythiaParameters;
    QPair<double, double> scanTimeMinMax;
};

Err TargetDecoyCandidatePairScoretron::init(
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanTimeMS1,
        MsReaderPointerAcc *msReaderPointerAcc,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        TurboXIC *turboXICMS1
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isTrue(turboXICMS1->isInit()); ree;

    m_pythiaParameters = pythiaParameters;
    m_msReaderPointerAcc = msReaderPointerAcc;
    m_turboXICMS1 = turboXICMS1;

    if (!diaTargetFrames->isEmpty()) {
        e = ErrorUtils::isNotEmpty(scanNumberVsScanTimeMS1); ree;
        m_diaTargetFrames = diaTargetFrames;
        m_ms1Frame = scanNumberVsScanTimeMS1;
    }

    ERR_RETURN
}

namespace {

    Err buildMzHashedVsCount(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyPointers,
            QHash<MzHashed, int> *mzHashedVsCount
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyPointers); eee_absorb;

        for (TargetDecoyCandidatePair* tdcp : targetDecoyPointers) {
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsTarget()) {
                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                (*mzHashedVsCount)[mzHashed]++;
            }
            for (const MS2Ion &ms2Ion : tdcp->ms2IonsDecoy()) {
                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                (*mzHashedVsCount)[mzHashed]++;
            }
        }

        ERR_RETURN
    }

    QVector<QPair<Err, QVector<CandidateScores>>> parallelScoreLogic(
            const QVector<TargetDecoyPairParallelInput> &inputs,
            const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
            TurboXIC *turboXICMS1
            ) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        e = ErrorUtils::isNotEmpty(inputs); rree;
        e = ErrorUtils::isTrue(turboXICMS1->isInit()); rree;

        QVector<QPair<Err, QVector<CandidateScores>>> outputs;

        for (const TargetDecoyPairParallelInput &pi : inputs) {

            QVector<CandidateScores> allCandidateScores;

            MsCalibratomatic msCalibratomatic = pi.msCalibratomatic;

            QHash<MzHashed, int> mzHashedVsCount;
            e = buildMzHashedVsCount(pi.targetDecoyPointers, &mzHashedVsCount); rree;

            CandidateScorertron candidateScorertron;
            e = candidateScorertron.init(
                    pi.diaTargetFrame,
                    scanNumberVsScanTime,
                    pi.pythiaParameters,
                    pi.targetKey,
                    pi.scanTimeMinMax,
                    pi.topNMs2Ions,
                    mzHashedVsCount,
                    &msCalibratomatic,
                    turboXICMS1
            ); rree;

            for (TargetDecoyCandidatePair* tdcp : pi.targetDecoyPointers) {

                CandidateScores candidateScoresTarget;
                e = candidateScorertron.calculateScores(
                        tdcp->ms2IonsTarget(),
                        tdcp,
                        &candidateScoresTarget
                ); rree;
                candidateScoresTarget.isDecoy = false;
                allCandidateScores.push_back(candidateScoresTarget);

                CandidateScores candidateScoresDecoy;
                e = candidateScorertron.calculateScores(
                        tdcp->ms2IonsDecoy(),
                        tdcp,
                        &candidateScoresDecoy
                ); rree;
                candidateScoresDecoy.isDecoy = true;
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
Err TargetDecoyCandidatePairScoretron::scoreTargetDecoyPairs(
        int topNMS2Ions,
        const QPair<ScanTime , ScanTime > &scanTimeMinMax,
        const MsCalibratomatic &msCalibratomatic,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
        QVector<CandidateScores> *candidateScoresVec
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isFalse(m_diaTargetFrames->isEmpty()); ree;

    candidateScoresVec->clear();

    QVector<TargetDecoyPairParallelInput> parallelInputs;
    e = buildParallelInput(
            topNMS2Ions,
            scanTimeMinMax,
            msCalibratomatic,
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

    const auto refineHillsLogicBinder = std::bind(
            parallelScoreLogic,
            std::placeholders::_1,
            m_msReaderPointerAcc->ptr->getScanNumberVsScanTime(),
            m_turboXICMS1
    );

    QFuture<QVector<QPair<Err, QVector<CandidateScores>>>> futures = QtConcurrent::mapped(
            parallelInputsTranched,
            refineHillsLogicBinder
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

    //TODO write here or pass.

#else
    for(const QVector<TargetDecoyPairParallelInput> &tdppis : parallelInputsTranched) {

        const QVector<QPair<Err, QVector<CandidateScores>>> results = parallelScoreLogic(
                tdppis,
                m_ms1Frame,
                m_msReaderPointerAcc->ptr->getScanNumberVsScanTime()
                ); ree;

        for (const QPair<Err, QVector<CandidateScores>> &res : results) {
            e = res.first; ree;
            candidateScoresVec->append(res.second);
        }
    }
#endif

    ERR_RETURN
}

bool TargetDecoyCandidatePairScoretron::isInit() {

    return m_pythiaParameters.isValid()
       && !m_diaTargetFrames->isEmpty()
       && !m_ms1Frame.isEmpty();
}

Err TargetDecoyCandidatePairScoretron::buildParallelInput(
        int topNMS2Ions,
        const QPair<double, double> &scanTimeMinMax,
        const MsCalibratomatic &msCalibratomatic,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers,
        QVector<TargetDecoyPairParallelInput> *input
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(mzTargetKeyVsTargetDecoyCandidatePointers->isEmpty()); ree;
    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(m_msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isNotEmpty(m_ms1Frame); ree;


    for (const MsScanInfo &msScanInfo : m_msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos()) {

        TargetDecoyPairParallelInput tdppi;
        tdppi.topNMs2Ions = topNMS2Ions;
        tdppi.targetKey = msScanInfo.targetKey();
        tdppi.msCalibratomatic = msCalibratomatic;
        tdppi.pythiaParameters = m_pythiaParameters;
        tdppi.targetDecoyPointers = mzTargetKeyVsTargetDecoyCandidatePointers->value(tdppi.targetKey);
        tdppi.scanTimeMinMax = scanTimeMinMax;
        tdppi.diaTargetFrame = m_diaTargetFrames->value(msScanInfo.targetKey());

        input->push_back(tdppi);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron::setPythiaParameters(const PythiaParameters &pythiaParameters) {
    ERR_INIT
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_pythiaParameters = pythiaParameters;
    ERR_RETURN
}
