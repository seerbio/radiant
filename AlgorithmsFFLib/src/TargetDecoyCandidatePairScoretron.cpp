//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "CandidateScores.h"
#include "CandidateScorertron.h"
#include "FeatureFinderHillBuilder.h"
#include "MsCalibratomatic.h"

#include <QtConcurrent/QtConcurrent>


TargetDecoyCandidatePairScoretron::TargetDecoyCandidatePairScoretron()
: m_diaTargetFrames(nullptr)
, m_msReaderPointerAcc(nullptr)
{}


class TargetDecoyPairParallelInput {

public:
    MzTargetKey targetKey;
    MsCalibratomatic msCalibratomatic;
    QMap<ScanNumber, ScanPoints*> *diaTargetFrame = nullptr;
    QMap<ScanNumber, ScanPoints> ms1Frame;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    int topNMs2Ions = -1.0;
    PythiaParameters pythiaParameters;
    QPair<double, double> scanTimeMinMax;
};

Err TargetDecoyCandidatePairScoretron::init(
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanTimeMS1,
        MsReaderPointerAcc *msReaderPointerAcc,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;

    m_pythiaParameters = pythiaParameters;
    m_msReaderPointerAcc = msReaderPointerAcc;

    if (!diaTargetFrames->isEmpty()) {
        e = ErrorUtils::isNotEmpty(scanNumberVsScanTimeMS1); ree;
        m_diaTargetFrames = diaTargetFrames;
        m_ms1Frame = scanNumberVsScanTimeMS1;
    }

    ERR_RETURN
}

namespace {

    QPair<Err, QVector<CandidateScores>> parallelScoreLogic(const TargetDecoyPairParallelInput &pi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        e = ErrorUtils::isFalse(pi.diaTargetFrame->isEmpty()); rree;
        e = ErrorUtils::isNotEmpty(pi.ms1Frame); rree;

        QVector<CandidateScores> allCandidateScores;

        MsCalibratomatic msCalibratomatic = pi.msCalibratomatic;

        FeatureFinderParameters featureFinderParameters(pi.pythiaParameters);

        QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsMS1;
        QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsMS2;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;

        QMap<ScanNumber, ScanPoints*> scanNumberVsScanPointsMS2Pntrs;
        for (auto it = scanNumberVsScanPointsMS2.begin(); it != scanNumberVsScanPointsMS2.end(); it++) {
            scanNumberVsScanPointsMS2Pntrs.insert(it.key(), &it.value());
        }

        FeatureFinderHillBuilder featureFinderHillBuilderMS2;
        e = featureFinderHillBuilderMS2.init(featureFinderParameters); rree;
        e = featureFinderHillBuilderMS2.buildHills(*pi.diaTargetFrame); rree;

        CandidateScorertron candidateScorertron;
        e = candidateScorertron.init(
                pi.scanNumberVsScanTime,
                pi.ms1Frame,
                pi.pythiaParameters,
                pi.targetKey,
                pi.scanTimeMinMax,
                pi.topNMs2Ions,
                &msCalibratomatic,
                &featureFinderHillBuilderMS2
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
            qDebug() << "Target key processed in" << pi.targetKey << et.elapsed() << "mSec";
        }

        return {e, allCandidateScores};
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

#define PARALLEL_SCORE
#ifdef PARALLEL_SCORE
    QFuture<QPair<Err, QVector<CandidateScores>>> futures = QtConcurrent::mapped(
            parallelInputs,
            parallelScoreLogic
    );
    futures.waitForFinished();

    for (const QPair<Err, QVector<CandidateScores>> &res : futures) {
        e = res.first; ree;
        const QVector<CandidateScores> &candidateScoresTargetMz = res.second;
        candidateScoresVec->append(candidateScoresTargetMz);
    }

    //TODO write here or pass.

#else
    for(const TargetDecoyPairParallelInput &tdppi : parallelInputs) {
        const QPair<Err, QVector<CandidateScores>> res = parallelScoreLogic(tdppi); ree;
        e = res.first; ree;
        candidateScoresVec->append(res.second);
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
        tdppi.ms1Frame = m_ms1Frame;
        tdppi.diaTargetFrame = &(*m_diaTargetFrames)[msScanInfo.targetKey()];
        tdppi.scanNumberVsScanTime = m_msReaderPointerAcc->ptr->getScanNumberVsScanTime();

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
