//
// Created by anichols on 10/18/23.
//

#include "TargetDecoyCandidatePairScoretron.h"

#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>

class TargetDecoyPairParallelInput {

public:
    UniqueMsInfoScanKey msInfoScanKey;
    QMap<ScanNumber, ScanPoints>* diaTargetFrame;
    QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
    double ppmTol = -1.0;
    int topNMs2Ions = -1.0;
};


Err TargetDecoyCandidatePairScoretron::init(
        const PythiaParameters &pythiaParameters,
        const QVector<MsScanInfo> &msScanInfos,
        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
        TargetDecoyCandidatePairManager *targetDecoyCandidatePairManager
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(msScanInfos); ree;
    e = ErrorUtils::isNotEmpty(*diaTargetFrames); ree;
    e = ErrorUtils::isTrue(targetDecoyCandidatePairManager->isInit()); ree;

    m_pythiaParameters = pythiaParameters;
    m_msScanInfos = msScanInfos;
    m_diaTargetFrames = diaTargetFrames;
    m_targetDecoyCandidatePairManager = targetDecoyCandidatePairManager;

    ERR_RETURN
}

namespace {

    Err extractMS2Ions(
            const QVector<MS2Ion> &ms2Ions,
            int scanNumberMin,
            int scanNumberMax,
            double ppmTol,
            TurboXIC *turboXic
            ) {

        ERR_INIT

        QMap<MzHashed, XICPoints> cachedPoints;

        for (const MS2Ion &ms2Ion : ms2Ions) {

            const double mzTol = MathUtils::calculatePPM(ms2Ion.mz, ppmTol);
            const double mzMin = ms2Ion.mz - mzTol;
            const double mzMax = ms2Ion.mz + mzTol;

            XICPoints xicPoints;

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            if (cachedPoints.contains(mzHashed)) {
                xicPoints = cachedPoints.value(mzHashed);
                continue;
            }

            xicPoints = turboXic->extractPointsXIC(
                    mzMin,
                    mzMax,
                    scanNumberMin,
                    scanNumberMax
            );

            cachedPoints.insert(mzHashed, xicPoints);
        }

        ERR_RETURN
    }

    Err parallelScoreLogic(const TargetDecoyPairParallelInput &pi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        TurboXIC turboXic;
        e = turboXic.init(pi.diaTargetFrame); ree;

        double scanNumberRtreeMin;
        double scanNumberRtreeMax;
        double mzRtreeMin;
        double mzRtreeMax;
        e = turboXic.getRTreeLimits(
                &scanNumberRtreeMin,
                &scanNumberRtreeMax,
                &mzRtreeMin,
                &mzRtreeMax
        ); ree;

        for (const TargetDecoyCandidatePair* ptr : pi.targetDecoyPointers) {

            QVector<MS2Ion> ms2IonsTarget = ptr->ms2IonsTarget();
            const int topNTarget = std::min(pi.topNMs2Ions, ms2IonsTarget.size());
            ms2IonsTarget.resize(topNTarget);

            QVector<MS2Ion> ms2IonsDecoy = ptr->ms2IonsDecoy();
            const int topNDecoy = std::min(pi.topNMs2Ions, ms2IonsDecoy.size());
            ms2IonsDecoy.resize(topNDecoy);
        }

        qDebug() << pi.msInfoScanKey << et.elapsed() << "mSec";

        ERR_RETURN
    }

}//namespace
Err TargetDecoyCandidatePairScoretron::scoreTargetDecoyPairs(
        double randomSelectionFraction,
        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(m_msScanInfos); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;

    QVector<TargetDecoyPairParallelInput> parallelInputs;
    e = buildParallelInput(
            randomSelectionFraction,
            &parallelInputs
            ); ree;

//#define PARALLEL_SCORE
#ifdef PARALLEL_SCORE
    QFuture<Err> futures = QtConcurrent::mapped(
            parallelInputs,
            parallelScoreLogic
    );
    futures.waitForFinished();

    for (Err res : futures) {
        e = res; ree;
    }
#else
    for(const TargetDecoyPairParallelInput &tdppi : parallelInputs) {
        e = parallelScoreLogic(tdppi); ree;
    }
#endif

    for (const TargetDecoyPairParallelInput &tdppi : parallelInputs) {
        scoredTargetDecoyPointers->append(tdppi.targetDecoyPointers);
    }

    ERR_RETURN
}

Err TargetDecoyCandidatePairScoretron::buildParallelInput(
        double randomSelectionFraction,
        QVector<TargetDecoyPairParallelInput> *input
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(m_msScanInfos); ree;
    e = ErrorUtils::isNotEmpty(*m_diaTargetFrames); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager->isInit()); ree;

    for (const MsScanInfo &msScanInfo : m_msScanInfos) {

        const double mzMin = msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower;
        const double mzMax = msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper;

        TargetDecoyPairParallelInput tdppi;
        tdppi.topNMs2Ions = m_pythiaParameters.topNMs2Ions;
        tdppi.ppmTol = m_pythiaParameters.ms2ExtractionWidthPPM;
        tdppi.diaTargetFrame = &(*m_diaTargetFrames)[msScanInfo.targetScanKey()];
        tdppi.msInfoScanKey = msScanInfo.targetScanKey();

        if (randomSelectionFraction < 0) {
            e = m_targetDecoyCandidatePairManager->getTargetDecoyCandidatePairPointers(
                    mzMin,
                    mzMax,
                    &tdppi.targetDecoyPointers
            ); ree;

        } else {
            e = m_targetDecoyCandidatePairManager->getTargetDecoyCandidatePairPointers(
                    mzMin,
                    mzMax,
                    randomSelectionFraction,
                    &tdppi.targetDecoyPointers
            ); ree;

        }

        input->push_back(tdppi);

    }

    ERR_RETURN
}
