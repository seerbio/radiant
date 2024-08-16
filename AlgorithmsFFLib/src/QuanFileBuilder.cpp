//
// Created by andrewnichols on 8/15/24.
//

#include "QuanFileBuilder.h"

#include "CandidateScores.h"
#include "MsFrame.h"
#include "ParquetReader.h"
#include "TurboXIC.h"
#include "TurboXIC2D.h"
#include "QuanReader.h"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>

namespace {

    QMap<MzTargetKey, QVector<CandidateScores*>> buildMzTargetKeyVsCandidateScoresPntrs(
        const QVector<CandidateScores*> &candidateScores
        ) {

        QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs;

        for (CandidateScores *cs : candidateScores) {
            mzTargetKeyVsCandidateScoresPntrs[cs->targetKey].push_back(cs);
        }

        return mzTargetKeyVsCandidateScoresPntrs;
    }

    struct ParallelInput {
        MzTargetKey targetKey;
        MsFrame *msFrame = nullptr;
        QVector<CandidateScores*> candidateScores;
        float ppmTol = -1.0;
    };

    Err buildParallelInputs(
        const QMap<MzTargetKey, QVector<CandidateScores*>> &mzTargetKeyVsCandidateScoresPntrs,
        const QMap<MzTargetKey, MsFrame*> &mzTargetKeyVsMsFrame,
        float ppmTol,
        QVector<ParallelInput> *parallelInputs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzTargetKeyVsCandidateScoresPntrs); ree;
        e = ErrorUtils::isNotEmpty(mzTargetKeyVsMsFrame); ree;

        parallelInputs->clear();

        for (const MzTargetKey &mzk : mzTargetKeyVsCandidateScoresPntrs.keys()) {

            e = ErrorUtils::contains(mzk, mzTargetKeyVsMsFrame); ree;

            ParallelInput pi;
            pi.targetKey = mzk;
            pi.candidateScores = mzTargetKeyVsCandidateScoresPntrs.value(pi.targetKey);
            pi.msFrame = mzTargetKeyVsMsFrame.value(pi.targetKey);
            pi.ppmTol = ppmTol;

            parallelInputs->push_back(pi);
        }

        ERR_RETURN
    }

    Err findBestXICPointsKey(
        float mzTransition,
        const QVector<MzHashed> &mzHahsedVec,
        MzHashed *bestKey
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzHahsedVec); ree;

        const MzHashed mzTransitionHashed
                = MathUtils::hashDecimal(mzTransition, S_GLOBAL_SETTINGS.HASHING_PRECISION);

        const int closestIndex = MathUtils::closest(mzHahsedVec, mzTransitionHashed);

        *bestKey = mzHahsedVec.at(closestIndex);

        constexpr int errorCheckDistance = 100;
        e = ErrorUtils::isTrue(std::abs(*bestKey - mzTransitionHashed) < errorCheckDistance); ree;

        ERR_RETURN
    }

    Err convertXICPointsToVectors(
        const XICPoints &xicPoints,
        MsFrame *msFrame,
        QVector<float> *intensityVec,
        QVector<float> *scanTimeVec
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msFrame->isValid()); ree;

        for (const XICPoint& xic : xicPoints) {
            intensityVec->push_back(xic.intensity);
            scanTimeVec->push_back(msFrame->scanTimeFromFrameIndex(xic.scanNumber));
        }

        ERR_RETURN
    }

    Err buildQuanReaderRow(
        const QHash<MzHashed, XICPoints> &mzHashedVsXICPoints,
        MsFrame *msFrame,
        const CandidateScores *cs,
        QuanReaderRow *qr
        ) {

        ERR_INIT

        qr->peptideStringWithMods = cs->targetDecoyCandidatePair->peptideStringWithMods();
        qr->charge = cs->targetDecoyCandidatePair->charge();
        qr->targetKey = cs->targetKey;
        qr->mzSearched1 = cs->featuresArray[CandidateScores::Features::MzSearched1];
        qr->mzSearched2 = cs->featuresArray[CandidateScores::Features::MzSearched2];
        qr->mzSearched3 = cs->featuresArray[CandidateScores::Features::MzSearched3];
        qr->mzSearched4 = cs->featuresArray[CandidateScores::Features::MzSearched4];
        qr->mzSearched5 = cs->featuresArray[CandidateScores::Features::MzSearched5];
        qr->mzSearched6 = cs->featuresArray[CandidateScores::Features::MzSearched6];
        qr->classifierScore = static_cast<float>(cs->classifierScore);
        qr->discScore = static_cast<float>(cs->discriminantScore);
        qr->qValue = static_cast<float>(cs->qValue);
        qr->mzInterferences = cs->mzInterferences;

        constexpr int top6 = 6;
        const QVector<float> mzSearched = cs->featuresArray.mid(CandidateScores::Features::MzSearched1, top6);

        const QVector<MzHashed> mzHashedKeysVec = mzHashedVsXICPoints.keys().toVector();
        for (int i = 0; i < mzSearched.size(); i++) {

            float mzTransition = mzSearched.at(i);

            constexpr float mzMin = 10.0;
            if (mzTransition < mzMin) {
                continue;
            }

            MzHashed mzHashedKey;
            e = findBestXICPointsKey(mzTransition, mzHashedKeysVec, &mzHashedKey); ree;
            e = ErrorUtils::contains(mzHashedKey, mzHashedVsXICPoints); ree;

            const XICPoints &xicPoints = mzHashedVsXICPoints.value(mzHashedKey);

            QVector<float> intensityVec;
            QVector<float> scanTimeVec;
            e = convertXICPointsToVectors(
                xicPoints,
                msFrame,
                &intensityVec,
                &scanTimeVec
            ); ree;

            switch (i){
                case 0:
                    qr->intensityValsMz1 = intensityVec;
                    qr->scanTimesMz1 = scanTimeVec;
                    break;
                case 1:
                    qr->intensityValsMz2 = intensityVec;
                    qr->scanTimesMz2 = scanTimeVec;
                    break;
                case 2:
                    qr->intensityValsMz3 = intensityVec;
                    qr->scanTimesMz3 = scanTimeVec;
                    break;
                case 3:
                    qr->intensityValsMz4 = intensityVec;
                    qr->scanTimesMz4 = scanTimeVec;
                    break;
                case 4:
                    qr->intensityValsMz5 = intensityVec;
                    qr->scanTimesMz5 = scanTimeVec;
                    break;
                case 5:
                    qr->intensityValsMz6 = intensityVec;
                    qr->scanTimesMz6 = scanTimeVec;
                    break;
                default:
                    rrr(eValueError);
            }
        }

        ERR_RETURN
    }

    QPair<Err, QVector<QuanReaderRow>> parallelLogic(const ParallelInput &pi) {

        ERR_INIT
        TurboXIC2D turboXIC;
        e = turboXIC.init(pi.msFrame->frameIndexVsScanPoints());

        QVector<QuanReaderRow> quanReaderRows;

        for (const CandidateScores *cs : pi.candidateScores) {

            e = ErrorUtils::isAboveThreshold(
                cs->frameIndexEnd,
                cs->frameIndexStart,
                ErrorUtilsParam::ExcludeThreshold
                ); rree;

            const int peakLength = cs->frameIndexEnd - cs->frameIndexStart;
            const int bufferLength = std::max(1, static_cast<int>(std::round(peakLength / 2.0)));

            const FrameIndex frameIndexStart = cs->frameIndexStart - bufferLength;
            const FrameIndex frameIndexEnd = cs->frameIndexStart + bufferLength;

            QVector<MS2Ion> ms2Ions = cs->isDecoy
                            ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
                            : cs->targetDecoyCandidatePair->ms2IonsTarget();

            constexpr int top6 = 6;
            ms2Ions.resize(std::min(top6, ms2Ions.size()));

            QHash<MzHashed, XICPoints> mzHashedVsXICPoints;

            for (const MS2Ion &ms2Ion : ms2Ions) {

                const float mzTol = MathUtils::calculatePPM(ms2Ion.mz, pi.ppmTol);
                const float mzStart = ms2Ion.mz - mzTol;
                const float mzEnd = ms2Ion.mz + mzTol;

                const XICPoints xicPoints = turboXIC.extractPointsXIC(
                    mzStart,
                    mzEnd,
                    static_cast<float>(frameIndexStart),
                    static_cast<float>(frameIndexEnd)
                    );

                const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                mzHashedVsXICPoints.insert(mzHashed, xicPoints);
            }

            QuanReaderRow quanReaderRow;
            e = buildQuanReaderRow(
                mzHashedVsXICPoints,
                pi.msFrame,
                cs,
                &quanReaderRow
            ); rree;

            quanReaderRows.push_back(quanReaderRow);
        }

        return {e, quanReaderRows};
    }

}//namespace
Err QuanFileBuilder::buildQuanFile(
    const QVector<CandidateScores*> &candidateScores,
    const QMap<MzTargetKey, MsFrame*> &mzTargetKeyVsMsFrame,
    const QString& outputFilePath,
    float ppmTolerance
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;
    e = ErrorUtils::isNotEmpty(outputFilePath); ree;
    e = ErrorUtils::isNotEmpty(mzTargetKeyVsMsFrame); ree;

    QElapsedTimer et;
    et.start();

    const QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs
                                        = buildMzTargetKeyVsCandidateScoresPntrs(candidateScores);

    QVector<ParallelInput> parallelInputs;
    e = buildParallelInputs(
        mzTargetKeyVsCandidateScoresPntrs,
        mzTargetKeyVsMsFrame,
        ppmTolerance,
        &parallelInputs
        ); ree;

    QVector<QuanReaderRow> quanReaderRows;

#define RUN_PARALLEL_BUILDER
#ifdef RUN_PARALLEL_BUILDER
    QFuture<QPair<Err, QVector<QuanReaderRow>>> futures = QtConcurrent::mapped(
        parallelInputs,
        parallelLogic
        );
    futures.waitForFinished();

    for (const QPair<Err, QVector<QuanReaderRow>> &pr : futures) {
        e = pr.first; ree;
        quanReaderRows.append(pr.second);
    }
#else
    for (const ParallelInput &pi : parallelInputs) {
        QPair<Err, QVector<QuanReaderRow>> errVsQuanReaderRows = parallelLogic(pi);
        e = errVsQuanReaderRows.first;
        quanReaderRows.append(errVsQuanReaderRows.second);
    }
#endif

    e = ParquetReader::write(quanReaderRows, outputFilePath); ree;

    qDebug() << "Quan file built in" << et.elapsed() << "mSec";

    ERR_RETURN
}
