//
// Created by andrewnichols on 8/12/24.
//

#include "QuanTransitionRefinertron.h"

#include "CandidateScores.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"

#include <QtConcurrent/QtConcurrent>

QuanTransitionRefinertron::QuanTransitionRefinertron(
    double ppm,
    int frameIndexBuffer
    )
: m_ppmThreshold(ppm)
, m_frameIndexBuffer(frameIndexBuffer)
{}

namespace {

    QMap<MzTargetKey, QVector<CandidateScores*>> groupCandidateScoresByTargetKey(const QVector<CandidateScores*> &candidateScoresesPntrs) {

        QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs;

        for (CandidateScores *cs : candidateScoresesPntrs) {
            mzTargetKeyVsCandidateScoresPntrs[cs->targetKey].push_back(cs);
        }

        return mzTargetKeyVsCandidateScoresPntrs;
    }

    void filterNonInterferingTransitions(QVector<QPair<MzHashed, int>> * mzHashedVsCountPairs) {

        const auto terminatorLogic = [](const QPair<MzHashed, int> &pr) {
            return pr.second <= 1;
        };

        const auto terminator = std::remove_if(
            mzHashedVsCountPairs->begin(),
            mzHashedVsCountPairs->end(),
            terminatorLogic
            );

        mzHashedVsCountPairs->erase(terminator, mzHashedVsCountPairs->end());

    }

    Err refineTransitionsParallelLogic(
        const QVector<CandidateScores*> &candidateScoresPntrs,
        double ppmThreshold,
        int frameIndexBuffer
        ) {

        ERR_INIT

        QMap<FrameIndex, QVector<CandidateScores*>> frameIndexVsCandidateScoresPntrs;
        for (CandidateScores *cs : candidateScoresPntrs) {
            frameIndexVsCandidateScoresPntrs[cs->frameIndex].push_back(cs);
        }

        for (FrameIndex frameIndex : frameIndexVsCandidateScoresPntrs.keys()) {

            const int frameIndexMin = frameIndex - frameIndexBuffer;
            const int frameIndexMax = frameIndex + frameIndexBuffer;

            QVector<float> mzVals;
            for (int fi = frameIndexMin; fi <= frameIndexMax; fi++) {
                const QVector<CandidateScores*> &candidateScoreses = frameIndexVsCandidateScoresPntrs.value(fi);

                for (CandidateScores *cs : candidateScoreses) {

                    const QVector<MS2Ion> &ms2Ions = cs->isDecoy
                                                   ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
                                                   : cs->targetDecoyCandidatePair->ms2IonsTarget();

                    std::transform(
                        ms2Ions.begin(),
                        ms2Ions.end(),
                        std::back_inserter(mzVals),
                        [](const MS2Ion &ms2Ion){return ms2Ion.mz;}
                        );

                }
            }

            QHash<MzHashed, int> mzHashedVsCount;
            e = QuanTransitionRefinertron::groupTransitionsByMz(
                mzVals,
                ppmThreshold,
                &mzHashedVsCount
                ); ree;

            QVector<QPair<MzHashed, int>> mzHashedVsCountPairs
                    = ParallelUtils::convertHashToVectorPairs(mzHashedVsCount);

            filterNonInterferingTransitions(&mzHashedVsCountPairs);

            if (mzHashedVsCountPairs.isEmpty()) {
                continue;
            }

            QVector<float> mzIntererencesFrame;
            std::transform(
                mzHashedVsCountPairs.begin(),
                mzHashedVsCountPairs.end(),
                std::back_inserter(mzIntererencesFrame),
                [](const QPair<MzHashed, int> &pr) {
                    return MathUtils::unHashDecimal<float>(pr.first, S_GLOBAL_SETTINGS.HASHING_PRECISION);
                });

            const QVector<CandidateScores*> candidateScoresFrameIndex
                        = frameIndexVsCandidateScoresPntrs.value(frameIndex);

            for (CandidateScores *cs : candidateScoresFrameIndex) {

                const QVector<MS2Ion> &ms2Ions = cs->isDecoy
                                               ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
                                               : cs->targetDecoyCandidatePair->ms2IonsTarget();

                for (const MS2Ion &ms2Ion : ms2Ions) {
                    const int closestMzHashedVsCountPairsIndex = MathUtils::closest(mzIntererencesFrame, ms2Ion.mz);
                    const float mzInterenceFound = mzIntererencesFrame.at(closestMzHashedVsCountPairsIndex);

                    const float ppmAccuracy = MathUtils::calculateMassAccuracyPPM(ms2Ion.mz, mzInterenceFound);

                    if (std::abs(ppmAccuracy) > ppmThreshold) {
                        continue;
                    }

                    cs->mzInterferences.push_back(ms2Ion.mz);
                }

            }

        }

        ERR_RETURN
    }

}//namespace
Err QuanTransitionRefinertron::refineTransitions(const QVector<CandidateScores*> &candidateScoresPntrs) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScoresPntrs); ree;

    const QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs
                                                                = groupCandidateScoresByTargetKey(candidateScoresPntrs);

    const auto transitionRefinerBinder = std::bind(
        refineTransitionsParallelLogic,
        std::placeholders::_1,
        m_ppmThreshold,
        m_frameIndexBuffer
        );

    QFuture<Err> futures = QtConcurrent::mapped(
        mzTargetKeyVsCandidateScoresPntrs,
        transitionRefinerBinder
        );
    futures.waitForFinished();


    ERR_RETURN
}

Err QuanTransitionRefinertron::groupTransitionsByMz(
    const QVector<float> &mzVals,
    double ppmThreshold,
    QHash<MzHashed, int> *mzHashedVsCount
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(mzVals); ree;
    e = ErrorUtils::isAboveThreshold(ppmThreshold, 0.0, ErrorUtilsParam::ExcludeThreshold); ree;

    QVector<float> mzValsSorted = mzVals;
    std::sort(mzValsSorted.begin(), mzValsSorted.end());

    QVector<float> currentGrouping;
    for (float mzVal : mzValsSorted) {

        const auto currentMzValMean = static_cast<float>(MathUtils::mean(currentGrouping));
        const float ppmAccuracy = MathUtils::calculateMassAccuracyPPM(currentMzValMean, mzVal);

        if (ppmAccuracy > ppmThreshold && !currentGrouping.isEmpty()) {
            const MzHashed mzHashedMean = MathUtils::hashDecimal(currentMzValMean, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            mzHashedVsCount->insert(mzHashedMean, currentGrouping.size());
            currentGrouping.clear();
        }

        currentGrouping.push_back(mzVal);
    }

    if (!currentGrouping.isEmpty()) {
        const auto currentMzValMean = static_cast<float>(MathUtils::mean(currentGrouping));
        const MzHashed mzHashedMean = MathUtils::hashDecimal(currentMzValMean, S_GLOBAL_SETTINGS.HASHING_PRECISION);
        mzHashedVsCount->insert(mzHashedMean, currentGrouping.size());
    }

    ERR_RETURN
}


