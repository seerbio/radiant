//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron.h"

#include "ErrorUtils.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

#include "EigenUtils.h"

#include <Eigen/Dense>


CandidateScorertron::CandidateScorertron()
: m_topNMS2Ions(-1)
{}

namespace{

    PeakIntegratomaticParameters buildPeakIntegratomaticParams(const PythiaParameters  &pythiaParameters) {

        PeakIntegratomaticParameters params;
        params.smoothCount = pythiaParameters.smoothCount;
        params.filterLength = pythiaParameters.filterLength;
        params.sigma = pythiaParameters.sigma;
        params.signalToNoiseRatio = pythiaParameters.signalToNoiseRatio;

        return params;
    }

}//namespace
Err CandidateScorertron::init(const PythiaParameters &pythiaParameters, int topNMS2Ions) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(topNMS2Ions > 0); ree;

    m_pythiaParameters = pythiaParameters;
    m_topNMS2Ions = topNMS2Ions;

    const PeakIntegratomaticParameters ffParams = buildPeakIntegratomaticParams(m_pythiaParameters);
    e = m_peakIntegratomatic.init(ffParams); ree;

    ERR_RETURN
}

namespace {

    FrameIndex getMaxFrameIndexFromMzHashedVsXICPoints(const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints) {

        FrameIndex frameIndex = -1;

        const QList<XICPoints> &xicPoints = mzHashedVsXICPoints.values();
        for (const XICPoints &xp : xicPoints) {

            if (xp.scanNumbersVsIntensityVals.isEmpty()) {
                continue;
            }

            frameIndex = std::max(frameIndex, xp.scanNumbersVsIntensityVals.lastKey());
        }

        return frameIndex;
    }

    Err buildMzHashedVsIonPresence(
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            QMap<MzHashed, QVector<double>> *mzHashedVsIonPresence
            ) {

        ERR_INIT

        mzHashedVsIonPresence->clear();

        e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;

        const FrameIndex frameIndexMax = getMaxFrameIndexFromMzHashedVsXICPoints(mzHashedVsXICPoints);
        if (frameIndexMax < 0) {
            ERR_RETURN
        }

        for (auto it = mzHashedVsXICPoints.begin(); it != mzHashedVsXICPoints.end(); it++) {

            const MzHashed mzHashed = it.key();
            const XICPoints &xicPoints = it.value();

            if (xicPoints.scanNumbersVsIntensityVals.isEmpty()) {
                mzHashedVsIonPresence->insert(mzHashed, {});
                continue;
            }

            Eigen::VectorX<double> vecNormalized = EigenUtils::convertQMapToEigenVector(
                    xicPoints.scanNumbersVsIntensityVals,
                    frameIndexMax
            );
            const double denom = vecNormalized.maxCoeff();
            vecNormalized = vecNormalized.array() / denom;

            Eigen::VectorX<double> vecPresence = vecNormalized.array() / vecNormalized.array();
            EigenUtils::replaceNaN(0.0, &vecPresence);
            mzHashedVsIonPresence->insert(mzHashed, EigenUtils::convertEigenVectorToQVector(vecPresence));
        }

        ERR_RETURN
    }

    int getMaxRowCount(
            const QVector<MS2Ion> &ms2Ions,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence
    ) {

        int maxRowCount = 0;
        for (const MS2Ion &ms2Ion : ms2Ions) {

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const QVector<double> &ionPresenceVec = mzHashedVsIonPresence.value(mzHashed);

            maxRowCount = std::max(maxRowCount, ionPresenceVec.size());
        }

        return maxRowCount;
    }

    Eigen::MatrixX<double> buildSummingMatrix(
            const QVector<MS2Ion> &ms2Ions,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
            int topNMs2Ions
    ) {

        const int cols = topNMs2Ions;
        const int rows = getMaxRowCount(
                ms2Ions,
                mzHashedVsIonPresence
        );

        if (rows == 0) {
            return {};
        }

        Eigen::MatrixX<double> mat(rows, cols);
        mat.setZero();

        for (int i = 0; i < cols; i++) {

            if (i >= ms2Ions.size()) {
                break;
            }

            const MS2Ion &ms2Ion = ms2Ions.at(i);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const QVector<double> &ionPresenceVec = mzHashedVsIonPresence.value(mzHashed);

            if (ionPresenceVec.isEmpty()) {
                continue;
            }

            const Eigen::VectorX<double> eVec = EigenUtils::convertQVectorToEigenVector(ionPresenceVec);
            mat.col(i) = eVec;
        }

        return mat;
    }

}//namespace
Err CandidateScorertron::calculateScores(
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
        const QVector<MS2Ion> &ms2Ions,
        ScanTime scanTimePredicted
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;

    QMap<MzHashed, QVector<double>> mzHashedVsIonPresence;
    e = buildMzHashedVsIonPresence(
            mzHashedVsXICPoints,
            &mzHashedVsIonPresence
            ); ree;

    if (mzHashedVsIonPresence.isEmpty()) {
        ERR_RETURN
    }

    const Eigen::MatrixX<double> presenceMatrix = buildSummingMatrix(
            ms2Ions,
            mzHashedVsIonPresence,
            m_topNMS2Ions
    );

    if (presenceMatrix.rows() == 0) {
        ERR_RETURN
    }

    const Eigen::VectorX<double> summedMatVec = presenceMatrix.rowwise().sum();
    const QVector<double> summedMatVecToVec = EigenUtils::convertEigenVectorToQVector(summedMatVec);

    QVector<PeakIntegrationIndexes> peakIntegrationIndexes;

    if (scanTimePredicted > 0) {

        e = findCandidateIntegrations(
                summedMatVecToVec,
                &peakIntegrationIndexes
        ); ree;

    }

    else {

//        e = findCandidateIntegrations(
//                summedMatVecToVec,
//                m_pythiaParameters.minFoundMzPeaks,
//                scanTimePredicted,
//                m_pythiaParameters.scanTimeWindowMinutes,
//                &peakIntegrationIndexes
//        ); ree;

    }

    if (peakIntegrationIndexes.isEmpty()) {
        ERR_RETURN
    }

////    e = buildScores(
////            candidatePeptideTarget,
////            peakIntegrationIndexes,
////            summedMatVecToVec,
////            scoredCandidate
////    ); ree;

    ERR_RETURN
}

namespace {

    void sortPeakIntegrationsDescMaxSumFound(
            const QVector<double> &summedMatToVec,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
    ) {

        const auto sortLogic = [summedMatToVec](const PeakIntegrationIndexes &l, const PeakIntegrationIndexes &r){

            const int peakWidthL = l.second - l.first;
            const int peakWidthR = r.second - r.first;

            const QVector<double> summedMatVecMaxL = summedMatToVec.mid(l.first, peakWidthL);
            const QVector<double> summedMatVecMaxR = summedMatToVec.mid(r.first, peakWidthR);

            const double summedPresenceIntegrationMaxL = *std::max_element(summedMatVecMaxL.begin(), summedMatVecMaxL.end());
            const double summedPresenceIntegrationMaxR = *std::max_element(summedMatVecMaxR.begin(), summedMatVecMaxR.end());

            if (MathUtils::tSame(summedPresenceIntegrationMaxL, summedPresenceIntegrationMaxR)) {
                return std::accumulate(summedMatVecMaxL.begin(), summedMatVecMaxL.end(), 0.0)
                       < std::accumulate(summedMatVecMaxR.begin(), summedMatVecMaxR.end(), 0.0);
            }

            return summedPresenceIntegrationMaxL < summedPresenceIntegrationMaxR;
        };

        std::sort(peakIntegrationIndexes->rbegin(), peakIntegrationIndexes->rend(), sortLogic);
    }

}//namespace
Err CandidateScorertron::findCandidateIntegrations(
        const QVector<double> &summedMatToVec,
        QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
) {

    ERR_INIT

    QVector<double> summedPresenceVecSmoothed;
    e = m_peakIntegratomatic.findAllPeaksLimitsInXIC(
            summedMatToVec,
            peakIntegrationIndexes,
            &summedPresenceVecSmoothed
    ); ree;

//    const int minPeakWidth = 3;
//    filterSummedVecPeakIntegrationsByPeakWidth(
//            summedMatToVec,
//            m_pythiaParameters.minFoundMzPeaks,
//            minPeakWidth,
//            peakIntegrationIndexes
//    );
//
//    sortPeakIntegrationsDescMaxSumFound(
//            summedMatToVec,
//            peakIntegrationIndexes
//    );
//
//    const int topNPeakIntegrations = 2;
//    const int peakIntegrationsMaxSize
//            = std::min(topNPeakIntegrations, peakIntegrationIndexes->size());
//
//    peakIntegrationIndexes->resize(peakIntegrationsMaxSize);

    ERR_RETURN
}


