//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron.h"

#include "CandidateScores.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrame.h"
#include "ScoreOverseer.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"


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
Err CandidateScorertron::init(
        const QMap<ScanNumber, ScanPoints> &scanNumberVsScanPointsMS1,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const PythiaParameters &pythiaParameters,
        int topNMS2Ions
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanNumberVsScanPointsMS1); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(topNMS2Ions > 0); ree;

    MsFrame ms1Frame;
    e = ms1Frame.init(
            scanNumberVsScanPointsMS1,
            scanNumberVsScanTime
            ); ree;

    e = m_turboXICMS1.init(ms1Frame.frameIndexVsScanPoints()); ree;

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

            const QMap<ScanNumber, double> &scanNumbersVsIntensityVals = xp.scanNumbersVsIntensityVals;

            if (scanNumbersVsIntensityVals.isEmpty()) {
                continue;
            }

            frameIndex = std::max(frameIndex, scanNumbersVsIntensityVals.lastKey());
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

            const QMap<ScanNumber, double> &scanNumbersVsIntensityVals = xicPoints.scanNumbersVsIntensityVals;


            if (scanNumbersVsIntensityVals.isEmpty()) {
                mzHashedVsIonPresence->insert(mzHashed, {});
                continue;
            }

            Eigen::VectorX<double> vecNormalized = EigenUtils::convertQMapToEigenVector(
                    scanNumbersVsIntensityVals,
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
        const TargetDecoyCandidatePair* targetDecoyCandidatePair,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPointsIsotopeShadows,
        double scanTimePredicted,
        MsFrame *msFrame,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;
    e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;

    QMap<MzHashed, QVector<double>> mzHashedVsIonPresence;
    e = buildMzHashedVsIonPresence(
            mzHashedVsXICPoints,
            &mzHashedVsIonPresence
            ); ree;

    if (mzHashedVsIonPresence.isEmpty()) {
        ERR_RETURN
    }

    const Eigen::MatrixX<double> presenceMatrix = buildSummingMatrix(
            ms2IonsTheoretical,
            mzHashedVsIonPresence,
            m_topNMS2Ions
    );

    const Eigen::VectorX<double> summedPresenceMatrixVec = presenceMatrix.rowwise().sum();
    const QVector<double> summedMatVecToVec = EigenUtils::convertEigenVectorToQVector(summedPresenceMatrixVec);

    QVector<PeakIntegrationIndexes> peakIntegrationIndexes;

    e = findCandidateIntegrations(
            summedMatVecToVec,
            &peakIntegrationIndexes
    ); ree;

    ScoreOverseer scoreOverseer(
            m_topNMS2Ions,
            m_pythiaParameters.cosineSimToAnchorThreshold,
            m_pythiaParameters.ms2ExtractionWidthPPM,
            summedMatVecToVec,
            &m_turboXICMS1
            );

    e = scoreOverseer.buildScores(
            targetDecoyCandidatePair,
            peakIntegrationIndexes,
            ms2IonsTheoretical,
            mzHashedVsXICPoints,
            ms2IonsTheoreticalIsotopeShadows,
            mzHashedVsXICPointsIsotopeShadows,
            m_pythiaParameters.ms2ExtractionWidthPPM,
            scanTimePredicted,
            msFrame,
            candidateScores
            ); ree;

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


    void filterSummedVecPeakIntegrationsByPeakWidth(
            const QVector<double> &summedMatToVec,
            int summedMzPresenceMin,
            int peakWidthMin,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
            ) {

        const auto terminatorLogic = [summedMatToVec, summedMzPresenceMin, peakWidthMin](const PeakIntegrationIndexes &pii){

            const int peakWidth = pii.second - pii.first;
            if (peakWidth < peakWidthMin) {
                return true;
            }

            const QVector<double> summedMatVecMax = summedMatToVec.mid(pii.first, peakWidth);
            const double summedPresenceIntegrationMax = *std::max_element(summedMatVecMax.begin(), summedMatVecMax.end());

            return summedPresenceIntegrationMax < summedMzPresenceMin;
        };

        const auto terminator = std::remove_if(peakIntegrationIndexes->begin(), peakIntegrationIndexes->end(), terminatorLogic);
        peakIntegrationIndexes->erase(terminator, peakIntegrationIndexes->end());
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

    const int minPeakWidth = 2;
    filterSummedVecPeakIntegrationsByPeakWidth(
            summedMatToVec,
            m_pythiaParameters.minFoundMzPeaks,
            minPeakWidth,
            peakIntegrationIndexes
            );

    sortPeakIntegrationsDescMaxSumFound(
            summedMatToVec,
            peakIntegrationIndexes
            );

    const int topNPeakIntegrations = 2;
    const int peakIntegrationsMaxSize
            = std::min(topNPeakIntegrations, peakIntegrationIndexes->size());

    peakIntegrationIndexes->resize(peakIntegrationsMaxSize);

    ERR_RETURN
}
