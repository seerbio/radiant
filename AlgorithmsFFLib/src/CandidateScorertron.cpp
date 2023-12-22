//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron.h"

#include "CandidateScores.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "ScoreOverseer.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

class Q_DECL_HIDDEN CandidateScorertron::Private {

public:

    Private();
    ~Private();

    Err init(
            const PythiaParameters &pythiaParameters,
            int topNMS2Ions,
            const MzTargetKey &targetKey,
            MsFrame *ms1Frame
            );

    Eigen::VectorX<float> m_gaussKernel;
    ScoreOverseer m_scoreOverseer;

};

CandidateScorertron::Private::Private() {}

CandidateScorertron::Private::~Private() {}

Err CandidateScorertron::Private::init(
        const PythiaParameters &pythiaParameters,
        int topNMS2Ions,
        const MzTargetKey &targetKey,
        MsFrame *ms1Frame
        ) {

    ERR_INIT

    m_gaussKernel = EigenKernelUtils::buildGaussianFilter1D<float>(
            pythiaParameters.filterLength,
            pythiaParameters.sigma
    );

    e = ErrorUtils::isFalse(MathUtils::tZero(m_gaussKernel.maxCoeff())); ree;
    m_gaussKernel /= m_gaussKernel.maxCoeff();

    e = m_scoreOverseer.init(
            pythiaParameters,
            targetKey,
            ms1Frame
            ); ree;

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

CandidateScorertron::CandidateScorertron()
: m_topNMS2Ions(-1)
, m_msCalibratomatic(nullptr)
, d_ptr(new Private())
{}

CandidateScorertron::~CandidateScorertron() {}

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
        const QMap<ScanNumber, ScanPoints*> &diaTargetFrame,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const QMap<ScanNumber, ScanPoints> &scanPointsMS1,
        const PythiaParameters &pythiaParameters,
        const MzTargetKey &targetKey,
        const QPair<double, double> &scanTimeMinMax,
        int topNMS2Ions,
        const QHash<MzHashed, int> &mzHashedVsCount,
        MsCalibratomatic *msCalibratomatic
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(topNMS2Ions > 0); ree;

    m_pythiaParameters = pythiaParameters;
    m_topNMS2Ions = topNMS2Ions;
    m_scanTimeMinMax = scanTimeMinMax;
    m_mzHashedVsCount = mzHashedVsCount;

    const PeakIntegratomaticParameters ffParams = buildPeakIntegratomaticParams(m_pythiaParameters);
    e = m_peakIntegratomatic.init(ffParams); ree;

    m_msCalibratomatic = msCalibratomatic;

    MsFrame msFrame;
    e = msFrame.init(diaTargetFrame, scanNumberVsScanTime); ree;
    e = m_turboXIC.init(msFrame.frameIndexVsScanPoints()); ree;

    e = m_turboXIC.getRTreeLimits(
            &m_mzMin,
            &m_mzMax
            );

    m_scanNumberVsScanPointsMS1 = scanPointsMS1;
    for (auto it = m_scanNumberVsScanPointsMS1.begin(); it != m_scanNumberVsScanPointsMS1.end(); it++) {
        m_scanNumberVsScanPointsMS1Pntrs.insert(it.key(), &it.value());
    }

    e = m_ms1Frame.init(m_scanNumberVsScanPointsMS1Pntrs, scanNumberVsScanTime); ree;

    e = d_ptr->init(
            pythiaParameters,
            m_topNMS2Ions,
            targetKey,
            &m_ms1Frame
    ); ree;

    m_targetKey = targetKey;

    ERR_RETURN
}

namespace {

    void buildMS2TheoreticalIsotopeShadows(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            QVector<MS2Ion> *ms2IonsTheoreticalIsotopeShadows
            ) {

        std::transform(
        ms2IonsTheoretical.begin(),
        ms2IonsTheoretical.end(),
        std::back_inserter(*ms2IonsTheoreticalIsotopeShadows),
        [](const MS2Ion &ms2Ion){
            const double isoChargeDiff = S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge;
            MS2Ion ms2IonNew = ms2Ion;
            ms2IonNew.mz -= static_cast<float>(isoChargeDiff);
            return ms2IonNew;
        }
        );
    }

    Err buildIntegrationVector(
            const QHash<MzHashed , XICPoints> &mzHashedVsXICPoints,
            QVector<float> *integrationVector
            ) {

        ERR_INIT

        FrameIndex frameIndexAllXICsMax = -1;
        for (const XICPoints &xicPoints : mzHashedVsXICPoints) {

            if (xicPoints.empty()) {
                continue;
            }

            const FrameIndex frameIndexXICMax = std::max_element(
                    xicPoints.begin(),
                    xicPoints.end(),
                    [](const XICPoint &l, const XICPoint &r){return l.scanNumber < r.scanNumber;})->scanNumber;

            frameIndexAllXICsMax = std::max(frameIndexAllXICsMax, frameIndexXICMax);
        }

        const int buffer = 2;
        const int vecSize = frameIndexAllXICsMax + buffer;

        QVector<float> eigenLoadingVector(vecSize, 0.0f);
        eigenLoadingVector.reserve(vecSize);

        for (const XICPoints &xic : mzHashedVsXICPoints) {

            for (const XICPoint &xp : xic) {
                const FrameIndex frameIndex = xp.scanNumber;

                eigenLoadingVector[frameIndex] += 1.0;
            }
        }

        *integrationVector = eigenLoadingVector;

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::calculateScores(
        const QVector<MS2Ion> &_ms2IonsTheoretical,
        TargetDecoyCandidatePair* targetDecoyCandidatePair,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(_ms2IonsTheoretical); ree;

    candidateScores->targetKey = m_targetKey;
    candidateScores->targetDecoyCandidatePair = targetDecoyCandidatePair;
    candidateScores->initFeaturesArray();

    const int topNMS2Ions = std::min(m_topNMS2Ions, _ms2IonsTheoretical.size());
    QVector<MS2Ion> ms2IonsTheoretical = _ms2IonsTheoretical;
    ms2IonsTheoretical.resize(topNMS2Ions);

    const int nStdDevsScanTimeWindow = 3;

    FrameIndex frameIndexPredictedMin = -1;
    FrameIndex frameIndexPredictedMax = -1;
    //TODO test whether it is better to filter out the points before integration, i.e. in extractHills(),
    // or filter out the integrations that don't intersect the predicted time range afterward afterwards.
    // I'm leaning towards the latter as this will not cut off a peak in the middle.
    if (m_msCalibratomatic->isInit()) {

        const float scanTimeWindow = m_msCalibratomatic->scanTimeStDev(nStdDevsScanTimeWindow);

        e = m_msCalibratomatic->predictScanTime(
                targetDecoyCandidatePair->iRt(),
                &candidateScores->scanTimePredicted
        ); ree;

        e = m_ms1Frame.frameIndexFromScanTime(
                candidateScores->scanTimePredicted - scanTimeWindow,
                &frameIndexPredictedMin
        ); ree;

        e = m_ms1Frame.frameIndexFromScanTime(
                candidateScores->scanTimePredicted + scanTimeWindow,
                &frameIndexPredictedMax
        ); ree;

    }

    QHash<MzHashed , XICPoints> mzHashedVsXICPoints;
    e = extractXICs(
            ms2IonsTheoretical,
            &mzHashedVsXICPoints
            ); ree;

    const QList<XICPoints> &xicPoints = mzHashedVsXICPoints.values();
    const int mzIonsFound = static_cast<int>(std::count_if(
            xicPoints.begin(),
            xicPoints.end(),
            [](const XICPoints &xic){return !xic.empty();}
            ));

    if (mzIonsFound < m_pythiaParameters.minFoundMzPeaks) {
        ERR_RETURN
    }

    QVector<float> integrationVector;
    e = buildIntegrationVector(mzHashedVsXICPoints, &integrationVector); ree;

    const double filterDeltaScoreValue = 2.0; //TODO consider making this settable

    QVector<QPair<PeakIntegrationIndexes, float>> peakIntegrationIndexesVsIntensity;
    e = findCandidateIntegrations(
            integrationVector,
            filterDeltaScoreValue,
            frameIndexPredictedMin,
            frameIndexPredictedMax,
            &peakIntegrationIndexesVsIntensity
            ); ree;

    if (peakIntegrationIndexesVsIntensity.isEmpty()) {
        ERR_RETURN
    }

    QVector<MS2Ion> ms2IonsTheoreticalIsotopeShadows;
    buildMS2TheoreticalIsotopeShadows(
            ms2IonsTheoretical,
            &ms2IonsTheoreticalIsotopeShadows
    );

    e = processPeakIntegrationIndexes(
            peakIntegrationIndexesVsIntensity,
            ms2IonsTheoretical,
            mzHashedVsXICPoints,
            ms2IonsTheoreticalIsotopeShadows,
            targetDecoyCandidatePair,
            candidateScores
            ); ree;

    ERR_RETURN
}

Err CandidateScorertron::extractXICs(
        const QVector<MS2Ion> &ms2IonsTheoretical,
        QHash<MzHashed , XICPoints> *mzHashedVsXICPoints
        ) {
    ERR_INIT

    for (const MS2Ion &ms2Ion : ms2IonsTheoretical) {

        const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

        if (ms2Ion.mz < m_mzMin || ms2Ion.mz > m_mzMax) {
            continue;
        }

        if (m_mzHashedVsXICPointsCached.contains(mzHashed)) {
            mzHashedVsXICPoints->insert(mzHashed, m_mzHashedVsXICPointsCached.value(mzHashed));

            if (m_mzHashedVsCount.value(mzHashed) == 2) {
                m_mzHashedVsXICPointsCached.remove(mzHashed);
            }

            continue;
        }

        const float mzTol = MathUtils::calculatePPM(
                ms2Ion.mz,
                static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
        );
        const float mzMin = ms2Ion.mz - mzTol;
        const float mzMax = ms2Ion.mz + mzTol;

        const XICPoints xicPoints = m_turboXIC.extractPointsXIC(
                mzMin,
                mzMax
                );

        mzHashedVsXICPoints->insert(mzHashed, xicPoints);

        const int cacheSizeMax = 2e4;
        if (m_mzHashedVsXICPointsCached.size() > cacheSizeMax) {
            m_mzHashedVsXICPointsCached.clear();
        }

        const int minCountToCache = 2;
        if (m_mzHashedVsCount.value(mzHashed) >= minCountToCache) {
            m_mzHashedVsXICPointsCached.insert(mzHashed, xicPoints);
        }
    }

    ERR_RETURN
}

namespace {

    void filterByFrameIndexes(
            FrameIndex frameIndexPredictedMin,
            FrameIndex frameIndexPredictedMax,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
            ) {

        if (frameIndexPredictedMin < 0 || frameIndexPredictedMax < 0) {
            return;
        }

        const auto terminatorLogic = [frameIndexPredictedMin, frameIndexPredictedMax](const QPair<PeakIntegrationIndexes, float> &pr){
            return pr.first.first > frameIndexPredictedMax || pr.first.second < frameIndexPredictedMin;
        };

        const auto terminator = std::remove_if(
                peakIntegrationIndexesVsIntensity->begin(),
                peakIntegrationIndexesVsIntensity->end(),
                terminatorLogic
                );

        peakIntegrationIndexesVsIntensity->erase(terminator, peakIntegrationIndexesVsIntensity->end());
    }

    void filterSummedVecPeakIntegrationsByPeakWidth(
            const QVector<float> &summedMatToVec,
            int summedMzPresenceMin,
            int peakWidthMin,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
            ) {

        const auto terminatorLogic
            = [summedMatToVec, summedMzPresenceMin, peakWidthMin](const QPair<PeakIntegrationIndexes, float> &pii){

            const int peakWidth = pii.first.second - pii.first.first;
            if (peakWidth < peakWidthMin) {
                return true;
            }

            const QVector<float> summedMatVecMax = summedMatToVec.mid(pii.first.first, peakWidth);
            const float summedPresenceIntegrationMax = *std::max_element(summedMatVecMax.begin(), summedMatVecMax.end());

            return (summedPresenceIntegrationMax < static_cast<float>(summedMzPresenceMin));
        };

        const auto terminator = std::remove_if(
                peakIntegrationIndexesVsIntensity->begin(),
                peakIntegrationIndexesVsIntensity->end(),
                terminatorLogic
                );

        peakIntegrationIndexesVsIntensity->erase(terminator, peakIntegrationIndexesVsIntensity->end());
    }

    Err filterByDeltaScore(
            double filterDeltaScoreValue,
            QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
            ) {

        ERR_INIT

        if (peakIntegrationIndexesVsIntensity->size() < 2) {
            ERR_RETURN;
        }

        using PR = QPair<PeakIntegrationIndexes, float>;

        std::sort(
                peakIntegrationIndexesVsIntensity->rbegin(),
                peakIntegrationIndexesVsIntensity->rend(),
                [](const PR &l, const PR &r){return l.second < r.second;}
                );

        QVector<QPair<PeakIntegrationIndexes, float>> peakIntegrationIndexesVsIntensityFiltered;
        for (int i = 0; i < peakIntegrationIndexesVsIntensity->size(); i++) {

            if (i == peakIntegrationIndexesVsIntensity->size() - 1) {
                peakIntegrationIndexesVsIntensityFiltered.push_back(peakIntegrationIndexesVsIntensity->at(i));
                break;
            }

            const float &topIntensity = peakIntegrationIndexesVsIntensity->at(i).second;
            const float &bottomIntensity = peakIntegrationIndexesVsIntensity->at(i+1).second;

            e = ErrorUtils::isTrue(topIntensity >= bottomIntensity); ree;

            peakIntegrationIndexesVsIntensityFiltered.push_back(peakIntegrationIndexesVsIntensity->at(i));

            if (topIntensity - bottomIntensity > filterDeltaScoreValue) {
                break;
            }
        }

        *peakIntegrationIndexesVsIntensity = peakIntegrationIndexesVsIntensityFiltered;

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::findCandidateIntegrations(
        const QVector<float> &integrationVector,
        double filterDeltaScoreValue,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
) {

    ERR_INIT

    e = simpleIntegrator(
            integrationVector,
            peakIntegrationIndexesVsIntensity
            ); ree;

    filterByFrameIndexes(
            frameIndexPredictedMin,
            frameIndexPredictedMax,
            peakIntegrationIndexesVsIntensity
            );

    const int minPeakWidth = 1;
    filterSummedVecPeakIntegrationsByPeakWidth(
            integrationVector,
            m_pythiaParameters.minFoundMzPeaks,
            minPeakWidth,
            peakIntegrationIndexesVsIntensity
            );

    const int maxPeakIntegrationIndexesReturned = 2; //TODO consider making this settable
    peakIntegrationIndexesVsIntensity->resize(
            std::min(maxPeakIntegrationIndexesReturned, peakIntegrationIndexesVsIntensity->size())
            );

    e = filterByDeltaScore(filterDeltaScoreValue, peakIntegrationIndexesVsIntensity); ree;

    ERR_RETURN
}

Err CandidateScorertron::simpleIntegrator(
        const QVector<float> &vec,
        QVector<QPair<PeakIntegrationIndexes, float>> *peakIntegrationIndexesVsIntensity
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(vec); ree;

    Eigen::VectorX<float> eVec = EigenUtils::convertQVectorToEigenVector(vec);
    EigenUtils::thresholdVector(static_cast<float>(1.01), &eVec);
    for (int smoothCount = 0; smoothCount < m_pythiaParameters.smoothCount; smoothCount++) {
        eVec = EigenKernelUtils::convolveVectorWithKernel(eVec, d_ptr->m_gaussKernel);
    }

    const int topNApexes = 10;
    const QMap<int, float> vecApexs = EigenUtils::apexes(eVec);

    if (vecApexs.isEmpty()) {
        ERR_RETURN
    }

    Eigen::VectorX<float> apexes =EigenUtils::convertQMapToEigenVector(vecApexs, vecApexs.lastKey() + 1);
    QVector<QPair<int, float>> apexPairs = EigenUtils::returnTopXIndexAndValues(apexes, topNApexes);

    const int maxPeakIntegrations = 5;
    apexPairs.resize(std::min(maxPeakIntegrations, apexPairs.size()));
    for (const QPair<int, float> &pr : apexPairs) {

        const int apexIndex = pr.first;
        const float apexValue = pr.second;

        if (MathUtils::tZero(apexValue)) {
            continue;
        }

        const float stopThreshold = apexValue * m_pythiaParameters.stopThresholdFraction;

        float rightStopVal = apexValue;
        int rightStopIndex = apexIndex;

        int rightCurrentIndex = apexIndex;
        while (rightCurrentIndex < eVec.size()) {

            const float currentValue = eVec(rightCurrentIndex);
            if (currentValue < stopThreshold) {
                rightStopIndex = rightCurrentIndex;
                break;
            }

            if (currentValue <= rightStopVal) {
                rightStopVal = currentValue;
                rightStopIndex = rightCurrentIndex;
                rightCurrentIndex++;
                continue;
            }

            break;
        }

        float leftStopVal = apexValue;
        int leftStopIndex = apexIndex;

        int leftCurrentIndex = apexIndex;
        while (leftCurrentIndex < eVec.size()) {

            const float currentValue = eVec(leftCurrentIndex);
            if (currentValue < stopThreshold) {
                leftStopIndex = leftCurrentIndex;
                break;
            }

            if (currentValue <= leftStopVal) {
                leftStopVal = currentValue;
                leftStopIndex = leftCurrentIndex;
                leftCurrentIndex--;
                continue;
            }

            break;
        }

        peakIntegrationIndexesVsIntensity->push_back({
            {std::max(leftStopIndex, 0), std::min(rightStopIndex, vec.size() -1)},
            apexValue}
            );

        for (int i = leftStopIndex; i <= rightStopIndex; i++) {
            eVec.coeffRef(i) = 0.0;
        }
    }

    ERR_RETURN
}

Err CandidateScorertron::processPeakIntegrationIndexes(
        const QVector<QPair<PeakIntegrationIndexes, float>> &peakIntegrationIndexesVsIntensity,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QHash<MzHashed , XICPoints> &mzHashedVsXICPoints,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
        TargetDecoyCandidatePair* targetDecoyCandidatePair,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    QHash<MzHashed , XICPoints> mzHashedVsXICPointsShadows;
    e = extractXICs(
            ms2IonsTheoreticalIsotopeShadows,
            &mzHashedVsXICPointsShadows
            ); ree;

    const auto unfragPrecursorMass = static_cast<float>(BiophysicalCalcs::calculateThomsonFromMass(
            targetDecoyCandidatePair->mass(),
            targetDecoyCandidatePair->charge(),
            0
    ));

    QHash<MzHashed , XICPoints> mzHashedVsXICPointsUnfragPrecursor;
    MS2Ion ms2IonUnfragPrecursor;
    ms2IonUnfragPrecursor.mz = unfragPrecursorMass;
    e = extractXICs(
            {ms2IonUnfragPrecursor},
            &mzHashedVsXICPointsUnfragPrecursor
    ); ree;

    //TODO figure out a way to empirially pick the best integration index instead
    // of wasting cycles analyzing the top 2; put break at the bottom of loop as temp fix.
    for (const QPair<PeakIntegrationIndexes, float> &pii : peakIntegrationIndexesVsIntensity) {

        e = ErrorUtils::isTrue(pii.first.second > pii.first.first); ree;

        CandidateScores candidateScoresPII;
        e = d_ptr->m_scoreOverseer.buildScores(
                targetDecoyCandidatePair,
                m_scanTimeMinMax,
                pii.first,
                ms2IonsTheoretical,
                mzHashedVsXICPoints,
                ms2IonsTheoreticalIsotopeShadows,
                mzHashedVsXICPointsShadows,
                ms2IonUnfragPrecursor,
                mzHashedVsXICPointsUnfragPrecursor,
                &candidateScoresPII
                ); ree;

        if (candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100]
            < candidateScoresPII.featuresArray[CandidateScores::Features::CosineSimSum100]) {
            *candidateScores = candidateScoresPII;
        }

        break;
    }

    ERR_RETURN
}
