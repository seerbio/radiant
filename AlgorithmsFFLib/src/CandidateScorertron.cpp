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
, m_featureFinderHillsBuilderMS2(nullptr)
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
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const QMap<ScanNumber, ScanPoints> &scanPointsMS1,
        const PythiaParameters &pythiaParameters,
        const MzTargetKey &targetKey,
        int topNMS2Ions,
        MsCalibratomatic *msCalibratomatic,
        FeatureFinderHillBuilder *featureFinderHillsBuilderMS2
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(topNMS2Ions > 0); ree;

    m_pythiaParameters = pythiaParameters;
    m_topNMS2Ions = topNMS2Ions;

    const PeakIntegratomaticParameters ffParams = buildPeakIntegratomaticParams(m_pythiaParameters);
    e = m_peakIntegratomatic.init(ffParams); ree;

    m_msCalibratomatic = msCalibratomatic;
    m_featureFinderHillsBuilderMS2 = featureFinderHillsBuilderMS2;

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
            ms2IonNew.mz -= isoChargeDiff;
            return ms2IonNew;
        }
        );
    }

    Err buildIntegrationVector(
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            QVector<float> *integrationVector
            ) {

        ERR_INIT

        const QPair<FrameIndex, FrameIndex> minMaxFrameIndexes = ScoreOverseer::getMinMaxFrameIndexes(mzHashedVsfeatureFinderHills);
        const int buffer = 2;
        const int vecSize = minMaxFrameIndexes.second + buffer;

        QVector<float> eigenLoadingVector(vecSize, static_cast<float>(0.0));
        eigenLoadingVector.reserve(vecSize);

        for (const QVector<FeatureFinderHill*> &ffhs : mzHashedVsfeatureFinderHills) {
            for (FeatureFinderHill *ffh : ffhs) {
                const QVector<FrameIndex> frameIndexes = ffh->scanNumberIndexes();
                for (int i = 0; i < frameIndexes.size(); i++) {
                    e = ErrorUtils::isTrue(frameIndexes[i] < vecSize); ree;
                    eigenLoadingVector[frameIndexes[i]] += 1.0;
                }
            }
        }

        *integrationVector = eigenLoadingVector;

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::calculateScores(
        const TargetDecoyCandidatePair* targetDecoyCandidatePair,
        const QVector<MS2Ion> &_ms2IonsTheoretical,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(_ms2IonsTheoretical); ree;

    candidateScores->peptideStringWithMods = targetDecoyCandidatePair->peptideStringWithMods();
    candidateScores->charge = targetDecoyCandidatePair->charge();
    candidateScores->targetKey = m_targetKey;

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

        const double scanTimeWindow = m_msCalibratomatic->scanTimeStDev(nStdDevsScanTimeWindow);

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

        //TODO write code to filter peakIntegration Indexes by FrameIndex
    }

    QHash<MzHashed , QVector<FeatureFinderHill*>> mzHashedVsfeatureFinderHills;
    e = extractHills(
            ms2IonsTheoretical,
            &mzHashedVsfeatureFinderHills
            ); ree;

    const QList<QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsVals = mzHashedVsfeatureFinderHills.values();
    const int mzIonsFound = static_cast<int>(std::count_if(
            mzHashedVsfeatureFinderHillsVals.begin(),
            mzHashedVsfeatureFinderHillsVals.end(),
            [](const QVector<FeatureFinderHill*> &ffhs){return !ffhs.isEmpty();}
            ));

    if (mzIonsFound < m_pythiaParameters.minFoundMzPeaks) {
        ERR_RETURN
    }

    QVector<float> integrationVector;
    e = buildIntegrationVector(mzHashedVsfeatureFinderHills, &integrationVector); ree;

    QVector<PeakIntegrationIndexes> peakIntegrationIndexes;
    e = findCandidateIntegrations(
            integrationVector,
            &peakIntegrationIndexes
            ); ree;

    if (peakIntegrationIndexes.isEmpty()) {
        ERR_RETURN
    }

    QVector<MS2Ion> ms2IonsTheoreticalIsotopeShadows;
    buildMS2TheoreticalIsotopeShadows(
            ms2IonsTheoretical,
            &ms2IonsTheoreticalIsotopeShadows
    );

    e = processPeakIntegrationIndexes(
            targetDecoyCandidatePair,
            peakIntegrationIndexes,
            ms2IonsTheoretical,
            mzHashedVsfeatureFinderHills,
            ms2IonsTheoreticalIsotopeShadows,
            candidateScores
            ); ree;

//    Eigen::MatrixX<float> mat;
//    Eigen::MatrixX<float> matShadows;
//    e = buildAlignmentMatricies(
//            ms2IonsTheoretical,
//            ms2IonsTheoreticalIsotopeShadows,
//            m_pythiaParameters.filterLength,
//            mzHashedVsfeatureFinderHills,
//            mzHashedVsfeatureFinderHillsShadows,
//            &mat,
//            &matShadows
//            ); ree;

//    const Eigen::MatrixX<float> matSmoothed = applyGaussSmoothRowWiseToMatrix(
//            mat,
//            m_pythiaParameters,
//            d_ptr->gaussKernel
//            );
//    const Eigen::MatrixX<float> matShadowsSmoothed = applyGaussSmoothRowWiseToMatrix(
//            matShadows,
//            m_pythiaParameters,
//            d_ptr->gaussKernel
//    );
//
//    Eigen::MatrixX<float> matIsotopesSubtracted = matSmoothed.array() - matShadowsSmoothed.array();
//    EigenUtils::thresholdMatrix(static_cast<float>(0.0), &matIsotopesSubtracted);

    ERR_RETURN
}

Err CandidateScorertron::extractHills(
        const QVector<MS2Ion> &ms2IonsTheoretical,
        QHash<MzHashed , QVector<FeatureFinderHill*>> *mzHashedVsfeatureFinderHills
        ) {
    ERR_INIT

    QHash<MzHashed , QVector<FeatureFinderHill*>> mzHashedVsfeatureFinderHillsTemp;

    for (const MS2Ion &ms2Ion : ms2IonsTheoretical) {

        const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

        if (m_mzHashedVsFeatureFinderHillsCached.contains(mzHashed)) {
            mzHashedVsfeatureFinderHillsTemp.insert(mzHashed, m_mzHashedVsFeatureFinderHillsCached.value(mzHashed));
            continue;
        }

        const double mzTol = MathUtils::calculatePPM(
                ms2Ion.mz,
                m_pythiaParameters.ms2ExtractionWidthPPM
        );
        const double mzMin = ms2Ion.mz - mzTol;
        const double mzMax = ms2Ion.mz + mzTol;

        QVector<FeatureFinderHill*> mzFeatureFinderHills;
        e = m_featureFinderHillsBuilderMS2->getHills(
                mzMin,
                mzMax,
                &mzFeatureFinderHills
        ); ree;

        mzHashedVsfeatureFinderHillsTemp[mzHashed].append(mzFeatureFinderHills);
        m_mzHashedVsFeatureFinderHillsCached.insert(mzHashed, mzHashedVsfeatureFinderHillsTemp.value(mzHashed));
    }

    *mzHashedVsfeatureFinderHills = mzHashedVsfeatureFinderHillsTemp;
    ERR_RETURN
}

namespace {

    void filterSummedVecPeakIntegrationsByPeakWidth(
            const QVector<float> &summedMatToVec,
            int summedMzPresenceMin,
            int peakWidthMin,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
            ) {

        const auto terminatorLogic = [summedMatToVec, summedMzPresenceMin, peakWidthMin](const PeakIntegrationIndexes &pii){

            const int peakWidth = pii.second - pii.first;
            if (peakWidth < peakWidthMin) {
                return true;
            }

            const QVector<float> summedMatVecMax = summedMatToVec.mid(pii.first, peakWidth);
            const float summedPresenceIntegrationMax = *std::max_element(summedMatVecMax.begin(), summedMatVecMax.end());

            return (summedPresenceIntegrationMax < static_cast<float>(summedMzPresenceMin));
        };

        const auto terminator = std::remove_if(peakIntegrationIndexes->begin(), peakIntegrationIndexes->end(), terminatorLogic);
        peakIntegrationIndexes->erase(terminator, peakIntegrationIndexes->end());
    }

}//namespace
Err CandidateScorertron::findCandidateIntegrations(
        const QVector<float> &integrationVector,
        QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
) {

    ERR_INIT

    e = simpleIntegrator(
            integrationVector,
            peakIntegrationIndexes
            ); ree;

    const int minPeakWidth = 2;
    filterSummedVecPeakIntegrationsByPeakWidth(
            integrationVector,
            m_pythiaParameters.minFoundMzPeaks,
            minPeakWidth,
            peakIntegrationIndexes
            );

    const int maxPeakIntegrationIndexesReturned = 2;
    peakIntegrationIndexes->resize(std::min(maxPeakIntegrationIndexesReturned, peakIntegrationIndexes->size()));

    ERR_RETURN
}

Err CandidateScorertron::simpleIntegrator(
        const QVector<float> &vec,
        QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
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

        peakIntegrationIndexes->push_back({leftStopIndex, rightStopIndex});
        for (int i = leftStopIndex; i <= rightStopIndex; i++) {
            eVec.coeffRef(i) = 0.0;
        }
    }

    ERR_RETURN
}

namespace {

    void filterFeatureFinderHillsByFrameIndex(
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            FrameIndex frameIndexMin,
            FrameIndex frameIndexMax,
            QHash<MzHashed , QVector<FeatureFinderHill*>> *mzHashedVsfeatureFinderHillsFiltered
    ) {
        for (auto it = mzHashedVsfeatureFinderHills.begin(); it != mzHashedVsfeatureFinderHills.end(); it++) {
            const MzHashed mzHashed = it.key();
            const QVector<FeatureFinderHill*> &ffhs = it.value();
            for (FeatureFinderHill *ffh : ffhs) {
                const QPair<FrameIndex, FrameIndex> ffhFrameIndexMinMax = ffh->scanNumberIndexMinMax();

                if(ffhFrameIndexMinMax.second < frameIndexMin || ffhFrameIndexMinMax.first > frameIndexMax) {
                    continue;
                }
                (*mzHashedVsfeatureFinderHillsFiltered)[mzHashed].push_back(ffh);
            }
        }
    }

}//namespace
Err CandidateScorertron::processPeakIntegrationIndexes(
        const TargetDecoyCandidatePair* targetDecoyCandidatePair,
        QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    QHash<MzHashed , QVector<FeatureFinderHill*>> mzHashedVsfeatureFinderHillsShadows;
    e = extractHills(
            ms2IonsTheoreticalIsotopeShadows,
            &mzHashedVsfeatureFinderHillsShadows
            ); ree;

    //TODO figure out a way to empirially pick the best integration index instead
    // of wasting cycles analyzing the top 2;
    for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {

        e = ErrorUtils::isTrue(pii.second > pii.first); ree;

        QHash<MzHashed , QVector<FeatureFinderHill*>> mzHashedVsfeatureFinderHillsFiltered;
        filterFeatureFinderHillsByFrameIndex(
                mzHashedVsfeatureFinderHills,
                pii.first,
                pii.second,
                &mzHashedVsfeatureFinderHillsFiltered
                );

        QHash<MzHashed , QVector<FeatureFinderHill*>> mzHashedVsfeatureFinderHillsShadowsFiltered;
        filterFeatureFinderHillsByFrameIndex(
                mzHashedVsfeatureFinderHillsShadows,
                pii.first,
                pii.second,
                &mzHashedVsfeatureFinderHillsShadowsFiltered
                );

        CandidateScores candidateScoresPII;
        e = d_ptr->m_scoreOverseer.buildScores(
                targetDecoyCandidatePair,
                pii,
                ms2IonsTheoretical,
                mzHashedVsfeatureFinderHillsFiltered,
                ms2IonsTheoreticalIsotopeShadows,
                mzHashedVsfeatureFinderHillsShadowsFiltered,
                &candidateScoresPII
                ); ree;

        if (candidateScores->cosineSimSum100 < candidateScoresPII.cosineSimSum100) {
            *candidateScores = candidateScoresPII;
        }

    }

    ERR_RETURN
}
