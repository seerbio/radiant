//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron.h"

#include "CandidateScores.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "MsFrame.h"
//#include "ScoreOverseer.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

class Q_DECL_HIDDEN CandidateScorertron::Private {

public:

    Private();
    ~Private();

    void initGaussKernel(const PythiaParameters &pythiaParameters);

    Eigen::VectorX<float> gaussKernel;

};

CandidateScorertron::Private::Private() {}

CandidateScorertron::Private::~Private() {}

void CandidateScorertron::Private::initGaussKernel(const PythiaParameters &pythiaParameters) {

    gaussKernel = EigenKernelUtils::buildGaussianFilter1D<float>(
            pythiaParameters.filterLength,
            pythiaParameters.sigma
    );

    gaussKernel /= gaussKernel.maxCoeff();
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

CandidateScorertron::CandidateScorertron()
: m_topNMS2Ions(-1)
, m_msCalibratomatic(nullptr)
, m_featureFinderHillsBuilderMS1(nullptr)
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
        const PythiaParameters &pythiaParameters,
        int topNMS2Ions,
        MsCalibratomatic *msCalibratomatic,
        FeatureFinderHillBuilder *featureFinderHillsBuilderMS1,
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
    m_featureFinderHillsBuilderMS1 = featureFinderHillsBuilderMS1;
    m_featureFinderHillsBuilderMS2 = featureFinderHillsBuilderMS2;

    d_ptr->initGaussKernel(pythiaParameters);

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

    QPair<FrameIndex, FrameIndex> getMinMaxFrameIndexes(
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills
            ) {

        QPair<FrameIndex, FrameIndex> minMaxFrameIndexes = {
                std::numeric_limits<FrameIndex>::max(),
                std::numeric_limits<FrameIndex>::min()
                };

        for (const QVector<FeatureFinderHill*> &ffhs : mzHashedVsfeatureFinderHills.values()) {
            for (FeatureFinderHill* ffh : ffhs) {
                const QPair<FrameIndex, FrameIndex> &ffhMinMaxFrameIndex = ffh->scanNumberIndexMinMax();
                minMaxFrameIndexes.first = std::min(minMaxFrameIndexes.first, ffhMinMaxFrameIndex.first);
                minMaxFrameIndexes.second = std::max(minMaxFrameIndexes.second, ffhMinMaxFrameIndex.second);
            }
        }

        return minMaxFrameIndexes;
    }

    Err buildIntegrationVector(
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            QVector<float> *integrationVector
            ) {

        ERR_INIT

        const QPair<FrameIndex, FrameIndex> minMaxFrameIndexes = getMinMaxFrameIndexes(mzHashedVsfeatureFinderHills);
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

    Err loadMzHashedVsFeatureFinderHillsToMatrix(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            int frameIndexMin,
            Eigen::MatrixX<float> *mat
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(frameIndexMin >= 0); ree;
        e = ErrorUtils::isTrue(mat->rows() > 0); ree;

        const int rowCount = static_cast<int>(mat->rows());
        mat->setZero();

        for (int col = 0; col < ms2IonsTheoretical.size(); col++) {

            const MS2Ion &ms2Ion = ms2IonsTheoretical.at(col);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            QVector<float> eigenVectorLoader(rowCount);
            eigenVectorLoader.reserve(rowCount);

            const QVector<FeatureFinderHill*> &ffhs = mzHashedVsfeatureFinderHills.value(mzHashed);
            for (FeatureFinderHill *ffh : ffhs) {
                const QVector<FrameIndex> frameIndexes = ffh->scanNumberIndexes();
                const QVector<float> intensities = ffh->intensities();

                e = ErrorUtils::isEqual(frameIndexes.size(), intensities.size()); ree;

                for (int i = 0; i < frameIndexes.size(); i++) {

                    const int frameIndexAdjusted = frameIndexes.at(i) - frameIndexMin;

                    if (frameIndexAdjusted < 0 || frameIndexAdjusted >= rowCount) {
                        continue;
                    }

                    eigenVectorLoader[frameIndexAdjusted] = intensities.at(i);
                }
            }

            mat->col(col) = EigenUtils::convertQVectorToEigenVector(eigenVectorLoader);
        }

        ERR_RETURN
    }

    Err buildAlignmentMatricies(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QVector<MS2Ion> &ms2IonsTheoreticalShadows,
            int smoothFilterLength,
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsShadows,
            Eigen::MatrixX<float> *mat,
            Eigen::MatrixX<float> *matShadows
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;
        e = ErrorUtils::isNotEmpty(ms2IonsTheoreticalShadows); ree;
        e = ErrorUtils::isEqual(ms2IonsTheoretical.size(), ms2IonsTheoreticalShadows.size()); ree;

        const int cols = ms2IonsTheoretical.size();

        const QPair<FrameIndex, FrameIndex> minMaxFrameIndex = getMinMaxFrameIndexes(mzHashedVsfeatureFinderHills);
        e = ErrorUtils::isTrue(minMaxFrameIndex.second >= minMaxFrameIndex.first); ree;
        const int rows = std::max(minMaxFrameIndex.second - minMaxFrameIndex.first + 1, smoothFilterLength);

        mat->resize(rows, cols);
        e = loadMzHashedVsFeatureFinderHillsToMatrix(
                ms2IonsTheoretical,
                mzHashedVsfeatureFinderHills,
                minMaxFrameIndex.first,
                mat
                ); ree;

        matShadows->resize(rows, cols);
        e = loadMzHashedVsFeatureFinderHillsToMatrix(
                ms2IonsTheoreticalShadows,
                mzHashedVsfeatureFinderHillsShadows,
                minMaxFrameIndex.first,
                matShadows
        ); ree;

        ERR_RETURN
    }

    Eigen::MatrixX<float> applyGaussSmoothRowWiseToMatrix(
            const Eigen::MatrixX<float> &mat,
            const PythiaParameters &pythiaParameters,
            const Eigen::VectorX<float> &gaussKernel
            ) {

        Eigen::MatrixX<float> matSmoothed = mat;
        for (int smoothIter = 0; smoothIter < pythiaParameters.smoothCount; smoothIter++) {
            matSmoothed = EigenKernelUtils::applyKernelToEachColumnInMatrix(matSmoothed, gaussKernel);
        }

        return matSmoothed;
    }

//    Err calculateMS1Corr(
//            const Eigen::VectorX<double> &bestAnchorColumn,
//            const PeakIntegrationIndexes &peakIntegrationIndexes,
//            double mzTarget,
//            double ppmTol,
//            TurboXIC *turboXic,
//            double *cosineSimMS1
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isTrue(bestAnchorColumn.size() > 0); ree;
//
//        const double massTol = MathUtils::calculatePPM(mzTarget, ppmTol);
//        const double mzStart = mzTarget - massTol;
//        const double mzEnd = mzTarget + massTol;
//
//        const XICPoints xicPoints = turboXic->extractPointsXIC(
//                mzStart,
//                mzEnd,
//                peakIntegrationIndexes.first,
//                peakIntegrationIndexes.second
//        );
//
//        Eigen::VectorX<double> ms1Vec(static_cast<int>(bestAnchorColumn.size()));
//        ms1Vec.setZero();
//
//        const QMap<ScanNumber, double> &scanNumbersVsIntensityVals = xicPoints.scanNumbersVsIntensity;
//
//        for (auto it = scanNumbersVsIntensityVals.begin(); it != scanNumbersVsIntensityVals.end(); it++) {
//            const FrameIndex frameIndex = it.key() - peakIntegrationIndexes.first;
//
//            if (frameIndex >= ms1Vec.size()) {
//                continue;
//            }
//
//            ms1Vec.coeffRef(frameIndex) = it.value();
//        }
//
//        ms1Vec = applyGaussSmoothRowWiseToMatrix(ms1Vec);
//
//        *cosineSimMS1 = EigenUtils::cosineSimilarity(bestAnchorColumn, ms1Vec);
//
//        ERR_RETURN
//    }

}//namespace
Err CandidateScorertron::calculateScores(
        const TargetDecoyCandidatePair* targetDecoyCandidatePair,
        const QVector<MS2Ion> &_ms2IonsTheoretical,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(_ms2IonsTheoretical); ree;

    const int topNMS2Ions = std::min(m_topNMS2Ions, _ms2IonsTheoretical.size());
    QVector<MS2Ion> ms2IonsTheoretical = _ms2IonsTheoretical;
    ms2IonsTheoretical.resize(topNMS2Ions);

    QVector<MS2Ion> ms2IonsTheoreticalIsotopeShadows;
    buildMS2TheoreticalIsotopeShadows(
            ms2IonsTheoretical,
            &ms2IonsTheoreticalIsotopeShadows
            );

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

    e = processPeakIntegrationIndexes(
            targetDecoyCandidatePair,
            peakIntegrationIndexes,
            ms2IonsTheoretical,
            mzHashedVsfeatureFinderHills,
            ms2IonsTheoreticalIsotopeShadows
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




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

//    const Eigen::MatrixX<float>

//    //LAST BIGGEST
//    QMap<MzHashed, QVector<double>> mzHashedVsIonPresence;
//    e = buildMzHashedVsIonPresence(
//            mzHashedVsXICPoints,
//            &mzHashedVsIonPresence
//            ); ree;
//
//    if (mzHashedVsIonPresence.isEmpty()) {
//        ERR_RETURN
//    }
//
//    const Eigen::MatrixX<double> presenceMatrix = buildSummingMatrix(
//            ms2IonsTheoretical,
//            mzHashedVsIonPresence,
//            m_topNMS2Ions
//            );
//
//    const Eigen::VectorX<double> summedPresenceMatrixVec = presenceMatrix.rowwise().sum();
//    const QVector<double> summedMatVecToVec = EigenUtils::convertEigenVectorToQVector(summedPresenceMatrixVec);
//
//    QVector<PeakIntegrationIndexes> peakIntegrationIndexes;
//
//    //NEXT BIGGEST TIME
//    e = findCandidateIntegrations(
//            summedMatVecToVec,
//            &peakIntegrationIndexes
//    ); ree;
//
////    ScoreOverseer scoreOverseer(
////            m_topNMS2Ions,
////            m_pythiaParameters.cosineSimToAnchorThreshold,
////            m_pythiaParameters.ms2ExtractionWidthPPM,
////            summedMatVecToVec,
////            &m_turboXICMS1
////            );
////
////    //BIGGEST TIME
////    e = scoreOverseer.buildScores(
////            targetDecoyCandidatePair,
////            peakIntegrationIndexes,
////            ms2IonsTheoretical,
////            mzHashedVsXICPoints,
////            ms2IonsTheoreticalIsotopeShadows,
////            mzHashedVsXICPointsIsotopeShadows,
////            scanTimePredicted,
////            m_pythiaParameters.subtractShadows,
////            msFrame,
////            candidateScores
////            ); ree;

    ERR_RETURN
}

namespace {

    void filterFeatureFinderHills(
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
        eVec = EigenKernelUtils::convolveVectorWithKernel(eVec, d_ptr->gaussKernel);
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

Err CandidateScorertron::processPeakIntegrationIndexes(
        const TargetDecoyCandidatePair* targetDecoyCandidatePair,
        QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows
        ) {

    ERR_INIT

    QHash<MzHashed , QVector<FeatureFinderHill*>> mzHashedVsfeatureFinderHillsShadows;
    e = extractHills(
            ms2IonsTheoreticalIsotopeShadows,
            &mzHashedVsfeatureFinderHillsShadows
            ); ree;

    for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {

        e = ErrorUtils::isTrue(pii.second > pii.first); ree;

        QHash<MzHashed , QVector<FeatureFinderHill*>> mzHashedVsfeatureFinderHillsFiltered;
        filterFeatureFinderHills(
                mzHashedVsfeatureFinderHills,
                pii.first,
                pii.second,
                &mzHashedVsfeatureFinderHillsFiltered
                );

        QHash<MzHashed , QVector<FeatureFinderHill*>> mzHashedVsfeatureFinderHillsShadowsFiltered;
        filterFeatureFinderHills(
                mzHashedVsfeatureFinderHillsShadows,
                pii.first,
                pii.second,
                &mzHashedVsfeatureFinderHillsShadowsFiltered
                );
    }

    ERR_RETURN
}

