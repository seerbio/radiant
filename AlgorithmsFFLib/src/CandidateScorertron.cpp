//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron.h"

#include "CandidateScores.h"
#include "DiscriminantScoretron.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "IsotopicDistributionBuilder.h"
#include "ObjectCSVWriters.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"
#include "XICPeakManager.h"


class Q_DECL_HIDDEN CandidateScorertron::Private {
public:

    using NominalMass = int;
    using AveragineVec = QVector<float>;

    Private();
    ~Private();

    Err init(const PythiaParameters &pythiaParameters);

    Eigen::VectorX<float> m_kernelIntegration;
    Eigen::VectorX<float> m_kernelMs2;

private:
    PythiaParameters m_pythiaParameters;

};

CandidateScorertron::Private::Private() {}
CandidateScorertron::Private::~Private() {}

Err CandidateScorertron::Private::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_pythiaParameters = pythiaParameters;

    constexpr int order = 1;
    constexpr int derivative = 0;
    constexpr int rate = 1;

    //TODO consider writing algo to optimizing filter length parameters as well as smoothCount later in the code.

    Eigen::MatrixX<float> kernel;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        m_pythiaParameters.filterLengthMS2,
        order,
        derivative,
        rate,
        &kernel
        ); ree;
    const Eigen::VectorX<float> kernelVec(kernel);
    m_kernelMs2 = kernelVec;

    Eigen::MatrixX<float> kernelIntegration;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        m_pythiaParameters.filterLengthIntegration,
        order,
        derivative,
        rate,
        &kernelIntegration
        ); ree;
    const Eigen::VectorX<float> kernelIntegrationVec(kernelIntegration);
    m_kernelIntegration = kernelIntegrationVec;

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

CandidateScorertron::CandidateScorertron()
: m_topNMS2Ions(-1)
, m_xicPeakManager(nullptr)
, m_msFrameMzTarget(nullptr)
, m_msFrameMS1(nullptr)
, m_turboXicMS1(nullptr)
, d_ptr(QScopedPointer<Private>(new Private))
, m_minPeakCount(3.9)
, m_scanTimeRange(0)
, m_useExtendedScores(false)
, m_useNeuralNetworkScores(false)
, m_useTopNIntegrationsParam(false)
{}

CandidateScorertron::~CandidateScorertron() {}

Err CandidateScorertron::init(
    const PythiaParameters &pythiaParameters,
    const MsCalibratomatic &msCalibratomatic,
    const MzTargetKey &mzTargetKey,
    int topNMS2Ions,
    float minPeakCount,
    float scanTimeRange,
    const QMap<NominalMzMass, QVector<float>> &averagineTable,
    bool useExtendedScores,
    bool useNeuralNetworkScores,
    bool useTopNIntegrationsParameter,
    XICPeakManager *xicPeakManager,
    MsFrame *msFrameMzTarget,
    TurboXIC *turboXicMS1,
    MsFrame *msFrameMS1
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;
    e = ErrorUtils::isTrue(msFrameMzTarget->isValid()); ree;
    e = ErrorUtils::isTrue(turboXicMS1->isInit()); ree;
    e = ErrorUtils::isTrue(msFrameMS1->isValid()); ree;
    e = ErrorUtils::isNotEmpty(mzTargetKey); ree;
    e = ErrorUtils::isNotEmpty(averagineTable); ree;
    e = ErrorUtils::isFalse(MathUtils::tZero(scanTimeRange)); ree;
    e = ErrorUtils::isAboveThreshold(
        topNMS2Ions,
        S_GLOBAL_SETTINGS.MIN_MS2_IONS,
        ErrorUtilsParam::IncludeThreshold)
    ; ree;

    e = ErrorUtils::isAboveThreshold(
        minPeakCount,
        1.0f,
        ErrorUtilsParam::ExcludeThreshold
        ); ree;

    m_pythiaParameters = pythiaParameters;
    m_topNMS2Ions = topNMS2Ions;
    m_mzTargetKey = mzTargetKey;
    m_xicPeakManager = xicPeakManager;
    m_msFrameMzTarget = msFrameMzTarget;
    m_turboXicMS1 = turboXicMS1;
    m_scanTimeRange = scanTimeRange;
    m_averagineTable = averagineTable;
    m_msFrameMS1 = msFrameMS1;
    m_useExtendedScores = useExtendedScores;
    m_useNeuralNetworkScores = useNeuralNetworkScores;
    m_minPeakCount = minPeakCount;
    m_useTopNIntegrationsParam = useTopNIntegrationsParameter;

    if (msCalibratomatic.isInitRT()) {
        m_msCalibratomatic = msCalibratomatic;
    }

    e = d_ptr->init(m_pythiaParameters); ree;

    ERR_RETURN
}

class MatriciesAndVecs {

public:

    Eigen::MatrixX<float> intensityMatrix100;
    Eigen::MatrixX<float> intensityMatrix100Shadow;
    Eigen::MatrixX<float> intensityMatrix45;

    Eigen::MatrixX<float> mzMatrix100;

    Eigen::VectorX<float> intensityVec;
    Eigen::VectorX<float> ionCountVec;
    Eigen::VectorX<float> integrationVecCosineSim;
    Eigen::VectorX<float> productVec;

    [[nodiscard]] bool intensityMatriciesAreValid() const {
        return intensityMatrix100.size() > 0 && mzMatrix100.size() > 0;
    }

    [[nodiscard]] bool integrationVecIsValid() const {
        return ionCountVec.size() > 0;
    }
};

class BestCorrelationResult {

public:
    QVector<float> peakCorrelations;
    QVector<float> peakCorrelations45;
    QVector<float> peakCorrelationsWindow1p5X;
    QVector<float> peakCorrelationsWindow2X;
    float peakCorrelationsSum = -1.0;
    Eigen::MatrixX<float> matBlockTrimmedIntensity;
    Eigen::MatrixX<float> matBlockTrimmedIntensityWindow1p5X;
    Eigen::MatrixX<float> matBlockTrimmedIntensityWindow2X;
    Eigen::MatrixX<float> matBlockTrimmedIntensity45;
    Eigen::MatrixX<float> matBlockTrimmedIntensityShadows;
    Eigen::MatrixX<float> matBlockTrimmedMz;
    int bestAnchorColumnIndex = -1;
    int bestAnchorRowIndex = -1;
    QVector<int> apexStarts;
    PeakIntegrationIndexes peakIntegrationIndexes;

    [[nodiscard]] QVector<int> columnNonZeroPeakLengths() const {

        const int cols = static_cast<int>(matBlockTrimmedIntensity.cols());
        QVector<int> peakLengths(cols);

        for (int col = 0; col < cols; col++) {
            peakLengths[col] = EigenUtils::nonZeros(matBlockTrimmedIntensity.col(col));
        }

        return peakLengths;
    }
};

namespace {

    Err sortBestCorrelationResult(QVector<BestCorrelationResult> *bestCorrelationResults) {

        ERR_INIT

        e = ErrorUtils::isFalse(bestCorrelationResults->isEmpty()); ree;

        const float bestScore = std::max_element(
            bestCorrelationResults->begin(),
            bestCorrelationResults->end(),
            [](const BestCorrelationResult &l, const BestCorrelationResult &r){return l.peakCorrelationsSum < r.peakCorrelationsSum;}
            )->peakCorrelationsSum;

        const auto terminatorLogic = [bestScore](const BestCorrelationResult &r) {
            constexpr float tolerance = 1.0;
            return (bestScore - r.peakCorrelationsSum) > tolerance;
        };
        const auto terminator = std::remove_if(
            bestCorrelationResults->begin(),
            bestCorrelationResults->end(),
            terminatorLogic
            );

        bestCorrelationResults->erase(terminator, bestCorrelationResults->end());

        std::sort(
            bestCorrelationResults->begin(),
            bestCorrelationResults->end(),
            [](const BestCorrelationResult &l, const BestCorrelationResult &r) {

// #define USE_PEAK_LENGTH
#ifdef USE_PEAK_LENGTH
                const int lLen = l.peakIntegrationIndexes.second - l.peakIntegrationIndexes.first;
                const int rLen = r.peakIntegrationIndexes.second - r.peakIntegrationIndexes.first;

                if (lLen == rLen) {
                    return l.matBlockTrimmedIntensity.coeff(l.bestAnchorRowIndex, l.bestAnchorColumnIndex)
                         > r.matBlockTrimmedIntensity.coeff(r.bestAnchorRowIndex, r.bestAnchorColumnIndex);
                }

                return lLen > rLen;
#else
                return l.matBlockTrimmedIntensity.sum() > r.matBlockTrimmedIntensity.sum();
#endif

            });

        ERR_RETURN;
    }

}//namespace
Err CandidateScorertron::calculateScores(
    const QVector<MS2Ion> &ms2Ions,
    const QVector<float> &weights,
    TargetDecoyCandidatePair* targetDecoyCandidatePair,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2Ions); ree;

    candidateScores->targetDecoyCandidatePair = targetDecoyCandidatePair;
    candidateScores->initFeaturesArray();
    candidateScores->targetKey = m_mzTargetKey;

    FrameIndex frameIndexPredictedMin;
    FrameIndex frameIndexPredictedMax;
    e = setPredictedFrameIndexes(
        targetDecoyCandidatePair->iRt(),
        candidateScores,
        &frameIndexPredictedMin,
        &frameIndexPredictedMax
        );

    MatriciesAndVecs matriciesAndVecs;
    e = initMatricesdAndVecs(
        ms2Ions,
        frameIndexPredictedMin,
        frameIndexPredictedMax,
        &matriciesAndVecs
        ); ree;

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
    e = EigenUtils::simpleIntegrator(
        matriciesAndVecs.productVec,
        m_pythiaParameters.stopThresholdFractionMS2,
        &peakIntegrationsVsIntensities
        ); ree;

    if (peakIntegrationsVsIntensities.isEmpty()) {
        ERR_RETURN
    }

    QVector<BestCorrelationResult> bestCorrelationResults;
    e = processIntegrationVectorPeakIntegrations(
        matriciesAndVecs,
        peakIntegrationsVsIntensities,
        &bestCorrelationResults
        ); ree;

    constexpr int multiplierForKeySettingByTen = 10;
    const int nominalMass
        = static_cast<int>((std::round(targetDecoyCandidatePair->mass() / multiplierForKeySettingByTen) * multiplierForKeySettingByTen));
    e = ErrorUtils::isTrue(m_averagineTable.contains(nominalMass)); ;
    const QVector<float> ms1Averagine = m_averagineTable.value(nominalMass);

// #define TURBO_SCORE
#ifdef TURBO_SCORE
    e = sortBestCorrelationResult(&bestCorrelationResults); ree;

    e = setCandidateScores(
        targetDecoyCandidatePair,
        bestCorrelationResults,
        ms1Averagine,
        candidateScores
        ); ree;
#else
    QVector<FeaturesArray> candidateScoresFeatureArrays;
    QVector<CandidateScores> candidateScoresFeatures;
    for (const BestCorrelationResult &bcr : bestCorrelationResults) {
        CandidateScores cs = *candidateScores;
        e = setCandidateScores(
            targetDecoyCandidatePair,
            {bcr},
            ms1Averagine,
            &cs
            ); ree;

        const FeaturesArray fa = DiscriminantScoretron::scoreVectorLogic(
            m_useExtendedScores,
            m_useNeuralNetworkScores,
            &cs);

        candidateScoresFeatures.push_back(cs);
        candidateScoresFeatureArrays.push_back(fa);
    }

    QVector<FeaturesArray*> featuresArraysPntrs;
    std::transform(
        candidateScoresFeatureArrays.begin(),
        candidateScoresFeatureArrays.end(),
        std::back_inserter(featuresArraysPntrs),
        [](FeaturesArray &f){return &f;}
        );

    constexpr int threadCount = 1;
    QVector<float> discScores;
    e = DiscriminantScoretron::applyWeights(
        weights,
        threadCount,
        featuresArraysPntrs,
        &discScores
        );

    QVector<QPair<float, CandidateScores>> candidateScoresPairs;
    for (int i = 0; i < discScores.size(); ++i) {
        candidateScoresPairs.push_back({discScores[i], candidateScoresFeatures[i]});
    }

    const float maxDisScore =  *std::max_element(discScores.begin(), discScores.end());
    const int bestIndex = MathUtils::closest(discScores, maxDisScore);

    *candidateScores = candidateScoresFeatures[bestIndex];
#endif

// #define TROUBLE_SHOOT_INTEGRATION
#ifdef TROUBLE_SHOOT_INTEGRATION
    if (targetDecoyCandidatePair->peptideStringWithMods() == "TVC(UniMod:4)LPDGSFPSGSEC(UniMod:4)HISGWGVTETGK"
        && targetDecoyCandidatePair->charge() == 3
        && !candidateScores->isDecoy
        ) {

        QMap<FrameIndex, BestCorrelationResult> bestCorrelationResultsInds;
        for (const BestCorrelationResult &bcr : bestCorrelationResults) {
            bestCorrelationResultsInds.insert(bcr.peakIntegrationIndexes.first, bcr);
        }

        std::sort(
            candidateScoresPairs.rbegin(),
            candidateScoresPairs.rend(),
            [](const QPair<float, CandidateScores> &l, const QPair<float, CandidateScores> &r){return l.first < r.first;}
            );

        for (QPair<float, CandidateScores> &b : candidateScoresPairs) {
            qDebug()
            << b.second.frameIndexStart << ","
            << b.second.frameIndexEnd
            << b.first
            << DiscriminantScoretron::scoreVectorLogic(true, false, &b.second)
            << bestCorrelationResultsInds.value(b.second.frameIndexStart).apexStarts
            << "SDLKFJSDL";

            // std::cout << bestCorrelationResultsInds.value(b.second.frameIndexStart).matBlockTrimmedIntensity << std::endl;;
            // std::cout << " **************** " << std::endl;

        }

        qDebug() << weights;

        const QString &intensityVecPath = QStringLiteral("/home/andrewnichols/Repos/Graphing/intensity.csv");
        const QString &prodVecPath = QStringLiteral("/home/andrewnichols/Repos/Graphing/prod.csv");

        e = ObjectCSVWriters::writeVectorToFile(
            EigenUtils::convertEigenVectorToQVector(matriciesAndVecs.intensityVec),
            intensityVecPath
            ); ree;

        e = ObjectCSVWriters::writeVectorToFile(
            EigenUtils::convertEigenVectorToQVector(matriciesAndVecs.productVec),
            prodVecPath
            ); ree;
    }
#endif

    ERR_RETURN
}

namespace {

    void filterXICPointsByAccuracyPPM(
        float mzVal,
        float ppmTol,
        XICPoints *xicPoints
        ) {

        const float massTol = MathUtils::calculatePPM(mzVal, ppmTol);

        const auto terminatorLogic = [mzVal, massTol](const XICPoint &p) {
            return !(mzVal - massTol < p.mz && p.mz < mzVal + massTol);
        };

        const auto terminator = std::remove_if(
            xicPoints->begin(),
            xicPoints->end(),
            terminatorLogic
            );

        xicPoints->erase(terminator, xicPoints->end());
    }

    void filterXICPointsByFrameIndex(
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        XICPoints *xicPoints
        ) {

        const auto terminatorLogic = [frameIndexPredictedMin, frameIndexPredictedMax](const XICPoint &p) {
            return !(frameIndexPredictedMin < p.scanNumber && p.scanNumber < frameIndexPredictedMax);
        };

        const auto terminator = std::remove_if(
            xicPoints->begin(),
            xicPoints->end(),
            terminatorLogic
            );

        xicPoints->erase(terminator, xicPoints->end());
    }

    Err getXICs(
        const QVector<MS2Ion> &ms2Ions,
        float ppmTol,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        XICPeakManager *xicPeakManager,
        QVector<XICPoints> *xicPointsVec100,
        QVector<XICPoints> *xicPointsVec100Shadows,
        QVector<XICPoints> *xicPointsVec45
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;

        xicPointsVec100->reserve(ms2Ions.size());
        xicPointsVec100Shadows->reserve(ms2Ions.size());
        xicPointsVec45->reserve(ms2Ions.size());

        for (int i = 0; i < ms2Ions.size(); i++) {

            const MS2Ion &ms2Ion = ms2Ions.at(i);

            XICPoints xicPoints;
            e = xicPeakManager->getXIC(ms2Ion.mz, &xicPoints); ree;

            if (frameIndexPredictedMax > 0) {
                filterXICPointsByFrameIndex(
                frameIndexPredictedMin,
                frameIndexPredictedMax,
                &xicPoints
                );
            }

            if (xicPoints.empty()) {
                xicPointsVec100->push_back({});
                xicPointsVec45->push_back({});
                xicPointsVec100Shadows->push_back({});
                continue;
            }

            XICPoints xicPointsShadows;
            const float isotopeDistanceThomsons = S_GLOBAL_SETTINGS.ISO_DIFF / ms2Ion.charge;
            e = xicPeakManager->getXIC(ms2Ion.mz - isotopeDistanceThomsons, &xicPointsShadows); ree;
            if (xicPointsShadows.empty()) {
                xicPointsVec100Shadows->push_back({});
            }
            else {
                if (frameIndexPredictedMax > 0) {
                    filterXICPointsByFrameIndex(
                        frameIndexPredictedMin,
                        frameIndexPredictedMax,
                        &xicPointsShadows
                        );
                }
                xicPointsVec100Shadows->push_back(xicPointsShadows);
            }

            xicPointsVec100->push_back(xicPoints);

            filterXICPointsByAccuracyPPM(
                ms2Ion.mz,
                ppmTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
                &xicPoints
                );
            xicPointsVec45->push_back(xicPoints);

        }

        ERR_RETURN
    }

    FrameIndex findFrameIndexMaxXICPointsVec (const QVector<XICPoints> &xicPointsVec100) {

        FrameIndex frameIndexMax = -1;

        for (const XICPoints &xicPoints : xicPointsVec100) {

            if (xicPoints.empty()) {
                continue;
            }

            const FrameIndex frameIndexMaxXICPoints = std::max_element(
                xicPoints.begin(),
                xicPoints.end(),
                [](const XICPoint &l, const XICPoint &r){return l.scanNumber < r.scanNumber;}
                )->scanNumber;
            frameIndexMax = std::max(frameIndexMax, frameIndexMaxXICPoints);
        }

        return frameIndexMax;
    }

    Err buildEigenMatrix(
        const QVector<XICPoints> &xicPointsVec,
        const Eigen::VectorX<float> &kernelMs2,
        FrameIndex frameIndexMax,
        bool buildMzMatrix,
        int smoothCount,
        Eigen::MatrixX<float> *matIntensity,
        Eigen::MatrixX<float> *matMz
        ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(xicPointsVec); ree;

        const FrameIndex frameIndexBuffer = 2;

        const int rows = frameIndexMax + frameIndexBuffer;

        matIntensity->resize(rows, xicPointsVec.size());
        matIntensity->setZero();

        matMz->resize(rows, xicPointsVec.size());
        matMz->setZero();

        for (int col = 0; col < xicPointsVec.size(); col++) {

            const XICPoints &xicPointsCol = xicPointsVec.at(col);
            for (const XICPoint &p : xicPointsCol) {

                if (p.scanNumber >= rows) {
                    break;
                }

                matIntensity->coeffRef(p.scanNumber, col) = p.intensity;
                if (buildMzMatrix) {
                    matMz->coeffRef(p.scanNumber, col) = p.mz;
                }
            }
        }

        for (int smooth = 0; smooth < smoothCount; smooth++) {
            *matIntensity= EigenKernelUtils::applyKernelToEachColumnInMatrix(*matIntensity, kernelMs2);
        }

        ERR_RETURN
    }

    Err buildIntegrationVector(
        const MatriciesAndVecs &matriciesAndVecs,
        const Eigen::VectorX<float> &kernelIntegration,
        float minPeakCount,
        int maxAnchorColumnIndex,
        Eigen::VectorX<float> *ionCountVec
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;

        Eigen::MatrixX<float> matCount = matriciesAndVecs.intensityMatrix100.leftCols(
                std::min(maxAnchorColumnIndex,
                static_cast<int>(matriciesAndVecs.intensityMatrix100.cols())
                ));

        constexpr float intensityThresholdVal = 0.1;
        constexpr float countValue = 1.0;

        matCount = (matCount.array() > intensityThresholdVal).select(countValue, matCount);
        EigenUtils::thresholdMatrix(0.0f, &matCount);

        Eigen::VectorX<float> integrationVecLocal = matCount.rowwise().sum();
        EigenUtils::thresholdVector(minPeakCount, &integrationVecLocal);

        *ionCountVec = EigenKernelUtils::convolveVectorWithKernel(
        integrationVecLocal,
        kernelIntegration
        );

        ERR_RETURN
    }

    Err buildIntegrationVectorCosineSim(
        const MatriciesAndVecs &matriciesAndVecs,
        const QVector<MS2Ion> &ms2IonsTheo,
        const Eigen::VectorX<float> &kernelIntegration,
        int maxAnchorColumnIndex,
        Eigen::VectorX<float> *integrationVecCosineSim
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;
        e = ErrorUtils::isNotEmpty(ms2IonsTheo); ree;
        e = ErrorUtils::isEqual(ms2IonsTheo.size(), static_cast<int>(matriciesAndVecs.intensityMatrix100.cols())); ree;

        const int columnCount = std::min(maxAnchorColumnIndex, static_cast<int>(matriciesAndVecs.intensityMatrix100.cols()));

        Eigen::MatrixX<float> matIntensityVals = matriciesAndVecs.intensityMatrix100.leftCols(columnCount);

        const int matRows = static_cast<int>(matIntensityVals.rows());
        Eigen::MatrixX<Intensity> matMs2IonsIntensityVals(matRows, matIntensityVals.cols());

        for (int i = 0; i < columnCount; i++) {
            const MS2Ion &ms2Ion = ms2IonsTheo.at(i);
            const QVector<Intensity> ms2IonIntensityCol(matRows, ms2Ion.intensity);

            const Eigen::VectorX<Intensity> v = EigenUtils::convertQVectorToEigenVector(ms2IonIntensityCol);
            matMs2IonsIntensityVals.col(i) = v;
        }

        Eigen::VectorX<float> integrationVecCosineSimLocal;
        e = EigenUtils::rowWiseCosineSimilarOfMatrices(
            matIntensityVals,
            matMs2IonsIntensityVals,
            &integrationVecCosineSimLocal
            ); ree;

        integrationVecCosineSimLocal = integrationVecCosineSimLocal.array().pow(1.0f/3.0f);
        *integrationVecCosineSim = EigenKernelUtils::convolveVectorWithKernel(
            integrationVecCosineSimLocal,
            kernelIntegration
            );

        ERR_RETURN

    }

}//namespace
Err CandidateScorertron::initMatricesdAndVecs(
        const QVector<MS2Ion> &ms2Ions,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        MatriciesAndVecs *matriciesAndVecs
        ) const {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;
        e = ErrorUtils::isAboveThreshold(m_topNMS2Ions, 0, ErrorUtilsParam::ExcludeThreshold); ree;

        QVector<MS2Ion> ms2IonsResized = ms2Ions;
        ms2IonsResized.resize(std::min(m_topNMS2Ions, ms2Ions.size()));

        const bool ms2IonsAreSorted = std::is_sorted(
            ms2IonsResized.rbegin(),
            ms2IonsResized.rend(),
            [](const MS2Ion &l, const MS2Ion &r){return l.intensity < r.intensity;}
            );
        e = ErrorUtils::isTrue(ms2IonsAreSorted); ree;

        QVector<XICPoints> xicPointsVec100;
        QVector<XICPoints> xicPointsVec100Shadow;
        QVector<XICPoints> xicPointsVec45;
        e = getXICs(
            ms2IonsResized,
            static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM),
            frameIndexPredictedMin,
            frameIndexPredictedMax,
            m_xicPeakManager,
            &xicPointsVec100,
            &xicPointsVec100Shadow,
            &xicPointsVec45
            ); ree;

        e = ErrorUtils::isEqual(xicPointsVec100.size(), xicPointsVec45.size()); ree;
        e = ErrorUtils::isEqual(xicPointsVec100.size(), xicPointsVec100Shadow.size()); ree;

        const FrameIndex frameIndexMax = findFrameIndexMaxXICPointsVec(xicPointsVec100);

        constexpr int smoothCount = 1;
        e = buildEigenMatrix(
            xicPointsVec100,
            d_ptr->m_kernelMs2,
            frameIndexMax,
            true,
            smoothCount,
            &matriciesAndVecs->intensityMatrix100,
            &matriciesAndVecs->mzMatrix100
            ); ree;
        matriciesAndVecs->intensityVec = matriciesAndVecs->intensityMatrix100.rowwise().sum();

        Eigen::MatrixX<float> unused;
        e = buildEigenMatrix(
            xicPointsVec100Shadow,
            d_ptr->m_kernelMs2,
            frameIndexMax,
            false,
            smoothCount,
            &matriciesAndVecs->intensityMatrix100Shadow,
            &unused
            ); ree;

        matriciesAndVecs->intensityMatrix100 = m_pythiaParameters.subtractShadows
                          ? matriciesAndVecs->intensityMatrix100 - matriciesAndVecs->intensityMatrix100Shadow
                          : matriciesAndVecs->intensityMatrix100;
        EigenUtils::thresholdMatrix(0.0f, &matriciesAndVecs->intensityMatrix100);

        e = buildIntegrationVector(
            *matriciesAndVecs,
            d_ptr->m_kernelIntegration,
            m_minPeakCount,
            m_pythiaParameters.maxAnchorColumnIndex,
            &matriciesAndVecs->ionCountVec
            ); ree;

        e = buildIntegrationVectorCosineSim(
            *matriciesAndVecs,
            ms2IonsResized,
            d_ptr->m_kernelIntegration,
            m_pythiaParameters.maxAnchorColumnIndex,
            &matriciesAndVecs->integrationVecCosineSim
            ); ree;

        matriciesAndVecs->productVec = matriciesAndVecs->ionCountVec.array()
                                     * matriciesAndVecs->integrationVecCosineSim.array();

        constexpr int noSmooths = 0;
        e = buildEigenMatrix(
            xicPointsVec45,
            d_ptr->m_kernelMs2,
            frameIndexMax,
            false,
            noSmooths,
            &matriciesAndVecs->intensityMatrix45,
            &unused
            ); ree;

        ERR_RETURN
    }

Err CandidateScorertron::setPredictedFrameIndexes(
    float iRT,
    CandidateScores *candidateScores,
    FrameIndex *frameIndexPredictedMin,
    FrameIndex *frameIndexPredictedMax
    ) const {

    ERR_INIT

    if (m_msCalibratomatic.isInitRT()) {

        const float scanTimeWindow
            = m_msCalibratomatic.scanTimeStDev(static_cast<float>(m_pythiaParameters.scanTimeWindowStDevs));

        e = m_msCalibratomatic.predictScanTime(
                iRT,
                &candidateScores->scanTimePredicted
        ); ree;

        e = m_msFrameMzTarget->frameIndexFromScanTime(
                candidateScores->scanTimePredicted - scanTimeWindow,
                frameIndexPredictedMin
        ); ree;

        e = m_msFrameMzTarget->frameIndexFromScanTime(
                candidateScores->scanTimePredicted + scanTimeWindow,
                frameIndexPredictedMax
        ); ree;

        ERR_RETURN
    }

    *frameIndexPredictedMin = 0;
    *frameIndexPredictedMax = 0;

    ERR_RETURN
}

namespace {

    QPair<PeakIntegrationIndexes, Intensity> correctPeakIntegrationForSingleRow(
        const QPair<PeakIntegrationIndexes, Intensity> &pii,
        const Eigen::VectorX<float> &ionCountVec
        ) {
        QPair<PeakIntegrationIndexes, Intensity> piiWorking = pii;
        if (piiWorking.first.first == piiWorking.first.second) {
            constexpr int bufferDistance = 1;
            piiWorking.first.first = std::max(0, piiWorking.first.first - bufferDistance);
            piiWorking.first.second = std::min(
                static_cast<int>(ionCountVec.size()) - 1,
                piiWorking.first.second + bufferDistance
                );
        }
        return piiWorking;
    }

    QVector<QVector<int>> getMatrxColumnApexes(const Eigen::MatrixX<float> &matBlock) {

        QVector<QVector<int>> apexIndexesByColumn;
        apexIndexesByColumn.reserve(matBlock.cols());
        for (int col = 0; col < matBlock.cols(); col++) {
            const Eigen::VectorX<float> &column = matBlock.col(col);
            apexIndexesByColumn.push_back(EigenUtils::apexesIndexesOnly(column));
        }

        return apexIndexesByColumn;
    }

    Eigen::MatrixX<float> buildMatBlockApexes(
        const Eigen::MatrixX<float> &matBlock,
        const QVector<QVector<int>> &apexIndexesByColumn
        ) {

        Eigen::MatrixX<float> matBlockApexes(matBlock.rows(), matBlock.cols());
        matBlockApexes.setZero();
        for (int col = 0; col < apexIndexesByColumn.size(); col++) {
            const QVector<int> &columnApexes = apexIndexesByColumn.at(col);
            for (int row : columnApexes) {
                matBlockApexes.coeffRef(row, col) = matBlock.coeff(row, col);
            }
        }
        return matBlockApexes;
    }

    QVector<int> findStartApexes(
            const Eigen::MatrixX<float> &matBlockApexes,
            int apexIndex
            ) {

            QVector<int> apexesStartsToUse(static_cast<int>(matBlockApexes.cols()), -1);
            for (int col = 0; col < matBlockApexes.cols(); col++) {

                const Eigen::VectorX<float> &apexColumn = matBlockApexes.col(col);

                if (apexColumn.coeff(apexIndex) > 0) {
                    apexesStartsToUse[col] = apexIndex;
                    continue;
                }

                for (int rowFromCenter = 1; rowFromCenter < apexColumn.size(); rowFromCenter++) {

                    const int rowLeftIndex = std::max(apexIndex - rowFromCenter, 0);
                    const int rowRightIndex = std::min(apexIndex + rowFromCenter, static_cast<int>(apexColumn.size() - 1));

                    int rowLeftIndexValue = -1;
                    int rowRightIndexValue = -1;

                    if (rowLeftIndex >= 0) {
                        rowLeftIndexValue = static_cast<int>(std::ceil(apexColumn.coeff(rowLeftIndex)));
                    }

                    if (rowRightIndex < apexColumn.size()) {
                        rowRightIndexValue = static_cast<int>(std::ceil(apexColumn.coeff(rowRightIndex)));
                    }

                    if (rowLeftIndexValue > 0 && rowRightIndexValue > 0) {
                        const int higherIndex = rowLeftIndexValue >= rowRightIndexValue ? rowLeftIndex : rowRightIndex;
                        apexesStartsToUse[col] = higherIndex;
                        break;
                    }
                    else if (rowLeftIndexValue > 0) {
                        apexesStartsToUse[col] = rowLeftIndex;
                        break;
                    }
                    else if (rowRightIndexValue > 0) {
                        apexesStartsToUse[col] = rowRightIndex;
                        break;
                    }

                    if (rowLeftIndex <= 0 && rowRightIndex >= apexColumn.size() - 1) {
                        apexesStartsToUse[col] = -1;
                        break;
                    }
                }
            }

            return apexesStartsToUse;
        }

    QPair<int, int> simpleIntegratorTrimmer(
        const Eigen::VectorX<float> &vec,
        int apexRowIndex,
        float stopThresholdFraction
        ) {

        const float nearZero = 0.001;
        const float apexIndexValue = vec.coeff(apexRowIndex);
        const float stopThresholdValue = apexIndexValue * stopThresholdFraction;

        int rightStopIndex = apexRowIndex;
        int rightCurrentIndex = apexRowIndex;
        while (rightCurrentIndex < vec.size()) {

            const float currentValue = vec.coeff(rightCurrentIndex);
            if (currentValue < nearZero || currentValue <= stopThresholdValue) {
                rightStopIndex = rightCurrentIndex;
                break;
            }

            if (currentValue > nearZero) {
                rightStopIndex = rightCurrentIndex;
                rightCurrentIndex++;
                continue;
            }

            break;
        }

        int leftStopIndex = apexRowIndex;
        int leftCurrentIndex = apexRowIndex;
        while (leftCurrentIndex >= 0) {

            const float currentValue = vec(leftCurrentIndex);
            if (currentValue < nearZero || currentValue <= stopThresholdValue) {
                leftStopIndex = leftCurrentIndex;
                break;
            }

            if (currentValue > nearZero) {
                leftStopIndex = leftCurrentIndex;
                leftCurrentIndex--;
                continue;
            }

            break;
        }

        return {std::max(leftStopIndex, 0), std::min(rightStopIndex, static_cast<int>(vec.size()) - 1)};
    }

    Eigen::MatrixX<float> trimMatrixBlock(
        const Eigen::MatrixX<float> &matBlock,
        const QVector<int> &apexStarts,
        float stopThresholdFraction
        ) {

        Eigen::MatrixX<float> matBlockTrimmedColumns(matBlock.rows(), matBlock.cols());
        matBlockTrimmedColumns.setZero();

        for (int col = 0; col < matBlock.cols(); col++) {

            const Eigen::VectorX<float> &colVec = matBlock.col(col);
            const QPair<int, int> peakLimits = simpleIntegratorTrimmer(
                colVec,
                apexStarts.at(col),
                stopThresholdFraction
                );

            for(int row = peakLimits.first; row <= peakLimits.second; row++) {
                matBlockTrimmedColumns.coeffRef(row, col) = colVec.coeff(row);
            }
        }
        return matBlockTrimmedColumns;
    }

    Err findBestAnchorColumn(
        const Eigen::MatrixX<float> &matBlockTrimmed,
        const QVector<int> &apexStarts,
        int maxAnchorColumnIndex,
        QVector<float> *peakCorrelations,
        int *bestAnchorColumnIndex,
        int *bestAnchorRowIndex
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matBlockTrimmed.size() > 0); ree;
        e = ErrorUtils::isNotEmpty(apexStarts); ree;
        e = ErrorUtils::isEqual(apexStarts.size(), static_cast<int>(matBlockTrimmed.cols())); ree;

        const int colCount = std::min(static_cast<int>(matBlockTrimmed.cols()), maxAnchorColumnIndex);

        float bestCosineSimSum = 0.0;
        for (int anchorCol = 0; anchorCol < colCount; anchorCol++) {

            const Eigen::VectorX<float> &anchor = matBlockTrimmed.col(anchorCol);

            QVector<float> corrs;
            corrs.reserve(colCount);
            float cosineSimSum = 0.0;
            for (int col = 0; col < colCount; col++) {
                const Eigen::VectorX<float> &c = matBlockTrimmed.col(col);

                float cosineSim;
                e = EigenUtils::cosineSimilarity(anchor, c, &cosineSim); ree;
                cosineSimSum += cosineSim / (col + 1.0);
                corrs.push_back(cosineSim);
            }

            if (cosineSimSum > bestCosineSimSum) {
                *bestAnchorColumnIndex = anchorCol;
                bestCosineSimSum = cosineSimSum;
            }
        }

        const Eigen::VectorX<float> &anchor = matBlockTrimmed.col(*bestAnchorColumnIndex);
        const QVector<float> anchorVec = EigenUtils::convertEigenVectorToQVector(anchor);

        QVector<float> corrs;

        const int colCountFull = matBlockTrimmed.cols();
        corrs.reserve(colCountFull);
        for (int col = 0; col < colCountFull; col++) {
            const Eigen::VectorX<float> &c = matBlockTrimmed.col(col);

            float cosineSim;
            e = EigenUtils::cosineSimilarity(anchor, c, &cosineSim); ree;

            corrs.push_back(cosineSim);
        }
        *peakCorrelations = corrs;

        *bestAnchorRowIndex = MathUtils::closest(anchorVec, anchor.maxCoeff());

        ERR_RETURN
    }

    Err calculatePeakCorrelations(
        const Eigen::MatrixX<float> &matBlockTrimmed,
        int bestAnchorColumnIndex,
        QVector<float> *peakCorrelations
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matBlockTrimmed.size() > 0); ree;
        e = ErrorUtils::isAboveThreshold(bestAnchorColumnIndex, 0, ErrorUtilsParam::IncludeThreshold); ree;

        peakCorrelations->clear();

        const Eigen::VectorX<float> &anchorColumn = matBlockTrimmed.col(bestAnchorColumnIndex);

        peakCorrelations->reserve(static_cast<int>(matBlockTrimmed.cols()));
        for (int col = 0; col < matBlockTrimmed.cols(); col++) {

            if (col == bestAnchorColumnIndex) {
                peakCorrelations->push_back(1.0);
                continue;
            }

            const Eigen::VectorX<float> &matCol = matBlockTrimmed.col(col);

            float cosineSim;
            e = EigenUtils::cosineSimilarity(matCol, anchorColumn, &cosineSim); ree;
            peakCorrelations->push_back(cosineSim);
        }

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::processIntegrationVectorPeakIntegrations(
    const MatriciesAndVecs &matriciesAndVecs,
    const QVector<QPair<PeakIntegrationIndexes, Intensity>> &peakIntegrationsVsIntensity,
    QVector<BestCorrelationResult> *bestCorrelationResults
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;
    e = ErrorUtils::isTrue(matriciesAndVecs.integrationVecIsValid()); ree;

    const int maxRows = static_cast<int>(matriciesAndVecs.intensityMatrix100.rows());
    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensityResized = peakIntegrationsVsIntensity;

    if (m_useTopNIntegrationsParam) {
        peakIntegrationsVsIntensityResized.resize(std::min(
            m_pythiaParameters.topNIntegrations,
            peakIntegrationsVsIntensityResized.size()
            ));
    }

    for (const QPair<PeakIntegrationIndexes, Intensity> &pii : peakIntegrationsVsIntensityResized) {

        const QPair<PeakIntegrationIndexes, Intensity> piiWorking = correctPeakIntegrationForSingleRow(
            pii,
            matriciesAndVecs.ionCountVec
            );

        const int ogPeakLength = piiWorking.first.second - piiWorking.first.first + 1;

        Eigen::MatrixX<float> matBlock = matriciesAndVecs.intensityMatrix100.block(
              piiWorking.first.first,
              0,
              ogPeakLength,
              matriciesAndVecs.intensityMatrix100.cols()
              ).eval();

        const QVector<QVector<int>> apexIndexesByColumn = getMatrxColumnApexes(matBlock);
        const Eigen::MatrixX<float> matBlockApexes = buildMatBlockApexes(
            matBlock,
            apexIndexesByColumn
            );

        const Eigen::VectorX<float> integrationVecSegment = matriciesAndVecs.ionCountVec.segment(
            piiWorking.first.first,
            ogPeakLength
            ).eval();

        const QPair<int, float> apexIndex = EigenUtils::returnTopIndexAndValue(integrationVecSegment);

        const QVector<int> apexStarts = findStartApexes(matBlockApexes, apexIndex.first);
        e = ErrorUtils::isEqual(apexStarts.size(), static_cast<int>(matBlockApexes.cols())); ree;

        constexpr float stopThresholdFraction = 0.0;
        const Eigen::MatrixX<float> matBlockTrimmed = trimMatrixBlock(
            matBlock,
            apexStarts,
            stopThresholdFraction
            );

        QVector<float> peakCorrelations;
        int bestAnchorColumnIndex = -1;
        int bestAnchorRowIndex = -1;
        e = findBestAnchorColumn(
            matBlockTrimmed,
            apexStarts,
            m_pythiaParameters.maxAnchorColumnIndex,
            &peakCorrelations,
            &bestAnchorColumnIndex,
            &bestAnchorRowIndex
            ); ree;

        const int topN = std::min(6, peakCorrelations.size());
        const float peakCorrelationsSum
                = std::accumulate(peakCorrelations.begin(), peakCorrelations.begin() + topN, 0.0f);

        BestCorrelationResult bestCorrelationResultPii;

        bestCorrelationResultPii.peakCorrelations = peakCorrelations;
        bestCorrelationResultPii.peakCorrelationsSum = peakCorrelationsSum;
        bestCorrelationResultPii.matBlockTrimmedIntensity = matBlockTrimmed;
        bestCorrelationResultPii.bestAnchorColumnIndex = bestAnchorColumnIndex;
        bestCorrelationResultPii.bestAnchorRowIndex = bestAnchorRowIndex;
        bestCorrelationResultPii.apexStarts = apexStarts;
        bestCorrelationResultPii.peakIntegrationIndexes = piiWorking.first;

        constexpr float windowMultiplier1p5X = 0.5;
        const auto frameIndex1p5XMin
            = std::max(static_cast<FrameIndex>(std::round(piiWorking.first.first - (windowMultiplier1p5X * ogPeakLength))), 0);
        const auto frameIndex1p5XMax
            = std::min(static_cast<FrameIndex>(std::round(piiWorking.first.first + (windowMultiplier1p5X * ogPeakLength))), maxRows);
        const int peakLength1p5X = frameIndex1p5XMax - frameIndex1p5XMin;
        Eigen::MatrixX<float> matBlock1p5X = matriciesAndVecs.intensityMatrix100.block(
              frameIndex1p5XMin,
              0,
              peakLength1p5X,
              matriciesAndVecs.intensityMatrix100.cols()
              ).eval();

        constexpr float windowMultiplier2X = 1.0;
        const auto frameIndex2XMin
            = std::max(static_cast<FrameIndex>(std::round(piiWorking.first.first - (windowMultiplier2X * ogPeakLength))), 0);
        const auto frameIndex2XMax
            = std::min(static_cast<FrameIndex>(std::round(piiWorking.first.first + (windowMultiplier2X * ogPeakLength))), maxRows);
        const int peakLength2X = frameIndex2XMax - frameIndex2XMin;
        Eigen::MatrixX<float> matBlock2X = matriciesAndVecs.intensityMatrix100.block(
                      frameIndex2XMin,
                      0,
                      peakLength2X,
                      matriciesAndVecs.intensityMatrix100.cols()
                      ).eval();

        bestCorrelationResultPii.matBlockTrimmedIntensityWindow1p5X = trimMatrixBlock(
            matBlock1p5X,
            apexStarts,
            stopThresholdFraction
            );

        e = calculatePeakCorrelations(
            bestCorrelationResultPii.matBlockTrimmedIntensityWindow1p5X,
            bestAnchorColumnIndex,
            &bestCorrelationResultPii.peakCorrelationsWindow1p5X
            ); ree;

        bestCorrelationResultPii.matBlockTrimmedIntensityWindow2X = trimMatrixBlock(
                    matBlock2X,
                    apexStarts,
                    stopThresholdFraction
                    );

        e = calculatePeakCorrelations(
            bestCorrelationResultPii.matBlockTrimmedIntensityWindow2X,
            bestAnchorColumnIndex,
            &bestCorrelationResultPii.peakCorrelationsWindow2X
            ); ree;

        const PeakIntegrationIndexes &p = piiWorking.first;
        const int pSize = p.second - p.first + 1;
        bestCorrelationResultPii.matBlockTrimmedMz = matriciesAndVecs.mzMatrix100.block(
            p.first,
            0,
            pSize,
            matBlockTrimmed.cols()
            ).eval();

        bestCorrelationResultPii.matBlockTrimmedIntensityShadows = matriciesAndVecs.intensityMatrix100Shadow.block(
            p.first,
            0,
            pSize,
            matBlockTrimmed.cols()
            ).eval();

        bestCorrelationResultPii.matBlockTrimmedIntensity45 = matriciesAndVecs.intensityMatrix45.block(
            p.first,
            0,
            pSize,
            matBlockTrimmed.cols()
            ).eval();

        e = calculatePeakCorrelations(
            bestCorrelationResultPii.matBlockTrimmedIntensity45,
            bestAnchorColumnIndex,
            &bestCorrelationResultPii.peakCorrelations45
            ); ree;

        bestCorrelationResults->push_back(bestCorrelationResultPii);

// #define OUTPUT_MATS
#ifdef OUTPUT_MATS
        qDebug() << peakCorrelations << std::accumulate(peakCorrelations.begin(), peakCorrelations.end(), 0.0f) << "SLDFJSDLFLJ";
        for (int as : apexStarts) {
            std::cout << as << " " << apexIndex.first << std::endl;
        }
        std::cout << piiWorking.first.first << " " << piiWorking.first.second << std::endl;
        std::cout << matBlock << std::endl;
        std::cout << bestAnchorColumnIndex << " **********" << std::endl;
#endif
    }

    ERR_RETURN
}

namespace {

    float calculatedCosineSimSumGreaterThan80(const QVector<float> &peakCorrelations) {

        const int top6 = std::min(6, peakCorrelations.size());
        constexpr float cosineSimToAnchorThreshold = 0.8;
        const float scoreGreater = std::accumulate(
            peakCorrelations.begin(),
            peakCorrelations.begin() + top6,
            0.0f,
            [cosineSimToAnchorThreshold](float sum, float f){return f > cosineSimToAnchorThreshold ? f + sum : 0.0f + sum;}
            );

        return scoreGreater;
    }

    Err setCosineSimilarityMetrics(
        const TargetDecoyCandidatePair *targetDecoyCandidatePair,
        const BestCorrelationResult& bestCorrelationResult,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        //TODO add kldiv and pearsons corr.

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.size() > 0); ree;
        QVector<MS2Ion> ms2IonsTheo = candidateScores->isDecoy
                                    ? targetDecoyCandidatePair->ms2IonsDecoy()
                                    : targetDecoyCandidatePair->ms2IonsTarget();
        ms2IonsTheo.resize(static_cast<int>(bestCorrelationResult.matBlockTrimmedIntensity.cols()));

        QVector<float> intensitiesTheo;
        std::transform(
            ms2IonsTheo.begin(),
            ms2IonsTheo.end(),
            std::back_inserter(intensitiesTheo),
            [](const MS2Ion &ms2Ion){return ms2Ion.intensity;}
            );

        const Eigen::VectorX<float> intensitiesTheoVec = EigenUtils::convertQVectorToEigenVector(intensitiesTheo);

        Eigen::MatrixX<float> intensitiesTheoMat(
            bestCorrelationResult.matBlockTrimmedIntensity.rows(),
            bestCorrelationResult.matBlockTrimmedIntensity.cols()
        );

        for (int row = 0; row < intensitiesTheoMat.rows(); row++) {
            intensitiesTheoMat.row(row) = intensitiesTheoVec;
        }

        Eigen::VectorX<float> cosineSimsByRow;
         e = EigenUtils::rowWiseCosineSimilarOfMatrices(
             bestCorrelationResult.matBlockTrimmedIntensity,
             intensitiesTheoMat,
             &cosineSimsByRow
             ); ree;

        Eigen::MatrixX<float> matSummed = bestCorrelationResult.matBlockTrimmedIntensity;
        for (int row = 0; row < matSummed.rows(); row++) {
            Eigen::VectorX<float> r = matSummed.row(row).array() / std::max(matSummed.row(row).sum(), 0.001f);
            matSummed.row(row) = r;
        }

        for (int row = 0; row < intensitiesTheoMat.rows(); row++) {
            Eigen::VectorX<float> r = intensitiesTheoMat.row(row).array() / std::max(intensitiesTheoMat.row(row).sum(), 0.001f);
            intensitiesTheoMat.row(row) = r;
        }

        Eigen::VectorX<float> klDivByRow(intensitiesTheoMat.rows());
        klDivByRow.setZero();
        for (int row = 0; row < intensitiesTheoMat.rows(); row++) {
            float klDiv;
            e = EigenUtils::klDivergence(
                Eigen::VectorX<float>(matSummed.row(row)),
                Eigen::VectorX<float>(intensitiesTheoMat.row(row)),
                &klDiv);
            klDivByRow.coeffRef(row) = klDiv;
        }

        for (int i = 0; i < bestCorrelationResult.peakCorrelations.size(); i++) {
            candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor1 + i] = bestCorrelationResult.peakCorrelations.at(i);
        }

        const float cosineSimMax = cosineSimsByRow.coeff(bestCorrelationResult.bestAnchorRowIndex);
        candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrum] = cosineSimMax;
        candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumCubed]
                                                    = static_cast<float>(std::pow(cosineSimMax, 3));

        const float klDivMax = klDivByRow.coeff(bestCorrelationResult.bestAnchorRowIndex);
        candidateScores->featuresArray[CandidateScores::Features::KlDivSpectrum] = klDivMax;
        candidateScores->featuresArray[CandidateScores::Features::KlDivSpectrumCubeRoot]
                                                    = static_cast<float>(std::pow(cosineSimMax, 1.0f/3.0f));

        const float cosineSimRowsSummed = cosineSimsByRow.sum();
        if (MathUtils::tZero(cosineSimRowsSummed)) {
            candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumOverTime] = 0.0f;
            candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumOverTimeCubed] = 0.0f;

            ERR_RETURN
        }

        const int nonZeroRows = EigenUtils::nonZeros(cosineSimsByRow);
        const float cosineSimSpectrumOverTime = cosineSimRowsSummed / static_cast<float>(nonZeroRows);

        candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumOverTime] = cosineSimSpectrumOverTime;
        candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumOverTimeCubed]
                                                                = static_cast<float>(std::pow(cosineSimSpectrumOverTime, 3));

        QVector<float> cosineSimsByRowSansZeros;
        cosineSimsByRowSansZeros.reserve(nonZeroRows);
        for (int i = 0; i < cosineSimsByRow.size(); i++) {
            const float f = cosineSimsByRow.coeff(i);
            if (MathUtils::tZero(f)) {
                continue;
            }
            cosineSimsByRowSansZeros.push_back(f);
        }
        candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumStDev]
                                                        = static_cast<float>(MathUtils::stDev(cosineSimsByRowSansZeros));

        const QVector<float> &cosineSimToAnchorVec = bestCorrelationResult.peakCorrelations;

        const int top6 = std::min(6, cosineSimToAnchorVec.size());
        const float cosineSimSumTop = std::accumulate(
            cosineSimToAnchorVec.begin(),
            cosineSimToAnchorVec.begin() + top6,
            0.0f
    );
        candidateScores->featuresArray[CandidateScores::Features::CosineSimSumTop] = cosineSimSumTop;

        float cosineSimSumBottom = 0.1;
        if (cosineSimToAnchorVec.size() > top6) {
            cosineSimSumBottom = std::accumulate(
                    cosineSimToAnchorVec.begin() + top6 + 1,
                    cosineSimToAnchorVec.end(),
                    std::numeric_limits<float>::min()
            );
        }

        candidateScores->featuresArray[CandidateScores::Features::CosineSimSumBottom] = cosineSimSumBottom;

        candidateScores->featuresArray[CandidateScores::Features::TopBottomRatio]
                = std::log(std::max(1.0f, cosineSimSumTop) / (cosineSimSumTop + cosineSimSumBottom + 1.0f));

        candidateScores->featuresArray[CandidateScores::Features::TopBottomRatioNorm]
                = cosineSimSumBottom / static_cast<float>(candidateScores->targetDecoyCandidatePair->totalFragmentCount());

        ERR_RETURN
    }

    constexpr int arraySizeMax = 12;

    Err buildScanTimeVectorFromFrameIndexIntegrationLimits(
        const PeakIntegrationIndexes &peakIntegrationIndexes,
        MsFrame *msFrameMzTarget,
        QVector<ScanTime> *scanTimes
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msFrameMzTarget->isValid()); ree;
        e = ErrorUtils::isTrue(peakIntegrationIndexes.second > peakIntegrationIndexes.first); ree;

        scanTimes->reserve(peakIntegrationIndexes.second - peakIntegrationIndexes.first + 1);
        for (int frameIndex = peakIntegrationIndexes.first; frameIndex <= peakIntegrationIndexes.second; ++frameIndex) {
            scanTimes->push_back(msFrameMzTarget->scanTimeFromFrameIndex(frameIndex));
        }

        ERR_RETURN
    }

    Err calculateTrapezoidalArea(
        const Eigen::MatrixX<float> &xicMat,
        const QVector<float> &scanTimes,
        Eigen::VectorX<float> *areas
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scanTimes); ree;
        e = ErrorUtils::isEqual(scanTimes.size(), static_cast<int>(xicMat.rows())); ree;

        const Eigen::VectorX<float> scanTimesVec = EigenUtils::convertQVectorToEigenVector(scanTimes);

        const int vSize = scanTimes.size() - 1;
        const Eigen::VectorX<float> dScanTimes = scanTimesVec.tail(vSize) - scanTimesVec.head(vSize);
        const Eigen::MatrixX<float> xicMatSums = xicMat.bottomRows(vSize) + xicMat.topRows(vSize);

        *areas = (dScanTimes.transpose() * xicMatSums) * 0.5;

        ERR_RETURN
    }

    Err setFoundMs2Ions(
        const QVector<BestCorrelationResult> &bestCorrelationResults,
        MsFrame *msFrameMzTarget,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        const BestCorrelationResult &bestCorrelationResult = bestCorrelationResults.front();

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedMz.size() > 0); ree;

        const Eigen::MatrixX<float> &matMz = bestCorrelationResult.matBlockTrimmedMz;

        QVector<float> mzMeanValsFound;
        mzMeanValsFound.reserve(static_cast<int>(matMz.cols()));
        QVector<float> stdMeanValsFound;
        stdMeanValsFound.reserve(static_cast<int>(matMz.cols()));

        for (int col = 0; col < matMz.cols(); col++) {

            const Eigen::VectorX<float> &colMz = matMz.col(col);

            if(MathUtils::tZero(colMz.maxCoeff())) {
                mzMeanValsFound.push_back(-1.0f);
                stdMeanValsFound.push_back(-1.0f);
                continue;
            }

            const Eigen::VectorX<float> &colMzNonZero = EigenUtils::removeZeroElements(colMz);

            mzMeanValsFound.push_back(colMzNonZero.mean());
            stdMeanValsFound.push_back(static_cast<float>(EigenUtils::calculateStDevOfVector(colMzNonZero)));
        }

        for (int i = 0; i < std::min(mzMeanValsFound.size(), arraySizeMax); i++) {
            candidateScores->featuresArray[CandidateScores::Features::MzFoundMean1 + i] = mzMeanValsFound.at(i);
        }

        for (int i = 0; i < std::min(stdMeanValsFound.size(), arraySizeMax); i++) {
            candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev1 + i] = stdMeanValsFound.at(i);
        }

        QVector<ScanTime> scanTimes;
        e = buildScanTimeVectorFromFrameIndexIntegrationLimits(
            bestCorrelationResult.peakIntegrationIndexes,
            msFrameMzTarget,
            &scanTimes
            );

        const Eigen::VectorX<float> intensitySums = bestCorrelationResult.matBlockTrimmedIntensity.colwise().sum();

        const Eigen::VectorX<float> intensitySumsNormalized
                = intensitySums / std::max(static_cast<float>(bestCorrelationResult.matBlockTrimmedIntensity.maxCoeff()), 1.0f);

        Eigen::VectorX<float> trapAreas;
        e = calculateTrapezoidalArea(
            bestCorrelationResult.matBlockTrimmedIntensity,
            scanTimes,
            &trapAreas
            ); ree;
        const QVector<float> trapAreasQVec = EigenUtils::convertEigenVectorToQVector(trapAreas);
        candidateScores->trapAreas = QVector<float>(arraySizeMax, -1.0);
        std::copy(trapAreasQVec.begin(), trapAreasQVec.end(), candidateScores->trapAreas.begin());

        for (int i = 0; i < std::min(static_cast<int>(intensitySums.size()), arraySizeMax); i++) {
            candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax1 + i] = intensitySums.coeff(i);
        }

        for (int i = 0; i < std::min(static_cast<int>(intensitySumsNormalized.size()), arraySizeMax); i++) {
            candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMaxNorm1 + i] = intensitySumsNormalized.coeff(i);
        }

        if (bestCorrelationResults.size() > 1) {
            const BestCorrelationResult &bestCorrelationResultAlt = bestCorrelationResults.at(1);
            const Eigen::VectorX<float> intensitySumsAlt = bestCorrelationResultAlt.matBlockTrimmedIntensity.colwise().sum();

            candidateScores->intensityValsAlt.resize(static_cast<int>(intensitySumsAlt.size()));
            for (int i = 0; i < std::min(static_cast<int>(intensitySumsAlt.size()), arraySizeMax); i++) {
                candidateScores->intensityValsAlt[i] = intensitySumsAlt.coeff(i);
            }
        }

        ERR_RETURN
    }

    Err setPeakShapeRatios(
        const BestCorrelationResult &bestCorrelationResult,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.size() > 0);

        constexpr double chunkDivision = 3.0;

        const int bestColumnIndex = bestCorrelationResult.bestAnchorColumnIndex;
        const Eigen::VectorX<float> &bestAnchorColumn
                                        = bestCorrelationResult.matBlockTrimmedIntensity.col(bestColumnIndex);

        QVector<float> bestAnchorColumnVec = EigenUtils::convertEigenVectorToQVector(bestAnchorColumn);

        const auto terminatorLogic = [](double d){return d < 1.0;};
        const auto terminator = std::remove_if(bestAnchorColumnVec.begin(), bestAnchorColumnVec.end(), terminatorLogic);
        bestAnchorColumnVec.erase(terminator, bestAnchorColumnVec.end());

        const int chunkSize = std::max(1, static_cast<int>(std::round(bestAnchorColumnVec.size() / chunkDivision)));
        const double bestAnchorColumnVecSum = std::accumulate(bestAnchorColumnVec.begin(), bestAnchorColumnVec.end(), 0.0001);

        if (bestAnchorColumnVec.size() < chunkDivision) {
            candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio1] = std::numeric_limits<float>::min();
            candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio2] = 1.0;
            candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio3] = std::numeric_limits<float>::min();
        }
        else {

            candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio1] = std::accumulate(
                    bestAnchorColumnVec.begin(),
                    bestAnchorColumnVec.begin() + chunkSize,
                    std::numeric_limits<float>::min()
            ) / bestAnchorColumnVecSum;

            candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio2] = std::accumulate(
                    bestAnchorColumnVec.begin() + chunkSize,
                    bestAnchorColumnVec.begin() + (chunkSize * 2),
                    std::numeric_limits<float>::min()
            ) / bestAnchorColumnVecSum;

            candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio3] = std::accumulate(
                    bestAnchorColumnVec.begin() + (chunkSize * 2),
                    bestAnchorColumnVec.end(),
                    std::numeric_limits<float>::min()
            ) / bestAnchorColumnVecSum;

        }

        ERR_RETURN
    }

    Err setAminoAcidFrequencies(CandidateScores *candidateScores) {

        ERR_INIT

        QMap<QChar, int> aminoAcidCounts = {
            {'A', 0},
            {'C', 0},
            {'D', 0},
            {'E', 0},
            {'F', 0},
            {'G', 0},
            {'H', 0},
            {'I', 0},
            {'K', 0},
            {'L', 0},
            {'M', 0},
            {'N', 0},
            {'P', 0},
            {'Q', 0},
            {'R', 0},
            {'S', 0},
            {'T', 0},
            {'V', 0},
            {'W', 0},
            {'Y', 0},
            {'B', 0},
            {'J', 0},
            {'O', 0},
            {'U', 0},
            {'X', 0},
            {'Z', 0}
        };

        const QString peptideString = candidateScores->targetDecoyCandidatePair->peptideString();

        for (const QChar aminoAcid : peptideString) {

            if (!aminoAcidCounts.contains(aminoAcid)) {
                qDebug() << peptideString << "missing amino acid" << aminoAcid;
                continue;
            }

            e = ErrorUtils::isTrue(aminoAcidCounts.contains(aminoAcid)); ree;
            aminoAcidCounts[aminoAcid]++;
        }

        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountA] = static_cast<float>(aminoAcidCounts['A']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountC] = static_cast<float>(aminoAcidCounts['C']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountD] = static_cast<float>(aminoAcidCounts['D']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountE] = static_cast<float>(aminoAcidCounts['E']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountF] = static_cast<float>(aminoAcidCounts['F']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountG] = static_cast<float>(aminoAcidCounts['G']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountH] = static_cast<float>(aminoAcidCounts['H']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountI] = static_cast<float>(aminoAcidCounts['I']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountK] = static_cast<float>(aminoAcidCounts['K']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountL] = static_cast<float>(aminoAcidCounts['L']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountM] = static_cast<float>(aminoAcidCounts['M']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountN] = static_cast<float>(aminoAcidCounts['N']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountP] = static_cast<float>(aminoAcidCounts['P']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountQ] = static_cast<float>(aminoAcidCounts['Q']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountR] = static_cast<float>(aminoAcidCounts['R']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountS] = static_cast<float>(aminoAcidCounts['S']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountT] = static_cast<float>(aminoAcidCounts['T']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountV] = static_cast<float>(aminoAcidCounts['V']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountW] = static_cast<float>(aminoAcidCounts['W']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountY] = static_cast<float>(aminoAcidCounts['Y']);

        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountB] = static_cast<float>(aminoAcidCounts['B']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountJ] = static_cast<float>(aminoAcidCounts['J']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountO] = static_cast<float>(aminoAcidCounts['O']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountU] = static_cast<float>(aminoAcidCounts['U']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountX] = static_cast<float>(aminoAcidCounts['X']);
        candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountZ] = static_cast<float>(aminoAcidCounts['Z']);

        ERR_RETURN
    }

    Err setShadowCorrelations(
        const BestCorrelationResult &bestCorrelationResult,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.sum() > 0); ree;

        for (int col = 0; col < bestCorrelationResult.matBlockTrimmedIntensity.cols(); col++) {

            const Eigen::VectorX<float> &v1 = bestCorrelationResult.matBlockTrimmedIntensity.col(col);
            const Eigen::VectorX<float> &v2 = bestCorrelationResult.matBlockTrimmedIntensityShadows.col(col);
            const float v1Max = v1.maxCoeff();
            const float v2Max = v2.maxCoeff();

            float cosineSim;
            e = EigenUtils::cosineSimilarity(v1, v2, &cosineSim); ree;

            candidateScores->featuresArray[CandidateScores::Features::ShadowsCosineSimSum] += cosineSim;
            candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor1 + col] = std::max(cosineSim, 0.0f);

        }

        ERR_RETURN
    }

    Err setMzPeakLengthRelatedScores(
        const BestCorrelationResult &bestCorrelationResult,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.sum() > 0); ree;

        const QVector<int> mzPeakLengthsVec = bestCorrelationResult.columnNonZeroPeakLengths();

        const int mzPeakLengthsSum = std::accumulate(
            mzPeakLengthsVec.begin(),
            mzPeakLengthsVec.end(),
            0
            );

        QVector<float> mzPeakLengthsNormalized;
        if (mzPeakLengthsSum != 0) {
            std::transform(
                    mzPeakLengthsVec.begin(),
                    mzPeakLengthsVec.end(),
                    std::back_inserter(mzPeakLengthsNormalized),
                    [mzPeakLengthsSum](int i){return i / static_cast<double>(mzPeakLengthsSum);}
            );
        }
        for (int i = 0; i < std::min(mzPeakLengthsNormalized.size(), arraySizeMax); i++) {
            candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm1 + i] = mzPeakLengthsNormalized.at(i);
        }

        ERR_RETURN
    }

    Err setMs1Averagine(
        const QVector<float> &ms1Averagine,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        const Eigen::VectorX<float> averagineVec = EigenUtils::convertQVectorToEigenVector(ms1Averagine);

        const QVector<float> ms1IsoDisActualVec = {
            candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1PreMono],
            candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1],
            candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso1],
            candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso2]
            };
        const Eigen::VectorX<float> ms1IsoDistActual = EigenUtils::convertQVectorToEigenVector(ms1IsoDisActualVec);

        float cosineSimAveragine;
        e = EigenUtils::cosineSimilarity(
            averagineVec,
            ms1IsoDistActual,
            &cosineSimAveragine
            ); ree;

        candidateScores->featuresArray[CandidateScores::Features::MS1Averagine] = cosineSimAveragine;

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::setCandidateScores(
    const TargetDecoyCandidatePair *targetDecoyCandidatePair,
    const QVector<BestCorrelationResult> &bestCorrelationResults,
    const QVector<float> &ms1Averagine,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    const BestCorrelationResult &bestCorrelationResult = bestCorrelationResults.front();

    e = ErrorUtils::isNotEmpty(bestCorrelationResult.peakCorrelations); ree;
    e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.size() > 0); ree;

    candidateScores->initFeaturesArray();

    candidateScores->frameIndex = bestCorrelationResult.peakIntegrationIndexes.first
                                + bestCorrelationResult.apexStarts.at(bestCorrelationResult.bestAnchorColumnIndex);

    candidateScores->frameIndexStart = bestCorrelationResult.peakIntegrationIndexes.first;
    candidateScores->frameIndexEnd = bestCorrelationResult.peakIntegrationIndexes.second;

    candidateScores->scanNumber = m_msFrameMzTarget->scanNumberFromFrameIndex(candidateScores->frameIndex);
    candidateScores->scanNumberStart = m_msFrameMzTarget->scanNumberFromFrameIndex(candidateScores->frameIndexStart);
    candidateScores->scanNumberEnd = m_msFrameMzTarget->scanNumberFromFrameIndex(candidateScores->frameIndexEnd);

    if(bestCorrelationResults.size() > 1) {

        const BestCorrelationResult &bestCorrelationResultAlt = bestCorrelationResults.at(1);

        candidateScores->scanTimeStartAlt
            = m_msFrameMzTarget->scanTimeFromFrameIndex(bestCorrelationResultAlt.peakIntegrationIndexes.first);

        candidateScores->scanTimeEndAlt
            = m_msFrameMzTarget->scanTimeFromFrameIndex(bestCorrelationResultAlt.peakIntegrationIndexes.second);

        candidateScores->scanTimeAlt
                = m_msFrameMzTarget->scanTimeFromFrameIndex(
                    bestCorrelationResultAlt.peakIntegrationIndexes.first + bestCorrelationResult.bestAnchorRowIndex
                    );

        candidateScores->cosineSimSum100Alt = bestCorrelationResultAlt.peakCorrelationsSum;
    }

    candidateScores->scanTime = m_msFrameMzTarget->scanTimeFromScanNumber(candidateScores->scanNumber);
    candidateScores->scanTimeStart = m_msFrameMzTarget->scanTimeFromScanNumber(candidateScores->scanNumberStart);
    candidateScores->scanTimeEnd = m_msFrameMzTarget->scanTimeFromScanNumber(candidateScores->scanNumberEnd);

    const int top6 = std::min(6, bestCorrelationResult.peakCorrelations.size());
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100] = std::accumulate(
            bestCorrelationResult.peakCorrelations.begin(),
            bestCorrelationResult.peakCorrelations.begin() + top6,
            std::numeric_limits<float>::min()
            );

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100Top12] = std::accumulate(
            bestCorrelationResult.peakCorrelations.begin(),
            bestCorrelationResult.peakCorrelations.end(),
            std::numeric_limits<float>::min()
            );

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100Window1p5X] = std::accumulate(
        bestCorrelationResult.peakCorrelationsWindow1p5X.begin(),
        bestCorrelationResult.peakCorrelationsWindow1p5X.begin() + top6,
        std::numeric_limits<float>::min()
        );

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100Window2X] = std::accumulate(
        bestCorrelationResult.peakCorrelationsWindow2X.begin(),
        bestCorrelationResult.peakCorrelationsWindow2X.begin() + top6,
        std::numeric_limits<float>::min()
        );

    candidateScores->featuresArray[CandidateScores::Features::ScanIonCount]
        = static_cast<float>(m_msFrameMzTarget->getScanPointsByScanNumber(candidateScores->scanNumber)->size());

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100GreaterThan80]
                                            = calculatedCosineSimSumGreaterThan80(bestCorrelationResult.peakCorrelations);

    candidateScores->featuresArray[CandidateScores::Features::TheoFragmentCount]
                                                            = static_cast<float>(targetDecoyCandidatePair->totalFragmentCount());

    candidateScores->featuresArray[CandidateScores::Features::TotalIntensityLog]
                                        = std::log(std::max(bestCorrelationResult.matBlockTrimmedIntensity.sum(),
                                                    std::numeric_limits<float>::min()));

     candidateScores->featuresArray[CandidateScores::Features::TotalIntensityPeakHeights]
                                = bestCorrelationResult.matBlockTrimmedIntensity.colwise().maxCoeff().sum();

    candidateScores->featuresArray[CandidateScores::Features::TotalIntensityRaw]
                                            = std::exp(candidateScores->featuresArray[CandidateScores::Features::TotalIntensityLog]);

    candidateScores->featuresArray[CandidateScores::Features::Charge]
                                                        = static_cast<float>(candidateScores->targetDecoyCandidatePair->charge());

    const float scanTimeDelta = candidateScores->scanTime - candidateScores->scanTimePredicted;
    candidateScores->featuresArray[CandidateScores::Features::ScanTimeDelta] = scanTimeDelta;
    candidateScores->featuresArray[CandidateScores::Features::ScanTimeDeltaAbs] = std::abs(scanTimeDelta);

    candidateScores->featuresArray[CandidateScores::Features::ScanTimePredicted] = candidateScores->scanTimePredicted;

    const double pdScanTime = std::sqrt(std::min(std::abs(scanTimeDelta), m_scanTimeRange) / m_scanTimeRange);
    candidateScores->featuresArray[CandidateScores::Features::ScanTimePdAbs] = static_cast<float>(pdScanTime);
    candidateScores->featuresArray[CandidateScores::Features::ScanTimePd] = scanTimeDelta < 0
                                                                             ? -static_cast<float>(std::abs(pdScanTime))
                                                                             : static_cast<float>(std::abs(pdScanTime));

    const double pepLength = (-10.0 + candidateScores->targetDecoyCandidatePair->peptideString().size()) / 10.0;
    candidateScores->featuresArray[CandidateScores::Features::PeptideLengthNorm] = static_cast<float>(pepLength);

    const auto mz = candidateScores->targetDecoyCandidatePair->mz();
    candidateScores->featuresArray[CandidateScores::Features::MzNorm] = (mz - 600.0f) * 0.002f;
    candidateScores->featuresArray[CandidateScores::Features::IRTPredicted] = candidateScores->targetDecoyCandidatePair->iRt();
    candidateScores->featuresArray[CandidateScores::Features::Mass] = candidateScores->targetDecoyCandidatePair->mass();

    const float mzTargetKey = MathUtils::unHashDecimal<float>(m_mzTargetKey.toInt(), S_GLOBAL_SETTINGS.HASHING_PRECISION);
    candidateScores->featuresArray[CandidateScores::Features::TargetWindowLocation] = mzTargetKey - mz;
    candidateScores->featuresArray[CandidateScores::Features::TargetWindowLocationAbs] = std::abs(mzTargetKey - mz);

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum45]
        = std::max(std::accumulate(bestCorrelationResult.peakCorrelations45.begin(), bestCorrelationResult.peakCorrelations45.begin() + top6, 0.0f), std::numeric_limits<float>::min());

    int bestAlignmentMatrixRowIndex = bestCorrelationResult.bestAnchorRowIndex;
    candidateScores->featuresArray[CandidateScores::Features::AllignedMaxIndexesCount] = static_cast<float>(std::count_if(
        bestCorrelationResult.apexStarts.begin(),
        bestCorrelationResult.apexStarts.end(),
        [bestAlignmentMatrixRowIndex](int i){return i == bestAlignmentMatrixRowIndex;}
        ));

    candidateScores->featuresArray[CandidateScores::Features::AlignmentIndexMean] = MathUtils::mean(bestCorrelationResult.apexStarts);
    candidateScores->featuresArray[CandidateScores::Features::AlignmentIndexStDev] = MathUtils::stDev(bestCorrelationResult.apexStarts);
    candidateScores->featuresArray[CandidateScores::Features::AlignmentCombinedScore]
                        = candidateScores->featuresArray[CandidateScores::Features::AlignmentIndexMean]
                        * candidateScores->featuresArray[CandidateScores::Features::AllignedMaxIndexesCount];

    candidateScores->featuresArray[CandidateScores::Features::MatrixZeroPercentage] = static_cast<float>((bestCorrelationResult.matBlockTrimmedIntensity.array() < 1).count())
                                                                                    / bestCorrelationResult.matBlockTrimmedIntensity.size();
    e = setAminoAcidFrequencies(candidateScores); ree;

    e = setCosineSimilarityMetrics(
        targetDecoyCandidatePair,
        bestCorrelationResult,
        candidateScores
        );ree

    e = setMs1RelatedScores(
        targetDecoyCandidatePair,
        bestCorrelationResult,
        static_cast<float>(m_pythiaParameters.ms1ExtractionWidthPPM),
        candidateScores
        ); ree;

    const QVector<MS2Ion> ms2IonsTheoritical = candidateScores->isDecoy
                                             ? targetDecoyCandidatePair->ms2IonsDecoy()
                                             : targetDecoyCandidatePair->ms2IonsTarget();

    e = setFoundMs2Ions(
        bestCorrelationResults,
        m_msFrameMzTarget,
        candidateScores
        ); ree;

    e = setPeakShapeRatios(bestCorrelationResult, candidateScores); ree;

    e = setShadowCorrelations(bestCorrelationResult, candidateScores); ree;

    e = setMzPeakLengthRelatedScores(bestCorrelationResult, candidateScores); ree;

    e = setMs1Averagine(ms1Averagine, candidateScores); ree;

    ERR_RETURN
}

namespace {

    Eigen::VectorX<float> buildXicPointsVector(
        const XICPoints &xicPoints,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax
        ) {

        Eigen::VectorX<float> vec(frameIndexMax + 1);
        vec.setZero();

        for (const auto &[mz, intensity, scanNumber] : xicPoints) {

            if (!(frameIndexMin <= scanNumber && scanNumber <= frameIndexMax)) {
                continue;
            }
            vec.coeffRef(scanNumber) = intensity;
        }

        return vec;
    }

    Err calculateMs1Scores(
        const Eigen::VectorX<float> &kernel,
        const Eigen::VectorX<float> &_anchorColumn,
        float mzToExtract,
        float massTol,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        TurboXIC *turboXicMS1,
        float *cosineSimMS1,
        float *mzMS1Mean,
        float *mzMS1StDev,
        float *ppmDiff,
        float *ms1IntensitySum
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(turboXicMS1->isInit()); ree;
        e = ErrorUtils::isTrue(_anchorColumn.size() > 0); ree;

        Eigen::VectorX<float> anchorColumn = _anchorColumn;

        *cosineSimMS1 = 0.01;
        *mzMS1Mean = 0.01;
        *mzMS1StDev = 10.0;
        *ppmDiff = 100.0;
        *ms1IntensitySum = 0.0;

        XICPoints xicPoints = turboXicMS1->extractPointsXIC(
            mzToExtract - massTol,
            mzToExtract + massTol
            );

        TurboXIC::filterXICPointsByScanNumber(
            frameIndexMin,
            frameIndexMax,
            &xicPoints
            );

        if (xicPoints.empty()) {
            ERR_RETURN
        }

        Eigen::VectorX<float> xicVec = buildXicPointsVector(
            xicPoints,
            frameIndexMin,
            frameIndexMax
            ).segment(frameIndexMin, frameIndexMax - frameIndexMin + 1).eval();

        if (xicVec.size() < anchorColumn.size()) {
            Eigen::VectorX<float> xicVecResized(anchorColumn.size());
            xicVecResized.setZero();
            const int stepSize = static_cast<int>(std::round(anchorColumn.size() / xicVec.size()));

            int ogVecIndex = 0;
            for(int i = 0; i >= anchorColumn.size(); i += stepSize) {

                if (ogVecIndex < xicVec.size()) {
                    break;
                }

                xicVecResized.coeffRef(i) = xicVec.coeff(ogVecIndex++);
            }

            xicVec = xicVecResized;
            if (stepSize >= 1) {
                xicVec = EigenKernelUtils::convolveVectorWithKernel(xicVec, kernel);
                anchorColumn = EigenKernelUtils::convolveVectorWithKernel(anchorColumn, kernel);
            }
        }
        else if (xicVec.size() > anchorColumn.size()) {

            Eigen::VectorX<float> anchorColumnResized(xicVec.size());
            anchorColumnResized.setZero();
            const int stepSize = static_cast<int>(std::round(xicVec.size() / anchorColumn.size()));

            int ogVecIndex = 0;
            for(int i = 0; i < xicVec.size(); i += stepSize) {

                if (ogVecIndex >= anchorColumn.size()) {
                    break;
                }

                anchorColumnResized.coeffRef(i) = anchorColumn.coeff(ogVecIndex++);
            }

            anchorColumn = anchorColumnResized;
            if (stepSize > 1) {
                xicVec = EigenKernelUtils::convolveVectorWithKernel(xicVec, kernel);
                anchorColumn = EigenKernelUtils::convolveVectorWithKernel(anchorColumn, kernel);
            }
        }

        e = EigenUtils::cosineSimilarity(
            anchorColumn,
            xicVec,
            cosineSimMS1
            ); ree;

        const int nonZeros = EigenUtils::nonZeros(xicVec);
        if (nonZeros < 1) {
            ERR_RETURN
        }

        const auto meanLogic = [](float sum, const XICPoint &p){return sum + p.mz;};
        *mzMS1Mean = std::accumulate(xicPoints.begin(), xicPoints.end(), 0.0f, meanLogic) / static_cast<float>(nonZeros);

        const auto stDevLogic = [mzMS1Mean](float sum, const XICPoint &p){return sum + std::pow((p.mz - *mzMS1Mean) , 2);};
        *mzMS1StDev = std::sqrt(std::accumulate(xicPoints.begin(), xicPoints.end(), 0.0f, stDevLogic)
                    / static_cast<float>(nonZeros));

        *ppmDiff = MathUtils::calculateMassAccuracyPPM(mzToExtract, *mzMS1Mean);
        *ms1IntensitySum = xicVec.sum();

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::setMs1RelatedScores(
    const TargetDecoyCandidatePair *targetDecoyCandidatePair,
    const BestCorrelationResult &bestCorrelationResult,
    float ppmTol,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(m_turboXicMS1->isInit()); ree;
    const FrameIndex frameIndexMin = bestCorrelationResult.peakIntegrationIndexes.first;
    const FrameIndex frameIndexMax = bestCorrelationResult.peakIntegrationIndexes.second;

    FrameIndex frameIndexMinMS1;
    e = m_msFrameMS1->frameIndexFromScanTime(candidateScores->scanTimeStart, &frameIndexMinMS1); ree;

    FrameIndex frameIndexMaxMS1;
    e = m_msFrameMS1->frameIndexFromScanTime(candidateScores->scanTimeEnd, &frameIndexMaxMS1); ree;

    // qDebug() << frameIndexMin << frameIndexMinMS1 << frameIndexMax << frameIndexMaxMS1 << "SDLFDKJsl";

    const float isotopeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / targetDecoyCandidatePair->charge();

    const float monoIsotopeMz = targetDecoyCandidatePair->mz();
    const float monoIsotopeShadowMz = monoIsotopeMz - isotopeDistance;
    const float c13isotopeMz1 = monoIsotopeMz + isotopeDistance;
    const float c13isotopeMz2 = monoIsotopeMz + (isotopeDistance * 2);

    const float massTol = MathUtils::calculatePPM(monoIsotopeMz, ppmTol);

    const Eigen::VectorX<float> anchorColumn
            = bestCorrelationResult.matBlockTrimmedIntensity.col(bestCorrelationResult.bestAnchorColumnIndex);

    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        monoIsotopeMz,
        massTol,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound100],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound100],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound100PPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound100]
        ); ree;

    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        monoIsotopeMz,
        massTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim45MS1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound45],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound45],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound45PPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound45]
        ); ree;

    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        monoIsotopeShadowMz,
        massTol,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1PreMono],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMono],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundPreMono],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMonoPPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundPreMono]
        ); ree;

    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        c13isotopeMz1,
        massTol,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1PPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundIso1]
        ); ree;
    
    e = calculateMs1Scores(
        d_ptr->m_kernelMs2,
        anchorColumn,
        c13isotopeMz2,
        massTol,
        frameIndexMinMS1,
        frameIndexMaxMS1,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso2],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso2],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso2],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso2PPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundIso2]
        ); ree;

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100MS1]
        = candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1]
        + candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso1]
        + candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso2]
        - candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1PreMono];

    ERR_RETURN
}

