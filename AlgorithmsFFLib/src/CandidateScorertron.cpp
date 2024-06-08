//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron.h"

#include "CandidateScores.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "IsotopicDistributionBuilder.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"
#include "XICPeakManager.h"


class Q_DECL_HIDDEN CandidateScorertron::Private {
public:

    using NominalMass = int;
    using AveragineVec = QVector<float>;

    Private();
    ~Private();

    Err init();

    Eigen::VectorX<float> m_kernelIntegration;
    Eigen::VectorX<float> m_kernelMs2;

};

CandidateScorertron::Private::Private() {}
CandidateScorertron::Private::~Private() {}

Err CandidateScorertron::Private::init() {

    ERR_INIT

    //TODO consider writing algo to optimizing filter length parameters as well as smoothCount later in the code.
    constexpr int filterLengthMs2 = 3;
    constexpr int filterLengthIntegration = 5;
    constexpr int order = 1;
    constexpr int derivative = 0;
    constexpr int rate = 1;

    Eigen::MatrixX<float> kernel;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        filterLengthMs2,
        order,
        derivative,
        rate,
        &kernel
        ); ree;
    const Eigen::VectorX<float> kernelVec(kernel);
    m_kernelMs2 = kernelVec;

    // m_kernelMs2 = EigenKernelUtils::buildGaussianFilter1D<float>(filterLengthMs2, 0.75);

    Eigen::MatrixX<float> kernelIntegration;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
        filterLengthIntegration,
        order,
        derivative,
        rate,
        &kernelIntegration
        ); ree;
    const Eigen::VectorX<float> kernelIntegrationVec(kernelIntegration);
    m_kernelIntegration = kernelIntegrationVec;

    // m_kernelIntegration = EigenKernelUtils::buildGaussianFilter1D<float>(filterLengthIntegration, 0.75);

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

CandidateScorertron::CandidateScorertron()
: m_topNMS2Ions(-1)
, m_xicPeakManager(nullptr)
, m_msFrameMzTarget(nullptr)
, m_turboXicMS1(nullptr)
, d_ptr(QScopedPointer<Private>(new Private))
, m_minPeakCount(3.9)
, m_scanTimeRange(0)
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
    XICPeakManager *xicPeakManager,
    MsFrame *msFrameMzTarget,
    TurboXIC *turboXicMS1
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;
    e = ErrorUtils::isTrue(msFrameMzTarget->isValid()); ree;
    e = ErrorUtils::isTrue(turboXicMS1->isInit()); ree;
    e = ErrorUtils::isNotEmpty(mzTargetKey); ree;
    e = ErrorUtils::isNotEmpty(averagineTable); ree;
    e = ErrorUtils::isFalse(MathUtils::tZero(scanTimeRange)); ree;
    e = ErrorUtils::isAboveThreshold(
        topNMS2Ions,
        S_GLOBAL_SETTINGS.MIN_MS2_IONS,
        ErrorUtilsParam::IncludeThreshold); ree;

    e = ErrorUtils::isAboveThreshold(
        minPeakCount,
        1.0f,
        ErrorUtilsParam::ExcludeThreshold); ree;

    m_pythiaParameters = pythiaParameters;
    m_topNMS2Ions = topNMS2Ions;
    m_mzTargetKey = mzTargetKey;
    m_xicPeakManager = xicPeakManager;
    m_msFrameMzTarget = msFrameMzTarget;
    m_turboXicMS1 = turboXicMS1;
    m_scanTimeRange = scanTimeRange;
    m_averagineTable = averagineTable;

    if (msCalibratomatic.isInitRT()) {
        m_msCalibratomatic = msCalibratomatic;
    }

    e = d_ptr->init(); ree;

    ERR_RETURN
}

class MatriciesAndVecs {

public:

    Eigen::MatrixX<float> intensityMatrix100;
    Eigen::MatrixX<float> intensityMatrix100Shadow;
    Eigen::MatrixX<float> intensityMatrix45;
    Eigen::MatrixX<float> intensityMatrix20;

    Eigen::MatrixX<float> mzMatrix100;

    Eigen::VectorX<float> integrationVec;
    Eigen::VectorX<float> integrationVecCosineSim;
    Eigen::VectorX<float> productVec;

    [[nodiscard]] bool intensityMatriciesAreValid() const {
        return intensityMatrix100.size() > 0 && mzMatrix100.size() > 0;
    }

    [[nodiscard]] bool integrationVecIsValid() const {
        return integrationVec.size() > 0;
    }
};

class BestCorrelationResult {

public:
    QVector<float> peakCorrelations;
    QVector<float> peakCorrelations45;
    QVector<float> peakCorrelations20;
    QVector<float> peakCorrelationsWindow1p5X;
    QVector<float> peakCorrelationsWindow2X;
    float peakCorrelationsSum = -1.0;
    Eigen::MatrixX<float> matBlockTrimmedIntensity;
    Eigen::MatrixX<float> matBlockTrimmedIntensityWindow1p5X;
    Eigen::MatrixX<float> matBlockTrimmedIntensityWindow2X;
    Eigen::MatrixX<float> matBlockTrimmedIntensity45;
    Eigen::MatrixX<float> matBlockTrimmedIntensity20;
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
        QVector<XICPoints> *xicPointsVec45,
        QVector<XICPoints> *xicPointsVec20
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;

        xicPointsVec100->reserve(ms2Ions.size());
        xicPointsVec100Shadows->reserve(ms2Ions.size());
        xicPointsVec45->reserve(ms2Ions.size());
        xicPointsVec20->reserve(ms2Ions.size());

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
                xicPointsVec20->push_back({});
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

            filterXICPointsByAccuracyPPM(
                ms2Ion.mz,
                ppmTol * S_GLOBAL_SETTINGS.TIGHT_2_FRACTION,
                &xicPoints
                );
            xicPointsVec20->push_back(xicPoints);
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
        Eigen::VectorX<float> *integrationVec
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;

        constexpr float intensityThresholdVal = 0.1;
        constexpr float countValue = 1.0;

        Eigen::MatrixX<float> matCount = matriciesAndVecs.intensityMatrix100;
        matCount = (matCount.array() > intensityThresholdVal).select(countValue, matCount);
        EigenUtils::thresholdMatrix(0.0f, &matCount);

        Eigen::VectorX<float> integrationVecLocal = matCount.rowwise().sum();
        EigenUtils::thresholdVector(minPeakCount, &integrationVecLocal);

        *integrationVec = EigenKernelUtils::convolveVectorWithKernel(
        integrationVecLocal,
        kernelIntegration
        );

        ERR_RETURN
    }

    Err buildIntegrationVectorCosineSim(
        const MatriciesAndVecs &matriciesAndVecs,
        const QVector<MS2Ion> &ms2IonsTheo,
        const Eigen::VectorX<float> &kernelIntegration,
        Eigen::VectorX<float> *integrationVecCosineSim
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;
        e = ErrorUtils::isNotEmpty(ms2IonsTheo); ree;
        e = ErrorUtils::isEqual(ms2IonsTheo.size(), static_cast<int>(matriciesAndVecs.intensityMatrix100.cols())); ree;

        Eigen::MatrixX<float> matIntensityVals = matriciesAndVecs.intensityMatrix100;

        const int matRows = static_cast<int>(matIntensityVals.rows());
        Eigen::MatrixX<Intensity> matMs2IonsIntensityVals(matRows, matIntensityVals.cols());


        for (int i = 0; i < ms2IonsTheo.size(); i++) {
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

    Err initMatricesdAndVecs(
        const QVector<MS2Ion> &ms2Ions,
        const Eigen::VectorX<float> &kernelMs2,
        const Eigen::VectorX<float> &kernelIntegration,
        FrameIndex frameIndexPredictedMin,
        FrameIndex frameIndexPredictedMax,
        int topNMS2Ions,
        float ppmTol,
        float minPeakCount,
        bool subtractShadows,
        XICPeakManager *xicPeakManager,
        MatriciesAndVecs *matriciesAndVecs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;
        e = ErrorUtils::isAboveThreshold(topNMS2Ions, 0, ErrorUtilsParam::ExcludeThreshold); ree;

        QVector<MS2Ion> ms2IonsResized = ms2Ions;
        ms2IonsResized.resize(std::min(topNMS2Ions, ms2Ions.size()));

        const bool ms2IonsAreSorted = std::is_sorted(
            ms2IonsResized.rbegin(),
            ms2IonsResized.rend(),
            [](const MS2Ion &l, const MS2Ion &r){return l.intensity < r.intensity;}
            );
        e = ErrorUtils::isTrue(ms2IonsAreSorted); ree;

        QVector<XICPoints> xicPointsVec100;
        QVector<XICPoints> xicPointsVec100Shadow;
        QVector<XICPoints> xicPointsVec45;
        QVector<XICPoints> xicPointsVec20;
        e = getXICs(
            ms2IonsResized,
            ppmTol,
            frameIndexPredictedMin,
            frameIndexPredictedMax,
            xicPeakManager,
            &xicPointsVec100,
            &xicPointsVec100Shadow,
            &xicPointsVec45,
            &xicPointsVec20
            ); ree;

        e = ErrorUtils::isEqual(xicPointsVec100.size(), xicPointsVec45.size()); ree;
        e = ErrorUtils::isEqual(xicPointsVec100.size(), xicPointsVec20.size()); ree;
        e = ErrorUtils::isEqual(xicPointsVec100.size(), xicPointsVec100Shadow.size()); ree;

        const FrameIndex frameIndexMax = findFrameIndexMaxXICPointsVec(xicPointsVec100);

        constexpr int smoothCount = 2;
        e = buildEigenMatrix(
            xicPointsVec100,
            kernelMs2,
            frameIndexMax,
            true,
            smoothCount,
            &matriciesAndVecs->intensityMatrix100,
            &matriciesAndVecs->mzMatrix100
            ); ree;

        Eigen::MatrixX<float> unused;
        e = buildEigenMatrix(
            xicPointsVec100Shadow,
            kernelMs2,
            frameIndexMax,
            false,
            smoothCount,
            &matriciesAndVecs->intensityMatrix100Shadow,
            &unused
            ); ree;

        matriciesAndVecs->intensityMatrix100 = subtractShadows
                          ? matriciesAndVecs->intensityMatrix100 - matriciesAndVecs->intensityMatrix100Shadow
                          : matriciesAndVecs->intensityMatrix100;
        EigenUtils::thresholdMatrix(0.0f, &matriciesAndVecs->intensityMatrix100);

        e = buildIntegrationVector(
            *matriciesAndVecs,
            kernelIntegration,
            minPeakCount,
            &matriciesAndVecs->integrationVec
            ); ree;

        e = buildIntegrationVectorCosineSim(
            *matriciesAndVecs,
            ms2IonsResized,
            kernelIntegration,
            &matriciesAndVecs->integrationVecCosineSim
            ); ree;

        matriciesAndVecs->productVec = matriciesAndVecs->integrationVec.array() * matriciesAndVecs->integrationVecCosineSim.array();

        // for (int i = 0; i < matriciesAndVecs->integrationVecCosineSim.size(); i++) {
        //     std::cout << matriciesAndVecs->integrationVecCosineSim.coeff(i) << ", ";
        // }
        // std::cout << std::endl;
        //
        // for (int i = 0; i < matriciesAndVecs->integrationVec.size(); i++) {
        //     std::cout << matriciesAndVecs->integrationVec.coeff(i) << ", ";
        // }
        // for (int i = 0; i < matriciesAndVecs->productVec.size(); i++) {
        //     std::cout << matriciesAndVecs->productVec.coeff(i) << ", ";
        // }
        // std::cout << std::endl;
        // std::cout << "************" << std::endl;

        constexpr int noSmooths = 0;
        e = buildEigenMatrix(
            xicPointsVec45,
            kernelMs2,
            frameIndexMax,
            false,
            noSmooths,
            &matriciesAndVecs->intensityMatrix45,
            &unused
            ); ree;

        e = buildEigenMatrix(
            xicPointsVec20,
            kernelMs2,
            frameIndexMax,
            false,
            noSmooths,
            &matriciesAndVecs->intensityMatrix20,
            &unused
            ); ree;

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::calculateScores(
    const QVector<MS2Ion> &ms2Ions,
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
        d_ptr->m_kernelMs2,
        d_ptr->m_kernelIntegration,
        frameIndexPredictedMin,
        frameIndexPredictedMax,
        m_topNMS2Ions,
        static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM),
        m_minPeakCount,
        m_pythiaParameters.subtractShadows,
        m_xicPeakManager,
        &matriciesAndVecs
        ); ree;

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
    e = EigenUtils::simpleIntegrator(
        matriciesAndVecs.productVec,
        m_pythiaParameters.stopThresholdFraction,
        &peakIntegrationsVsIntensities
        ); ree;

    if (peakIntegrationsVsIntensities.isEmpty()) {
        ERR_RETURN
    }

    BestCorrelationResult bestCorrelationResult;
    e = processIntegrationVectorPeakIntegrations(
        matriciesAndVecs,
        peakIntegrationsVsIntensities,
        &bestCorrelationResult
        ); ree;


    const int nominalMass = static_cast<int>((std::round(targetDecoyCandidatePair->mass() / 10) * 10));
    e = ErrorUtils::isTrue(m_averagineTable.contains(nominalMass)); ;
    const QVector<float> ms1Averagine = m_averagineTable.value(nominalMass);

    e = setCandidateScores(
        targetDecoyCandidatePair,
        bestCorrelationResult,
        ms1Averagine,
        candidateScores
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
            = m_msCalibratomatic.scanTimeStDev(m_pythiaParameters.scanTimeWindowStDevs);

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
        const Eigen::VectorX<float> &integrationVec
        ) {
        QPair<PeakIntegrationIndexes, Intensity> piiWorking = pii;
        if (piiWorking.first.first == piiWorking.first.second) {
            piiWorking.first.first = std::max(0, piiWorking.first.first - 1);
            piiWorking.first.second = std::min(
                static_cast<int>(integrationVec.size()) - 1,
                piiWorking.first.second + 1
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

            QVector<int> apexesStartsToUse;
            apexesStartsToUse.reserve(static_cast<int>(matBlockApexes.cols()));
            for (int col = 0; col < matBlockApexes.cols(); col++) {

                const Eigen::VectorX<float> &apexColumn = matBlockApexes.col(col);

                if (apexColumn.coeff(apexIndex) > 0) {
                    apexesStartsToUse.push_back(apexIndex);
                    continue;
                }

                for (int rowFromCenter = 1; rowFromCenter < apexColumn.size(); rowFromCenter++) {

                    const int rowLeftIndex = apexIndex - rowFromCenter;
                    const int rowRightIndex = apexIndex + rowFromCenter;

                    int rowLeftIndexValue = -1;
                    int rowRightIndexValue = -1;

                    if (rowLeftIndex >= 0) {
                        rowLeftIndexValue = static_cast<int>(apexColumn.coeff(rowLeftIndex));
                    }

                    if (rowRightIndex < apexColumn.size()) {
                        rowRightIndexValue = static_cast<int>(apexColumn.coeff(rowRightIndex));
                    }

                    if (rowLeftIndexValue > 0 && rowRightIndexValue > 0) {
                        const int higherIndex = rowLeftIndexValue >= rowRightIndexValue ? rowLeftIndex : rowRightIndex;
                        apexesStartsToUse.push_back(higherIndex);
                        break;
                    }
                    else if (rowLeftIndexValue > 0) {
                        apexesStartsToUse.push_back(rowLeftIndex);
                        break;
                    }
                    else if (rowRightIndexValue > 0) {
                        apexesStartsToUse.push_back(rowRightIndex);
                        break;
                    }

                    if (rowLeftIndex <= 0 && rowRightIndex >= apexColumn.size() - 1) {
                        apexesStartsToUse.push_back(-1);
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
        QVector<float> *peakCorrelations,
        int *bestAnchorColumnIndex,
        int *bestAnchorRowIndex
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matBlockTrimmed.size() > 0); ree;
        e = ErrorUtils::isNotEmpty(apexStarts); ree;
        e = ErrorUtils::isEqual(apexStarts.size(), static_cast<int>(matBlockTrimmed.cols())); ree;

        const int colCount = static_cast<int>(matBlockTrimmed.cols());

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
                cosineSimSum += cosineSim;
                corrs.push_back(cosineSim);
            }

            if (cosineSimSum > bestCosineSimSum) {
                *bestAnchorColumnIndex = anchorCol;
                *peakCorrelations = corrs;
                bestCosineSimSum = cosineSimSum;
            }
        }

        const Eigen::VectorX<float> &anchor = matBlockTrimmed.col(*bestAnchorColumnIndex);
        const QVector<float> anchorVec = EigenUtils::convertEigenVectorToQVector(anchor);

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
    BestCorrelationResult *bestCorrelationResult
    ) const {

    ERR_INIT

    e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;
    e = ErrorUtils::isTrue(matriciesAndVecs.integrationVecIsValid()); ree;

    const int maxRows = matriciesAndVecs.intensityMatrix100.rows();
    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensityResized = peakIntegrationsVsIntensity;

    constexpr int topNIntegrations = 15; //TODO make this settable
    peakIntegrationsVsIntensityResized.resize(std::min(topNIntegrations, peakIntegrationsVsIntensityResized.size()));

    for (const QPair<PeakIntegrationIndexes, Intensity> &pii : peakIntegrationsVsIntensityResized) {

        const QPair<PeakIntegrationIndexes, Intensity> piiWorking = correctPeakIntegrationForSingleRow(
            pii,
            matriciesAndVecs.integrationVec
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

        const Eigen::VectorX<float> integrationVecSegment = matriciesAndVecs.integrationVec.segment(
            piiWorking.first.first,
            piiWorking.first.second - piiWorking.first.first + 1
            ).eval();

        const QPair<int, float> apexIndex = EigenUtils::returnTopIndexAndValue(integrationVecSegment);

        const QVector<int> apexStarts = findStartApexes(matBlockApexes, apexIndex.first);
        e = ErrorUtils::isEqual(apexStarts.size(), static_cast<int>(matBlockApexes.cols()));

        const Eigen::MatrixX<float> matBlockTrimmed = trimMatrixBlock(
            matBlock,
            apexStarts,
            m_pythiaParameters.stopThresholdFraction
            );

        QVector<float> peakCorrelations;
        int bestAnchorColumnIndex = -1;
        int bestAnchorRowIndex = -1;
        e = findBestAnchorColumn(
            matBlockTrimmed,
            apexStarts,
            &peakCorrelations,
            &bestAnchorColumnIndex,
            &bestAnchorRowIndex
            ); ree;

        const float peakCorrelationsSum
                = std::accumulate(peakCorrelations.begin(), peakCorrelations.end(), 0.0f);

        if (peakCorrelationsSum > bestCorrelationResult->peakCorrelationsSum) {

            bestCorrelationResult->peakCorrelations = peakCorrelations;
            bestCorrelationResult->peakCorrelationsSum = peakCorrelationsSum;
            bestCorrelationResult->matBlockTrimmedIntensity = matBlockTrimmed;
            bestCorrelationResult->bestAnchorColumnIndex = bestAnchorColumnIndex;
            bestCorrelationResult->bestAnchorRowIndex = bestAnchorRowIndex;
            bestCorrelationResult->apexStarts = apexStarts;
            bestCorrelationResult->peakIntegrationIndexes = piiWorking.first;

            constexpr float windowMultiplier1p5X = 1.5;
            const auto frameIndex1p5XMin
                = std::max(static_cast<FrameIndex>(std::round(piiWorking.first.first - (windowMultiplier1p5X * ogPeakLength))), 0);
            const auto frameIndex1p5XMax
                = std::min(static_cast<FrameIndex>(std::round(piiWorking.first.first + (windowMultiplier1p5X * ogPeakLength))), maxRows);
            const int peakLength1p5X = frameIndex1p5XMax - frameIndex1p5XMin + 1;
            Eigen::MatrixX<float> matBlock1p5X = matriciesAndVecs.intensityMatrix100.block(
                  frameIndex1p5XMin,
                  0,
                  peakLength1p5X,
                  matriciesAndVecs.intensityMatrix100.cols()
                  ).eval();

            constexpr float windowMultiplier2X = 2.0;
            const auto frameIndex2XMin
                = std::max(static_cast<FrameIndex>(std::round(piiWorking.first.first - (windowMultiplier2X * ogPeakLength))), 0);
            const auto frameIndex2XMax
                = std::min(static_cast<FrameIndex>(std::round(piiWorking.first.first + (windowMultiplier2X * ogPeakLength))), maxRows);
            const int peakLength2X = frameIndex2XMax - frameIndex2XMin + 1;
            Eigen::MatrixX<float> matBlock2X = matriciesAndVecs.intensityMatrix100.block(
                          frameIndex2XMin,
                          0,
                          peakLength2X,
                          matriciesAndVecs.intensityMatrix100.cols()
                          ).eval();

            bestCorrelationResult->matBlockTrimmedIntensityWindow1p5X = trimMatrixBlock(
                matBlock1p5X,
                apexStarts,
                m_pythiaParameters.stopThresholdFraction
                );
            e = calculatePeakCorrelations(
                bestCorrelationResult->matBlockTrimmedIntensityWindow1p5X,
                bestAnchorColumnIndex,
                &bestCorrelationResult->peakCorrelationsWindow1p5X
                ); ree;

            bestCorrelationResult->matBlockTrimmedIntensityWindow2X = trimMatrixBlock(
                        matBlock2X,
                        apexStarts,
                        m_pythiaParameters.stopThresholdFraction
                        );

            e = calculatePeakCorrelations(
                bestCorrelationResult->matBlockTrimmedIntensityWindow2X,
                bestAnchorColumnIndex,
                &bestCorrelationResult->peakCorrelationsWindow2X
                ); ree;

            const PeakIntegrationIndexes &p = piiWorking.first;
            const int pSize = p.second - p.first + 1;
            bestCorrelationResult->matBlockTrimmedMz = matriciesAndVecs.mzMatrix100.block(
                p.first,
                0,
                pSize,
                matBlockTrimmed.cols()
                ).eval();

            bestCorrelationResult->matBlockTrimmedIntensityShadows = matriciesAndVecs.intensityMatrix100Shadow.block(
                p.first,
                0,
                pSize,
                matBlockTrimmed.cols()
                ).eval();

            bestCorrelationResult->matBlockTrimmedIntensity45 = matriciesAndVecs.intensityMatrix45.block(
                p.first,
                0,
                pSize,
                matBlockTrimmed.cols()
                ).eval();

            bestCorrelationResult->matBlockTrimmedIntensity20 = matriciesAndVecs.intensityMatrix20.block(
                p.first,
                0,
                pSize,
                matBlockTrimmed.cols()
                ).eval();

            e = calculatePeakCorrelations(
                bestCorrelationResult->matBlockTrimmedIntensity45,
                bestAnchorColumnIndex,
                &bestCorrelationResult->peakCorrelations45
                ); ree;

            e = calculatePeakCorrelations(
                bestCorrelationResult->matBlockTrimmedIntensity20,
                bestAnchorColumnIndex,
                &bestCorrelationResult->peakCorrelations20
                ); ree;

        }
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

    Err setTheoreticalMs2Ions(
        const QVector<MS2Ion> &ms2IonsTheoretical,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(ms2IonsTheoretical.size() <= arraySizeMax); ree;
        for (int i = 0; i < ms2IonsTheoretical.size(); i++) {
            const MS2Ion &ms2Ion = ms2IonsTheoretical.at(i);
            candidateScores->featuresArray[CandidateScores::Features::MzSearched1 + i] = static_cast<float>(ms2Ion.mz);
            candidateScores->featuresArray[CandidateScores::Features::TheoIntensity1 + i] = ms2Ion.intensity;
        }

        ERR_RETURN
    }

    Err setFoundMs2Ions(
        const BestCorrelationResult &bestCorrelationResult,
        CandidateScores *candidateScores
        ) {

        ERR_INIT

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

        e = ErrorUtils::isTrue(mzMeanValsFound.size() <= arraySizeMax); ree;
        for (int i = 0; i < mzMeanValsFound.size(); i++) {
            candidateScores->featuresArray[CandidateScores::Features::MzFoundMean1 + i] = mzMeanValsFound.at(i);
        }

        e = ErrorUtils::isTrue(stdMeanValsFound.size() <= arraySizeMax); ree;
        for (int i = 0; i < stdMeanValsFound.size(); i++) {
            candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev1 + i] = stdMeanValsFound.at(i);
        }

        Eigen::VectorX<float> intensitySums = bestCorrelationResult.matBlockTrimmedIntensity.colwise().sum();
        intensitySums /= intensitySums.sum();

        e = ErrorUtils::isTrue(intensitySums.size() <= arraySizeMax); ree;
        for (int i = 0; i < intensitySums.size(); i++) {
            candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax1 + i] = intensitySums.coeff(i);
        }

        ERR_RETURN
    }

    Err setMassAccuracies(CandidateScores *candidateScores) {

        ERR_INIT

        constexpr int top12 = 12;

        for (int i = 0; i < top12; i++) {

            const float mzSearched = candidateScores->featuresArray[CandidateScores::Features::MzSearched1 + i];
            const float mzFound = candidateScores->featuresArray[CandidateScores::Features::MzFoundMean1 + i];

            if (MathUtils::tZero(mzSearched)) {
                continue;
            }

            const float absAccuracy = std::abs(MathUtils::calculateMassAccuracyPPM(mzSearched, mzFound));

            candidateScores->featuresArray[CandidateScores::Features::MzAccuracy1 + i] = 1.0f / std::min(std::max(absAccuracy, 0.001f), 200.0f);

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

            float intensityRatio = MathUtils::tZero(v1Max) ? 1.0f : v2Max / v1Max;
            intensityRatio = MathUtils::tZero(v2Max) ? 0.0f : intensityRatio;
            candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio1 + col] = intensityRatio;
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
        e = ErrorUtils::isTrue(mzPeakLengthsNormalized.size() <= arraySizeMax); ree;
        for (int i = 0; i < mzPeakLengthsNormalized.size(); i++) {
            candidateScores->featuresArray[CandidateScores::Features::MzPeakLengthsNorm1 + i] = mzPeakLengthsNormalized.at(i);
        }

        const QVector<int> &columnApexIndexes = bestCorrelationResult.apexStarts;
        const double columnApexIndexesMean = MathUtils::mean(columnApexIndexes);
        const double columnApexIndexesSize = columnApexIndexes.size();
        QVector<float> columnApexIndexRatiosToAnchor;
        std::transform(
                columnApexIndexes.begin(),
                columnApexIndexes.end(),
                std::back_inserter(columnApexIndexRatiosToAnchor),
                [columnApexIndexesMean, columnApexIndexesSize](int i){return (i - columnApexIndexesMean) / columnApexIndexesSize;}
                );

        e = ErrorUtils::isTrue(columnApexIndexRatiosToAnchor.size() <= arraySizeMax); ree;
        for (int i = 0; i < columnApexIndexRatiosToAnchor.size(); i++) {
            candidateScores->featuresArray[CandidateScores::Features::ColumnApexIndexRatiosToAnchor1 + i] = columnApexIndexRatiosToAnchor.at(i);
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
    const BestCorrelationResult& bestCorrelationResult,
    const QVector<float> &ms1Averagine,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(bestCorrelationResult.peakCorrelations); ree;
    e = ErrorUtils::isTrue(bestCorrelationResult.matBlockTrimmedIntensity.size() > 0); ree;

    candidateScores->initFeaturesArray();

    candidateScores->scanNumber = m_msFrameMzTarget->scanNumberFromFrameIndex(
        bestCorrelationResult.peakIntegrationIndexes.first
        + bestCorrelationResult.apexStarts.at(bestCorrelationResult.bestAnchorColumnIndex)
        );
    candidateScores->scanNumberStart = m_msFrameMzTarget->scanNumberFromFrameIndex(bestCorrelationResult.peakIntegrationIndexes.first);
    candidateScores->scanNumberEnd = m_msFrameMzTarget->scanNumberFromFrameIndex(bestCorrelationResult.peakIntegrationIndexes.second);

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

    candidateScores->featuresArray[CandidateScores::Features::ScanIonCount] = static_cast<float>(m_msFrameMzTarget->getScanPointsByScanNumber(candidateScores->scanNumber)->size());

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

    const float scanTimeDelta = std::abs(candidateScores->scanTime - candidateScores->scanTimePredicted);
    candidateScores->featuresArray[CandidateScores::Features::ScanTimeDelta] = scanTimeDelta;

    candidateScores->featuresArray[CandidateScores::Features::ScanTimePredicted] = candidateScores->scanTimePredicted;

    const double pdScanTime = std::sqrt(std::min(std::abs(scanTimeDelta), m_scanTimeRange) / m_scanTimeRange);
    candidateScores->featuresArray[CandidateScores::Features::ScanTimePd] = static_cast<float>(pdScanTime);

    const double pepLength = (-10.0 + candidateScores->targetDecoyCandidatePair->peptideString().size()) / 10.0;
    candidateScores->featuresArray[CandidateScores::Features::PeptideLengthNorm] = static_cast<float>(pepLength);

    const auto mz = candidateScores->targetDecoyCandidatePair->mz();
    candidateScores->featuresArray[CandidateScores::Features::MzNorm] = (mz - 600.0f) * 0.002f;
    candidateScores->featuresArray[CandidateScores::Features::IRTPredicted] = candidateScores->targetDecoyCandidatePair->iRt();
    candidateScores->featuresArray[CandidateScores::Features::Mass] = candidateScores->targetDecoyCandidatePair->mass();

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum45]
        = std::max(std::accumulate(bestCorrelationResult.peakCorrelations45.begin(), bestCorrelationResult.peakCorrelations45.begin() + top6, 0.0f), std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum20]
        = std::max(std::accumulate(bestCorrelationResult.peakCorrelations20.begin(), bestCorrelationResult.peakCorrelations20.begin() + top6, 0.0f), std::numeric_limits<float>::min());

    int bestAlignmentMatrixRowIndex = bestCorrelationResult.bestAnchorRowIndex;
    candidateScores->featuresArray[CandidateScores::Features::AllignedMaxIndexesCount] = static_cast<float>(std::count_if(
        bestCorrelationResult.apexStarts.begin(),
        bestCorrelationResult.apexStarts.end(),
        [bestAlignmentMatrixRowIndex](int i){return i == bestAlignmentMatrixRowIndex;}
        ));

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

    e = setTheoreticalMs2Ions(ms2IonsTheoritical, candidateScores); ree;

    e = setFoundMs2Ions(bestCorrelationResult, candidateScores); ree;

    e = setPeakShapeRatios(bestCorrelationResult, candidateScores); ree;

    e = setMassAccuracies(candidateScores); ree;

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
        const Eigen::VectorX<float> &anchorColumn,
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
        e = ErrorUtils::isTrue(anchorColumn.size() > 0); ree;

        *cosineSimMS1 = 0.01;
        *mzMS1Mean = 0.01;
        *mzMS1StDev = 10.0;
        *ppmDiff = 100.0;
        *ms1IntensitySum = 0.0;

        XICPoints xicPoints = turboXicMS1->extractPointsXIC(
            mzToExtract - massTol,
            mzToExtract + massTol
            );

        TurboXIC::filterXICPointsByScanNumber(frameIndexMin, frameIndexMax, &xicPoints);

        if (xicPoints.empty()) {
            ERR_RETURN
        }

        const Eigen::VectorX<float> xicVec = buildXicPointsVector(
            xicPoints,
            frameIndexMin,
            frameIndexMax
            ).segment(frameIndexMin, frameIndexMax - frameIndexMin + 1);

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

    const float isotopeDistance = S_GLOBAL_SETTINGS.ISO_DIFF / targetDecoyCandidatePair->charge();

    const float monoIsotopeMz = targetDecoyCandidatePair->mz();
    const float monoIsotopeShadowMz = monoIsotopeMz - isotopeDistance;
    const float c13isotopeMz1 = monoIsotopeMz + isotopeDistance;
    const float c13isotopeMz2 = monoIsotopeMz + (isotopeDistance * 2);

    const float massTol = MathUtils::calculatePPM(monoIsotopeMz, ppmTol);

    const Eigen::VectorX<float> anchorColumn
            = bestCorrelationResult.matBlockTrimmedIntensity.col(bestCorrelationResult.bestAnchorColumnIndex);

    e = calculateMs1Scores(
        anchorColumn,
        monoIsotopeMz,
        massTol,
        frameIndexMin,
        frameIndexMax,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound100],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound100],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound100PPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound100]
        ); ree;

    e = calculateMs1Scores(
        anchorColumn,
        monoIsotopeMz,
        massTol * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
        frameIndexMin,
        frameIndexMax,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim45MS1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound45],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound45],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound45PPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound45]
        ); ree;

    e = calculateMs1Scores(
        anchorColumn,
        monoIsotopeMz,
        massTol * S_GLOBAL_SETTINGS.TIGHT_2_FRACTION,
        frameIndexMin,
        frameIndexMax,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim20MS1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound20],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound20],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound20PPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound20]
        ); ree;

    e = calculateMs1Scores(
        anchorColumn,
        monoIsotopeShadowMz,
        massTol,
        frameIndexMin,
        frameIndexMax,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1PreMono],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMono],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundPreMono],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMonoPPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundPreMono]
        ); ree;

    e = calculateMs1Scores(
        anchorColumn,
        c13isotopeMz1,
        massTol,
        frameIndexMin,
        frameIndexMax,
        m_turboXicMS1,
        &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso1],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1PPM],
        &candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundIso1]
        ); ree;


    e = calculateMs1Scores(
        anchorColumn,
        c13isotopeMz2,
        massTol,
        frameIndexMin,
        frameIndexMax,
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

