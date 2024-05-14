//
// Created by anichols on 10/19/23.
//

#include "CandidateScorertron2.h"

#include "CandidateScores.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHillBuilder.h"
#include "IsotopicDistributionBuilder.h"
#include "ScoreOverseer.h"
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

    constexpr int filterLengthMs2 = 5;
    constexpr int filterLengthIntegration = 5;
    constexpr int order = 2;
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

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

CandidateScorertron::CandidateScorertron()
: m_topNMS2Ions(-1)
, m_xicPeakManager(nullptr)
, m_msFrameMzTarget(nullptr)
, d_ptr(new Private())
{}

CandidateScorertron::~CandidateScorertron() {}

Err CandidateScorertron::init(
    const PythiaParameters &pythiaParameters,
    const MsCalibratomatic &msCalibratomatic,
    int topNMS2Ions,
    XICPeakManager *xicPeakManager,
    MsFrame *msFrameMzTarget
    ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;
    e = ErrorUtils::isTrue(msFrameMzTarget->isValid()); ree;
    e = ErrorUtils::isAboveThreshold(
        topNMS2Ions,
        S_GLOBAL_SETTINGS.MIN_MS2_IONS,
        ErrorUtilsParam::IncludeThreshold)
    ; ree;

    m_pythiaParameters = pythiaParameters;
    m_topNMS2Ions = topNMS2Ions;
    m_xicPeakManager = xicPeakManager;
    m_msFrameMzTarget = msFrameMzTarget;

    if (msCalibratomatic.isInit()) {
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
    Eigen::MatrixX<float> mzMatrix100Shadow;
    Eigen::MatrixX<float> mzMatrix45;
    Eigen::MatrixX<float> mzMatrix20;

    Eigen::VectorX<float> integrationVec;

    [[nodiscard]] bool intensityMatriciesAreValid() const {
        return intensityMatrix100.size() > 0;
    }

    [[nodiscard]] bool integrationVecIsValid() const {
        return integrationVec.size() > 0;
    }
};

class BestCorrelationResult {

public:
    QVector<float> peakCorrelations;
    float peakCorrelationsSum = -1.0;
    Eigen::MatrixX<float> matBlockTrimmed;
    int bestAnchorColumnIndex = -1;
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

    Err getXICs(
        const QVector<MS2Ion> &ms2Ions,
        float ppmTol,
        XICPeakManager *xicPeakManager,
        QVector<XICPoints> *xicPointsVec100,
        QVector<XICPoints> *xicPointsVec100Shadows,
        QVector<XICPoints> *xicPointsVec45,
        QVector<XICPoints> *xicPointsVec20
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;
        e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;

        xicPointsVec100->reserve(ms2Ions.size());
        xicPointsVec100Shadows->reserve(ms2Ions.size());
        xicPointsVec45->reserve(ms2Ions.size());
        xicPointsVec20->reserve(ms2Ions.size());

        for (int i = 0; i < ms2Ions.size(); i++) {

            const MS2Ion &ms2Ion = ms2Ions.at(i);

            XICPoints xicPoints;
            e = xicPeakManager->getXIC(ms2Ion.mz, &xicPoints); ree;

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
                matMz->coeffRef(p.scanNumber, col) = p.mz;
            }
        }

        *matIntensity= EigenKernelUtils::applyKernelToEachColumnInMatrix(*matIntensity, kernelMs2);

        ERR_RETURN
    }

    Err buildIntegrationVector(
        const MatriciesAndVecs &matriciesAndVecs,
        const Eigen::VectorX<float> &kernelIntegration,
        Eigen::VectorX<float> *integrationVec
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matriciesAndVecs.intensityMatriciesAreValid()); ree;

        constexpr float intensityThresholdVal = 0.1;
        constexpr float countValue = 1.0;
        constexpr float minPeakCount = 3.9;

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

    Err initMatricesdAndVecs(
        const QVector<MS2Ion> &ms2Ions,
        const Eigen::VectorX<float> &kernelMs2,
        const Eigen::VectorX<float> &kernelIntegration,
        int topNMS2Ions,
        float ppmTol,
        XICPeakManager *xicPeakManager,
        MatriciesAndVecs *matriciesAndVecs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2Ions); ree;
        e = ErrorUtils::isTrue(xicPeakManager->isValid()); ree;
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

        e = buildEigenMatrix(
            xicPointsVec100,
            kernelMs2,
            frameIndexMax,
            &matriciesAndVecs->intensityMatrix100,
            &matriciesAndVecs->mzMatrix100
            ); ree;

        e = buildEigenMatrix(
            xicPointsVec45,
            kernelMs2,
            frameIndexMax,
            &matriciesAndVecs->intensityMatrix45,
            &matriciesAndVecs->mzMatrix45
            ); ree;

        e = buildEigenMatrix(
            xicPointsVec20,
            kernelMs2,
            frameIndexMax,
            &matriciesAndVecs->intensityMatrix20,
            &matriciesAndVecs->mzMatrix20
            ); ree;

        e = buildEigenMatrix(
            xicPointsVec100Shadow,
            kernelMs2,
            frameIndexMax,
            &matriciesAndVecs->intensityMatrix100Shadow,
            &matriciesAndVecs->mzMatrix100Shadow
            ); ree;

        e = buildIntegrationVector(
            *matriciesAndVecs,
            kernelIntegration,
            &matriciesAndVecs->integrationVec
            ); ree;

        ERR_RETURN
    }

}//namespace
Err CandidateScorertron::calculateScores(
    const QVector<MS2Ion> &ms2Ions,
    TargetDecoyCandidatePair* tdcp,
    CandidateScores *candidateScores
    ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2Ions); ree;

    MatriciesAndVecs matriciesAndVecs;
    e = initMatricesdAndVecs(
        ms2Ions,
        d_ptr->m_kernelMs2,
        d_ptr->m_kernelIntegration,
        m_topNMS2Ions,
        static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM),
        m_xicPeakManager,
        &matriciesAndVecs
        ); ree;

    if (m_msCalibratomatic.isInit()) {
        //TODO write code to truncate vector based on predicted scan time/frame index
    }

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensities;
    e = EigenUtils::simpleIntegrator(
        matriciesAndVecs.integrationVec,
        m_pythiaParameters.stopThresholdFraction,
        &peakIntegrationsVsIntensities
        ); ree;

    BestCorrelationResult bestCorrelationResult;
    e = processIntegrationVectorPeakIntegrations(
        matriciesAndVecs,
        peakIntegrationsVsIntensities,
        &bestCorrelationResult
        ); ree;

    qDebug() << bestCorrelationResult.peakCorrelations << bestCorrelationResult.peakCorrelationsSum << "SDLFJSDLJL";

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
        int *bestAnchorColumnIndex
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(matBlockTrimmed.size() > 0); ree;
        e = ErrorUtils::isNotEmpty(apexStarts); ree;
        e = ErrorUtils::isEqual(apexStarts.size(), static_cast<int>(matBlockTrimmed.cols())); ree;

        Eigen::MatrixX<float> matCount = matBlockTrimmed;
        matCount = (matBlockTrimmed.array() > 0.0f).select(1.0, matCount);
        const Eigen::VectorX<float> colSizes = matCount.colwise().sum();

        struct Temp {
            int fragIonIndex = -1;
            int apexIndex = -1;
            int peakLength = -1;
            float intensity = -1.0;
        };

        QVector<Temp> temps;
        for (int i = 0; i < apexStarts.size(); i++) {
            Temp t;
            t.fragIonIndex = i;
            t.apexIndex = apexStarts.at(i);
            t.peakLength = colSizes.coeff(i);
            t.intensity = matBlockTrimmed.coeff(t.apexIndex, i);
            temps.push_back(t);
        }

        const auto sortLogic = [](const Temp &l, const Temp &r) {
            if (l.apexIndex == r.apexIndex) {
                if (l.peakLength == r.peakLength) {
                    return l.intensity > r.intensity;
                }
                return l.peakLength > r.peakLength;
            }
            return l.apexIndex < r.apexIndex;
        };
        std::sort(temps.begin(), temps.end(), sortLogic);

        QHash<int, int> apexIndexVsCount;
        for (const Temp &t : temps) {
            apexIndexVsCount[t.apexIndex]++;
        }

        int bestApexIndex = -1;
        int bestApexIndexCnt = -1;
        for (auto it = apexIndexVsCount.begin(); it != apexIndexVsCount.end(); ++it) {
            if (it.value() > bestApexIndexCnt) {
                bestApexIndex = it.key();
                bestApexIndexCnt = it.value();
            }
        }

        for (const Temp &t : temps) {
            if (t.apexIndex == bestApexIndex) {
                *bestAnchorColumnIndex = t.fragIonIndex;
                break;
            }
        }

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

    QVector<QPair<PeakIntegrationIndexes, Intensity>> peakIntegrationsVsIntensityResized = peakIntegrationsVsIntensity;

    constexpr int topNIntegrations = 10; //TODO make this settable
    peakIntegrationsVsIntensityResized.resize(topNIntegrations);

    for (const QPair<PeakIntegrationIndexes, Intensity> &pii : peakIntegrationsVsIntensity) {

        const QPair<PeakIntegrationIndexes, Intensity> piiWorking = correctPeakIntegrationForSingleRow(
            pii,
            matriciesAndVecs.integrationVec
            );

        Eigen::MatrixX<float> mat = m_pythiaParameters.subtractShadows
                                  ? matriciesAndVecs.intensityMatrix100 - matriciesAndVecs.intensityMatrix100Shadow
                                  : matriciesAndVecs.intensityMatrix100;
        EigenUtils::thresholdMatrix(0.0f, &mat);

        Eigen::MatrixX<float> matBlock = mat.block(
              piiWorking.first.first,
              0,
              piiWorking.first.second - piiWorking.first.first + 1,
              mat.cols()
              ).eval();

        const QVector<QVector<int>> apexIndexesByColumn = getMatrxColumnApexes(matBlock);
        const Eigen::MatrixX<float> matBlockApexes = buildMatBlockApexes(
            matBlock,
            apexIndexesByColumn
            );

        const Eigen::VectorX<float> integrationVecSegment = matriciesAndVecs.integrationVec.segment(
            piiWorking.first.first,
            piiWorking.first.second - piiWorking.first.first + 1
            );

        const QPair<int, float> apexIndex = EigenUtils::returnTopIndexAndValue(integrationVecSegment);

        const QVector<int> apexStarts = findStartApexes(matBlockApexes, apexIndex.first);
        e = ErrorUtils::isEqual(apexStarts.size(), static_cast<int>(matBlockApexes.cols()));

        const Eigen::MatrixX<float> matBlockTrimmed = trimMatrixBlock(
            matBlock,
            apexStarts,
            m_pythiaParameters.stopThresholdFraction
            );

        int bestAnchorColumnIndex = -1;
        e = findBestAnchorColumn(
            matBlockTrimmed,
            apexStarts,
            &bestAnchorColumnIndex
            ); ree;

        QVector<float> peakCorrelations;
        e = calculatePeakCorrelations(
            matBlockTrimmed,
            bestAnchorColumnIndex,
            &peakCorrelations
            ); ree;

        const float peakCorrelationsSum
                = std::accumulate(peakCorrelations.begin(), peakCorrelations.end(), 0.0f);

        if (peakCorrelationsSum > bestCorrelationResult->peakCorrelationsSum) {
            bestCorrelationResult->peakCorrelations = peakCorrelations;
            bestCorrelationResult->peakCorrelationsSum = peakCorrelationsSum;
            bestCorrelationResult->matBlockTrimmed = matBlockTrimmed;
            bestCorrelationResult->bestAnchorColumnIndex = bestAnchorColumnIndex;
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
