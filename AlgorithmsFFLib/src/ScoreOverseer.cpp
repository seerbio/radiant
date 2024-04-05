//
// Created by anichols on 10/20/23.
//

#include "ScoreOverseer.h"

#include "BiophysicalCalcs.h"
#include "CandidateScores.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MS2Ion.h"
#include "ParallelUtils.h"
#include "TargetDecoyCandidatePair.h"


class Q_DECL_HIDDEN ScoreOverseer::Private
{

public:

    Private();
    ~Private();

    Err init(const PythiaParameters &pythiaParameters);

    Err buildAlignmentMatricies(
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QVector<MS2Ion> &ms2IonsTheoreticalShadows,
            const QHash<MzHashed , XICPoints> &mzHashedVsXICPoints,
            const QHash<MzHashed , XICPoints> &mzHashedVsXICPointsShadows,
            int *bestAlignmentMatrixRow,
            QVector<int> *columnApexIndexes,
            QVector<float> *intensityFoundMaxVec,
            QVector<int> *columnPeakLengths,
            QPair<int, int> *alignmentMatrixLimits
    );

    Err calculateCandidateAllignmentMetrics(
            QVector<float> *cosineSimsIndividual,
            QVector<float> *cosineSimsShadowIndividual,
            QVector<float> *shadowsIntensityRatioVec,
            Eigen::VectorX<float> *bestAnchorColumn,
            int *bestAnchorColumnIndex
    );

    Err calculateCosineSimSumTight(
            int anchorColumnIndex,
            float *cosineSimSum45,
            float *cosineSimSum20
    );

    Eigen::VectorX<float> m_gaussKernel;
    Eigen::MatrixX<float> m_intensityMatrix100;
    Eigen::MatrixX<float> m_intensityMatrix100Shadow;
    Eigen::MatrixX<float> m_intensityMatrix45;
    Eigen::MatrixX<float> m_intensityMatrix20;

    Eigen::VectorX<float> m_intensityMatrix100ApexRow;

    QVector<float> m_mzMeanValsFound;
    QVector<float> m_stdMeanValsFound;

    PythiaParameters m_pythiaParams;

};

ScoreOverseer::Private::Private() {}

ScoreOverseer::Private::~Private() {}

Err ScoreOverseer::Private::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

    m_gaussKernel = EigenKernelUtils::buildGaussianFilter1D<float>(
            pythiaParameters.filterLength,
            pythiaParameters.sigma
    );

    m_gaussKernel /= m_gaussKernel.maxCoeff();
    m_pythiaParams = pythiaParameters;

    ERR_RETURN
}

namespace {

    Err loadMzHashedVsFeatureFinderHillsToMatrix(
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed, XICPoints> &mzHashedVsXICPoints,
            Eigen::MatrixX<float> *mat,
            QVector<float> *mzMeanValsFound,
            QVector<float> *stdMeanValsFound
            ) {

        ERR_INIT

        e = ErrorUtils::isAboveThreshold(
                peakIntegrationIndexes.first,
                0,
                ErrorUtilsParam::IncludeThreshold
                ); ree;

        e = ErrorUtils::isTrue(mat->rows() > 0); ree;

        mat->setZero();

        const int rowCount = static_cast<int>(mat->rows());

        mzMeanValsFound->resize(ms2IonsTheoretical.size());
        stdMeanValsFound->resize(ms2IonsTheoretical.size());
        mzMeanValsFound->reserve(ms2IonsTheoretical.size());
        stdMeanValsFound->reserve(ms2IonsTheoretical.size());

        for (int col = 0; col < ms2IonsTheoretical.size(); col++) {

            if (col >= mat->cols()) {
                break;
            }

            const MS2Ion &ms2Ion = ms2IonsTheoretical.at(col);
            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const XICPoints &xicPoints = mzHashedVsXICPoints.value(mzHashed);

            QVector<float> eigenVectorLoaderIntensity(rowCount, 0.0f);
            eigenVectorLoaderIntensity.reserve(rowCount);

            QVector<float> mzVals;
            for (int i = 0; i < xicPoints.size(); i++) {

                if (xicPoints[i].scanNumber < peakIntegrationIndexes.first
                 || xicPoints[i].scanNumber > peakIntegrationIndexes.second) {
                    continue;
                }

                const int frameIndexAdjusted = xicPoints[i].scanNumber - peakIntegrationIndexes.first;

                if (frameIndexAdjusted < 0 || frameIndexAdjusted >= rowCount) {
                    continue;
                }

                eigenVectorLoaderIntensity[frameIndexAdjusted] = xicPoints[i].intensity;
                mzVals.push_back(xicPoints[i].mz);
            }

            (*mzMeanValsFound)[col] = mzVals.isEmpty() ? -1.0f : static_cast<float>(MathUtils::mean(mzVals));
            (*stdMeanValsFound)[col] = mzVals.isEmpty() ? -1.0f : static_cast<float>(MathUtils::stDev(mzVals));

            mat->col(col) = EigenUtils::convertQVectorToEigenVector(eigenVectorLoaderIntensity);
        }

        ERR_RETURN
    }

    Err loadMzHashedVsFeatureFinderHillsToMatrixTight(
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed, XICPoints> &mzHashedVsXICPoints,
            float ppmTol,
            const QPair<int, int> &alignmentMatrixLimits,
            Eigen::MatrixX<float> *matTight1,
            Eigen::MatrixX<float> *matTight2
            ) {

        ERR_INIT

        e = ErrorUtils::isAboveThreshold(peakIntegrationIndexes.first, 0, ErrorUtilsParam::IncludeThreshold); ree;
        e = ErrorUtils::isBelowThreshold(
                peakIntegrationIndexes.first,
                peakIntegrationIndexes.second,
                ErrorUtilsParam::ExcludeThreshold
                ); ree;

        e = ErrorUtils::isTrue(matTight1->rows() > 0); ree;
        e = ErrorUtils::isTrue(matTight2->rows() > 0); ree;

        const int rowCount = static_cast<int>(matTight1->rows());
        matTight1->setZero();
        matTight2->setZero();

        const int top6OrLess = std::min(6, ms2IonsTheoretical.size());

        for (int col = 0; col < top6OrLess; col++) {

            const MS2Ion &ms2Ion = ms2IonsTheoretical.at(col);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

            QVector<float> eigenVectorLoaderTight1(rowCount);
            eigenVectorLoaderTight1.reserve(rowCount);

            QVector<float> eigenVectorLoaderTight2(rowCount);
            eigenVectorLoaderTight2.reserve(rowCount);

            const XICPoints &xicPoints = mzHashedVsXICPoints.value(mzHashed);

            for (int i = 0; i < xicPoints.size(); i++) {

                if (xicPoints[i].scanNumber < peakIntegrationIndexes.first
                 || xicPoints[i].scanNumber > peakIntegrationIndexes.second) {
                    continue;
                }

                const int frameIndexAdjusted = xicPoints[i].scanNumber - peakIntegrationIndexes.first;
                if (frameIndexAdjusted < 0 || frameIndexAdjusted >= rowCount) {
                    continue;
                }

                const float mzTolTight1 = MathUtils::calculatePPM(ms2Ion.mz, S_GLOBAL_SETTINGS.TIGHT_1_FRACTION * ppmTol);
                const float mzMinTight1 = ms2Ion.mz - mzTolTight1;
                const float mzMaxTight1 = ms2Ion.mz + mzTolTight1;

                const double mz = xicPoints[i].mz;

                if (mz < mzMinTight1 || mz > mzMaxTight1) {
                    continue;
                }

                eigenVectorLoaderTight1[frameIndexAdjusted] = xicPoints[i].intensity;

                const double mzTolTight2 = MathUtils::calculatePPM(ms2Ion.mz, S_GLOBAL_SETTINGS.TIGHT_2_FRACTION * ppmTol);
                const double mzMinTight2 = ms2Ion.mz - mzTolTight2;
                const double mzMaxTight2 = ms2Ion.mz + mzTolTight2;

                if (mz < mzMinTight2 || mz > mzMaxTight2) {
                    continue;
                }

                eigenVectorLoaderTight2[frameIndexAdjusted] = xicPoints[i].intensity;

            }

            matTight1->col(col) = EigenUtils::convertQVectorToEigenVector(eigenVectorLoaderTight1);
            matTight2->col(col) = EigenUtils::convertQVectorToEigenVector(eigenVectorLoaderTight2);
        }

        *matTight1 = matTight1->block(
                alignmentMatrixLimits.first,
                0,
                alignmentMatrixLimits.second - alignmentMatrixLimits.first + 1,
                matTight1->cols()
                ).eval();

        *matTight2 = matTight2->block(
                alignmentMatrixLimits.first,
                0,
                alignmentMatrixLimits.second - alignmentMatrixLimits.first + 1,
                matTight2->cols()
        ).eval();

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

    Err findBestApexRowInMatrix(
            const Eigen::MatrixX<float> &intensityMatrix,
            int *bestApexRowIndex
            ) {

        ERR_INIT

        using RowMax = float;
        using RowSum = float;

        struct SearchObject {
            RowMax rowMax = -1.0;
            RowSum rowSum = -1.0;
            int index = -1;
        };

        const Eigen::VectorX<float> rowSums = intensityMatrix.rowwise().sum();
        Eigen::MatrixX<float> rowCountMat = (intensityMatrix.array() / intensityMatrix.array());
        EigenUtils::replaceNaN(0.0f, &rowCountMat);
        const Eigen::VectorX<float> rowCount = rowCountMat.rowwise().sum();

        const QVector<float> rowSumsVec = EigenUtils::convertEigenVectorToQVector(rowSums);
        const QVector<float> rowCountVec = EigenUtils::convertEigenVectorToQVector(rowCount);

        QVector<QPair<RowMax, RowSum>> rowMaxVsRowSum;
        e = ParallelUtils::zip(rowCountVec, rowSumsVec, &rowMaxVsRowSum); ree;

        int counter = 0;
        const auto insertLogic = [&counter](const QPair<RowMax, RowSum> &item){
            SearchObject so;
            so.rowMax = item.first;
            so.rowSum = item.second;
            so.index = counter++;
            return so;
        };

        QVector<SearchObject> searchObjects;
        std::transform(
                rowMaxVsRowSum.begin(),
                rowMaxVsRowSum.end(),
                std::back_inserter(searchObjects),
                insertLogic
        );

        const auto sortRowMaxThenRowSumDescLogic = [](const SearchObject &l, const SearchObject &r){
            if (static_cast<int>(l.rowMax) == static_cast<int>(r.rowMax)) {
                return l.rowSum < r.rowSum;
            }
            return static_cast<int>(l.rowMax) < static_cast<int>(r.rowMax);
        };

        std::sort(searchObjects.rbegin(), searchObjects.rend(), sortRowMaxThenRowSumDescLogic);

        *bestApexRowIndex = searchObjects.front().index;

        ERR_RETURN
    }

    Err removeInterferingPeaksInColumn(
            int bestApexRowIndex,
            Eigen::VectorX<float> *vec
            ) {

        ERR_INIT

        e = ErrorUtils::isBelowThreshold(
                bestApexRowIndex,
                static_cast<int>(vec->rows()),
                ErrorUtilsParam::ExcludeThreshold
        ); ree;

        const double kindaNearZero = 0.0001;

        int leftStartIndex = bestApexRowIndex;
        bool startZeroingLeft = false;
        double previousValLeftStartIndex = vec->coeff(leftStartIndex);
        while (--leftStartIndex >= 0) {

            const double currentVal = vec->coeff(leftStartIndex);
            if (currentVal > kindaNearZero && !startZeroingLeft && currentVal < previousValLeftStartIndex) {
                previousValLeftStartIndex = currentVal;
                continue;
            }

            startZeroingLeft = true;
            vec->coeffRef(leftStartIndex) = 0.0;
        }

        int rightStartIndex = bestApexRowIndex;
        bool startZeroingRight = false;
        double previousValRightStartIndex = vec->coeff(rightStartIndex);
        while (++rightStartIndex < vec->rows()) {

            const double currentVal = vec->coeff(rightStartIndex);
            if (currentVal > kindaNearZero && !startZeroingRight && currentVal < previousValRightStartIndex) {
                previousValRightStartIndex = currentVal;
                continue;
            }

            startZeroingRight = true;
            vec->coeffRef(rightStartIndex) = 0.0;
        }

        ERR_RETURN
    }

    Err removeInterferingPeaksInMatrix(
            const Eigen::MatrixX<float> &_intensityMatrix,
            Eigen::MatrixX<float> *intensityMatrixNoInterference,
            int *bestApexRowIndex,
            QVector<int> *columnApexIndexes,
            QVector<float> *intensityFoundMaxVec,
            QVector<int> *columnPeakLengths
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(_intensityMatrix.maxCoeff() > 0); ree;

        Eigen::MatrixX<float> intensityMatrixCopy = _intensityMatrix;

        e = findBestApexRowInMatrix(
                intensityMatrixCopy,
                bestApexRowIndex
                ); ree;

        for (int col = 0; col < intensityMatrixCopy.cols(); col++) {

            Eigen::VectorX<float> colVec = intensityMatrixCopy.col(col);
            const QVector<int> columnApexes = EigenUtils::apexesIndexesOnly(colVec);

            if (columnApexes.isEmpty()) {
                continue;
            }

            const int bestApexIndexColumn = columnApexes.at(MathUtils::closest(columnApexes, *bestApexRowIndex));

            e = removeInterferingPeaksInColumn(bestApexIndexColumn, &colVec); ree;
            const float colApexIntensity = colVec.coeff(bestApexIndexColumn);

            intensityMatrixCopy.col(col) = colVec;
            columnApexIndexes->push_back(bestApexIndexColumn);
            intensityFoundMaxVec->push_back(colApexIntensity);

            Eigen::VectorX<float> zerosVec = colVec.array() / colVec.array();
            EigenUtils::replaceNaN(0.0f, &zerosVec);

            const int nonZeros = static_cast<int>(std::round(zerosVec.sum()));
            columnPeakLengths->push_back(static_cast<int>(nonZeros));

        }

        *intensityMatrixNoInterference = intensityMatrixCopy;

        ERR_RETURN
    }

    QPair<int, int> simpleIntegrator(
            const Eigen::MatrixX<float> &mat,
            int apexRowIndex
            ) {

        const float nearZero = 0.001;

        const Eigen::VectorX<float> vec = mat.rowwise().sum();

        int rightStopIndex = apexRowIndex;
        int rightCurrentIndex = apexRowIndex;
        while (rightCurrentIndex < mat.rows()) {

            const float currentValue = vec(rightCurrentIndex);
            if (currentValue < nearZero) {
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
        while (leftCurrentIndex < vec.size()) {

            const float currentValue = vec(leftCurrentIndex);
            if (currentValue < nearZero) {
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

    Eigen::MatrixX<float> removeEmptyRowsFromMatrix(
            const Eigen::MatrixX<float> &mat,
            int apexRowIndex,
            QPair<int, int> *anchorColumnPeakLimits
            ) {

        *anchorColumnPeakLimits = simpleIntegrator(
                mat,
                apexRowIndex
                );

        const Eigen::MatrixX<float> newMatrix = mat.block(
                anchorColumnPeakLimits->first,
                0,
                anchorColumnPeakLimits->second - anchorColumnPeakLimits->first + 1,
                mat.cols()
                ).eval();

        return newMatrix;
    }

}// namespace
Err ScoreOverseer::Private::buildAlignmentMatricies(
        const PeakIntegrationIndexes &peakIntegrationIndexes,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QVector<MS2Ion> &ms2IonsTheoreticalShadows,
        const QHash<MzHashed , XICPoints> &mzHashedVsXICPoints,
        const QHash<MzHashed , XICPoints> &mzHashedVsXICPointsShadows,
        int *bestAlignmentMatrixRowIndex,
        QVector<int> *columnApexIndexes,
        QVector<float> *intensityFoundMaxVec,
        QVector<int> *columnPeakLengths,
        QPair<int, int> *alignmentMatrixLimits
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;
    e = ErrorUtils::isNotEmpty(ms2IonsTheoreticalShadows); ree;
    e = ErrorUtils::isEqual(ms2IonsTheoretical.size(), ms2IonsTheoreticalShadows.size()); ree;
    e = ErrorUtils::isBelowThreshold(
            peakIntegrationIndexes.first,
            peakIntegrationIndexes.second,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;

    m_mzMeanValsFound.clear();
    m_stdMeanValsFound.clear();

    const int cols = ms2IonsTheoretical.size();
    const int rows = std::max(peakIntegrationIndexes.second - peakIntegrationIndexes.first + 1, m_pythiaParams.filterLength);

    m_intensityMatrix100.resize(rows, cols);
    e = loadMzHashedVsFeatureFinderHillsToMatrix(
            peakIntegrationIndexes,
            ms2IonsTheoretical,
            mzHashedVsXICPoints,
            &m_intensityMatrix100,
            &m_mzMeanValsFound,
            &m_stdMeanValsFound
            ); ree;

    Eigen::MatrixX<float> intensityMatrixNoInterference;
    e = removeInterferingPeaksInMatrix(
            m_intensityMatrix100,
            &intensityMatrixNoInterference,
            bestAlignmentMatrixRowIndex,
            columnApexIndexes,
            intensityFoundMaxVec,
            columnPeakLengths
            ); ree;

    m_intensityMatrix100ApexRow = intensityMatrixNoInterference.row(*bestAlignmentMatrixRowIndex);

    m_intensityMatrix100 = removeEmptyRowsFromMatrix(
            intensityMatrixNoInterference,
            *bestAlignmentMatrixRowIndex,
            alignmentMatrixLimits
            );

    m_intensityMatrix100 = applyGaussSmoothRowWiseToMatrix(
            m_intensityMatrix100,
            m_pythiaParams,
            m_gaussKernel
    );

//#define CHECK_ALIGNMENT_MATRIX
#ifdef CHECK_ALIGNMENT_MATRIX
    std::cout << m_intensityMatrix100 << std::endl;
    std::cout << "****" << std::endl;
#endif

    m_intensityMatrix100Shadow.resize(rows, cols);
    QVector<float> unused1;
    QVector<float> unused2;
    e = loadMzHashedVsFeatureFinderHillsToMatrix(
            peakIntegrationIndexes,
            ms2IonsTheoreticalShadows,
            mzHashedVsXICPointsShadows,
            &m_intensityMatrix100Shadow,
            &unused1,
            &unused2
    ); ree;

    const int rowSize = alignmentMatrixLimits->second - alignmentMatrixLimits->first + 1;
    if (rowSize < m_intensityMatrix100Shadow.rows()) {
        m_intensityMatrix100Shadow = m_intensityMatrix100Shadow.block(
                alignmentMatrixLimits->first,
                0,
                rowSize,
                m_intensityMatrix100Shadow.cols()
        ).eval();

        m_intensityMatrix100Shadow = applyGaussSmoothRowWiseToMatrix(
                m_intensityMatrix100Shadow,
                m_pythiaParams,
                m_gaussKernel
        );
    }

    const int tightColsMax = 6;
    m_intensityMatrix45.resize(rows, tightColsMax);
    m_intensityMatrix20.resize(rows, tightColsMax);
    e = loadMzHashedVsFeatureFinderHillsToMatrixTight(
            peakIntegrationIndexes,
            ms2IonsTheoretical,
            mzHashedVsXICPoints,
            static_cast<float>(m_pythiaParams.ms2ExtractionWidthPPM),
            *alignmentMatrixLimits,
            &m_intensityMatrix45,
            &m_intensityMatrix20
    ); ree;

    ERR_RETURN
}

namespace {

    Err calcBestCosineSimSumLogic(
            const Eigen::MatrixX<float> &intensityMatrixIntegratedLimitsSmoothed,
            double cosineSimToAnchorThreshold,
            QVector<float> *bestCosineSimsIndividual,
            Eigen::VectorX<float> *bestAnchorColumn,
            int *bestAnchorColumnIndex
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrixIntegratedLimitsSmoothed.rows() > 0); ree;

        float bestCosineSimSum = 0.0;
        *bestAnchorColumnIndex = -1;

        const int topCols6OrLess = std::min(6, static_cast<int>(intensityMatrixIntegratedLimitsSmoothed.cols()));

        for (int col = 0; col < topCols6OrLess; col++) {

            QVector<float> bestCosineSimsIndividualAnchor;
            const Eigen::VectorX<float> &anchorColumn = intensityMatrixIntegratedLimitsSmoothed.col(col);

            if (MathUtils::tZero(anchorColumn.maxCoeff())) {
                continue;
            }

            for (int j = 0; j < intensityMatrixIntegratedLimitsSmoothed.cols(); j++) {

                const Eigen::VectorX<float> &altColumn = intensityMatrixIntegratedLimitsSmoothed.col(j);
                if (MathUtils::tZero(altColumn.maxCoeff())) {
                    bestCosineSimsIndividualAnchor.push_back(0.0);
                    continue;
                }

                float cosineSimToAnchor;
                e = EigenUtils::cosineSimilarity(
                        anchorColumn,
                        altColumn,
                        &cosineSimToAnchor
                        ); ree;

                if (cosineSimToAnchor < cosineSimToAnchorThreshold) {
                    bestCosineSimsIndividualAnchor.push_back(0.0);
                    continue;
                }

                bestCosineSimsIndividualAnchor.push_back(cosineSimToAnchor);
            }

            const float cosineSimSumAnchor
                    = std::accumulate(bestCosineSimsIndividualAnchor.begin(), bestCosineSimsIndividualAnchor.end(), 0.0f);

            if (cosineSimSumAnchor > bestCosineSimSum) {
                bestCosineSimSum = cosineSimSumAnchor;
                *bestCosineSimsIndividual = bestCosineSimsIndividualAnchor;
                *bestAnchorColumn = anchorColumn;
                *bestAnchorColumnIndex = col;
            }
        }

        ERR_RETURN
    }

}//namespace
Err ScoreOverseer::Private::calculateCandidateAllignmentMetrics(
        QVector<float> *cosineSimsIndividual,
        QVector<float> *cosineSimsShadowIndividual,
        QVector<float> *shadowsIntensityRatioVec,
        Eigen::VectorX<float> *bestAnchorColumn,
        int *bestAnchorColumnIndex
        ) {

    ERR_INIT

    if (m_intensityMatrix100Shadow.cols() > 0) {

        e = ErrorUtils::isEqual(m_intensityMatrix100.cols(), m_intensityMatrix100Shadow.cols()); ree;
        e = ErrorUtils::isEqual(m_intensityMatrix100.rows(), m_intensityMatrix100Shadow.rows()); ree;

        if (m_pythiaParams.subtractShadows) {
            m_intensityMatrix100 -= m_intensityMatrix100Shadow;
            EigenUtils::thresholdMatrix(0.0f, &m_intensityMatrix100);

            if (MathUtils::tZero(m_intensityMatrix100.maxCoeff())) {
                ERR_RETURN
            }
        }

        for (int col = 0; col < m_intensityMatrix100.cols(); col++) {

            const Eigen::VectorX<float> &mzPeak = m_intensityMatrix100.col(col);
            const Eigen::VectorX<float> &mzPeakShadow = m_intensityMatrix100Shadow.col(col);

            float cosineSimToShadow;
            e = EigenUtils::cosineSimilarity(
                    mzPeak,
                    mzPeakShadow,
                    &cosineSimToShadow
                    ); ree;

            cosineSimsShadowIndividual->push_back(cosineSimToShadow);

            if (cosineSimToShadow > 0.0 && !MathUtils::tZero(cosineSimToShadow)) {
                const float ratio = mzPeak.maxCoeff() / mzPeakShadow.maxCoeff();
                shadowsIntensityRatioVec->push_back(ratio);
            }
            else{
                shadowsIntensityRatioVec->push_back(0.0);
            }
        }

    }

    if (MathUtils::tZero(m_intensityMatrix100.maxCoeff())) {
        QVector<float> v(static_cast<int>(m_intensityMatrix100.cols()) , std::numeric_limits<float>::min());
        *cosineSimsIndividual = v;
        ERR_RETURN
    }

    e = calcBestCosineSimSumLogic(
            m_intensityMatrix100,
            m_pythiaParams.cosineSimToAnchorThreshold,
            cosineSimsIndividual,
            bestAnchorColumn,
            bestAnchorColumnIndex
            ); ree;

    ERR_RETURN
}

Err ScoreOverseer::Private::calculateCosineSimSumTight(
        int anchorColumnIndex,
        float *cosineSimSum45,
        float *cosineSimSum20
        ) {

    ERR_INIT

    *cosineSimSum45 = 0.0;
    *cosineSimSum20 = 0.0;

    const int colCount100 = static_cast<int>(m_intensityMatrix100.cols());
    const int colCount45 = static_cast<int>(m_intensityMatrix45.cols());
    if (colCount45 == 0) {
        ERR_RETURN
    }

    Eigen::MatrixX<float> intensityMat100TopCols;
    Eigen::MatrixX<float> intensityMatrix45;
    Eigen::MatrixX<float> intensityMatrix20;

    if (colCount100 > colCount45) {
        intensityMat100TopCols = m_intensityMatrix100.block(0, 0, m_intensityMatrix100.rows(), colCount45).eval();
        intensityMatrix45 = m_intensityMatrix45;
        intensityMatrix20 = m_intensityMatrix20;
    }
    else{
        intensityMat100TopCols = m_intensityMatrix100;
        intensityMatrix45 = m_intensityMatrix45.block(0, 0, m_intensityMatrix45.rows(), colCount100).eval();
        intensityMatrix20 = m_intensityMatrix20.block(0, 0, m_intensityMatrix20.rows(), colCount100).eval();
    }

    Eigen::MatrixX<float> peakMaskMat = intensityMat100TopCols.array() / intensityMat100TopCols.array();
    EigenUtils::replaceNaN(std::numeric_limits<float>::min() , &peakMaskMat);

    e = ErrorUtils::isEqual(peakMaskMat.cols(), intensityMatrix45.cols()); ree;
    e = ErrorUtils::isEqual(peakMaskMat.cols(), intensityMatrix20.cols()); ree;

    intensityMatrix45 = intensityMatrix45.array() * peakMaskMat.array();
    intensityMatrix20 = intensityMatrix20.array() * peakMaskMat.array();

    const Eigen::VectorX<float> anchorCol100 = intensityMat100TopCols.col(anchorColumnIndex);
    const Eigen::VectorX<float> anchorCol45 = intensityMatrix45.col(anchorColumnIndex);
    const Eigen::VectorX<float> anchorCol20 = intensityMatrix20.col(anchorColumnIndex);

    for (int col = 0; col < intensityMatrix45.cols(); col++) {
        const Eigen::VectorX<float> col45 = intensityMatrix45.col(col);
        const Eigen::VectorX<float> col20 = intensityMatrix20.col(col);

        float cosineSim45;
        e = EigenUtils::cosineSimilarity(
                col45,
                anchorCol45,
                &cosineSim45
                ); ree;

        float cosineSim20;
        e = EigenUtils::cosineSimilarity(
                col20,
                anchorCol45, //NOTE: for some reason using anchorCol45 her insted of anchorCol20 yields
                &cosineSim20        //much better results. around 1.5K more PSMs for Astral data.
                ); ree;

        *cosineSimSum45 += cosineSim45;
        *cosineSimSum20 += cosineSim20;
    }

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


ScoreOverseer::ScoreOverseer()
: d_ptr(new Private())
{}

ScoreOverseer::~ScoreOverseer() {}

Err ScoreOverseer::init(
        const PythiaParameters &pythiaParameters,
        const MzTargetKey &targetKey,
        MsFrame *msFrameTandem,
        TurboXIC *turboXICMS1
) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = d_ptr->init(pythiaParameters); ree;
    e = ErrorUtils::isTrue(msFrameTandem->isValid()); ree;

    m_msFrameTandem = msFrameTandem;
    m_mzTargetKey = targetKey;
    m_turboXICMS1 = turboXICMS1;

    ERR_RETURN
}

namespace {

    Err calculateSpectrumMetrics(
            const Eigen::VectorX<float> &intensityApexVector,
            const QVector<MS2Ion> &ms2IonsCandidate,
            double *cosineSimSpectrum,
            double *klDivSpectrum
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2IonsCandidate); ree;
        e = ErrorUtils::isTrue(intensityApexVector.rows() > 0); ree;

        QVector<float> intensityVals;
        std::transform(
                ms2IonsCandidate.begin(),
                ms2IonsCandidate.end(),
                std::back_inserter(intensityVals),
                [](const MS2Ion &ms2Ion){return ms2Ion.intensity;}
        );

        const Eigen::VectorX<float> intensityValsTheo = EigenUtils::convertQVectorToEigenVector(intensityVals);

        Eigen::VectorX<float> intensityApexVectorAppended = Eigen::VectorX<float>::Zero(intensityVals.size());
        intensityApexVectorAppended.head(intensityApexVector.size()) = intensityApexVector;
        
        float cosineSimSpectrumF;
        e = EigenUtils::cosineSimilarity(
                intensityApexVectorAppended,
                intensityValsTheo,
                &cosineSimSpectrumF
                ); ree;

        *cosineSimSpectrum = cosineSimSpectrumF;

        float klDivSpectrumF;
        e = EigenUtils::klDivergence(
                intensityApexVector,
                intensityValsTheo,
                &klDivSpectrumF
        ); ree;
        *klDivSpectrum = klDivSpectrumF;

        ERR_RETURN
    }

    Err calculateSpectrumMetricsOverPeak(
            const Eigen::MatrixX<float> &intensityMatrix,
            const QVector<MS2Ion> &ms2IonsCandidate,
            int apexColumnIndex,
            float *cosineSimOverPeak
            ) {

        ERR_INIT

        const int columnSizeMax = 6;
        const int columnsSize = std::min(
                columnSizeMax,
                static_cast<int>(intensityMatrix.cols())
                );

        const Eigen::MatrixX<float> intensityMatrixTrunc = intensityMatrix.leftCols(columnsSize);

        QVector<float> ms2IonIntensities;
        std::transform(
                ms2IonsCandidate.begin(),
                ms2IonsCandidate.end(),
                std::back_inserter(ms2IonIntensities),
                [](const MS2Ion &ms2Ion){return ms2Ion.intensity;}
                );

        const Eigen::VectorX<float> ms2IonsVec
            = EigenUtils::convertQVectorToEigenVector(ms2IonIntensities.mid(0, columnsSize));

        Eigen::MatrixX<float> ms2IonMat(intensityMatrix.rows(), columnsSize);
        for (int row = 0; row < ms2IonMat.rows(); row++) {
            ms2IonMat.row(row) = ms2IonsVec;
        }

        Eigen::VectorX<float> cosineSimByRow;
        e = EigenUtils::rowWiseCosineSimilarOfMatrices(
                intensityMatrixTrunc,
                ms2IonMat,
                &cosineSimByRow
                ); ree;

        const Eigen::VectorX<float> &apexColumnVec = intensityMatrixTrunc.col(apexColumnIndex);
        const float apexColumnVecSum = apexColumnVec.sum();

        e = ErrorUtils::isFalse(MathUtils::tZero(apexColumnVecSum)); ree;
        *cosineSimOverPeak = cosineSimByRow.dot(apexColumnVec) / apexColumnVecSum;

        ERR_RETURN
    }

    Err calculateMS1Corr(
            const Eigen::VectorX<float> &bestAnchorColumn,
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            float mzTarget,
            float ppmTol,
            TurboXIC *turboXicMS1,
            Eigen::VectorX<float> *gaussKernel,
            float *cosineSimMS1,
            float *ms1MeanMzFound,
            float *ms1StDevMzFound,
            float *ms1IntensityFound
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(bestAnchorColumn.size() > 0); ree;

        const float massTol = MathUtils::calculatePPM(mzTarget, ppmTol);
        const float mzStart = mzTarget - massTol;
        const float mzEnd = mzTarget + massTol;

        const XICPoints xicPoints = turboXicMS1->extractPointsXIC(
                mzStart,
                mzEnd
                );

        if (xicPoints.empty()) {
            *cosineSimMS1 = std::numeric_limits<float>::min();
            ERR_RETURN
        }

        const int ms1VecSize = static_cast<int>(bestAnchorColumn.size());

        std::vector<float> ms1IntensityLoadingVec(ms1VecSize, 0.0f);
        ms1IntensityLoadingVec.reserve(ms1VecSize);

        std::vector<float> ms1MzValsLoadingVec;
        std::transform(
                xicPoints.begin(),
                xicPoints.end(),
                std::back_inserter(ms1MzValsLoadingVec),
                [](const XICPoint &xp){return xp.mz;}
                );

        for (int i = 0; i < xicPoints.size(); i++) {

            const XICPoint &xp = xicPoints.at(i);

            const FrameIndex frameIndex = xp.scanNumber - peakIntegrationIndexes.first;

            if (frameIndex >= ms1VecSize || frameIndex < 0) {
                continue;
            }

            ms1IntensityLoadingVec[frameIndex] = xp.intensity;

        }

        Eigen::VectorX<float> ms1Vec(ms1VecSize);
        ms1Vec.setZero();
        ms1Vec
            = EigenUtils::convertQVectorToEigenVector(QVector<float>(ms1IntensityLoadingVec.begin(), ms1IntensityLoadingVec.end()));

        if (MathUtils::tZero(ms1Vec.maxCoeff())) {
            *cosineSimMS1 = std::numeric_limits<float>::min();
            ERR_RETURN
        }

        ms1Vec = EigenKernelUtils::convolveVectorWithKernel(ms1Vec, *gaussKernel);
        EigenUtils::replaceNaN(std::numeric_limits<float>::min(), &ms1Vec);

        float cosineSimMS1Local;
        e = EigenUtils::cosineSimilarity(
                bestAnchorColumn,
                ms1Vec,
                &cosineSimMS1Local
                ); ree;

        if(std::isnan(cosineSimMS1Local)) {
            cosineSimMS1Local = std::numeric_limits<float>::min();
        }

        *cosineSimMS1 = cosineSimMS1Local;
        *ms1MeanMzFound = static_cast<float>(MathUtils::mean(ms1MzValsLoadingVec));
        *ms1StDevMzFound = static_cast<float>(MathUtils::stDev(ms1MzValsLoadingVec));
        *ms1IntensityFound = std::accumulate(ms1IntensityLoadingVec.begin(), ms1IntensityLoadingVec.end(), 0.0f);

        ERR_RETURN
    }

}//namespace
Err ScoreOverseer::buildScores(
        TargetDecoyCandidatePair* targetDecoyCandidatePair,
        const QPair<float, float> &scanTimeMinMax,
        const PeakIntegrationIndexes &peakIntegrationIndexes,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QHash<MzHashed , XICPoints> &mzHashedVsXICPoints,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
        const QHash<MzHashed , XICPoints> &mzHashedVsXICPointsShadows,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;
    e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;
    e = ErrorUtils::isNotEmpty(ms2IonsTheoreticalIsotopeShadows); ree;

    QElapsedTimer et;
    et.start();

    candidateScores->targetDecoyCandidatePair = targetDecoyCandidatePair;
    candidateScores->initFeaturesArray();

    int bestAlignmentMatrixRowIndex;
    QVector<int> columnApexIndexes;
    QVector<float> intensityFoundMaxVec;
    QVector<int> mzPeakLengthsVec;
    QPair<int, int> alignmentMatrixLimits;
    e = d_ptr->buildAlignmentMatricies(
            peakIntegrationIndexes,
            ms2IonsTheoretical,
            ms2IonsTheoreticalIsotopeShadows,
            mzHashedVsXICPoints,
            mzHashedVsXICPointsShadows,
            &bestAlignmentMatrixRowIndex,
            &columnApexIndexes,
            &intensityFoundMaxVec,
            &mzPeakLengthsVec,
            &alignmentMatrixLimits
            ); ree;

    candidateScores->targetKey = m_mzTargetKey;
    candidateScores->scanNumber = m_msFrameTandem->scanNumberFromFrameIndex(bestAlignmentMatrixRowIndex + peakIntegrationIndexes.first);
    candidateScores->scanTime = m_msFrameTandem->scanTimeFromScanNumber(candidateScores->scanNumber);

    const int arraySizeMax = 12;

    e = ErrorUtils::isTrue(d_ptr->m_mzMeanValsFound.size() <= arraySizeMax); ree;
    for (int i = 0; i < d_ptr->m_mzMeanValsFound.size(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::MzFoundMean1 + i] = d_ptr->m_mzMeanValsFound.at(i);
    }

    e = ErrorUtils::isTrue(d_ptr->m_stdMeanValsFound.size() <= arraySizeMax); ree;
    for (int i = 0; i < d_ptr->m_stdMeanValsFound.size(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::MzFoundStDev1 + i] = d_ptr->m_stdMeanValsFound.at(i);
    }

    const double totalIntensity = std::accumulate(intensityFoundMaxVec.begin(), intensityFoundMaxVec.end(), 0.0001);
    candidateScores->featuresArray[CandidateScores::Features::TotalIntensityLog] = std::log(std::max(totalIntensity, std::numeric_limits<double>::min()));

    const double totalArea = d_ptr->m_intensityMatrix100.sum();
    const Eigen::VectorX<float> matrixColumnSums = d_ptr->m_intensityMatrix100.colwise().sum(); //replace this w/ m_intensityMatrix100ApexRow if it doesn't work out
    e = ErrorUtils::isTrue(matrixColumnSums.rows() <= arraySizeMax); ree;
    for (int i = 0; i < matrixColumnSums.rows(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax1 + i] = static_cast<float>(matrixColumnSums.coeff(i) / totalArea);
    }

    e = ErrorUtils::isTrue(ms2IonsTheoretical.size() <= arraySizeMax); ree;
    for (int i = 0; i < ms2IonsTheoretical.size(); i++) {
        const MS2Ion &ms2Ion = ms2IonsTheoretical.at(i);
        candidateScores->featuresArray[CandidateScores::Features::MzSearched1 + i] = static_cast<float>(ms2Ion.mz);
        candidateScores->featuresArray[CandidateScores::Features::TheoIntensity1 + i] = ms2Ion.intensity;
    }

    candidateScores->featuresArray[CandidateScores::Features::Mass] = static_cast<float>(targetDecoyCandidatePair->mass());
    candidateScores->featuresArray[CandidateScores::Features::IRTPredicted] = static_cast<float>(targetDecoyCandidatePair->iRt());
    candidateScores->featuresArray[CandidateScores::Features::TheoFragmentCount] = static_cast<float>(targetDecoyCandidatePair->totalFragmentCount());
    candidateScores->featuresArray[CandidateScores::Features::ScanIonCount] = static_cast<float>(m_msFrameTandem->getScanPointsByScanNumber(candidateScores->scanNumber)->size());
    candidateScores->featuresArray[CandidateScores::Features::AllignedMaxIndexesCount] = static_cast<float>(std::count_if(
            columnApexIndexes.begin(),
            columnApexIndexes.end(),
            [bestAlignmentMatrixRowIndex](int i){return i == bestAlignmentMatrixRowIndex;}
            ));

    Eigen::VectorX<float> bestAnchorColumn;
    int bestAnchorColumnIndex;
    QVector<float> cosineSimToAnchorVec;
    QVector<float> cosineSimShadowsToAnchorVec;
    QVector<float> shadowsIntensityRatioVec;
    e = d_ptr->calculateCandidateAllignmentMetrics(
            &cosineSimToAnchorVec,
            &cosineSimShadowsToAnchorVec,
            &shadowsIntensityRatioVec,
            &bestAnchorColumn,
            &bestAnchorColumnIndex
    ); ree;

//#define CHECK_ALIGNMENT_MATRIX_BY_SEQUENCE
#ifdef CHECK_ALIGNMENT_MATRIX_BY_SEQUENCE
    if (targetDecoyCandidatePair->peptideStringWithMods() == "EAQGNSSAGVEAAEQRPVEDGER" && targetDecoyCandidatePair->charge() == 3) {
        std::cout << peakIntegrationIndexes.first << " " << peakIntegrationIndexes.second;
        std::cout << d_ptr->m_intensityMatrix100 << std::endl;
        for (float c : cosineSimToAnchorVec) {
            std::cout << c << std::endl;
        }
        std::cout << "**** " << bestAnchorColumnIndex << std::endl;
    }
#endif

    e = ErrorUtils::isTrue(cosineSimToAnchorVec.size() <= arraySizeMax); ree;

    for (int i = 0; i < cosineSimToAnchorVec.size(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor1 + i] = cosineSimToAnchorVec.at(i);
    }

    e = ErrorUtils::isTrue(cosineSimShadowsToAnchorVec.size() <= arraySizeMax); ree;
    for (int i = 0; i < cosineSimShadowsToAnchorVec.size(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor1 + i] = cosineSimShadowsToAnchorVec.at(i);
    }

    for (int i = 0; i < ms2IonsTheoretical.size(); i++) {
        const MS2Ion &ms2IonSearched = ms2IonsTheoretical.at(i);

        float mzFound = 0.0;
        if (i < d_ptr->m_mzMeanValsFound.size()) {
            mzFound = d_ptr->m_mzMeanValsFound.at(i);
        }

        float corr = 0.0;
        if (i < cosineSimToAnchorVec.size()) {
            corr = cosineSimToAnchorVec.at(i);
        }

        const float accuracy = corr * (std::abs(mzFound - ms2IonSearched.mz) / ms2IonSearched.mz);
        candidateScores->featuresArray[CandidateScores::Features::MzAccuracy1 + i] = accuracy;
    }

    e = ErrorUtils::isTrue(shadowsIntensityRatioVec.size() <= arraySizeMax); ree;
    for (int i = 0; i < shadowsIntensityRatioVec.size(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio1 + i] = shadowsIntensityRatioVec.at(i);
    }

    if (MathUtils::tZero(d_ptr->m_intensityMatrix100.maxCoeff()) || bestAnchorColumnIndex < 0) {
        ERR_RETURN
    }

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

    const float scanTimeDelta = std::abs(candidateScores->scanTime - candidateScores->scanTimePredicted);

    candidateScores->featuresArray[CandidateScores::Features::ScanTimeDelta] = scanTimeDelta;
    candidateScores->featuresArray[CandidateScores::Features::Charge] = static_cast<float>(candidateScores->targetDecoyCandidatePair->charge());

    const float scanTimeRange = std::max(
            std::numeric_limits<float>::min(),
            scanTimeMinMax.second - scanTimeMinMax.first
    );
    candidateScores->featuresArray[CandidateScores::Features::ScanTimeRange] = scanTimeRange;

    candidateScores->featuresArray[CandidateScores::Features::ScanTimePredicted] = candidateScores->scanTimePredicted;

    candidateScores->featuresArray[CandidateScores::Features::ShadowsCosineSimSum] = static_cast<float>(std::accumulate(
            cosineSimShadowsToAnchorVec.begin(),
            cosineSimShadowsToAnchorVec.end(),
            0.0
            ));

    const double pdScanTime = std::sqrt(std::min(std::abs(scanTimeDelta), scanTimeRange) / scanTimeRange);
    candidateScores->featuresArray[CandidateScores::Features::ScanTimePd] = static_cast<float>(pdScanTime);

    const double pepLength = (-10.0 + candidateScores->targetDecoyCandidatePair->peptideString().size()) / 10.0;
    candidateScores->featuresArray[CandidateScores::Features::PeptideLengthNorm] = static_cast<float>(pepLength);

    candidateScores->featuresArray[CandidateScores::Features::TheoFragmentCount]
                                            = static_cast<float>(candidateScores->targetDecoyCandidatePair->totalFragmentCount());

    if (bestAnchorColumn.rows() < 1) {
        ERR_RETURN
    }

    double cosineSimSpectrum;
    double klDivSpectrum;
    e = calculateSpectrumMetrics(
            d_ptr->m_intensityMatrix100ApexRow,
            ms2IonsTheoretical,
            &cosineSimSpectrum,
            &klDivSpectrum
    ); ree;

    float cosineSimOverPeak;
    e = calculateSpectrumMetricsOverPeak(
            d_ptr->m_intensityMatrix100,
            ms2IonsTheoretical,
            bestAnchorColumnIndex,
            &cosineSimOverPeak
            ); ree;

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrum] = static_cast<float>(cosineSimSpectrum);
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumOverTime] = static_cast<float>(cosineSimOverPeak);
    candidateScores->featuresArray[CandidateScores::Features::KlDivSpectrum] = static_cast<float>(klDivSpectrum);
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumCubed] = static_cast<float>(std::pow(std::max(0.0, cosineSimSpectrum), 3));
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumOverTimeCubed] = static_cast<float>(std::pow(std::max(0.0f, cosineSimOverPeak), 3));
    candidateScores->featuresArray[CandidateScores::Features::KlDivSpectrumCubeRoot] = static_cast<float>(std::pow(std::max(0.0, klDivSpectrum), 1 / 3.0));

    e = ErrorUtils::isTrue(m_turboXICMS1->isInit()); ree;

    float cosineSimMS1;
    float ms1MzMeanFound100;
    float ms1StDevMzFound100;
    float ms1IntensityFound100;
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            targetDecoyCandidatePair->mz(),
            static_cast<float>(d_ptr->m_pythiaParams.ms2ExtractionWidthPPM),
            m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &cosineSimMS1,
            &ms1MzMeanFound100,
            &ms1StDevMzFound100,
            &ms1IntensityFound100
            ); ree;

    candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1] = std::max(static_cast<float>(cosineSimMS1), std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound100] = std::max(ms1MzMeanFound100, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound100PPM]
            = MathUtils::calculateMassAccuracyPPM(targetDecoyCandidatePair->mz(), ms1MzMeanFound100);
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound100] = std::max(ms1StDevMzFound100, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound100] = std::max(ms1IntensityFound100, std::numeric_limits<float>::min());
    
    float cosineSim45MS1;
    float ms1MzMeanFound45;
    float ms1StDevMzFound45;
    float ms1IntensityFound45;
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            targetDecoyCandidatePair->mz(),
            static_cast<float>(d_ptr->m_pythiaParams.ms2ExtractionWidthPPM) * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
            m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &cosineSim45MS1,
            &ms1MzMeanFound45,
            &ms1StDevMzFound45,
            &ms1IntensityFound45
    ); ree;
    candidateScores->featuresArray[CandidateScores::Features::CosineSim45MS1] = std::max(static_cast<float>(cosineSim45MS1), std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound45] = std::max(ms1MzMeanFound45, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound45PPM]
            = MathUtils::calculateMassAccuracyPPM(targetDecoyCandidatePair->mz(), ms1MzMeanFound45);
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound45] = std::max(ms1StDevMzFound45, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound45] = std::max(ms1IntensityFound45, std::numeric_limits<float>::min());

    float cosineSim20MS1;
    float ms1MzMeanFound20;
    float ms1StDevMzFound20;
    float ms1IntensityFound20;
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            targetDecoyCandidatePair->mz(),
            static_cast<float>(d_ptr->m_pythiaParams.ms2ExtractionWidthPPM) * S_GLOBAL_SETTINGS.TIGHT_2_FRACTION,
            m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &cosineSim20MS1,
            &ms1MzMeanFound20,
            &ms1StDevMzFound20,
            &ms1IntensityFound20
    ); ree;
    candidateScores->featuresArray[CandidateScores::Features::CosineSim20MS1] = std::max(static_cast<float>(cosineSim20MS1), std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound20] = std::max(ms1MzMeanFound20, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFound20PPM]
            = MathUtils::calculateMassAccuracyPPM(targetDecoyCandidatePair->mz(), ms1MzMeanFound20);
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFound20] = std::max(ms1StDevMzFound20, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFound20] = std::max(ms1IntensityFound20, std::numeric_limits<float>::min());

    const float mzPreMono = BiophysicalCalcs::calculateThomsonFromMass(
            targetDecoyCandidatePair->mass(),
            targetDecoyCandidatePair->charge(),
            -1
    );
    float mzPreMonoCosineSimMS1;
    float ms1MzMeanFoundPreMono;
    float ms1StDevMzFoundPreMono;
    float ms1IntensityFoundPreMono;
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            mzPreMono,
            d_ptr->m_pythiaParams.ms2ExtractionWidthPPM,
            m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &mzPreMonoCosineSimMS1,
            &ms1MzMeanFoundPreMono,
            &ms1StDevMzFoundPreMono,
            &ms1IntensityFoundPreMono
    ); ree;
    candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1PreMono]
                    = std::max(static_cast<float>(mzPreMonoCosineSimMS1), std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMono]
                    = std::max(ms1MzMeanFoundPreMono, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundPreMonoPPM]
            = MathUtils::calculateMassAccuracyPPM(mzPreMono, ms1MzMeanFoundPreMono);
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundPreMono] = std::max(ms1StDevMzFoundPreMono, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundPreMono] = std::max(ms1IntensityFoundPreMono, std::numeric_limits<float>::min());

    const float mzIso1 = BiophysicalCalcs::calculateThomsonFromMass(
            targetDecoyCandidatePair->mass(),
            targetDecoyCandidatePair->charge(),
            1
    );
    float cosineSimMzIso1;
    float ms1MzMeanFoundIso1;
    float ms1StDevMzFoundIso1;
    float ms1IntensityFoundIso1;
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            mzIso1,
            d_ptr->m_pythiaParams.ms2ExtractionWidthPPM,
            m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &cosineSimMzIso1,
            &ms1MzMeanFoundIso1,
            &ms1StDevMzFoundIso1,
            &ms1IntensityFoundIso1
    ); ree;
    candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso1]
                = std::max(static_cast<float>(cosineSimMzIso1), std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1]
                = std::max(ms1MzMeanFoundIso1, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso1PPM]
            = MathUtils::calculateMassAccuracyPPM(mzIso1, ms1MzMeanFoundIso1);
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso1] = std::max(ms1StDevMzFoundIso1, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundIso1] = std::max(ms1IntensityFoundIso1, std::numeric_limits<float>::min());

    const float mzIso2 = BiophysicalCalcs::calculateThomsonFromMass(
            targetDecoyCandidatePair->mass(),
            targetDecoyCandidatePair->charge(),
            2
    );
    float cosineSimMzIso2;
    float ms1MzMeanFoundIso2;
    float ms1StDevMzFoundIso2;
    float ms1IntensityFoundIso2;
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            mzIso2,
            d_ptr->m_pythiaParams.ms2ExtractionWidthPPM,
            m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &cosineSimMzIso2,
            &ms1MzMeanFoundIso2,
            &ms1StDevMzFoundIso2,
            &ms1IntensityFoundIso2
    ); ree;
    candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso2]
                    = std::max(static_cast<float>(cosineSimMzIso2), std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso2]
            = std::max(ms1MzMeanFoundIso2, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzMeanFoundIso2PPM]
            = MathUtils::calculateMassAccuracyPPM(mzIso2, ms1MzMeanFoundIso2);
    candidateScores->featuresArray[CandidateScores::Features::Ms1MzStDevFoundIso2] = std::max(ms1StDevMzFoundIso2, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::Ms1IntensityFoundIso2] = std::max(ms1IntensityFoundIso2, std::numeric_limits<float>::min());

    float cosineSimSum45;
    float cosineSimSum20;
    e = d_ptr->calculateCosineSimSumTight(
            bestAnchorColumnIndex,
            &cosineSimSum45,
            &cosineSimSum20
    ); ree;
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum45] = std::max(cosineSimSum45, std::numeric_limits<float>::min());
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum20] = std::max(cosineSimSum20, std::numeric_limits<float>::min());

    const double chunkDivision = 3.0;
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
                0.0
        ) / bestAnchorColumnVecSum;

    }

    //TODO revisit if you want to use all 12 and have a separate one for top 6
    const int top6 = std::min(6, cosineSimToAnchorVec.size());
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100] = std::accumulate(
            cosineSimToAnchorVec.begin(),
            cosineSimToAnchorVec.begin() + top6,
            std::numeric_limits<float>::min()
            );

    const float cosineSimilarityMin = 0.8;
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum100GreaterThan80] = std::accumulate(
            cosineSimToAnchorVec.begin(),
            cosineSimToAnchorVec.begin() + top6,
            std::numeric_limits<float>::min(),
            [cosineSimilarityMin](float sum, float f){return f >= cosineSimilarityMin ? sum + f : sum;}
    );

    const float cosineSimSumTop6 = std::accumulate(
            cosineSimToAnchorVec.begin(),
            cosineSimToAnchorVec.begin() + top6,
            0.0
    );
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSumTop6] = cosineSimSumTop6;

    float cosineSimSumBottom6 = 0.1;
    if (cosineSimToAnchorVec.size() > top6) {
        cosineSimSumBottom6 = std::accumulate(
                cosineSimToAnchorVec.begin() + top6 + 1,
                cosineSimToAnchorVec.end(),
                std::numeric_limits<float>::min()
        );
    }

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSumBottom6] = cosineSimSumBottom6;

    candidateScores->featuresArray[CandidateScores::Features::TopBottomRatio]
            = std::log(std::max(1.0f, cosineSimSumTop6) / (cosineSimSumTop6 + cosineSimSumBottom6 + 1.0f));

    candidateScores->featuresArray[CandidateScores::Features::TopBottomRatioNorm]
            = cosineSimSumBottom6 / static_cast<float>(candidateScores->targetDecoyCandidatePair->totalFragmentCount());

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

    QString peptideString = candidateScores->targetDecoyCandidatePair->peptideString();

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

    const auto mz = static_cast<float>(BiophysicalCalcs::calculateThomsonFromMass(
            candidateScores->targetDecoyCandidatePair->mass(),
            candidateScores->targetDecoyCandidatePair->charge())
    );

    candidateScores->featuresArray[CandidateScores::Features::MzNorm] = (mz - 600.0f) * 0.002f;

    ERR_RETURN
}
