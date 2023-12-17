//
// Created by anichols on 10/20/23.
//

#include "ScoreOverseer.h"

#include "BiophysicalCalcs.h"
#include "CandidateScores.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FeatureFinderHill.h"
#include "FeatureFinderHillBuilder.h"
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
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QVector<MS2Ion> &ms2IonsTheoreticalShadows,
            const MS2Ion &ms2IonUnfragPrecursor,
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsShadows,
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &unfragPrecursorVsfeatureFinderHills,
            int *bestAlignmentMatrixRow,
            QVector<int> *columnApexIndexes,
            QVector<float> *intensityFoundMaxVec,
            QVector<int> *columnPeakLengths,
            QPair<FrameIndex, FrameIndex> *minMaxFrameIndex,
            QPair<int, int> *alignmentMatrixLimits
    );

    Err calculateCandidateAllignementMetrics(
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
    Eigen::MatrixX<float> m_unfragPrecursorVec;

    Eigen::VectorX<float> m_intensityMatrix100ApexRow;

    QVector<double> m_mzMeanValsFound;
    QVector<double> m_stdMeanValsFound;

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
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            int frameIndexMin,
            Eigen::MatrixX<float> *mat,
            QVector<double> *mzMeanValsFound,
            QVector<double> *stdMeanValsFound
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(frameIndexMin >= 0); ree;
        e = ErrorUtils::isTrue(mat->rows() > 0); ree;

        mat->setZero();

        const int rowCount = static_cast<int>(mat->rows());

        mzMeanValsFound->reserve(ms2IonsTheoretical.size());
        stdMeanValsFound->reserve(ms2IonsTheoretical.size());

        for (int col = 0; col < ms2IonsTheoretical.size(); col++) {

            if (col >= mat->cols()) {
                break;
            }

            const MS2Ion &ms2Ion = ms2IonsTheoretical.at(col);

            QVector<double> mzMeanValsFoundMS2Ion;
            QVector<double> stdMeanValsFoundMS2Ion;
            mzMeanValsFoundMS2Ion.reserve(ms2IonsTheoretical.size());
            stdMeanValsFoundMS2Ion.reserve(ms2IonsTheoretical.size());

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            QVector<float> eigenVectorLoader(rowCount, 0.0f);
            eigenVectorLoader.reserve(rowCount);

            const QVector<FeatureFinderHill*> &ffhs = mzHashedVsfeatureFinderHills.value(mzHashed);
            for (FeatureFinderHill *ffh : ffhs) {

                mzMeanValsFoundMS2Ion.push_back(ffh->mzMean());
                stdMeanValsFoundMS2Ion.push_back(ffh->mzStDev());

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

            mzMeanValsFound->push_back(MathUtils::mean(mzMeanValsFoundMS2Ion));
            stdMeanValsFound->push_back(MathUtils::mean(stdMeanValsFoundMS2Ion));
            mat->col(col) = EigenUtils::convertQVectorToEigenVector(eigenVectorLoader);
        }

        ERR_RETURN
    }

    Err loadMzHashedVsFeatureFinderHillsToMatrixTight(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            int frameIndexMin,
            double ppmTol,
            const QPair<int, int> &alignmentMatrixLimits,
            Eigen::MatrixX<float> *matTight1,
            Eigen::MatrixX<float> *matTight2
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(frameIndexMin >= 0); ree;
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

            const QVector<FeatureFinderHill*> &ffhs = mzHashedVsfeatureFinderHills.value(mzHashed);
            for (FeatureFinderHill *ffh : ffhs) {

                const QVector<FrameIndex> frameIndexes = ffh->scanNumberIndexes();
                const QVector<float> intensities = ffh->intensities();
                const QVector<double> mzVals = ffh->mzVals();

                e = ErrorUtils::isEqual(frameIndexes.size(), intensities.size()); ree;

                for (int i = 0; i < frameIndexes.size(); i++) {

                    const int frameIndexAdjusted = frameIndexes.at(i) - frameIndexMin;
                    if (frameIndexAdjusted < 0 || frameIndexAdjusted >= rowCount) {
                        continue;
                    }

                    const double mzTolTight1 = MathUtils::calculatePPM(ms2Ion.mz, S_GLOBAL_SETTINGS.TIGHT_1_FRACTION * ppmTol);
                    const double mzMinTight1 = ms2Ion.mz - mzTolTight1;
                    const double mzMaxTight1 = ms2Ion.mz + mzTolTight1;

                    const double mz = mzVals.at(i);

                    if (mz < mzMinTight1 || mz > mzMaxTight1) {
                        continue;
                    }

                    eigenVectorLoaderTight1[frameIndexAdjusted] = intensities.at(i);

                    const double mzTolTight2 = MathUtils::calculatePPM(ms2Ion.mz, S_GLOBAL_SETTINGS.TIGHT_2_FRACTION * ppmTol);
                    const double mzMinTight2 = ms2Ion.mz - mzTolTight2;
                    const double mzMaxTight2 = ms2Ion.mz + mzTolTight2;

                    if (mz < mzMinTight2 || mz > mzMaxTight2) {
                        continue;
                    }

                    eigenVectorLoaderTight2[frameIndexAdjusted] = intensities.at(i);

                }
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
        Eigen::VectorX<float> rowCount = rowCountMat.rowwise().sum();

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
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QVector<MS2Ion> &ms2IonsTheoreticalShadows,
        const MS2Ion &ms2IonUnfragPrecursor,
        const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsShadows,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &unfragPrecursorVsfeatureFinderHills,
        int *bestAlignmentMatrixRowIndex,
        QVector<int> *columnApexIndexes,
        QVector<float> *intensityFoundMaxVec,
        QVector<int> *columnPeakLengths,
        QPair<FrameIndex, FrameIndex> *minMaxFrameIndex,
        QPair<int, int> *alignmentMatrixLimits
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;
    e = ErrorUtils::isNotEmpty(ms2IonsTheoreticalShadows); ree;
    e = ErrorUtils::isEqual(ms2IonsTheoretical.size(), ms2IonsTheoreticalShadows.size()); ree;

    m_mzMeanValsFound.clear();
    m_stdMeanValsFound.clear();

    const int cols = ms2IonsTheoretical.size();

    *minMaxFrameIndex = ScoreOverseer::getMinMaxFrameIndexes(mzHashedVsfeatureFinderHills);

    e = ErrorUtils::isTrue(minMaxFrameIndex->second >= minMaxFrameIndex->first); ree;
    const int rows = std::max(minMaxFrameIndex->second - minMaxFrameIndex->first + 1, m_pythiaParams.filterLength);

    m_intensityMatrix100.resize(rows, cols);
    e = loadMzHashedVsFeatureFinderHillsToMatrix(
            ms2IonsTheoretical,
            mzHashedVsfeatureFinderHills,
            minMaxFrameIndex->first,
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

    m_intensityMatrix100Shadow.resize(rows, cols);
    QVector<double> unused1;
    QVector<double> unused2;
    e = loadMzHashedVsFeatureFinderHillsToMatrix(
            ms2IonsTheoreticalShadows,
            mzHashedVsfeatureFinderHillsShadows,
            minMaxFrameIndex->first,
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
    }

    if (m_intensityMatrix100Shadow.maxCoeff() > 0.0f) {

        m_intensityMatrix100Shadow = applyGaussSmoothRowWiseToMatrix(
                m_intensityMatrix100Shadow,
                m_pythiaParams,
                m_gaussKernel
        );

    }

    m_unfragPrecursorVec.resize(rows, 1);
    e = loadMzHashedVsFeatureFinderHillsToMatrix(
            {ms2IonUnfragPrecursor},
            unfragPrecursorVsfeatureFinderHills,
            minMaxFrameIndex->first,
            &m_unfragPrecursorVec,
            &unused1,
            &unused2
    ); ree;

    if (rowSize < m_unfragPrecursorVec.rows()) {
        m_unfragPrecursorVec = m_unfragPrecursorVec.block(
                alignmentMatrixLimits->first,
                0,
                rowSize,
                m_unfragPrecursorVec.cols()
        ).eval();
    }

    if (m_unfragPrecursorVec.maxCoeff() > 0.0f) {

        m_unfragPrecursorVec = applyGaussSmoothRowWiseToMatrix(
                m_unfragPrecursorVec,
                m_pythiaParams,
                m_gaussKernel
        );

    }

    const int tightColsMax = 6;
    m_intensityMatrix45.resize(rows, tightColsMax);
    m_intensityMatrix20.resize(rows, tightColsMax);
    e = loadMzHashedVsFeatureFinderHillsToMatrixTight(
            ms2IonsTheoretical,
            mzHashedVsfeatureFinderHills,
            minMaxFrameIndex->first,
            m_pythiaParams.ms2ExtractionWidthPPM,
            *alignmentMatrixLimits,
            &m_intensityMatrix45,
            &m_intensityMatrix20
    );
    ree;


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

                const float cosineSimToAnchor = EigenUtils::cosineSimilarity(anchorColumn, altColumn);

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
Err ScoreOverseer::Private::calculateCandidateAllignementMetrics(
        QVector<float> *cosineSimsIndividual,
        QVector<float> *cosineSimsShadowIndividual,
        QVector<float> *shadowsIntensityRatioVec,
        Eigen::VectorX<float> *bestAnchorColumn,
        int *bestAnchorColumnIndex
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_intensityMatrix100.nonZeros() > 1); ree;

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

            const float cosineSimToShadow = EigenUtils::cosineSimilarity(mzPeak, mzPeakShadow);
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
        QVector<float> v(static_cast<int>(m_intensityMatrix100.cols()) , 0.0f);
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
    EigenUtils::replaceNaN(0.0f , &peakMaskMat);

    e = ErrorUtils::isEqual(peakMaskMat.cols(), intensityMatrix45.cols()); ree;
    e = ErrorUtils::isEqual(peakMaskMat.cols(), intensityMatrix20.cols()); ree;

    intensityMatrix45 = intensityMatrix45.array() * peakMaskMat.array();
    intensityMatrix20 = intensityMatrix20.array() * peakMaskMat.array();

    intensityMatrix45 = applyGaussSmoothRowWiseToMatrix(intensityMatrix45, m_pythiaParams, m_gaussKernel);
    intensityMatrix20 = applyGaussSmoothRowWiseToMatrix(intensityMatrix20, m_pythiaParams, m_gaussKernel);

    const Eigen::VectorX<float> anchorCol45 = intensityMatrix45.col(anchorColumnIndex);
    const Eigen::VectorX<float> anchorCol20 = intensityMatrix20.col(anchorColumnIndex);

    for (int col = 0; col < intensityMatrix45.cols(); col++) {
        const Eigen::VectorX<float> col45 = intensityMatrix45.col(col);
        const Eigen::VectorX<float> col20 = intensityMatrix20.col(col);
        const float cosineSim45 = EigenUtils::cosineSimilarity(col45, anchorCol45);
        const float cosineSim20 = EigenUtils::cosineSimilarity(col45, anchorCol20);
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
        MsFrame *ms1Frame
) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = d_ptr->init(pythiaParameters); ree;
    e = ErrorUtils::isTrue(ms1Frame->isValid()); ree;

    m_ms1Frame = ms1Frame;
    m_mzTargetKey = targetKey;

    e = m_turboXICMS1.init(m_ms1Frame->frameIndexVsScanPoints()); ree;

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

        *cosineSimSpectrum = EigenUtils::cosineSimilarity(intensityApexVector, intensityValsTheo);
        *klDivSpectrum = EigenUtils::klDivergence(intensityApexVector, intensityValsTheo);

        ERR_RETURN
    }

    Err calculateMS1Corr(
            const Eigen::VectorX<float> &bestAnchorColumn,
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            double mzTarget,
            double ppmTol,
            TurboXIC *turboXic,
            Eigen::VectorX<float> *gaussKernel,
            double *cosineSimMS1
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(bestAnchorColumn.size() > 0); ree;

        const double massTol = MathUtils::calculatePPM(mzTarget, ppmTol);
        const double mzStart = mzTarget - massTol;
        const double mzEnd = mzTarget + massTol;

        const XICPoints xicPoints = turboXic->extractPointsXIC(
                mzStart,
                mzEnd,
                peakIntegrationIndexes.first,
                peakIntegrationIndexes.second
                );

        Eigen::VectorX<float> ms1Vec(static_cast<int>(bestAnchorColumn.size()));
        ms1Vec.setZero();

        const QMap<ScanNumber, double> &scanNumbersVsIntensityVals = xicPoints.scanNumbersVsIntensity;

        for (auto it = scanNumbersVsIntensityVals.begin(); it != scanNumbersVsIntensityVals.end(); it++) {
            const FrameIndex frameIndex = it.key() - peakIntegrationIndexes.first;

            if (frameIndex >= ms1Vec.size()) {
                continue;
            }

            ms1Vec.coeffRef(frameIndex) = static_cast<float>(it.value());
        }

        ms1Vec = EigenKernelUtils::convolveVectorWithKernel(ms1Vec, *gaussKernel);

        *cosineSimMS1 = EigenUtils::cosineSimilarity(bestAnchorColumn, ms1Vec);

        ERR_RETURN
    }

}//namespace
Err ScoreOverseer::buildScores(
        TargetDecoyCandidatePair* targetDecoyCandidatePair,
        const QPair<double, double> &scanTimeMinMax,
        const PeakIntegrationIndexes &peakIntegrationIndexes,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsShadows,
        const MS2Ion &ms2IonUnfragPrecursor,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &unfragPrecursorVsfeatureFinderHills,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;
    e = ErrorUtils::isNotEmpty(mzHashedVsfeatureFinderHills); ree;
    e = ErrorUtils::isNotEmpty(ms2IonsTheoreticalIsotopeShadows); ree;

    QElapsedTimer et;
    et.start();

    candidateScores->targetDecoyCandidatePair = targetDecoyCandidatePair;
    candidateScores->initFeaturesArray();

    int bestAlignmentMatrixRowIndex;
    QVector<int> columnApexIndexes;
    QVector<float> intensityFoundMaxVec;
    QVector<int> mzPeakLengthsVec;
    QPair<FrameIndex, FrameIndex> frameIndexMinMaxHills;
    QPair<int, int> alignmentMatrixLimits;
    e = d_ptr->buildAlignmentMatricies(
            ms2IonsTheoretical,
            ms2IonsTheoreticalIsotopeShadows,
            ms2IonUnfragPrecursor,
            mzHashedVsfeatureFinderHills,
            mzHashedVsfeatureFinderHillsShadows,
            unfragPrecursorVsfeatureFinderHills,
            &bestAlignmentMatrixRowIndex,
            &columnApexIndexes,
            &intensityFoundMaxVec,
            &mzPeakLengthsVec,
            &frameIndexMinMaxHills,
            &alignmentMatrixLimits
            );

    candidateScores->targetKey = m_mzTargetKey;
    candidateScores->scanNumber = m_ms1Frame->scanNumberFromFrameIndex(bestAlignmentMatrixRowIndex + frameIndexMinMaxHills.first);
    candidateScores->scanTime = m_ms1Frame->scanTimeFromScanNumber(candidateScores->scanNumber);

    //Figure out if you want to use this as well.
//    candidateScores->mzFoundStDevVec = d_ptr->m_stdMeanValsFound;

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
    const Eigen::VectorX<float> matrixColumnSums = d_ptr->m_intensityMatrix100.rowwise().sum(); //replace this w/ m_intensityMatrix100ApexRow if it doesn't work out
    e = ErrorUtils::isTrue(matrixColumnSums.cols() <= arraySizeMax); ree;
    for (int i = 0; i < matrixColumnSums.cols(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::IntensityFoundMax1 + i] = matrixColumnSums.coeff(i) / totalArea;
    }

    e = ErrorUtils::isTrue(ms2IonsTheoretical.size() <= arraySizeMax); ree;
    for (int i = 0; i < ms2IonsTheoretical.size(); i++) {
        const MS2Ion &ms2Ion = ms2IonsTheoretical.at(i);
        candidateScores->featuresArray[CandidateScores::Features::MzSearched1 + i] = ms2Ion.mz;
        candidateScores->featuresArray[CandidateScores::Features::TheoIntensity1 + i] = ms2Ion.intensity;
    }

    candidateScores->featuresArray[CandidateScores::Features::Mass] = targetDecoyCandidatePair->mass();
    candidateScores->featuresArray[CandidateScores::Features::IRTPredicted] = targetDecoyCandidatePair->iRt();
    candidateScores->featuresArray[CandidateScores::Features::TheoFragmentCount] = targetDecoyCandidatePair->totalFragmentCount();
    candidateScores->featuresArray[CandidateScores::Features::ScanIonCount] = m_ms1Frame->getScanPointsByScanNumber(candidateScores->scanNumber)->size();
    candidateScores->featuresArray[CandidateScores::Features::AllignedMaxIndexesCount] = static_cast<int>(std::count_if(
            columnApexIndexes.begin(),
            columnApexIndexes.end(),
            [bestAlignmentMatrixRowIndex](int i){return i == bestAlignmentMatrixRowIndex;}
            ));

    Eigen::VectorX<float> bestAnchorColumn;
    int bestAnchorColumnIndex;
    QVector<float> cosineSimToAnchorVec;
    QVector<float> cosineSimShadowsToAnchorVec;
    QVector<float> shadowsIntensityRatioVec;
    e = d_ptr->calculateCandidateAllignementMetrics(
            &cosineSimToAnchorVec,
            &cosineSimShadowsToAnchorVec,
            &shadowsIntensityRatioVec,
            &bestAnchorColumn,
            &bestAnchorColumnIndex
    ); ree;

    e = ErrorUtils::isTrue(cosineSimToAnchorVec.size() <= arraySizeMax); ree;
    for (int i = 0; i < cosineSimToAnchorVec.size(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::CosineSimToAnchor1 + i] = cosineSimToAnchorVec.at(i);
    }

    e = ErrorUtils::isTrue(cosineSimShadowsToAnchorVec.size() <= arraySizeMax); ree;
    for (int i = 0; i < cosineSimShadowsToAnchorVec.size(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::CosineSimShadowsToAnchor1 + i] = cosineSimShadowsToAnchorVec.at(i);
    }

    e = ErrorUtils::isTrue(shadowsIntensityRatioVec.size() <= arraySizeMax); ree;
    for (int i = 0; i < shadowsIntensityRatioVec.size(); i++) {
        candidateScores->featuresArray[CandidateScores::Features::ShadowsIntensityRatio1 + i] = shadowsIntensityRatioVec.at(i);
    }

    if (MathUtils::tZero(d_ptr->m_intensityMatrix100.maxCoeff()) || bestAnchorColumnIndex < 0) {
        ERR_RETURN
    }

    if (!MathUtils::tZero(d_ptr->m_unfragPrecursorVec.maxCoeff())) {

        const Eigen::VectorX<float> unfragPrecursorEigVed(d_ptr->m_unfragPrecursorVec);

        e = ErrorUtils::isEqual(bestAnchorColumn.size(), unfragPrecursorEigVed.rows()); ree;
        candidateScores->featuresArray[CandidateScores::Features::UnfragPrecursorCosineSim] = EigenUtils::cosineSimilarity(
                bestAnchorColumn,
                unfragPrecursorEigVed
                        );
    }

    const double columnApexIndexesMean = MathUtils::mean(columnApexIndexes);
    const double columnApexIndexesSize = columnApexIndexes.size();
    QVector<int> columnApexIndexRatiosToAnchor;
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
    QVector<double> mzPeakLengthsNormalized;
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

    const double scanTimeDelta = std::abs(candidateScores->scanTime - candidateScores->scanTimePredicted);
    candidateScores->featuresArray[CandidateScores::Features::ScanTimeDelta] = scanTimeDelta;
    candidateScores->featuresArray[CandidateScores::Features::ChargeNorm] = -2.0 + candidateScores->targetDecoyCandidatePair->charge();

    const double scanTimeRange = std::max(
            std::numeric_limits<double>::min(),
            scanTimeMinMax.second - scanTimeMinMax.first
    );
    candidateScores->featuresArray[CandidateScores::Features::ScanTimeRange] = scanTimeRange;

    candidateScores->featuresArray[CandidateScores::Features::ShadowsCosineSimSum] = std::accumulate(
            cosineSimShadowsToAnchorVec.begin(),
            cosineSimShadowsToAnchorVec.end(),
            0.0
            );

    const double pdScanTime = std::sqrt(std::min(std::abs(scanTimeDelta), scanTimeRange) / scanTimeRange);
    candidateScores->featuresArray[CandidateScores::Features::ScanTimePd] = pdScanTime;

    const double pepLength = (-10.0 + candidateScores->targetDecoyCandidatePair->peptideStringWithMods().size()) / 10.0;
    candidateScores->featuresArray[CandidateScores::Features::PeptideLengthNorm] = pepLength;

    candidateScores->featuresArray[CandidateScores::Features::TheoFragmentCount]
                                            = candidateScores->targetDecoyCandidatePair->totalFragmentCount();


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

    candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrum] = cosineSimSpectrum;
    candidateScores->featuresArray[CandidateScores::Features::KlDivSpectrum] = klDivSpectrum;
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSpectrumCubed] = std::pow(std::max(0.0, cosineSimSpectrum), 3);
    candidateScores->featuresArray[CandidateScores::Features::KlDivSpectrumCubeRoot] = std::pow(std::max(0.0, klDivSpectrum), 1 / 3.0);

    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            targetDecoyCandidatePair->mz(),
            d_ptr->m_pythiaParams.ms2ExtractionWidthPPM,
            &m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1]
            ); ree;



    const double mzPreMono = BiophysicalCalcs::calculateThomsonFromMass(
            targetDecoyCandidatePair->mass(),
            targetDecoyCandidatePair->charge(),
            -1
    );
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            mzPreMono,
            d_ptr->m_pythiaParams.ms2ExtractionWidthPPM,
            &m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1PreMono]
    ); ree;


    const double mzIso1 = BiophysicalCalcs::calculateThomsonFromMass(
            targetDecoyCandidatePair->mass(),
            targetDecoyCandidatePair->charge(),
            1
    );
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            mzIso1,
            d_ptr->m_pythiaParams.ms2ExtractionWidthPPM,
            &m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso1]
    ); ree;

    const double mzIso2 = BiophysicalCalcs::calculateThomsonFromMass(
            targetDecoyCandidatePair->mass(),
            targetDecoyCandidatePair->charge(),
            2
    );
    e = calculateMS1Corr(
            bestAnchorColumn,
            peakIntegrationIndexes,
            mzIso2,
            d_ptr->m_pythiaParams.ms2ExtractionWidthPPM,
            &m_turboXICMS1,
            &d_ptr->m_gaussKernel,
            &candidateScores->featuresArray[CandidateScores::Features::CosineSim100MS1Iso2]
    ); ree;

    float cosineSimSum45;
    float cosineSimSum20;
    e = d_ptr->calculateCosineSimSumTight(
            bestAnchorColumnIndex,
            &cosineSimSum45,
            &cosineSimSum20
    ); ree;
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum45] = cosineSimSum45;
    candidateScores->featuresArray[CandidateScores::Features::CosineSimSum20] = cosineSimSum20;

    const double chunkDivision = 3.0;
    QVector<float> bestAnchorColumnVec = EigenUtils::convertEigenVectorToQVector(bestAnchorColumn);
    const auto terminatorLogic = [](double d){return d < 1.0;};
    const auto terminator = std::remove_if(bestAnchorColumnVec.begin(), bestAnchorColumnVec.end(), terminatorLogic);
    bestAnchorColumnVec.erase(terminator, bestAnchorColumnVec.end());
    const int chunkSize = std::max(1, static_cast<int>(std::round(bestAnchorColumnVec.size() / chunkDivision)));
    const double bestAnchorColumnVecSum = std::accumulate(bestAnchorColumnVec.begin(), bestAnchorColumnVec.end(), 0.0001);
    if (bestAnchorColumnVec.size() < chunkDivision) {
        candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio1] = 0.0;
        candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio2] = 1.0;
        candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio3] = 0.0;
    }
    else {

        candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio1] = std::accumulate(
                bestAnchorColumnVec.begin(),
                bestAnchorColumnVec.begin() + chunkSize,
                0.0
        ) / bestAnchorColumnVecSum;

        candidateScores->featuresArray[CandidateScores::Features::PeakShapeRatio2] = std::accumulate(
                bestAnchorColumnVec.begin() + chunkSize,
                bestAnchorColumnVec.begin() + (chunkSize * 2),
                0.0
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
            0.0f
            );

    if (cosineSimToAnchorVec.size() > top6) {

        const float cosineSimSumTop6 = std::accumulate(
                cosineSimToAnchorVec.begin(),
                cosineSimToAnchorVec.begin() + top6,
                0.0f
        );
        candidateScores->featuresArray[CandidateScores::Features::CosineSimSumTop6] = cosineSimSumTop6;

        const float cosineSimSumBottom6 = std::accumulate(
                cosineSimToAnchorVec.begin() + top6 + 1,
                cosineSimToAnchorVec.end(),
                0.0f
        );
        candidateScores->featuresArray[CandidateScores::Features::CosineSimSumBottom6] = cosineSimSumBottom6;

        candidateScores->featuresArray[CandidateScores::Features::TopBottomRatio]
                = std::log(std::max(1.0f, cosineSimSumTop6) / (cosineSimSumTop6 + cosineSimSumBottom6 + 1.0));

        candidateScores->featuresArray[CandidateScores::Features::TopBottomRatioNorm]
                = cosineSimSumBottom6 / candidateScores->targetDecoyCandidatePair->totalFragmentCount();
    }

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
            {'Y', 0}
    };

    for (const QChar aminoAcid : candidateScores->targetDecoyCandidatePair->peptideStringWithMods()) {

        if (!aminoAcidCounts.contains(aminoAcid)) {
            qDebug() << candidateScores->targetDecoyCandidatePair->peptideStringWithMods() << "missing amino acid";
        }

        e = ErrorUtils::isTrue(aminoAcidCounts.contains(aminoAcid)); ree;
        aminoAcidCounts[aminoAcid]++;
    }

    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountA] = aminoAcidCounts['A'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountC] = aminoAcidCounts['C'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountD] = aminoAcidCounts['D'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountE] = aminoAcidCounts['E'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountF] = aminoAcidCounts['F'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountG] = aminoAcidCounts['G'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountH] = aminoAcidCounts['H'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountI] = aminoAcidCounts['I'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountK] = aminoAcidCounts['K'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountL] = aminoAcidCounts['L'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountM] = aminoAcidCounts['M'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountN] = aminoAcidCounts['N'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountP] = aminoAcidCounts['P'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountQ] = aminoAcidCounts['Q'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountR] = aminoAcidCounts['R'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountS] = aminoAcidCounts['S'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountT] = aminoAcidCounts['T'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountV] = aminoAcidCounts['V'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountW] = aminoAcidCounts['W'];
    candidateScores->featuresArray[CandidateScores::Features::AminoAcidCountY] = aminoAcidCounts['Y'];

    const double mz = BiophysicalCalcs::calculateThomsonFromMass(
            candidateScores->targetDecoyCandidatePair->mass(),
            candidateScores->targetDecoyCandidatePair->charge()
    );

    candidateScores->featuresArray[CandidateScores::Features::MzNorm] = (mz - 600.0) * 0.002;

    ERR_RETURN
}

QPair<FrameIndex, FrameIndex> ScoreOverseer::getMinMaxFrameIndexes(
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
