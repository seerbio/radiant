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
#include "MsFrame.h"
#include "ParallelUtils.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"


class Q_DECL_HIDDEN ScoreOverseer::Private
{

public:

    Private();
    ~Private();

    Err init(const PythiaParameters &pythiaParameters);

    Err buildAlignmentMatricies(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QVector<MS2Ion> &ms2IonsTheoreticalShadows,
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsShadows
    );

    Eigen::VectorX<float> m_gaussKernel;

    Eigen::MatrixX<float> m_intensityMatrix100;
    Eigen::MatrixX<float> m_intensityMatrix45;
    Eigen::MatrixX<float> m_intensityMatrix20;
    Eigen::MatrixX<float> m_intensityMatrix100Shadow;

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

    Err loadMzHashedVsFeatureFinderHillsToMatrixTight(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
            int frameIndexMin,
            double ppmTol,
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

}// namespace
Err ScoreOverseer::Private::buildAlignmentMatricies(
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QVector<MS2Ion> &ms2IonsTheoreticalShadows,
        const QHash<MzHashed, QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsShadows
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;
    e = ErrorUtils::isNotEmpty(ms2IonsTheoreticalShadows); ree;
    e = ErrorUtils::isEqual(ms2IonsTheoretical.size(), ms2IonsTheoreticalShadows.size()); ree;

    const int cols = ms2IonsTheoretical.size();

    const QPair<FrameIndex, FrameIndex> minMaxFrameIndex
            = ScoreOverseer::getMinMaxFrameIndexes(mzHashedVsfeatureFinderHills);

    e = ErrorUtils::isTrue(minMaxFrameIndex.second >= minMaxFrameIndex.first); ree;
    const int rows = std::max(minMaxFrameIndex.second - minMaxFrameIndex.first + 1, m_pythiaParams.filterLength);

    m_intensityMatrix100.resize(rows, cols);
    e = loadMzHashedVsFeatureFinderHillsToMatrix(
            ms2IonsTheoretical,
            mzHashedVsfeatureFinderHills,
            minMaxFrameIndex.first,
            &m_intensityMatrix100
    ); ree;
    m_intensityMatrix100 = applyGaussSmoothRowWiseToMatrix(
            m_intensityMatrix100,
            m_pythiaParams,
            m_gaussKernel
            );

    m_intensityMatrix100Shadow.resize(rows, cols);
    e = loadMzHashedVsFeatureFinderHillsToMatrix(
            ms2IonsTheoreticalShadows,
            mzHashedVsfeatureFinderHillsShadows,
            minMaxFrameIndex.first,
            &m_intensityMatrix100Shadow
    ); ree;
    m_intensityMatrix100Shadow = applyGaussSmoothRowWiseToMatrix(
            m_intensityMatrix100Shadow,
            m_pythiaParams,
            m_gaussKernel
            );

    const int tightColsMax = 6;
    m_intensityMatrix45.resize(rows, tightColsMax);
    m_intensityMatrix20.resize(rows, tightColsMax);
    e = loadMzHashedVsFeatureFinderHillsToMatrixTight(
            ms2IonsTheoretical,
            mzHashedVsfeatureFinderHills,
            minMaxFrameIndex.first,
            m_pythiaParams.ms2ExtractionWidthPPM,
            &m_intensityMatrix45,
            &m_intensityMatrix20
    ); ree;

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


ScoreOverseer::ScoreOverseer()
: d_ptr(new Private())
{}

ScoreOverseer::~ScoreOverseer() {}

namespace {


    Err findBestApexRowInMatrix(
            const Eigen::MatrixX<double> &intensityMatrixIntegratedLimits,
            const QVector<double> &peakCountsPerRow,
            int *bestApexRowIndex
    ) {

        ERR_INIT

        using RowMax = double;
        using RowSum = double;

        struct SearchObject {
            RowMax rowMax = -1.0;
            RowSum rowSum = -1.0;
            int index = -1;
        };

        const Eigen::VectorX<double> rowSums = intensityMatrixIntegratedLimits.rowwise().sum();
        const QVector<double> rowSumsVec = EigenUtils::convertEigenVectorToQVector(rowSums);

        QVector<QPair<RowMax, RowSum>> rowMaxVsRowSum;
        e = ParallelUtils::zip(peakCountsPerRow, rowSumsVec, &rowMaxVsRowSum); ree;

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
            Eigen::VectorX<double> *vec
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
            const Eigen::MatrixX<double> &_intensityMatrixIntegratedLimits,
            const QVector<double> &peakCountsPerRow,
            Eigen::MatrixX<double> *intensityMatrixIntegratedLimitsNoInterference,
            int *bestApexRowIndex,
            QVector<int> *apexIndexes
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peakCountsPerRow); ree;
        e = ErrorUtils::isTrue(_intensityMatrixIntegratedLimits.maxCoeff() > 0); ree;

        Eigen::MatrixX<double> intensityMatrixIntegratedLimits = _intensityMatrixIntegratedLimits;

        e = findBestApexRowInMatrix(intensityMatrixIntegratedLimits, peakCountsPerRow, bestApexRowIndex); ree;

        for (int col = 0; col < intensityMatrixIntegratedLimits.cols(); col++) {

            Eigen::VectorX<double> colVec = intensityMatrixIntegratedLimits.col(col);
            const QVector<int> columnApexes = EigenUtils::apexesIndexesOnly(colVec);

            if (columnApexes.isEmpty()) {
                continue;
            }

            const int bestApexIndexColumn = columnApexes.at(MathUtils::closest(columnApexes, *bestApexRowIndex));

            e = removeInterferingPeaksInColumn(bestApexIndexColumn, &colVec); ree;
            intensityMatrixIntegratedLimits.col(col) = colVec;
            apexIndexes->push_back(bestApexIndexColumn);

        }

        *intensityMatrixIntegratedLimitsNoInterference = intensityMatrixIntegratedLimits;

        ERR_RETURN
    }

    Err calcBestCosineSimSumLogic(
            const Eigen::MatrixX<double> &intensityMatrixIntegratedLimitsSmoothed,
            double cosineSimToAnchorThreshold,
            QVector<double> *bestCosineSimsIndividual,
            double *bestCosineSimSum,
            QVector<double> *bestIntensitiesIndividual,
            Eigen::VectorX<double> *bestAnchorColumn
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrixIntegratedLimitsSmoothed.rows() > 0); ree;

        *bestCosineSimSum = 0.0;

        const int topCols6OrLess = std::min(6, static_cast<int>(intensityMatrixIntegratedLimitsSmoothed.cols()));

        for (int i = 0; i < topCols6OrLess; i++) {

            QVector<double> bestCosineSimsIndividualAnchor;
            QVector<double> bestIndividualIntensities;
            const Eigen::VectorX<double> &anchorColumn = intensityMatrixIntegratedLimitsSmoothed.col(i);

            for (int j = 0; j < intensityMatrixIntegratedLimitsSmoothed.cols(); j++) {

                const Eigen::VectorX<double> &altColumn = intensityMatrixIntegratedLimitsSmoothed.col(j);

                const double cosineSimToAnchor = EigenUtils::cosineSimilarity(anchorColumn, altColumn);

                if (cosineSimToAnchor < cosineSimToAnchorThreshold) {
                    bestCosineSimsIndividualAnchor.push_back(0.0);
                    bestIndividualIntensities.push_back(0.0);
                    continue;
                }

                const double intensity = altColumn.maxCoeff();
                bestIndividualIntensities.push_back(intensity);
                bestCosineSimsIndividualAnchor.push_back(cosineSimToAnchor);
            }

            const double cosineSimSumAnchor
                    = std::accumulate(bestCosineSimsIndividualAnchor.begin(), bestCosineSimsIndividualAnchor.end(), 0.0);

            if (cosineSimSumAnchor > *bestCosineSimSum) {
                *bestIntensitiesIndividual = bestIndividualIntensities;
                *bestCosineSimSum = cosineSimSumAnchor;
                *bestCosineSimsIndividual = bestCosineSimsIndividualAnchor;
                *bestAnchorColumn = anchorColumn;
            }

        }

        ERR_RETURN
    }

    Err calcBestCosineSimSum(
            const Eigen::MatrixX<double> &intensityMatrix,
            const Eigen::MatrixX<double> &intensityMatrixShadows,
            const PeakIntegrationIndexes &pii,
            const QVector<double> &summedMatToVec,
            int topNMs2FragPeaks,
            double cosineSimToAnchorThreshold,
            bool subtractShadows,
            QVector<double> *cosineSimsIndividualIntegration,
            double *bestCosineSimSumIntegration,
            QVector<double> *bestIntensitiesIndividual,
            Eigen::VectorX<double> *bestAnchorColumn,
            FrameIndex *frameIndexIntensityApexIntegration,
            int *allignedMaxIndexesCount,
            QVector<double> *cosineSimsShadowIndividual,
            QVector<double> *shadowsIntensityRatioPeakIntVec
            ) {

        ERR_INIT

        Eigen::MatrixX<double> intensityMatrixIntegratedLimitsNoInterference;
//        e = buildAlignmentMatrix(
//                intensityMatrix,
//                pii,
//                summedMatToVec,
//                &intensityMatrixIntegratedLimitsNoInterference,
//                frameIndexIntensityApexIntegration,
//                allignedMaxIndexesCount
//                ); ree;
//
//        intensityMatrixIntegratedLimitsNoInterference
//                = applyGaussSmoothRowWiseToMatrix(intensityMatrixIntegratedLimitsNoInterference);

        if (intensityMatrixShadows.cols() > 0) {

            FrameIndex unused;
            Eigen::MatrixX<double> intensityMatrixIntegratedLimitsNoInterferenceShadows;
            int allignedMaxIndexesCountUnused;
//            e = buildAlignmentMatrix(
//                    intensityMatrixShadows,
//                    pii,
//                    summedMatToVec,
//                    &intensityMatrixIntegratedLimitsNoInterferenceShadows,
//                    &unused,
//                    &allignedMaxIndexesCountUnused
//            ); ree;

            e = ErrorUtils::isEqual(
                    intensityMatrixIntegratedLimitsNoInterference.cols(),
                    intensityMatrixIntegratedLimitsNoInterferenceShadows.cols()
                    ); ree;

            e = ErrorUtils::isEqual(
                    intensityMatrixIntegratedLimitsNoInterference.rows(),
                    intensityMatrixIntegratedLimitsNoInterferenceShadows.rows()
            ); ree;

//            intensityMatrixIntegratedLimitsNoInterferenceShadows
//                = applyGaussSmoothRowWiseToMatrix(intensityMatrixIntegratedLimitsNoInterferenceShadows);

            if (subtractShadows) {
                intensityMatrixIntegratedLimitsNoInterference -= intensityMatrixIntegratedLimitsNoInterferenceShadows;
                EigenUtils::thresholdMatrix(0.0, &intensityMatrixIntegratedLimitsNoInterference);
            }

            for (int col = 0; col < intensityMatrixIntegratedLimitsNoInterference.cols(); col++) {

                const Eigen::VectorX<double> &mzPeak = intensityMatrixIntegratedLimitsNoInterference.col(col);
                const Eigen::VectorX<double> &mzPeakShadow = intensityMatrixIntegratedLimitsNoInterferenceShadows.col(col);

                const double cosineSimToShadow = EigenUtils::cosineSimilarity(
                        mzPeak,
                        mzPeakShadow
                        );

                cosineSimsShadowIndividual->push_back(cosineSimToShadow);

                if (cosineSimToShadow > 0.0 && !MathUtils::tZero(cosineSimToShadow)) {
                    const double ratio = mzPeak.maxCoeff() / mzPeakShadow.maxCoeff();
                    shadowsIntensityRatioPeakIntVec->push_back(ratio);
                }
                else{
                    shadowsIntensityRatioPeakIntVec->push_back(0.0);
                }
            }

        }

        if (MathUtils::tZero(intensityMatrixIntegratedLimitsNoInterference.maxCoeff())) {
            const QVector<double> v(topNMs2FragPeaks, 0.0);
            *cosineSimsIndividualIntegration = v;
            *bestCosineSimSumIntegration = 0.0;
            *bestIntensitiesIndividual = v;
            ERR_RETURN
        }

        e = calcBestCosineSimSumLogic(
                intensityMatrixIntegratedLimitsNoInterference,
                cosineSimToAnchorThreshold,
                cosineSimsIndividualIntegration,
                bestCosineSimSumIntegration,
                bestIntensitiesIndividual,
                bestAnchorColumn
        ); ree;

        ERR_RETURN
    }

    Err calculateCandidateAllignementMetrics(
            const Eigen::MatrixX<double> &intensityMatrix,
            const Eigen::MatrixX<double> &intensityMatrixShadows,
            const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
            const QVector<double> &summedMatToVec,
            int topNMs2FragPeaks,
            double cosineSimToAnchorThreshold,
            bool subtractShadows,
            QVector<double> *cosineSimsIndividual,
            QVector<double> *intensityFoundMaxVec,
            QVector<double> *cosineSimsShadowIndividual,
            QVector<double> *shadowsIntensityRatioVec,
            FrameIndex *frameIndexIntensityApex,
            PeakIntegrationIndexes *bestPeakIntegrationIndexes,
            Eigen::VectorX<double> *bestAnchorColumn,
            int *allignedMaxIndexesCount
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrix.nonZeros() > 1); ree;
        e = ErrorUtils::isNotEmpty(peakIntegrationIndexes); ree;
        e = ErrorUtils::isNotEmpty(summedMatToVec); ree;

        double bestCosineSimSum = 0.0;

        for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {

            QVector<double> cosineSimsIndividualIntegration;
            QVector<double> intensitiesIndividual;
            QVector<double> cosineSimsShadowIndividualPeakIntVec;
            QVector<double> shadowsIntensityRatioPeakIntVec;
            double bestCosineSimSumIntegration;
            FrameIndex frameIndexIntensityApexIntegration;
            int allignedMaxIndexesCountTemp;
            e = calcBestCosineSimSum(
                    intensityMatrix,
                    intensityMatrixShadows,
                    pii,
                    summedMatToVec,
                    topNMs2FragPeaks,
                    cosineSimToAnchorThreshold,
                    subtractShadows,
                    &cosineSimsIndividualIntegration,
                    &bestCosineSimSumIntegration,
                    &intensitiesIndividual,
                    bestAnchorColumn,
                    &frameIndexIntensityApexIntegration,
                    &allignedMaxIndexesCountTemp,
                    &cosineSimsShadowIndividualPeakIntVec,
                    &shadowsIntensityRatioPeakIntVec
            ); ree;

            if (bestCosineSimSumIntegration > bestCosineSimSum) {
                bestCosineSimSum = bestCosineSimSumIntegration;
                *cosineSimsIndividual = cosineSimsIndividualIntegration;
                *frameIndexIntensityApex = frameIndexIntensityApexIntegration;
                *bestPeakIntegrationIndexes = pii;
                *intensityFoundMaxVec = intensitiesIndividual;
                *allignedMaxIndexesCount = allignedMaxIndexesCountTemp;
                *cosineSimsShadowIndividual = cosineSimsShadowIndividualPeakIntVec;
                *shadowsIntensityRatioVec = shadowsIntensityRatioPeakIntVec;
            }

        }

        ERR_RETURN
    }

    Err calculateSpectrumMetrics(
            const Eigen::VectorX<double> &intensityApexVector,
            const QVector<MS2Ion> &ms2IonsCandidate,
            double *cosineSimSpectrum,
            double *klDivSpectrum
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(ms2IonsCandidate); ree;
        e = ErrorUtils::isTrue(intensityApexVector.rows() > 0); ree;

        QVector<double> intensityVals;
        std::transform(
                ms2IonsCandidate.begin(),
                ms2IonsCandidate.end(),
                std::back_inserter(intensityVals),
                [](const MS2Ion &ms2Ion){return ms2Ion.intensity;}
        );

        const Eigen::VectorX<double> intensityValsTheo = EigenUtils::convertQVectorToEigenVector(intensityVals);

        *cosineSimSpectrum = EigenUtils::cosineSimilarity(intensityApexVector, intensityValsTheo);
        *klDivSpectrum = EigenUtils::klDivergence(intensityApexVector, intensityValsTheo);

        ERR_RETURN
    }

    Err calculateMassAccuracyMetrics(
            const QVector<double> &cosineSimToAnchorVec,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            QVector<double> *mzMeanValsFound,
            QVector<double> *stdMeanValsFound,
            QVector<double> *mzValsSearched
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;
        e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;
        e = ErrorUtils::isFalse(MathUtils::tSame(peakIntegrationIndexes.first, peakIntegrationIndexes.second)); ree
        e = ErrorUtils::isTrue(cosineSimToAnchorVec.size() >= ms2IonsTheoretical.size()); ree;

        for (int index = 0; index < ms2IonsTheoretical.size(); index++) {

            const MS2Ion &ms2Ion = ms2IonsTheoretical.at(index);
            const double cosineSimMs2Ion = cosineSimToAnchorVec.at(index);

            mzValsSearched->push_back(ms2Ion.mz);

            if (MathUtils::tZero(cosineSimMs2Ion)) {
                mzMeanValsFound->push_back(0.0);
                stdMeanValsFound->push_back(0.0);
                continue;
            }

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            e = ErrorUtils::isTrue(mzHashedVsXICPoints.contains(mzHashed)); ree;

            const XICPoints &xicPoints = mzHashedVsXICPoints.value(mzHashed);
            const QMap<ScanNumber, QVector<double>> scanNumberVsMzVals = xicPoints.scanNumberVsMzVals();

            QVector<double> mzValsForMz;
            for (int i = peakIntegrationIndexes.first; i < peakIntegrationIndexes.second; i++) {

                const QVector<double> &mzVals = scanNumberVsMzVals.value(i);

                if (mzVals.isEmpty()) {
                    continue;
                }

                mzValsForMz.append(mzVals);
            }

            if (mzValsForMz.isEmpty()) {
                mzMeanValsFound->push_back(0.0);
                stdMeanValsFound->push_back(0.0);
                continue;
            }

            const double mzMean = MathUtils::mean(mzValsForMz);
            const double mzStd = MathUtils::stDev(mzValsForMz);
            mzMeanValsFound->push_back(mzMean);
            stdMeanValsFound->push_back(mzStd);
        }

        ERR_RETURN
    }

    Err calculateMS1Corr(
            const Eigen::VectorX<double> &bestAnchorColumn,
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            double mzTarget,
            double ppmTol,
            TurboXIC *turboXic,
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

        Eigen::VectorX<double> ms1Vec(static_cast<int>(bestAnchorColumn.size()));
        ms1Vec.setZero();

        const QMap<ScanNumber, double> &scanNumbersVsIntensityVals = xicPoints.scanNumbersVsIntensity;

        for (auto it = scanNumbersVsIntensityVals.begin(); it != scanNumbersVsIntensityVals.end(); it++) {
            const FrameIndex frameIndex = it.key() - peakIntegrationIndexes.first;

            if (frameIndex >= ms1Vec.size()) {
                continue;
            }

            ms1Vec.coeffRef(frameIndex) = it.value();
        }

//        ms1Vec = applyGaussSmoothRowWiseToMatrix(ms1Vec);

        *cosineSimMS1 = EigenUtils::cosineSimilarity(bestAnchorColumn, ms1Vec);

        ERR_RETURN
    }

    Err calculateEmpiricalMzStats(
            const QVector<double> &cosineSimToAnchorVec,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints100,
            const QVector<MS2Ion> &ms2IonsTheoretical,
            PeakIntegrationIndexes &bestPeakIntegrationIndexes,
            QVector<double> *mzMeanValsFound,
            QVector<double> *stdMeanValsFound,
            QVector<double> *mzValsSearched,
            QVector<double> *theoApexIntensity
    ) {

        ERR_INIT

        e = ErrorUtils::isAboveThreshold(
                bestPeakIntegrationIndexes.second,
                bestPeakIntegrationIndexes.first,
                ErrorUtilsParam::ExcludeThreshold
        ); ree;

        e = calculateMassAccuracyMetrics(
                cosineSimToAnchorVec,
                mzHashedVsXICPoints100,
                ms2IonsTheoretical,
                bestPeakIntegrationIndexes,
                mzMeanValsFound,
                stdMeanValsFound,
                mzValsSearched
        ); ree;

        QVector<QPair<double, double>> mzFoundVsMzSearched;
        e = ParallelUtils::zip(
                *mzMeanValsFound,
                *mzValsSearched,
                &mzFoundVsMzSearched
        ); ree

        std::transform(
                ms2IonsTheoretical.begin(),
                ms2IonsTheoretical.end(),
                std::back_inserter(*theoApexIntensity),
                [](const MS2Ion &ms2Ion){return ms2Ion.intensity;}
        );

        ERR_RETURN
    }

    Err buildTightPPMToleranceData(
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints100,
            double ppmTol100,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints45,
            QMap<MzHashed, XICPoints> *mzHashedVsXICPoints20
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints100); ree;

        const double ppmTight1 = ppmTol100 * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION;
        const double ppmTight2 = ppmTol100 * S_GLOBAL_SETTINGS.TIGHT_2_FRACTION;

        for (auto it = mzHashedVsXICPoints100.begin(); it != mzHashedVsXICPoints100.end(); it++) {

            const MzHashed mzHashed = it.key();
            const XICPoints &xicPoints = it.value();
            const QMap<ScanNumber, ScanPoints> &scanNumbersVsScanPoints = xicPoints.scanNumbersVsScanPoints;

            const auto mz = MathUtils::unHashDecimal<double>(
                    mzHashed,
                    S_GLOBAL_SETTINGS.HASHING_PRECISION
                    );

            XICPoints xicPointsTight1New;
            XICPoints xicPointsTight2New;
            for (auto itt = scanNumbersVsScanPoints.begin(); itt != scanNumbersVsScanPoints.end(); itt++) {

                const ScanNumber scanNumber = itt.key();
                const ScanPoints &scanPoints = itt.value();

                for (const ScanPoint &sp : scanPoints) {
                    const double ppmDiff = std::abs(1e6 * (sp.x() - mz)) / mz;
                    if (ppmDiff <= ppmTight1) {
                        xicPointsTight1New.scanNumbersVsScanPoints[scanNumber].push_back(sp);
//                        xicPointsTight1New.scanNumbersVsIntensityVals[scanNumber] += sp.y();

                        if (ppmDiff <= ppmTight2) {
                            xicPointsTight2New.scanNumbersVsScanPoints[scanNumber].push_back(sp);
//                            xicPointsTight2New.scanNumbersVsIntensityVals[scanNumber] += sp.y();
                        }
                    }
                }
            }
            e = xicPointsTight1New.buildScanNumbersVsIntensityVals(); ree;
            e = xicPointsTight2New.buildScanNumbersVsIntensityVals(); ree;
            mzHashedVsXICPoints45->insert(mzHashed, xicPointsTight1New);
            mzHashedVsXICPoints20->insert(mzHashed, xicPointsTight2New);
        }

        ERR_RETURN
    }

}//namespace
Err ScoreOverseer::buildScores(
        const TargetDecoyCandidatePair* targetDecoyCandidatePair,
        QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHills,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
        const QHash<MzHashed , QVector<FeatureFinderHill*>> &mzHashedVsfeatureFinderHillsShadows,
        CandidateScores *candidateScores
        ) {

    ERR_INIT

    if (peakIntegrationIndexes.isEmpty()) {
        ERR_RETURN
    }

    e = ErrorUtils::isNotEmpty(peakIntegrationIndexes); ree;
    e = ErrorUtils::isNotEmpty(ms2IonsTheoretical); ree;
    e = ErrorUtils::isNotEmpty(mzHashedVsfeatureFinderHills); ree;
    e = ErrorUtils::isNotEmpty(ms2IonsTheoreticalIsotopeShadows); ree;

    QElapsedTimer et;
    et.start();

    e = d_ptr->buildAlignmentMatricies(
            ms2IonsTheoretical,
            ms2IonsTheoreticalIsotopeShadows,
            mzHashedVsfeatureFinderHills,
            mzHashedVsfeatureFinderHillsShadows
    );

//    qDebug() << "2" << et.nsecsElapsed() << et.restart();

    FrameIndex frameIndexIntensityApex;
    PeakIntegrationIndexes bestPeakIntegrationIndexes;
    Eigen::VectorX<double> bestAnchorColumn;
    int allignedMaxIndexesCount;
//    e = calculateCandidateAllignementMetrics(
//            d_ptr->m_intensityMatrix100,
//            d_ptr->m_intensityMatrix100Shadow,
//            peakIntegrationIndexes,
//            d_ptr->m_summedMatVecToVec,
//            static_cast<int>(d_ptr->m_intensityMatrix100.cols()),
//            d_ptr->m_cosineSimToAnchorThreshold,
//            subtractShadows,
//            &candidateScores->cosineSimToAnchorVec,
//            &candidateScores->intensityFoundMaxVec,
//            &candidateScores->cosineSimShadowsToAnchorVec,
//            &candidateScores->shadowsIntensityRatioVec,
//            &frameIndexIntensityApex,
//            &bestPeakIntegrationIndexes,
//            &bestAnchorColumn,
//            &allignedMaxIndexesCount
//    ); ree;

//    candidateScores->shadowsCosineSimSum = std::accumulate(
//            candidateScores->cosineSimShadowsToAnchorVec.begin(),
//            candidateScores->cosineSimShadowsToAnchorVec.end(),
//            0.0
//            );

//    qDebug() << "3" << et.nsecsElapsed() << et.restart();

//    if (bestAnchorColumn.rows() < 1) {
//        ERR_RETURN
//    }

//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            targetDecoyCandidatePair->mz(),
//            d_ptr->m_ppmTol,
//            m_turboXICMS1,
//            &candidateScores->cosineSim100MS1
//            ); ree;

//    qDebug() << "4" << et.nsecsElapsed() << et.restart();

//    e = calculateSpectrumMetrics(
//            d_ptr->m_intensityMatrix100.row(frameIndexIntensityApex),
//            ms2IonsTheoretical,
//            &candidateScores->cosineSimSpectrum,
//            &candidateScores->klDivSpectrum
//            ); ree;

//    qDebug() << "5" << et.nsecsElapsed() << et.restart();
//
//    e = calculateEmpiricalMzStats(
//            candidateScores->cosineSimToAnchorVec,
//            mzHashedVsXICPoints100,
//            ms2IonsTheoretical,
//            bestPeakIntegrationIndexes,
//            &candidateScores->mzFoundMeanVec,
//            &candidateScores->mzFoundStDevVec,
//            &candidateScores->mzSearchedVec,
//            &candidateScores->theoIntensityVec
//            );

//    qDebug() << "6" << et.nsecsElapsed() << et.restart();

    QVector<double> cosineSimsIndividual45;
    double unused1;
    Eigen::VectorX<double> unused2;
    QVector<double> unusedIntensity;
    int allignedMaxIndexesCountUnused;
    QVector<double> unusedShadows;
    QVector<double> unusedShadowsRatios;
    FrameIndex unused3;
//    e = calcBestCosineSimSum(
//            d_ptr->m_intensityMatrix45,
//            Eigen::MatrixX<double>(),
//            bestPeakIntegrationIndexes,
//            d_ptr->m_summedMatVecToVec,
//            static_cast<int>(d_ptr->m_intensityMatrix45.cols()),
//            d_ptr->m_cosineSimToAnchorThreshold,
//            subtractShadows,
//            &cosineSimsIndividual45,
//            &unused1,
//            &unusedIntensity,
//            &unused2,
//            &unused3,
//            &allignedMaxIndexesCountUnused,
//            &unusedShadows,
//            &unusedShadowsRatios
//            ); ree;

//    qDebug() << "7" << et.nsecsElapsed() << et.restart();

    QVector<double> cosineSimsIndividual20;
//    e = calcBestCosineSimSum(
//            d_ptr->m_intensityMatrix20,
//            Eigen::MatrixX<double>(),
//            bestPeakIntegrationIndexes,
//            d_ptr->m_summedMatVecToVec,
//            static_cast<int>(d_ptr->m_intensityMatrix20.cols()),
//            d_ptr->m_cosineSimToAnchorThreshold,
//            subtractShadows,
//            &cosineSimsIndividual20,
//            &unused1,
//            &unusedIntensity,
//            &unused2,
//            &unused3,
//            &allignedMaxIndexesCountUnused,
//            &unusedShadows,
//            &unusedShadowsRatios
//            ); ree;

//    qDebug() << "8" << et.nsecsElapsed() << et.restart();

//    candidateScores->cosineSimToAnchorVec.resize(candidateScores->mzSearchedVec.size());
//
//    candidateScores->intensityFoundMaxVec.resize(candidateScores->mzSearchedVec.size());
//
//    candidateScores->peptideStringWithMods = targetDecoyCandidatePair->peptideStringWithMods();
//
//    candidateScores->allignedMaxIndexesCount = allignedMaxIndexesCount;
//
//    const int top6 = 6;
//    candidateScores->cosineSimSum100 = std::accumulate(
//            candidateScores->cosineSimToAnchorVec.begin(),
//            candidateScores->cosineSimToAnchorVec.begin() + top6,
//            0.0
//            );
//
//    candidateScores->cosineSimSum45 = std::accumulate(
//            cosineSimsIndividual45.begin(),
//            cosineSimsIndividual45.begin() + top6,
//            0.0
//            );
//
//    candidateScores->cosineSimSum20 = std::accumulate(
//            cosineSimsIndividual20.begin(),
//            cosineSimsIndividual20.begin() + top6,
//            0.0
//            );
//
//    if (candidateScores->cosineSimToAnchorVec.size() > top6) {
//        candidateScores->cosineSimSumBottom6 = std::accumulate(
//                candidateScores->cosineSimToAnchorVec.begin() + top6 + 1,
//                candidateScores->cosineSimToAnchorVec.end(),
//                0.0
//        );
//    }
//
//    candidateScores->charge = targetDecoyCandidatePair->charge();
//    candidateScores->mass = targetDecoyCandidatePair->mass();
//    candidateScores->scanNumber = msFrame->scanNumberFromFrameIndex(frameIndexIntensityApex);
//    candidateScores->scanTime = msFrame->scanTimeFromScanNumber(candidateScores->scanNumber);
//    candidateScores->scanIonCount = msFrame->getScanPointsByScanNumber(candidateScores->scanNumber)->size();
//    candidateScores->iRTPredicted = targetDecoyCandidatePair->iRt();
//    candidateScores->scanTimePredicted = scanTimePredicted;
//    candidateScores->theoFragmentCount = targetDecoyCandidatePair->totalFragmentCount();
//    candidateScores->targetKey = targetDecoyCandidatePair->bestDiscriminateScoreKeyTarget();

//    const double precursorMz = BiophysicalCalcs::calculateThomsonFromMass(targetDecoyCandidatePair->mass(), targetDecoyCandidatePair->charge());
//    const double precursorMzIso1 = precursorMz + (S_GLOBAL_SETTINGS.ISO_DIFF / targetDecoyCandidatePair->charge());
//    const double precursorMzIso2 = precursorMz + ((2 * S_GLOBAL_SETTINGS.ISO_DIFF) / targetDecoyCandidatePair->charge());
//
//    const double chunkDivision = 3.0;
//
//    QVector<double> bestAnchorColumnVec = EigenUtils::convertEigenVectorToQVector(bestAnchorColumn);
//    const auto terminatorLogic = [](double d){return d < 1.0;};
//    const auto terminator = std::remove_if(bestAnchorColumnVec.begin(), bestAnchorColumnVec.end(), terminatorLogic);
//    bestAnchorColumnVec.erase(terminator, bestAnchorColumnVec.end());
//    const int chunkSize = std::max(1, static_cast<int>(std::round(bestAnchorColumnVec.size() / chunkDivision)));
//    const double bestAnchorColumnVecSum = std::accumulate(bestAnchorColumnVec.begin(), bestAnchorColumnVec.end(), 0.0001);
//
//    if (bestAnchorColumnVec.size() < chunkDivision) {
//        candidateScores->peakShapeRatio1 = 0.0;
//        candidateScores->peakShapeRatio2 = 1.0;
//        candidateScores->peakShapeRatio3 = 0.0;
//    }
//    else {
//
//        candidateScores->peakShapeRatio1 = std::accumulate(
//                bestAnchorColumnVec.begin(),
//                bestAnchorColumnVec.begin() + chunkSize,
//                0.0
//        ) / bestAnchorColumnVecSum;
//
//        candidateScores->peakShapeRatio2 = std::accumulate(
//                bestAnchorColumnVec.begin() + chunkSize,
//                bestAnchorColumnVec.begin() + (chunkSize * 2),
//                0.0
//        ) / bestAnchorColumnVecSum;
//
//        candidateScores->peakShapeRatio3 = std::accumulate(
//                bestAnchorColumnVec.begin() + (chunkSize * 2),
//                bestAnchorColumnVec.end(),
//                0.0
//        ) / bestAnchorColumnVecSum;
//
//    }

//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            precursorMz,
//            d_ptr->m_ppmTol,
//            m_turboXICMS1,
//            &candidateScores->cosineSim100MS1
//    ); ree;

//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            precursorMzIso1,
//            d_ptr->m_ppmTol,
//            m_turboXICMS1,
//            &candidateScores->cosineSim100MS1Iso1
//    ); ree;
//
//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            precursorMzIso2,
//            d_ptr->m_ppmTol,
//            m_turboXICMS1,
//            &candidateScores->cosineSim100MS1Iso2
//    ); ree;

    ERR_RETURN
}

Err ScoreOverseer::init(
        const PythiaParameters &pythiaParameters,
        FeatureFinderHillBuilder *featureFinderHillsBuilderMS1
        ) {
    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = d_ptr->init(pythiaParameters); ree;
    m_featureFinderHillsBuilderMS1 = featureFinderHillsBuilderMS1;
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
