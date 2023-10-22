//
// Created by anichols on 10/20/23.
//

#include "ScoreOverseer.h"

#include "CandidateScores.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MS2Ion.h"
#include "ParallelUtils.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"


class Q_DECL_HIDDEN ScoreOverseer::Private
{

public:

    Private(
            int topNMS2Ions,
            double cosineSimToAnchorThreshold,
            double ppmTol,
            const QVector<double> &summedMatVecToVec
            );
    ~Private();

    Err initMatricies(
            const QVector<MS2Ion> &ms2IonsTheoretical,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints100,
            const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPointsIsotopeShadows
            );

    Eigen::MatrixX<double> m_intensityMatrix100;
    Eigen::MatrixX<double> m_intensityMatrix100Shadow;

    int m_topNMS2Ions;
    double m_cosineSimToAnchorThreshold;
    double m_ppmTol;
    QVector<double> m_summedMatVecToVec;

};

ScoreOverseer::Private::Private(
        int topNMS2Ions,
        double cosineSimToAnchorThreshold,
        double ppmTol,
        const QVector<double> &summedMatVecToVec
        )
: m_topNMS2Ions(topNMS2Ions)
, m_cosineSimToAnchorThreshold(cosineSimToAnchorThreshold)
, m_ppmTol(ppmTol)
, m_summedMatVecToVec(summedMatVecToVec)
{}

ScoreOverseer::Private::~Private() {}

namespace {

    Eigen::MatrixX<double> buildIntensityVecMatrix(
            const QVector<MZION> &mzIons,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            int vecSize,
            int topNMzFrags
            ) {

        const int cols = topNMzFrags;
        const int rows = vecSize;

        Eigen::MatrixX<double> mat(rows, cols);
        mat.setZero();

        for (int i = 0; i < cols; i++) {

            if (i >= mzIons.size()) {
                break;
            }

            const MZION &mz = mzIons.at(i);

            const MzHashed mzHashed = MathUtils::hashDecimal(mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

            QVector<double> intensityVals = ParallelUtils::convertMapToVector(
                    mzHashedVsXICPoints.value(mzHashed).scanNumbersVsIntensityVals,
                    rows
            );

            if (intensityVals.isEmpty()) {
                intensityVals = QVector<double>(vecSize, 0.0);
            }

            const Eigen::VectorX<double> eVec = EigenUtils::convertQVectorToEigenVector(intensityVals);
            mat.col(i) = eVec;
        }

        return mat;
    }

}// namespace
Err ScoreOverseer::Private::initMatricies(
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints100,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPointsIsotopeShadows
        ) {

    ERR_INIT

    QVector<MZION> mzIonsTopN;
    std::transform(
            ms2IonsTheoretical.begin(),
            ms2IonsTheoretical.end(),
            std::back_inserter(mzIonsTopN),
            [](const MS2Ion &ms2Ion){return ms2Ion.mz;}
    );

    m_intensityMatrix100 = buildIntensityVecMatrix(
            mzIonsTopN,
            mzHashedVsXICPoints100,
            m_summedMatVecToVec.size(),
            m_topNMS2Ions
    );

    QVector<MZION> mzIonsShadow;
    std::transform(
            ms2IonsTheoreticalIsotopeShadows.begin(),
            ms2IonsTheoreticalIsotopeShadows.end(),
            std::back_inserter(mzIonsShadow),
            [](const MS2Ion &ms2Ion){return ms2Ion.mz;}
    );

    m_intensityMatrix100Shadow = buildIntensityVecMatrix(
            mzIonsShadow,
            mzHashedVsXICPointsIsotopeShadows,
            m_summedMatVecToVec.size(),
            mzIonsShadow.size()
    );

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


ScoreOverseer::ScoreOverseer(
        int topNMS2Ions,
        double cosineSimToAnchorThreshold,
        double ppmTol,
        const QVector<double> &summedMatVecToVec,
        TurboXIC *turboXICMS1
        )
: d_ptr(new Private(topNMS2Ions, cosineSimToAnchorThreshold, ppmTol, summedMatVecToVec))
, m_turboXICMS1(turboXICMS1)
{}

ScoreOverseer::~ScoreOverseer() {}

namespace {

    Eigen::MatrixX<double> applyGaussSmoothRowWiseToMatrix(const Eigen::MatrixX<double> &mat) {

        const int filterLength = 3;
        const double sigma = 1.0;
        const Eigen::VectorX<double> gaussKernel = EigenKernelUtils::buildGaussianFilter1D(filterLength, sigma);

        return EigenKernelUtils::applyKernelRowWiseToMatrix(mat, gaussKernel);
    }

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
            Eigen::VectorX<double> *vec) {

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
            Eigen::MatrixX<double> *intensityMatrixIntegratedLimitsNoInterference
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peakCountsPerRow); ree;
        e = ErrorUtils::isTrue(_intensityMatrixIntegratedLimits.maxCoeff() > 0); ree;

        Eigen::MatrixX<double> intensityMatrixIntegratedLimits = _intensityMatrixIntegratedLimits;

        int bestApexRowIndex;
        e = findBestApexRowInMatrix(intensityMatrixIntegratedLimits, peakCountsPerRow, &bestApexRowIndex); ree;

        for (int col = 0; col < intensityMatrixIntegratedLimits.cols(); col++) {

            Eigen::VectorX<double> colVec = intensityMatrixIntegratedLimits.col(col);
            const QMap<int, double> columnApexes = EigenUtils::apexes(colVec);

            const QList<int> &keys = columnApexes.keys();
            const int bestApexIndexColumn = keys.at(MathUtils::closest(keys.toVector(), bestApexRowIndex));

            e = removeInterferingPeaksInColumn(bestApexIndexColumn, &colVec); ree;
            intensityMatrixIntegratedLimits.col(col) = colVec;

        }

        *intensityMatrixIntegratedLimitsNoInterference = intensityMatrixIntegratedLimits;

        ERR_RETURN
    }

    Err calcBestCosineSimSumLogic(
            const Eigen::MatrixX<double> &intensityMatrixIntegratedLimits,
            double cosineSimToAnchorThreshold,
            QVector<double> *bestCosineSimsIndividual,
            double *bestCosineSimSum,
            Eigen::VectorX<double> *bestAnchorColumn
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrixIntegratedLimits.rows() > 0); ree;

        *bestCosineSimSum = 0.0;

        const Eigen::MatrixX<double> intensityMatrixIntegratedLimitsSmoothed
                = applyGaussSmoothRowWiseToMatrix(intensityMatrixIntegratedLimits);

        for (int i = 0; i < intensityMatrixIntegratedLimitsSmoothed.cols(); i++) {

            QVector<double> bestCosineSimsIndividualAnchor;
            const Eigen::VectorX<double> anchorColumn = intensityMatrixIntegratedLimitsSmoothed.col(i);

            const QPair<int, double> anchorFrameIndexMaxVsVal = EigenUtils::returnTopIndexAndValue(anchorColumn);

            for (int j = 0; j < intensityMatrixIntegratedLimitsSmoothed.cols(); j++) {

                const Eigen::VectorX<double> altColumn = intensityMatrixIntegratedLimitsSmoothed.col(j);
                const QPair<int, double> altColumnFrameIndexMaxVsVal = EigenUtils::returnTopIndexAndValue(altColumn);

                const double cosineSimToAnchor = EigenUtils::cosineSimilarity(anchorColumn, altColumn);

                if (cosineSimToAnchor < cosineSimToAnchorThreshold) {
                    bestCosineSimsIndividualAnchor.push_back(0.0);
                    continue;
                }

                bestCosineSimsIndividualAnchor.push_back(cosineSimToAnchor);
            }

            const double cosineSimSumAnchor
                    = std::accumulate(bestCosineSimsIndividualAnchor.begin(), bestCosineSimsIndividualAnchor.end(), 0.0);

            if (cosineSimSumAnchor > *bestCosineSimSum) {
                *bestCosineSimSum = cosineSimSumAnchor;
                *bestCosineSimsIndividual = bestCosineSimsIndividualAnchor;
                *bestAnchorColumn = anchorColumn;
            }

        }

        ERR_RETURN
    }

    Err buildAlignmentMatrix(
            const Eigen::MatrixX<double> &intensityMatrix,
            const PeakIntegrationIndexes &pii,
            const QVector<double> &summedMatToVec,
            int topNMs2FragPeaks,
            Eigen::MatrixX<double> *intensityMatrixIntegratedLimitsNoInterference,
            FrameIndex *frameIndexIntensityApexIntegration
    ) {

        ERR_INIT

        const int rowStart = pii.first;
        const int rowCount = pii.second < intensityMatrix.rows()
                             ? pii.second - pii.first + 1
                             : static_cast<int>(intensityMatrix.rows()) - pii.first;

        const int colStart = 0;
        const int colCount = topNMs2FragPeaks;

        const Eigen::MatrixX<double> intensityMatrixIntegratedLimits = intensityMatrix.block(
                rowStart,
                colStart,
                rowCount,
                colCount
        );

        if (MathUtils::tZero(intensityMatrixIntegratedLimits.maxCoeff())) {
            *intensityMatrixIntegratedLimitsNoInterference = intensityMatrixIntegratedLimits;
            ERR_RETURN
        }

        const QVector<double> peakCountsPerRow = summedMatToVec.mid(rowStart, rowCount);
        *frameIndexIntensityApexIntegration = MathUtils::findMaxIndexInVector(peakCountsPerRow) + rowStart;

        e = removeInterferingPeaksInMatrix(
                intensityMatrixIntegratedLimits,
                peakCountsPerRow,
                intensityMatrixIntegratedLimitsNoInterference
        ); ree;

//#define CHECK_ALIGNMENT
#ifdef CHECK_ALIGNMENT
        std::cout << "FrameIndex " << frameIndexIntensityApexIntegration << std::endl;
        std::cout << intensityMatrixIntegratedLimits << std::endl;
        std::cout << "no interference" << std::endl;
        std::cout << intensityMatrixIntegratedLimitsNoInterference << std::endl;
        std::cout << pii.first << ", " << pii.second << std::endl;
#endif

        ERR_RETURN
    }

    Err calcBestCosineSimSum(
            const Eigen::MatrixX<double> &intensityMatrix,
            const PeakIntegrationIndexes &pii,
            const QVector<double> &summedMatToVec,
            int topNMs2FragPeaks,
            double cosineSimToAnchorThreshold,
            QVector<double> *cosineSimsIndividualIntegration,
            double *bestCosineSimSumIntegration,
            Eigen::VectorX<double> *bestAnchorColumn,
            FrameIndex *frameIndexIntensityApexIntegration
    ) {

        ERR_INIT

        Eigen::MatrixX<double> intensityMatrixIntegratedLimitsNoInterference;

        e = buildAlignmentMatrix(
                intensityMatrix,
                pii,
                summedMatToVec,
                topNMs2FragPeaks,
                &intensityMatrixIntegratedLimitsNoInterference,
                frameIndexIntensityApexIntegration
        ); ree;

        if (MathUtils::tZero(intensityMatrixIntegratedLimitsNoInterference.maxCoeff())) {
            const QVector<double> v(topNMs2FragPeaks, 0.0);
            *cosineSimsIndividualIntegration = v;
            *bestCosineSimSumIntegration = 0.0;
            ERR_RETURN
        }

        e = calcBestCosineSimSumLogic(
                intensityMatrixIntegratedLimitsNoInterference,
                cosineSimToAnchorThreshold,
                cosineSimsIndividualIntegration,
                bestCosineSimSumIntegration,
                bestAnchorColumn
        ); ree;

        ERR_RETURN
    }

    Err calculateCandidateAllignementMetrics(
            const Eigen::MatrixX<double> &intensityMatrix,
            const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
            const QVector<double> &summedMatToVec,
            int topNMs2FragPeaks,
            double cosineSimToAnchorThreshold,
            QVector<double> *cosineSimsIndividual,
            FrameIndex *frameIndexIntensityApex,
            PeakIntegrationIndexes *bestPeakIntegrationIndexes,
            Eigen::VectorX<double> *bestAnchorColumn
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrix.nonZeros() > 1); ree;
        e = ErrorUtils::isNotEmpty(peakIntegrationIndexes); ree;
        e = ErrorUtils::isNotEmpty(summedMatToVec); ree;

        double bestCosineSimSum = 0.0;

        for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {

            QVector<double> cosineSimsIndividualIntegration;
            double bestCosineSimSumIntegration;
            FrameIndex frameIndexIntensityApexIntegration;
            e = calcBestCosineSimSum(
                    intensityMatrix,
                    pii,
                    summedMatToVec,
                    topNMs2FragPeaks,
                    cosineSimToAnchorThreshold,
                    &cosineSimsIndividualIntegration,
                    &bestCosineSimSumIntegration,
                    bestAnchorColumn,
                    &frameIndexIntensityApexIntegration
            ); ree;

            if (bestCosineSimSumIntegration > bestCosineSimSum) {
                bestCosineSimSum = bestCosineSimSumIntegration;
                *cosineSimsIndividual = cosineSimsIndividualIntegration;
                *frameIndexIntensityApex = frameIndexIntensityApexIntegration;
                *bestPeakIntegrationIndexes = pii;
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

        for (const MS2Ion &ms2Ion : ms2IonsTheoretical) {

            mzValsSearched->push_back(ms2Ion.mz);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

            if (!mzHashedVsXICPoints.contains(mzHashed)) {
                mzMeanValsFound->push_back(0.0);
                stdMeanValsFound->push_back(0.0);
                continue;
            }

            const XICPoints &xicPoints = mzHashedVsXICPoints.value(mzHashed);
            const QMap<ScanNumber, QVector<double>> scanNumberVsMzVals = xicPoints.scanNumberVsMzVals();

            QVector<double> mzValsForMz;
            for (int i = peakIntegrationIndexes.first; i < peakIntegrationIndexes.second; i++) {

                if (i >= peakIntegrationIndexes.second) {
                    const double mzMean = MathUtils::mean(mzValsForMz);
                    const double mzStd = MathUtils::stDev(mzValsForMz);
                    mzMeanValsFound->push_back(mzMean);
                    stdMeanValsFound->push_back(mzStd);
                    break;
                }

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

        const QMap<ScanNumber, double> &scanNumbersVsIntensityVals = xicPoints.scanNumbersVsIntensityVals;

        for (auto it = scanNumbersVsIntensityVals.begin(); it != scanNumbersVsIntensityVals.end(); it++) {
            const FrameIndex frameIndex = it.key() - peakIntegrationIndexes.first;

            if (frameIndex >= ms1Vec.size()) {
                continue;
            }

            ms1Vec.coeffRef(frameIndex) = it.value();
        }

        ms1Vec = applyGaussSmoothRowWiseToMatrix(ms1Vec);

        *cosineSimMS1 = EigenUtils::cosineSimilarity(bestAnchorColumn, ms1Vec);

        ERR_RETURN
    }

    QVector<double> calculatePeakShapeRatios(const Eigen::VectorX<double> &vec) {

        const Eigen::VectorX<double> vecTrimmed = EigenUtils::trimVector(vec);
        const int segmentSize = std::round(vecTrimmed.size() / 3.0);

        const double seg1Sum = vecTrimmed.segment(0, segmentSize).sum();
        const double seg2Sum = vecTrimmed.segment(segmentSize, segmentSize).sum();
        const double seg3Sum = vecTrimmed.segment(segmentSize * 2, segmentSize).sum();

        return {
                seg1Sum / seg2Sum,
                seg3Sum / seg2Sum,
                seg1Sum / std::max(seg3Sum, std::numeric_limits<double>::min())
        };
    }

}//namespace
Err ScoreOverseer::buildScores(
        const TargetDecoyCandidatePair* targetDecoyCandidatePair,
        const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
        const QVector<MS2Ion> &ms2IonsTheoretical,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints100,
        const QVector<MS2Ion> &ms2IonsTheoreticalIsotopeShadows,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPointsIsotopeShadows,
        CandidateScores *candidateScores
) {

    ERR_INIT

    QVector<MS2Ion> ms2IonsTheoreticalResized = ms2IonsTheoretical;
    ms2IonsTheoreticalResized.resize(d_ptr->m_topNMS2Ions);

    e = d_ptr->initMatricies(
            ms2IonsTheoretical,
            mzHashedVsXICPoints100,
            ms2IonsTheoreticalIsotopeShadows,
            mzHashedVsXICPointsIsotopeShadows
            ); ree;

    QVector<double> cosineSimsIndividual100;
    FrameIndex frameIndexIntensityApex;
    PeakIntegrationIndexes bestPeakIntegrationIndexes;
    Eigen::VectorX<double> bestAnchorColumn;
    e = calculateCandidateAllignementMetrics(
            d_ptr->m_intensityMatrix100,
            peakIntegrationIndexes,
            d_ptr->m_summedMatVecToVec,
            d_ptr->m_topNMS2Ions,
            d_ptr->m_cosineSimToAnchorThreshold,
            &cosineSimsIndividual100,
            &frameIndexIntensityApex,
            &bestPeakIntegrationIndexes,
            &bestAnchorColumn
    ); ree;


    e = calculateMS1Corr(
            bestAnchorColumn,
            bestPeakIntegrationIndexes,
            targetDecoyCandidatePair->mz(),
            d_ptr->m_ppmTol,
            m_turboXICMS1,
            &candidateScores->cosineSim100MS1
            ); ree;

    e = calculateSpectrumMetrics(
            d_ptr->m_intensityMatrix100.row(frameIndexIntensityApex),
            ms2IonsTheoreticalResized,
            &candidateScores->cosineSimSpectrum,
            &candidateScores->klDivSpectrum
    ); ree;

    candidateScores->cosineSimSum100
        = std::accumulate(cosineSimsIndividual100.begin(), cosineSimsIndividual100.end(), 0.0);


//    QVector<double> cosineSimsIndividualShadows;
//    double cosineSimSumShadows;
//    Eigen::VectorX<double> unused10;
//    FrameIndex unused11;
//    e = calcBestCosineSimSum(
//            d_ptr->m_intensityMatrix100Shadow,
//            bestPeakIntegrationIndexes,
//            d_ptr->m_summedMatVecToVec,
//            d_ptr->m_topNMS2Ions,
//            d_ptr->m_cosineSimToAnchorThreshold,
//            &cosineSimsIndividualShadows,
//            &cosineSimSumShadows,
//            &unused10,
//            &unused11
//    ); ree;

//    const QVector<double> peakShapeRatios = calculatePeakShapeRatios(bestAnchorColumn);
//
//    QVector<double> cosineSimsIndividual45;
//    double unused1;
//    Eigen::VectorX<double> unused2;
//    FrameIndex unused3;
//    e = calcBestCosineSimSum(
//            intensityMatrix45,
//            bestPeakIntegrationIndexes,
//            summedMatVecToVec,
//            m_topNMS2Ions,
//            m_pythiaParameters.cosineSimToAnchorThreshold,
//            &cosineSimsIndividual45,
//            &unused1,
//            &unused2,
//            &unused3
//            ); ree;
//
//    QVector<double> cosineSimsIndividual20;
//    double unused4;
//    Eigen::VectorX<double> unused5;
//    FrameIndex unused6;
//    e = calcBestCosineSimSum(
//            intensityMatrix20,
//            bestPeakIntegrationIndexes,
//            summedMatVecToVec,
//            m_topNMS2Ions,
//            m_pythiaParameters.cosineSimToAnchorThreshold,
//            &cosineSimsIndividual20,
//            &unused4,
//            &unused5,
//            &unused6
//    ); ree;
//
//    QVector<double> cosineSimsIndividualB2B3;
//    double cosineSimSumB2B3;
//    Eigen::VectorX<double> unused8;
//    FrameIndex unused9;
//    e = calcBestCosineSimSum(
//            intensityMatrixB2B3,
//            bestPeakIntegrationIndexes,
//            summedMatVecToVec,
//            candidatePeptide.ms2IonMzB2B3.size(),
//            m_pythiaParameters.cosineSimToAnchorThreshold,
//            &cosineSimsIndividualB2B3,
//            &cosineSimSumB2B3,
//            &unused8,
//            &unused9
//    ); ree;

//            d_ptr->m_intensityMatrix100,
//            peakIntegrationIndexes,
//            d_ptr->m_summedMatVecToVec,
//            d_ptr->m_topNMS2Ions,
//            d_ptr->m_cosineSimToAnchorThreshold,
//            &cosineSimsIndividual100,
//            &frameIndexIntensityApex,
//            &bestPeakIntegrationIndexes,
//            &bestAnchorColumn

//    QVector<double> intensityApexVals100 = EigenUtils::convertEigenVectorToQVector(
//            Eigen::VectorX<double>(intensityMatrix100.row(frameIndexIntensityApex))
//    );
//
//    double cosineSimSpectrum;
//    double klDivSpectrum;
//    e = calculateSpectrumMetrics(
//            intensityMatrix100.row(frameIndexIntensityApex),
//            candidatePeptide.ms2Ions,
//            &cosineSimSpectrum,
//            &klDivSpectrum
//    ); ree;
//
//    const double precursorMz = BiophysicalCalcs::calculateThomsonFromMass(candidatePeptide.mass, candidatePeptide.charge);
//    const double precursorMzIso1 = precursorMz + (S_GLOBAL_SETTINGS.ISO_DIFF / candidatePeptide.charge);
//    const double precursorMzIso2 = precursorMz + ((2 * S_GLOBAL_SETTINGS.ISO_DIFF) / candidatePeptide.charge);
//
//    double cosineSim100MS1;
//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            precursorMz,
//            m_pythiaParameters.ms2ExtractionWidthPPM,
//            &m_turboXICMS1,
//            &cosineSim100MS1
//    ); ree;
//
//    double cosineSim100MS1Iso1;
//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            precursorMzIso1,
//            m_pythiaParameters.ms2ExtractionWidthPPM,
//            &m_turboXICMS1,
//            &cosineSim100MS1Iso1
//    ); ree;
//
//    double cosineSim100MS1Iso2;
//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            precursorMzIso2,
//            m_pythiaParameters.ms2ExtractionWidthPPM,
//            &m_turboXICMS1,
//            &cosineSim100MS1Iso2
//    ); ree;
//
//    double cosineSim45MS1;
//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            BiophysicalCalcs::calculateThomsonFromMass(candidatePeptide.mass, candidatePeptide.charge),
//            m_pythiaParameters.ms2ExtractionWidthPPM * S_GLOBAL_SETTINGS.TIGHT_1_FRACTION,
//            &m_turboXICMS1,
//            &cosineSim45MS1
//    ); ree;
//
//    double cosineSim20MS1;
//    e = calculateMS1Corr(
//            bestAnchorColumn,
//            bestPeakIntegrationIndexes,
//            BiophysicalCalcs::calculateThomsonFromMass(candidatePeptide.mass, candidatePeptide.charge),
//            m_pythiaParameters.ms2ExtractionWidthPPM * S_GLOBAL_SETTINGS.TIGHT_2_FRACTION,
//            &m_turboXICMS1,
//            &cosineSim20MS1
//    ); ree;
//
//    QVector<double> mzMeanValsFound;
//    QVector<double> stdMeanValsFound;
//    QVector<double> mzValsSearched;
//    QVector<double> theoApexIntensity;
//    e = calculateEmpiricalMzStats(
//            candidatePeptide,
//            bestPeakIntegrationIndexes,
//            &mzMeanValsFound,
//            &stdMeanValsFound,
//            &mzValsSearched,
//            &theoApexIntensity
//            ); ree;
//
//    cosineSimsIndividual100.resize(mzValsSearched.size());
//    cosineSimsIndividual45.resize(mzValsSearched.size());
//    cosineSimsIndividual20.resize(mzValsSearched.size());
//    intensityApexVals100.resize(mzValsSearched.size());
//
//    const double bestCosineSimSum100 = std::accumulate(cosineSimsIndividual100.begin(), cosineSimsIndividual100.end(), 0.0);
//    const double bestCosineSimSum45 = std::accumulate(cosineSimsIndividual45.begin(), cosineSimsIndividual45.end(), 0.0);
//    const double bestCosineSimSum20 = std::accumulate(cosineSimsIndividual20.begin(), cosineSimsIndividual20.end(), 0.0);
//
//    scoredCandidate->peptideStringWithMods = candidatePeptide.peptideStringWithMods;
//    scoredCandidate->cosineSimSum100 = bestCosineSimSum100;
//    scoredCandidate->cosineSimSum45 = bestCosineSimSum45;
//    scoredCandidate->cosineSimSum20 = bestCosineSimSum20;
//    scoredCandidate->isDecoy = candidatePeptide.isDecoy;
//    scoredCandidate->charge = candidatePeptide.charge;
//    scoredCandidate->mass = candidatePeptide.mass;
//    scoredCandidate->scanNumber = m_msFrame.scanNumberFromFrameIndex(frameIndexIntensityApex);
//    scoredCandidate->scanTime = m_msFrame.scanTimeFromScanNumber(scoredCandidate->scanNumber);
//    scoredCandidate->scanIonCount = m_msFrame.getScanPointsByScanNumber(scoredCandidate->scanNumber).size();
//
//    scoredCandidate->scanTimePredicted = m_fragPredsPredictedScanTime.isEmpty()
//            ? -1.0
//            : m_fragPredsPredictedScanTime.value(candidatePeptide.peptideStringWithMods);
//
//    scoredCandidate->iRTPredicted = candidatePeptide.iRt;
//    scoredCandidate->mzSearchedVec = mzValsSearched;
//    scoredCandidate->theoIntensityVec = theoApexIntensity;
//    scoredCandidate->mzFoundMeanVec = mzMeanValsFound;
//    scoredCandidate->mzFoundStDevVec = stdMeanValsFound;
//    scoredCandidate->intensityFoundMaxVec = intensityApexVals100;
//    scoredCandidate->cosineSimToAnchorVec = cosineSimsIndividual100;
//    scoredCandidate->targetKey = m_uniqueMsInfoScanKey;
//    scoredCandidate->klDivSpectrum = klDivSpectrum;
//    scoredCandidate->cosineSimSpectrum = cosineSimSpectrum;
//    scoredCandidate->cosineSim100MS1 = cosineSim100MS1;
//    scoredCandidate->cosineSim100MS1Iso1 = cosineSim100MS1Iso1;
//    scoredCandidate->cosineSim100MS1Iso2 = cosineSim100MS1Iso2;
//    scoredCandidate->cosineSim45MS1 = cosineSim45MS1;
//    scoredCandidate->cosineSim20MS1 = cosineSim20MS1;
//    scoredCandidate->theoFragmentCount = candidatePeptide.totalFragmentCount;
//    scoredCandidate->peakShapeRatio1 = peakShapeRatios.at(0);
//    scoredCandidate->peakShapeRatio2 = peakShapeRatios.at(1);
//    scoredCandidate->peakShapeRatio3 = peakShapeRatios.at(2);
//    scoredCandidate->b2Corr = cosineSimsIndividualB2B3.at(0);
//    scoredCandidate->b3Corr = cosineSimsIndividualB2B3.at(1);
//    scoredCandidate->b2b3CosineSimSum = cosineSimSumB2B3;
//    scoredCandidate->shadowsCosineSimSum = cosineSimSumShadows;
//    scoredCandidate->cosineSimShadowsToAnchorVec = cosineSimsIndividualShadows;

    ERR_RETURN
}
