//
// Created by anichols on 9/5/23.
//

#include "CandidateProcessertron.h"

#include "BiophysicalCalcs.h"
#include "EigenKernelUtils.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsUtils.h"
#include "ParallelUtils.h"


CandidateProcessertron::CandidateProcessertron() : m_peakWidthMin(3) {}

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
Err CandidateProcessertron::init(
        const PythiaParameters &pythiaParameters,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
        const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
        const MsFrame &msFrame,
        const MsFrame &msFrameMS1,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const QMap<MzHashed, FrequencyPercent> &fragmentFrequencies,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;
    e = ErrorUtils::isNotEmpty(mzHashedVsIonPresence); ree;
    e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;
    e = ErrorUtils::isTrue(msFrame.isValid()); ree;
    e = ErrorUtils::isTrue(msFrameMS1.isValid()); ree;
    e = ErrorUtils::isNotEmpty(fragmentFrequencies); ree;
    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;

    m_pythiaParameters = pythiaParameters;
    m_mzHashedVsXICPoints = mzHashedVsXICPoints;
    m_mzHashedVsIonPresence = mzHashedVsIonPresence;
    m_scanNumberVsScanTime = scanNumberVsScanTime;
    m_msFrame = msFrame;
    m_msFrameMS1 = msFrameMS1;
    m_fragmentFrequencies = fragmentFrequencies;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;

    const PeakIntegratomaticParameters ffParams = buildPeakIntegratomaticParams(m_pythiaParameters);
    e = m_peakIntegratomatic.init(ffParams); ree;

    e = m_turboXICMS1.init(m_msFrameMS1.frameIndexVsScanPoints()); ree;

    ERR_RETURN
}

Err CandidateProcessertron::init(
        const PythiaParameters &pythiaParameters,
        const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
        const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
        const MsFrame &msFrame,
        const MsFrame &msFrameMS1,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        const QMap<MzHashed, FrequencyPercent> &fragmentFrequencies,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QMap<PeptideStringWithMods, ScanTime> &fragPredsPredictedScanTime
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(fragPredsPredictedScanTime);
    m_fragPredsPredictedScanTime = fragPredsPredictedScanTime;

    e = init(
            pythiaParameters,
            mzHashedVsXICPoints,
            mzHashedVsIonPresence,
            msFrame,
            msFrameMS1,
            scanNumberVsScanTime,
            fragmentFrequencies,
            uniqueMsInfoScanKey
            ); ree;

    qDebug() << "updating rt vals from iRT";

    ERR_RETURN
}

namespace {

    int getMaxRowCount(
            const CandidatePeptide &candidatePeptide,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence
    ) {

        int maxRowCount = 0;
        for (const MS2Ion &ms2Ion : candidatePeptide.ms2Ions) {

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const QVector<double> &ionPresenceVec = mzHashedVsIonPresence.value(mzHashed);

            maxRowCount = std::max(maxRowCount, ionPresenceVec.size());
        }

        return maxRowCount;
    }

    Eigen::MatrixX<double> buildSummingMatrix(
            const CandidatePeptide &candidatePeptide,
            const QMap<MzHashed, QVector<double>> &mzHashedVsIonPresence,
            int topNMs2Ions
    ) {

        const int cols = topNMs2Ions;
        const int rows = getMaxRowCount(
                candidatePeptide,
                mzHashedVsIonPresence
        );

        if (rows == 0) {
            return {};
        }

        Eigen::MatrixX<double> mat(rows, cols);
        mat.setZero();

        for (int i = 0; i < cols; i++) {

            if (i >= candidatePeptide.ms2Ions.size()) {
                break;
            }

            const MS2Ion &ms2Ion = candidatePeptide.ms2Ions.at(i);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);
            const QVector<double> &ionPresenceVec = mzHashedVsIonPresence.value(mzHashed);

            if (ionPresenceVec.isEmpty()) {
                continue;
            }

            const Eigen::VectorX<double> eVec = EigenUtils::convertQVectorToEigenVector(ionPresenceVec);
            mat.col(i) = eVec;
        }

        return mat;
    }

    void filterSummedVecPeakIntegrationsByPeakWidth(
            const QVector<double> &summedMatToVec,
            int summedMzPresenceMin,
            int peakWidthMin,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
    ) {

        const auto terminatorLogic = [summedMatToVec, summedMzPresenceMin, peakWidthMin](const PeakIntegrationIndexes &pii){

            const int peakWidth = pii.second - pii.first;
            if (peakWidth < peakWidthMin) {
                return true;
            }

            const QVector<double> summedMatVecMax = summedMatToVec.mid(pii.first, peakWidth);
            const double summedPresenceIntegrationMax = *std::max_element(summedMatVecMax.begin(), summedMatVecMax.end());

            return summedPresenceIntegrationMax < summedMzPresenceMin;
        };

        const auto terminator = std::remove_if(peakIntegrationIndexes->begin(), peakIntegrationIndexes->end(), terminatorLogic);
        peakIntegrationIndexes->erase(terminator, peakIntegrationIndexes->end());
    }

}//namespace
Err CandidateProcessertron::processCandidate(
        const CandidatePeptide &candidatePeptide,
        ScoredCandidate *scoredCandidate
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_mzHashedVsXICPoints); ree;
    e = ErrorUtils::isNotEmpty(m_mzHashedVsIonPresence); ree

    const Eigen::MatrixX<double> presenceMatrix = buildSummingMatrix(
            candidatePeptide,
            m_mzHashedVsIonPresence,
            m_pythiaParameters.topNMs2Ions
    );

    if (presenceMatrix.rows() == 0) {
        ERR_RETURN
    }

    const Eigen::VectorX<double> summedMatVec = presenceMatrix.rowwise().sum();
    const QVector<double> summedMatVecToVec = EigenUtils::convertEigenVectorToQVector(summedMatVec);

    QVector<PeakIntegrationIndexes> peakIntegrationIndexes;

    if (m_fragPredsPredictedScanTime.isEmpty()) {

        e = findCandidateIntegrations(
                summedMatVecToVec,
                &peakIntegrationIndexes
        ); ree;

    }
    else {

        const ScanTime scanTimeWindowTolerance = 5.0; // TODO auto set this in the calibration step.
        const ScanTime scanTimePredicted = m_fragPredsPredictedScanTime.value(candidatePeptide.peptideStringWithMods);
        e = findCandidateIntegrations(
                summedMatVecToVec,
                scanTimePredicted,
                scanTimeWindowTolerance,
                &peakIntegrationIndexes
        ); ree;

    }

    if (peakIntegrationIndexes.isEmpty()) {
        ERR_RETURN
    }

    e = buildScores(
            candidatePeptide,
            peakIntegrationIndexes,
            summedMatVecToVec,
            scoredCandidate
            ); ree;

    ERR_RETURN
}

namespace {

    void sortPeakIntegrationsDescMaxSumFound(
            const QVector<double> &summedMatToVec,
            QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
    ) {

        const auto sortLogic = [summedMatToVec](const PeakIntegrationIndexes &l, const PeakIntegrationIndexes &r){

            const int peakWidthL = l.second - l.first;
            const int peakWidthR = r.second - r.first;

            const QVector<double> summedMatVecMaxL = summedMatToVec.mid(l.first, peakWidthL);
            const QVector<double> summedMatVecMaxR = summedMatToVec.mid(r.first, peakWidthR);

            const double summedPresenceIntegrationMaxL = *std::max_element(summedMatVecMaxL.begin(), summedMatVecMaxL.end());
            const double summedPresenceIntegrationMaxR = *std::max_element(summedMatVecMaxR.begin(), summedMatVecMaxR.end());

            if (MathUtils::tSame(summedPresenceIntegrationMaxL, summedPresenceIntegrationMaxR)) {
                return std::accumulate(summedMatVecMaxL.begin(), summedMatVecMaxL.end(), 0.0)
                       < std::accumulate(summedMatVecMaxR.begin(), summedMatVecMaxR.end(), 0.0);
            }

            return summedPresenceIntegrationMaxL < summedPresenceIntegrationMaxR;
        };

        std::sort(peakIntegrationIndexes->rbegin(), peakIntegrationIndexes->rend(), sortLogic);
    }

}//namespace
Err CandidateProcessertron::findCandidateIntegrations(
        const QVector<double> &summedMatToVec,
        QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
        ) {

    ERR_INIT

    QVector<double> summedPresenceVecSmoothed;
    e = m_peakIntegratomatic.findAllPeaksLimitsInXIC(
            summedMatToVec,
            peakIntegrationIndexes,
            &summedPresenceVecSmoothed
    ); ree;

    filterSummedVecPeakIntegrationsByPeakWidth(
            summedMatToVec,
            m_pythiaParameters.minFoundMzPeaks,
            m_peakWidthMin,
            peakIntegrationIndexes
    );

    sortPeakIntegrationsDescMaxSumFound(
            summedMatToVec,
            peakIntegrationIndexes
    );

    ERR_RETURN
}

Err CandidateProcessertron::findCandidateIntegrations(
        const QVector<double> &summedMatToVec,
        double scanTime,
        double scanTimeTolerance,
        QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
) {

    ERR_INIT

    e = findCandidateIntegrations(
            summedMatToVec,
            peakIntegrationIndexes
            ); ree;

    e = filterSummedVecPeakIntegrationsByPredictedScanTime(
            scanTime,
            scanTimeTolerance,
            peakIntegrationIndexes
            ); ree;

    ERR_RETURN
}

Err CandidateProcessertron::filterSummedVecPeakIntegrationsByPredictedScanTime(
        ScanTime scanTime,
        double windowWidthMinutes,
        QVector<PeakIntegrationIndexes> *peakIntegrationIndexes
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_scanNumberVsScanTime); ree;
    e = ErrorUtils::isTrue(m_msFrame.isValid()); ree;

    const auto terminatorLogic = [&](const PeakIntegrationIndexes &pii){

        const ScanNumber scanNumberStart = m_msFrame.scanNumberFromFrameIndex(pii.first);
        const ScanNumber scanNumberEnd = m_msFrame.scanNumberFromFrameIndex(pii.second);
        const ScanTime scanTimeStart = m_scanNumberVsScanTime.value(scanNumberStart);
        const ScanTime scanTimeEnd = m_scanNumberVsScanTime.value(scanNumberEnd);
        const ScanTime meanScanTime = (scanTimeStart + scanTimeEnd) / 2.0;

        const ScanTime scanTimePredStart = scanTime - windowWidthMinutes;
        const ScanTime scanTimePredEnd = scanTime + windowWidthMinutes;

        return meanScanTime < scanTimePredStart || meanScanTime > scanTimePredEnd;
    };

    const auto terminator = std::remove_if(
            peakIntegrationIndexes->begin(),
            peakIntegrationIndexes->end(),
            terminatorLogic
    );

    peakIntegrationIndexes->erase(terminator, peakIntegrationIndexes->end());

    ERR_RETURN
}

namespace {

    Eigen::MatrixX<double> buildIntensityVecMatrix(
            const CandidatePeptide &candidatePeptide,
            const QMap<MzHashed, XICPoints> &mzHashedVsXICPoints,
            int vecSize,
            int topNMzFrags
    ) {

        const int cols = topNMzFrags;
        const int rows = vecSize;

        Eigen::MatrixX<double> mat(rows, cols);
        mat.setZero();

        for (int i = 0; i < cols; i++) {

            if (i >= candidatePeptide.ms2Ions.size()) {
                break;
            }

            const MS2Ion &ms2Ion = candidatePeptide.ms2Ions.at(i);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

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
        e = ErrorUtils::isTrue(_intensityMatrixIntegratedLimits.maxCoeff() > 0);

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

    Err calcBestCosineSimAndKLDivSum(
            const Eigen::MatrixX<double> &intensityMatrixIntegratedLimits,
            QVector<double> *bestCosineSimsIndividual,
            double *bestCosineSimSum,
            QVector<double> *bestKLDivIndividual,
            double *bestKLDivSum,
            QVector<int> *bestFrameIndexMaxDiffFromAnchorVec,
            Eigen::VectorX<double> *bestAnchorColumn
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrixIntegratedLimits.rows() > 0); ree;

        *bestCosineSimSum = 0.0;
        *bestKLDivSum = std::numeric_limits<double>::max();

        const Eigen::MatrixX<double> intensityMatrixIntegratedLimitsSmoothed
                = applyGaussSmoothRowWiseToMatrix(intensityMatrixIntegratedLimits);

        for (int i = 0; i < intensityMatrixIntegratedLimitsSmoothed.cols(); i++) {

            QVector<double> bestCosineSimsIndividualAnchor;
            QVector<double> bestKLDivIndividualAnchor;
            QVector<int> frameIndexMaxDiffFromAnchorVec;
            const Eigen::VectorX<double> anchorColumn = intensityMatrixIntegratedLimitsSmoothed.col(i);

            const QPair<int, double> anchorFrameIndexMaxVsVal = EigenUtils::returnTopIndexAndValue(anchorColumn);

            for (int j = 0; j < intensityMatrixIntegratedLimitsSmoothed.cols(); j++) {

                const Eigen::VectorX<double> altColumn = intensityMatrixIntegratedLimitsSmoothed.col(j);
                const QPair<int, double> altColumnFrameIndexMaxVsVal = EigenUtils::returnTopIndexAndValue(altColumn);

                frameIndexMaxDiffFromAnchorVec.push_back(anchorFrameIndexMaxVsVal.first - altColumnFrameIndexMaxVsVal.first);

                const double cosineSim = EigenUtils::cosineSimilarity(anchorColumn, altColumn);
                const double klDiv = EigenUtils::klDivergence(anchorColumn, altColumn);

                bestCosineSimsIndividualAnchor.push_back(cosineSim);
                bestKLDivIndividualAnchor.push_back(klDiv);
            }

            const double cosineSimSumAnchor
                    = std::accumulate(bestCosineSimsIndividualAnchor.begin(), bestCosineSimsIndividualAnchor.end(), 0.0);

            const double klDivSumAnchor
                    = std::accumulate(bestKLDivIndividualAnchor.begin(), bestKLDivIndividualAnchor.end(), 0.0);

            if (cosineSimSumAnchor > *bestCosineSimSum) {
                *bestCosineSimSum = cosineSimSumAnchor;
                *bestCosineSimsIndividual = bestCosineSimsIndividualAnchor;
                *bestKLDivSum = klDivSumAnchor;
                *bestKLDivIndividual = bestKLDivIndividualAnchor;
                *bestFrameIndexMaxDiffFromAnchorVec = frameIndexMaxDiffFromAnchorVec;
                *bestAnchorColumn = anchorColumn;
            }

        }

        ERR_RETURN
    }

    Err calculateIndividualPeakPointCount(
            const Eigen::MatrixX<double> &intensityMatrixIntegratedLimitsNoInterference,
            QVector<int> *individualPeakPointCount
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrixIntegratedLimitsNoInterference.size() > 0); ree;

        Eigen::MatrixX<double> mat = intensityMatrixIntegratedLimitsNoInterference;
        mat = (mat.array() > 0.0).select(1.0, mat);

        const Eigen::VectorX<int> colSumVec = mat.colwise().sum().cast<int>();
        *individualPeakPointCount = EigenUtils::convertEigenVectorToQVector(colSumVec);

        ERR_RETURN
    }

    Err calculateCandidateAllignementMetrics(
            const Eigen::MatrixX<double> &intensityMatrix,
            const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
            const QVector<double> &summedMatToVec,
            int topNMs2FragPeaks,
            QVector<double> *cosineSimsIndividual,
            QVector<double> *klDivsIndividual,
            FrameIndex *frameIndexIntensityApex,
            PeakIntegrationIndexes *bestPeakIntegrationIndexes,
            QVector<int> *bestFrameIndexMaxDiffFromAnchorVec,
            QVector<int> *individualPeakPointCount,
            Eigen::VectorX<double> *bestAnchorColumn
    ) {

        ERR_INIT

        e = ErrorUtils::isTrue(intensityMatrix.nonZeros() > 1); ree;
        e = ErrorUtils::isNotEmpty(peakIntegrationIndexes); ree;
        e = ErrorUtils::isNotEmpty(summedMatToVec); ree;

        double bestCosineSimSum = 0.0;

        for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {

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

            const QVector<double> peakCountsPerRow = summedMatToVec.mid(rowStart, rowCount);
            const FrameIndex frameIndexIntensityApexIntegration = MathUtils::findMaxIndexInVector(peakCountsPerRow) + rowStart;

            Eigen::MatrixX<double> intensityMatrixIntegratedLimitsNoInterference;
            e = removeInterferingPeaksInMatrix(
                    intensityMatrixIntegratedLimits,
                    peakCountsPerRow,
                    &intensityMatrixIntegratedLimitsNoInterference
            ); ree;

            e = calculateIndividualPeakPointCount(
                    intensityMatrixIntegratedLimitsNoInterference,
                    individualPeakPointCount
            ); ree;

            QVector<double> cosineSimsIndividualIntegration;
            double bestCosineSimSumIntegration;
            QVector<double> klDivsIndividualIntegration;
            double bestKlDivSimSumIntegration;
            QVector<int> bestFrameIndexMaxDiffFromAnchorVecIntegration;
            e = calcBestCosineSimAndKLDivSum(
                    intensityMatrixIntegratedLimitsNoInterference,
                    &cosineSimsIndividualIntegration,
                    &bestCosineSimSumIntegration,
                    &klDivsIndividualIntegration,
                    &bestKlDivSimSumIntegration,
                    &bestFrameIndexMaxDiffFromAnchorVecIntegration,
                    bestAnchorColumn
            ); ree;

//#define CHECK_ALIGNMENT
#ifdef CHECK_ALIGNMENT
            std::cout << "FrameIndex " << frameIndexIntensityApexIntegration << std::endl;
            std::cout << intensityMatrixIntegratedLimits << std::endl;
            std::cout << "no interference" << std::endl;
            std::cout << intensityMatrixIntegratedLimitsNoInterference << std::endl;
            for (const double v : cosineSimsIndividualIntegration) {
                std::cout << v << " ";
            }
            std::cout << std::endl;
            std::cout << bestCosineSimSumIntegration << std::endl;
            std::cout << pii.first << ", " << pii.second << std::endl;
#endif

            if (bestCosineSimSumIntegration > bestCosineSimSum) {
                bestCosineSimSum = bestCosineSimSumIntegration;
                *cosineSimsIndividual = cosineSimsIndividualIntegration;
                *klDivsIndividual = klDivsIndividualIntegration;
                *frameIndexIntensityApex = frameIndexIntensityApexIntegration;
                *bestPeakIntegrationIndexes = pii;
                *bestFrameIndexMaxDiffFromAnchorVec = bestFrameIndexMaxDiffFromAnchorVecIntegration;
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
            const CandidatePeptide &candidatePeptide,
            const PeakIntegrationIndexes &peakIntegrationIndexes,
            QVector<double> *mzMeanValsFound,
            QVector<double> *stdMeanValsFound,
            QVector<double> *mzValsSearched
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzHashedVsXICPoints); ree;
        e = ErrorUtils::isNotEmpty(candidatePeptide.ms2Ions); ree;
        e = ErrorUtils::isFalse(MathUtils::tSame(peakIntegrationIndexes.first, peakIntegrationIndexes.second)); ree

        for (const MS2Ion &ms2Ion : candidatePeptide.ms2Ions) {

            mzValsSearched->push_back(ms2Ion.mz);

            const MzHashed mzHashed = MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION);

            if (!mzHashedVsXICPoints.contains(mzHashed)) {
                mzMeanValsFound->push_back(0.0);
                stdMeanValsFound->push_back(0.0);
            }

            const XICPoints &xicPoints = mzHashedVsXICPoints.value(mzHashed);

            QVector<double> mzValsForMz;
            for (int i = peakIntegrationIndexes.first; i < peakIntegrationIndexes.second; i++) {

                if (i >= peakIntegrationIndexes.second) {
                    const double mzMean = MathUtils::mean(mzValsForMz);
                    const double mzStd = MathUtils::stDev(mzValsForMz);
                    mzMeanValsFound->push_back(mzMean);
                    stdMeanValsFound->push_back(mzStd);
                    break;
                }

                const QVector<double> &mzVals = xicPoints.scanNumberVsMzVals.value(i);
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

        for (auto it = xicPoints.scanNumbersVsIntensityVals.begin(); it != xicPoints.scanNumbersVsIntensityVals.end(); it++) {
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


}//namespace
Err CandidateProcessertron::buildScores(
        const CandidatePeptide &candidatePeptide,
        const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
        const QVector<double> &summedMatVecToVec,
        ScoredCandidate *scoredCandidate
        ) {
    ERR_INIT

    const Eigen::MatrixX<double> intensityMatrix = buildIntensityVecMatrix(
            candidatePeptide,
            m_mzHashedVsXICPoints,
            summedMatVecToVec.size(),
            m_pythiaParameters.topNMs2Ions
    );

    QVector<double> cosineSimsIndividual;
    QVector<double> klDivsIndividual;
    FrameIndex frameIndexIntensityApex;
    PeakIntegrationIndexes bestPeakIntegrationIndexes;
    QVector<int> bestFrameIndexMaxDiffFromAnchorVec;
    QVector<int> individualPeakPointCount;
    Eigen::VectorX<double> bestAnchorColumn;
    e = calculateCandidateAllignementMetrics(
            intensityMatrix,
            peakIntegrationIndexes,
            summedMatVecToVec,
            m_pythiaParameters.topNMs2Ions,
            &cosineSimsIndividual,
            &klDivsIndividual,
            &frameIndexIntensityApex,
            &bestPeakIntegrationIndexes,
            &bestFrameIndexMaxDiffFromAnchorVec,
            &individualPeakPointCount,
            &bestAnchorColumn
    ); ree;

    QVector<double> intensityApexVals = EigenUtils::convertEigenVectorToQVector(
            Eigen::VectorX<double>(intensityMatrix.row(frameIndexIntensityApex))
    );

    double cosineSimSpectrum;
    double klDivSpectrum;
    e = calculateSpectrumMetrics(
            intensityMatrix.row(frameIndexIntensityApex),
            candidatePeptide.ms2Ions,
            &cosineSimSpectrum,
            &klDivSpectrum
    ); ree;

    double cosineSimMS1;
    e = calculateMS1Corr(
            bestAnchorColumn,
            bestPeakIntegrationIndexes,
            BiophysicalCalcs::calculateThomsonFromMass(candidatePeptide.mass, candidatePeptide.charge),
            m_pythiaParameters.ms2ExtractionWidthPPM,
            &m_turboXICMS1,
            &cosineSimMS1
    ); ree;

    QVector<double> mzMeanValsFound;
    QVector<double> stdMeanValsFound;
    QVector<double> mzValsSearched;
    QVector<double> ppmDifference;
    QVector<double> theoApexIntensity;
    e = calculateEmpiricalMzStats(
            candidatePeptide,
            bestPeakIntegrationIndexes,
            &mzMeanValsFound,
            &stdMeanValsFound,
            &mzValsSearched,
            &ppmDifference,
            &theoApexIntensity
            ); ree;

    QVector<double> fragmentFrequencies;
    std::transform(
            candidatePeptide.ms2Ions.begin(),
            candidatePeptide.ms2Ions.end(),
            std::back_inserter(fragmentFrequencies),
            [&](const MS2Ion &ms2Ion){
                return m_fragmentFrequencies.value(MathUtils::hashDecimal(ms2Ion.mz, S_GLOBAL_SETTINGS.HASHING_PRECISION));
            }
    );

    cosineSimsIndividual.resize(mzValsSearched.size());
    klDivsIndividual.resize(mzValsSearched.size());
    intensityApexVals.resize(mzValsSearched.size());

    const double bestCosineSimSum = std::accumulate(cosineSimsIndividual.begin(), cosineSimsIndividual.end(), 0.0);
    const double bestKlDivSimSum = std::accumulate(klDivsIndividual.begin(), klDivsIndividual.end(), 0.0);

    double xCorr;
    e = MsUtils::calcXCorr(
            mzMeanValsFound,
            mzValsSearched,
            intensityApexVals,
            theoApexIntensity,
            m_pythiaParameters.ms2ExtractionWidthPPM,
            &xCorr
    ); ree;

    scoredCandidate->peptideStringWithMods = candidatePeptide.peptideStringWithMods;
    scoredCandidate->cosineSimSum = bestCosineSimSum;
    scoredCandidate->isDecoy = candidatePeptide.isDecoy;
    scoredCandidate->charge = candidatePeptide.charge;
    scoredCandidate->mass = candidatePeptide.mass;
    scoredCandidate->scanNumber = m_msFrame.scanNumberFromFrameIndex(frameIndexIntensityApex);
    scoredCandidate->scanTime = m_msFrame.scanTimeFromScanNumber(scoredCandidate->scanNumber);
    scoredCandidate->scanIonCount = m_msFrame.getScanPointsByScanNumber(scoredCandidate->scanNumber).size();
    scoredCandidate->xCorr = xCorr;

    scoredCandidate->scanTimePredicted = m_fragPredsPredictedScanTime.isEmpty()
            ? -1.0
            : m_fragPredsPredictedScanTime.value(candidatePeptide.peptideStringWithMods);

    scoredCandidate->iRTPredicted = candidatePeptide.iRt;
    scoredCandidate->mzSearchedVec = mzValsSearched;
    scoredCandidate->theoIntensityVec = theoApexIntensity;
    scoredCandidate->mzFoundMeanVec = mzMeanValsFound;
    scoredCandidate->mzFoundStDevVec = stdMeanValsFound;
    scoredCandidate->intensityFoundMaxVec = intensityApexVals;
    scoredCandidate->frameIndexMaxDiffFromAnchorVec = bestFrameIndexMaxDiffFromAnchorVec;
    scoredCandidate->cosineSimToAnchorVec = cosineSimsIndividual;
    scoredCandidate->klDivToAnchorVec = klDivsIndividual;
    scoredCandidate->peakPointCountFoundVec = individualPeakPointCount;
    scoredCandidate->fragmentFrequencyVec = fragmentFrequencies;
    scoredCandidate->targetKey = m_uniqueMsInfoScanKey;
    scoredCandidate->klDivSum = bestKlDivSimSum;
    scoredCandidate->klDivSpectrum = klDivSpectrum;
    scoredCandidate->cosineSimSpectrum = cosineSimSpectrum;
    scoredCandidate->cosineSimMS1 = cosineSimMS1;
    scoredCandidate->theoFragmentCount = candidatePeptide.totalFragmentCount;

    ERR_RETURN
}

Err CandidateProcessertron::calculateEmpiricalMzStats(
        const CandidatePeptide &candidatePeptide,
        PeakIntegrationIndexes &bestPeakIntegrationIndexes,
        QVector<double> *mzMeanValsFound,
        QVector<double> *stdMeanValsFound,
        QVector<double> *mzValsSearched,
        QVector<double> *ppmDifference,
        QVector<double> *theoApexIntensity
        ) {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(
            bestPeakIntegrationIndexes.second,
            bestPeakIntegrationIndexes.first,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;

    e = calculateMassAccuracyMetrics(
            m_mzHashedVsXICPoints,
            candidatePeptide,
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
            mzFoundVsMzSearched.begin(),
            mzFoundVsMzSearched.end(),
            std::back_inserter(*ppmDifference),
            [](const QPair<double, double> &p){return 1e6 * (p.first - p.second) / p.second;}
    );

    std::transform(
            candidatePeptide.ms2Ions.begin(),
            candidatePeptide.ms2Ions.end(),
            std::back_inserter(*theoApexIntensity),
            [](const MS2Ion &ms2Ion){return ms2Ion.intensity;}
    );

    ERR_RETURN
}


