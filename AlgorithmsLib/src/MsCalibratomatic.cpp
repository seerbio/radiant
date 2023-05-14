//
// Created by anichols on 4/16/23.
//

#include "MsCalibratomatic.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoretronProcessormatic.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"
#include "PSMsReader.h"
#include "TandemFragmentPredictotron.h"

#include <Eigen/Dense>

#include <QtConcurrent/QtConcurrent>


MsCalibratomatic::MsCalibratomatic()
: m_calPointK(3)
, m_stDevNew(-1.0)
{}

Err MsCalibratomatic::init(
        const PythiaParameters &pythiaParameters,
        const QString &firstPassSearchFilePath,
        int calPointK
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    e = ErrorUtils::isAboveThreshold(
            calPointK,
            0,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;
    m_calPointK = calPointK;

    e = ErrorUtils::fileExists(firstPassSearchFilePath); ree;
    m_firstPassSearchFilePath = firstPassSearchFilePath;

    e = buildCalibrator(); ree;

    qDebug() << "NearestNeighbors Treesize" << m_nnSearch.kdTreeSize();

    ERR_RETURN
}

namespace{

    void sortScoreDesc(QVector<PSMsReaderRow> *psmsReaderRows) {

        const auto sortLogic = [](const PSMsReaderRow &l, const PSMsReaderRow &r){
            return l.score < r.score;
        };

        std::sort(psmsReaderRows->rbegin(), psmsReaderRows->rend(), sortLogic);
    }

    Err buildNearestNeigbhorsDataMatrix(
            const QString &firstPassSearchFilePath,
            Eigen::MatrixX<double> *mat
            ) {

        ERR_INIT

        QVector<PSMsReaderRow> psmsReaderRows;
        e = ParquetReader::read(
                firstPassSearchFilePath,
                &psmsReaderRows
        ); ree;
        sortScoreDesc(&psmsReaderRows);

        QVector<double> scanNumbers;
        QVector<double> mzFounds;
        QVector<double> ppms;

        int decoyCount = 0;
        int rowCounter = 0;
        const double thresholdFDR = 0.01;
        for(const PSMsReaderRow &row : psmsReaderRows) {

            if (row.isDecoy) {
                decoyCount++;
            }

            if (static_cast<double>(decoyCount) / ++rowCounter >= thresholdFDR) {
                break;
            }

//            const QVector<double> &mzFound = row.mzFound;
//            const QVector<double> &mzSearched = row.mzSearched;
//
//            e = ErrorUtils::isEqual(mzFound.size(), mzSearched.size()); ree;
//
//            const auto scanNumberDouble = static_cast<double>(row.scanNumber);
//
//            for (int i = 0; i < mzFound.size(); i++) {
//
//                const double mzSrch = mzSearched.at(i);
//                const double mzFnd = mzFound.at(i);
//
//                e = ErrorUtils::isFalse(MathUtils::tZero(mzSrch)); ree;
//
//                const double ppm = 1e6 * ((mzFnd - mzSrch) / mzSrch);
//
//                scanNumbers.push_back(scanNumberDouble);
//                mzFounds.push_back(mzFnd);
//                ppms.push_back(ppm);
//            }

        }

        e = ErrorUtils::isEqual(scanNumbers.size(), mzFounds.size());
        e = ErrorUtils::isEqual(scanNumbers.size(), ppms.size());

        const int colCount = 3;
        mat->resize(mzFounds.size(), colCount);

        const Eigen::VectorX<double> scanNumbersVec = EigenUtils::convertQVectorToEigenVector(scanNumbers);
        const Eigen::VectorX<double> mzFoundsVec = EigenUtils::convertQVectorToEigenVector(mzFounds);
        const Eigen::VectorX<double> ppmsVec = EigenUtils::convertQVectorToEigenVector(ppms);

        mat->col(0) = scanNumbersVec;
        mat->col(1) = mzFoundsVec;
        mat->col(2) = ppmsVec;

        ERR_RETURN
    }

    Err filterDataMatrix(Eigen::MatrixX<double> *dataMatrix) {

        ERR_INIT

        const int expectedColSize = 3;

        e = ErrorUtils::isTrue(dataMatrix->rows()); ree;
        e = ErrorUtils::isEqual(
                static_cast<int>(dataMatrix->cols()),
                expectedColSize
                ); ree;

        const Eigen::VectorX<double> ppmColumn = dataMatrix->col(2);

        const double ppmMean = EigenUtils::calculateMeanOfVector(ppmColumn);
        const double ppmStDev = EigenUtils::calculateStDevOfVector(ppmColumn);

        qDebug() << "OG ppm:stDev" << ppmMean << ppmStDev;

        const double absThresholdPPM = std::abs(ppmMean * ppmStDev);

        const int thresholdColIdx = 2;
        EigenUtils::removeRowsAboveOrBelowThreshold(
                absThresholdPPM,
                thresholdColIdx,
                EigenUtils::ThresholderDirection::Above,
                dataMatrix
                );

        EigenUtils::removeRowsAboveOrBelowThreshold(
                -absThresholdPPM,
                thresholdColIdx,
                EigenUtils::ThresholderDirection::Below,
                dataMatrix
        );

        ERR_RETURN
    }


    Err buildNNInput(
            const Eigen::MatrixX<double> &mat,
            QVector<QPair<double, Coors>> *valuesVsTreePoints
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(mat.rows() > 0); ree;

        for (int row = 0; row < mat.rows(); row++) {

            const Eigen::VectorX<double> matRowQVec = mat.row(row);
            QVector<double> matRow = EigenUtils::convertEigenVectorToQVector(matRowQVec);

            const double ppm = matRow.back();
            matRow.pop_back();

            valuesVsTreePoints->push_back({ppm, matRow});
        }

        ERR_RETURN
    }

    Err calculateNewAccuracyMetrics(
            Eigen::MatrixX<double> dataMatrix,
            int calPointK,
            NearestNeighborsSearch *nnSearch,
            double *stDevNew
            ) {

        ERR_INIT

        e =  ErrorUtils::isTrue(dataMatrix.size() > 0); ree;

        QVector<QPair<DiffPPM, Coors>> valuesVsTreePoints;

        const int yValIdx = static_cast<int>(dataMatrix.cols()) - 1;
        const Eigen::VectorX<double> vec = dataMatrix.col(yValIdx);
        const QVector<double> ppmDiffVals = EigenUtils::convertEigenVectorToQVector(vec);

        buildNNInput(
                dataMatrix,
                &valuesVsTreePoints
        );

        const double ppmMean = MathUtils::mean(ppmDiffVals);
        const double stDev = MathUtils::stDev(ppmDiffVals);
        qDebug() << "OG Mean / StDev" << ppmMean << stDev;

        const QPair<QVector<DiffPPM>, QVector<QVector<double>>> unzippedTreePoints
                = ParallelUtils::unZip(valuesVsTreePoints);

        const int calPointsForMetricsBecauseFirstResultIsInTrainingSetWillBePopped = calPointK + 1;
        QVector<NNSearchResult> searchResults;
        e = nnSearch->kNearestNeighborsSearch(
                unzippedTreePoints.second,
                calPointsForMetricsBecauseFirstResultIsInTrainingSetWillBePopped,
                &searchResults
        ); ree;

        QVector<double> adjustedPPMs;
        for (int i = 0; i < searchResults.size(); i++) {

            const NNSearchResult &sr = searchResults.at(i);
            const double ogPPM = ppmDiffVals.at(i);

            const double ppmCorrectionMean = sr.values;

            adjustedPPMs.push_back(ogPPM - ppmCorrectionMean);
        }

        const double ppmMeanNew = MathUtils::mean(adjustedPPMs);
        *stDevNew = MathUtils::stDev(adjustedPPMs);
        qDebug() << "New Mean / StDev" << ppmMeanNew << *stDevNew;

        ERR_RETURN
    }

}//namespace
Err MsCalibratomatic::buildCalibrator() {

    ERR_INIT

    e = ErrorUtils::isTrue(m_params.isValid()); ree;
    e = ErrorUtils::fileExists(m_firstPassSearchFilePath); ree;

    Eigen::MatrixX<double> dataMatrix;
    e = buildNearestNeigbhorsDataMatrix(
            m_firstPassSearchFilePath,
            &dataMatrix
            ); ree;

    e = filterDataMatrix(&dataMatrix); ree;

    const double trainingFraction = 0.9;
    const int trainSize = static_cast<int>(dataMatrix.rows() * trainingFraction);
    Eigen::MatrixX<double> dataMatrixTrain = dataMatrix.topRows(trainSize);
    Eigen::MatrixX<double> dataMatrixTest = dataMatrix.bottomRows(dataMatrix.rows() - trainSize);
    qDebug() << "Datapoints:" << dataMatrixTrain.rows();

    QVector<QPair<double, Coors>> valuesVsTreePoints;
    buildNNInput(dataMatrixTrain, &valuesVsTreePoints);

    e = m_nnSearch.init(valuesVsTreePoints); ree;

    e = calculateNewAccuracyMetrics(
            dataMatrixTest,
            m_calPointK,
            &m_nnSearch,
            &m_stDevNew
            ); ree;

    ERR_RETURN
}

namespace {

    QVector<QVector<double>> buildNearestNeighborsInput(
            const QMap<FrameIndex, ScanPoints> &indexVsScanPoints
            ) {

        QVector<QVector<double>> nnInputs;

        for (auto it = indexVsScanPoints.begin(); it != indexVsScanPoints.end(); it++) {

            const FrameIndex frameIndex = it.key();
            const ScanPoints &scanPoints = it.value();

            for (const ScanPoint &sp : scanPoints) {
                nnInputs.push_back({static_cast<double>(frameIndex), sp.x()});
            }
        }

        return nnInputs;
    }

    QMap<FrameIndex, QVector<NNSearchResult>> groupNNSearchResultsByFrame(
            const QVector<NNSearchResult> &nnSearchResults
            ) {

        QMap<FrameIndex, QVector<NNSearchResult>> frameIndexVsNNSearchResults;

        for (const NNSearchResult &res : nnSearchResults) {

            const auto frameIndex = static_cast<int>(res.searchedCoor.at(0));
            frameIndexVsNNSearchResults[frameIndex].push_back(res);
        }

        return frameIndexVsNNSearchResults;
    }

    Err correctMzValues(
            const QMap<FrameIndex, ScanPoints> &indexVsScanPoints,
            const QVector<NNSearchResult> &nnSearchResults,
            QMap<FrameIndex, ScanPoints> *correctedIndexVsScanPoints
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(indexVsScanPoints); ree;
        e = ErrorUtils::isNotEmpty(nnSearchResults); ree;

        const QMap<FrameIndex, QVector<NNSearchResult>> frameIndexVsNNSearchResults
                =  groupNNSearchResultsByFrame(nnSearchResults);

        e = ErrorUtils::isEqual(
                indexVsScanPoints.size(),
                frameIndexVsNNSearchResults.size()
                ); ree;

        for (auto it = indexVsScanPoints.begin(); it != indexVsScanPoints.end(); it++) {

            const FrameIndex frameIndex = it.key();
            const ScanPoints &scanPoints = it.value();
            const QVector<NNSearchResult> &nnSearchResultsFrame = frameIndexVsNNSearchResults.value(frameIndex);

            e = ErrorUtils::isEqual(
                    scanPoints.size(),
                    nnSearchResultsFrame.size()
                    ); ree;

            for (int i = 0; i < scanPoints.size(); i++) {

                const ScanPoint &scanPoint = scanPoints.at(i);
                const NNSearchResult &nnSearchResult = nnSearchResultsFrame.at(i);

                const double ogMz = scanPoint.x();
                const double ogIntensity = scanPoint.y();
                const double searchedMz = nnSearchResult.searchedCoor.at(1);
                const int searchedFrameIndex = static_cast<int>(nnSearchResult.searchedCoor.at(0));

                e = ErrorUtils::isTrue(MathUtils::tSame(ogMz, searchedMz)); ree;
                e = ErrorUtils::isEqual(
                        frameIndex,
                        searchedFrameIndex
                        ); ree;

                const double meanDiffPPM = nnSearchResult.values;
                const double mzCorrectionAmount = (ogMz * meanDiffPPM) / 1e6;
                const double correctedMz = ogMz - mzCorrectionAmount;

                (*correctedIndexVsScanPoints)[frameIndex].push_back({correctedMz, ogIntensity});
            }

        }

        ERR_RETURN
    }

}//namespace
Err MsCalibratomatic::recalibratePoints(
        const QMap<FrameIndex, ScanPoints> &indexVsScanPoints,
        QMap<FrameIndex, ScanPoints> *recalIndexVsScanPoints
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(indexVsScanPoints); ree;

    const QVector<QVector<double>> nnInputs
            = buildNearestNeighborsInput(indexVsScanPoints);

    QVector<NNSearchResult> nnSearchResults;
    e = m_nnSearch.kNearestNeighborsSearch(
            nnInputs,
            m_calPointK,
            &nnSearchResults
            ); ree;

    e = correctMzValues(
            indexVsScanPoints,
            nnSearchResults,
            recalIndexVsScanPoints
            ); ree;

    ERR_RETURN
}

double MsCalibratomatic::newStDev() {
    return m_stDevNew;
}
