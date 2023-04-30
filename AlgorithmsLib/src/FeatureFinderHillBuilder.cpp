//
// Created by anichols on 12/10/22.
//

#include "FeatureFinderHillBuilder.h"

#include "CSVReader.h"
#include "ErrorUtils.h"
#include "MathUtils.h"
#include "MsReaderBase.h"
#include "ParallelUtils.h"
#include "PeakIntegratomatic.h"

#include <Eigen/Dense>

#include <QtConcurrent/QtConcurrent>

#include <iostream>


class Q_DECL_HIDDEN FeatureFinderHillBuilder::Private
{
public:

    Private();
    ~Private();

    Err init(const FeatureFinderParameters &featureFinderParameters);

    Err buildScanPointGroups(
            const QVector<ScanPoints> &allScanPoints,
            QVector<QVector<QVector<double>>> *groupedMzVals,
            QVector<QVector<QVector<double>>> *groupedIntensityVals
    );

    Err connectCentroidsInGroupedMzVals(
            const QVector<QVector<QVector<double>>> &groupedMzVals,
            double tolerancePPM,
            QVector<Eigen::MatrixX<int>> *connectedCentroids
    );

    Err buildHills(
            const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber,
            QVector<FeatureFinderHill> *featureFinderHills
    );

    Err refineHills(QVector<FeatureFinderHill> *featureFinderHills);

    void setRunParallel(bool runParallel);

private:

    FeatureFinderParameters m_featureFinderParameters;
    bool m_runParallel;

};


FeatureFinderHillBuilder::Private::Private() : m_runParallel(false) {}
FeatureFinderHillBuilder::Private::~Private() {}


Err FeatureFinderHillBuilder::Private::init(const FeatureFinderParameters &featureFinderParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(featureFinderParameters.isValid(), eValueError); ree;
    m_featureFinderParameters = featureFinderParameters;

    ERR_RETURN
}


namespace {

    struct ScanGroupingParallelInput {
        int skipScanCount = -1;
        bool isLastTranche = false;
        QVector<ScanPoints> allScanPointsTranch;
    };

    struct ScanGroupingParallelOutput {
        Err e = Error::eNoError;
        QVector<QVector<QVector<double>>> groupedMzVals;
        QVector<QVector<QVector<double>>> groupedIntensityVals;
    };

    QVector<ScanGroupingParallelInput> buildScanGroupingParallelInputs(
            const QVector<QVector<ScanPoints>> &tranchedScanPoints,
            int skipScanCount
    ) {

        QVector<ScanGroupingParallelInput> scanGroupingParallelInputs;

        const int minAllowableGroupCount = 2;

        for (const QVector<ScanPoints> &spVec : tranchedScanPoints) {

            if (spVec.size() < minAllowableGroupCount) {
                continue;
            }

            ScanGroupingParallelInput inp;
            inp.allScanPointsTranch = spVec;
            inp.skipScanCount = skipScanCount;
            inp.isLastTranche = spVec.size() < minAllowableGroupCount + skipScanCount;
            scanGroupingParallelInputs.push_back(inp);

            if (inp.isLastTranche) {
                break;
            }
        }

        return scanGroupingParallelInputs;
    }


    ScanGroupingParallelOutput scanGroupingParallelLogic(const ScanGroupingParallelInput &sgpi) {

        ScanGroupingParallelOutput output;

        const int defaultGroupPadding = 2;
        const int nonLastTrancheAdjustment = sgpi.isLastTranche ? 0 : sgpi.skipScanCount;

        for (
                ScanNumberIndex scanIndex = 0;
                scanIndex < sgpi.allScanPointsTranch.size() - (1 + nonLastTrancheAdjustment);
                scanIndex++
        ) {

            const ScanNumberIndex stopNextScanIndex
                = scanIndex + sgpi.skipScanCount + defaultGroupPadding;

            QVector<QVector<double>> mzValsGroup;
            QVector<QVector<double>> intensityValsGroup;

            for (ScanNumberIndex nextScanIndex = scanIndex; nextScanIndex < stopNextScanIndex; ++nextScanIndex) {

                if (nextScanIndex >= sgpi.allScanPointsTranch.size()) {
                    break;
                }

                const ScanPoints &scanPoints = sgpi.allScanPointsTranch.at(nextScanIndex);

                QVector<double> mzVals;
                QVector<double> intensityVals;
                output.e = MsReaderBase::splitScanPoints(
                        scanPoints,
                        &mzVals,
                        &intensityVals
                );

                if (output.e != eNoError) {
                    return output;
                }

                mzValsGroup.push_back(mzVals);
                intensityValsGroup.push_back(intensityVals);
            }

            output.groupedMzVals.push_back(mzValsGroup);
            output.groupedIntensityVals.push_back(intensityValsGroup);
        }

        return output;
    }


}//namespace
Err FeatureFinderHillBuilder::Private::buildScanPointGroups(
        const QVector<ScanPoints> &allScanPoints,
        QVector<QVector<QVector<double>>> *groupedMzVals,
        QVector<QVector<QVector<double>>> *groupedIntensityVals
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(allScanPoints); ree;
    groupedMzVals->clear();
    groupedIntensityVals->clear();

//  While this does speed things up a bit, it is not the most accurate.
//  Speeds things up by about 100 mSec.

//#define GROUP_SCANS_PARALLEL
#ifdef GROUP_SCANS_PARALLEL
    qDebug() << "Running buildScanPointGroups parallel";

    const int processorCount = ParallelUtils::numberOfAvailableSystemProcessors();
    const int defaultExtraGroup = 1;

    QVector<QVector<ScanPoints>> tranchedScanPoints;
    ParallelUtils::tranchVectorForParallelizationInOrder(
            allScanPoints,
            processorCount,
            m_featureFinderParameters.skipScanCount + defaultExtraGroup,
            &tranchedScanPoints
    );

    const QVector<ScanGroupingParallelInput> scanGroupingParallelInputs = buildScanGroupingParallelInputs(
            tranchedScanPoints,
            m_featureFinderParameters.skipScanCount
    );

    QFuture<ScanGroupingParallelOutput> futures = QtConcurrent::mapped(
            scanGroupingParallelInputs,
            scanGroupingParallelLogic
    );

    futures.waitForFinished();

    int counter = 0;
    for (const ScanGroupingParallelOutput &sgpi : futures) {
        e = sgpi.e; ree;
        groupedMzVals->append(sgpi.groupedMzVals);
        groupedIntensityVals->append(sgpi.groupedIntensityVals);
    }

#else
    qDebug() << "Running buildScanPointGroups serial";

    ScanGroupingParallelInput scanGroupingParallelInput;
    scanGroupingParallelInput.skipScanCount = m_featureFinderParameters.skipScanCount;
    scanGroupingParallelInput.allScanPointsTranch = allScanPoints;
    scanGroupingParallelInput.isLastTranche = true;

    const ScanGroupingParallelOutput output = scanGroupingParallelLogic(scanGroupingParallelInput);
    *groupedMzVals = output.groupedMzVals;
    *groupedIntensityVals = output.groupedIntensityVals;
    e = output.e; ree;
#endif

    ERR_RETURN
}


namespace {

    struct ConnectCentroidsInGroupedMzValsInput {
        QVector<QVector<double>> mzValsGroup;
        double tolerancePPM = -1.0;
    };


    QVector<ConnectCentroidsInGroupedMzValsInput> buildConnectCentroidsInGroupedMzValsInputVec(
            const QVector<QVector<QVector<double>>> &groupedMzVals,
            double tolerancePPM
    ) {

        QVector<ConnectCentroidsInGroupedMzValsInput> inputs;

        for (const QVector<QVector<double>> &mzValsGroup : groupedMzVals) {

            ConnectCentroidsInGroupedMzValsInput input;
            input.mzValsGroup = mzValsGroup;
            input.tolerancePPM = tolerancePPM;

            inputs.push_back(input);
        }

        return inputs;
    }


    Eigen::VectorX<int> findCentroidConnections(
            const QVector<double> &referenceScan,
            const QVector<double> &nextScan,
            double tolerancePPM
    ) {

        ScanNumberIndex currentNextScanIndex = 0;

        Eigen::VectorX<int> foundIndexes(referenceScan.size());
        foundIndexes.setOnes();
        foundIndexes *= -1;

        Eigen::VectorX<double> scoredIndexes(referenceScan.size());
        scoredIndexes.setOnes();
        scoredIndexes *= std::numeric_limits<double>::max();

        for (ScanNumberIndex refScanInd = 0; refScanInd < referenceScan.size(); refScanInd++) {

            const double mzRefScan = referenceScan.at(refScanInd);

            for (ScanNumberIndex nextScanInd = currentNextScanIndex; nextScanInd < nextScan.size(); nextScanInd++) {

                const double mzNextScan = nextScan.at(nextScanInd);
                const double mzTolerance = MathUtils::calculatePPM(mzRefScan, tolerancePPM);
                const double absoluteMzDistance = std::abs(mzNextScan - mzRefScan);

                if (absoluteMzDistance < mzTolerance) {

                    if (scoredIndexes.coeff(refScanInd) > absoluteMzDistance) {
                        foundIndexes.coeffRef(refScanInd) = nextScanInd;
                        scoredIndexes.coeffRef(refScanInd) = absoluteMzDistance;
                    }

                    continue;
                }

                else if (mzNextScan < mzRefScan) {
                    continue;
                }

                currentNextScanIndex = nextScanInd;
                break;
            }
        }

        return foundIndexes;
    }


    Eigen::MatrixX<int> connectCentroidsInGroupedMzValsParallelLogic(
            const ConnectCentroidsInGroupedMzValsInput &input
    ) {

        const QVector<double> &referenceScan = input.mzValsGroup.front();

        const int returnMatRows = input.mzValsGroup.size() - 1;
        const int returnMatCols = referenceScan.size();
        Eigen::MatrixX<int> returnMat(returnMatRows, returnMatCols);

        for (ScanNumberIndex nextScanIndex = 1; nextScanIndex < input.mzValsGroup.size(); nextScanIndex++) {

            const QVector<double> &nextScan = input.mzValsGroup.at(nextScanIndex);
            const Eigen::VectorX<int> connectedCentroidsNextScan = findCentroidConnections(
                    referenceScan,
                    nextScan,
                    input.tolerancePPM
            );

            returnMat.row(nextScanIndex - 1) = connectedCentroidsNextScan;
        }

        return returnMat;
    }


}//namespace
Err FeatureFinderHillBuilder::Private::connectCentroidsInGroupedMzVals(
        const QVector<QVector<QVector<double>>> &groupedMzVals,
        double tolerancePPM,
        QVector<Eigen::MatrixX<int>> *connectedCentroids
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(groupedMzVals); ree;

    QVector<ConnectCentroidsInGroupedMzValsInput> inputs = buildConnectCentroidsInGroupedMzValsInputVec(
            groupedMzVals,
            tolerancePPM
    );

    if (m_runParallel) {

        qDebug() << "Running connectCentroidsInGroupedMzVals parallel";

        QFuture<Eigen::MatrixX<int>> futures
                = QtConcurrent::mapped(inputs, connectCentroidsInGroupedMzValsParallelLogic);

        futures.waitForFinished();

        for (const Eigen::MatrixX<int> &mat: futures) {
            connectedCentroids->push_back(mat);
        }
    }

    else {

        qDebug() << "Running connectCentroidsInGroupedMzVals serial";

        for (const ConnectCentroidsInGroupedMzValsInput &input: inputs) {
            const Eigen::MatrixX<int> mat = connectCentroidsInGroupedMzValsParallelLogic(input);

            connectedCentroids->push_back(mat);
        }
    }

    ERR_RETURN
}


namespace {

    void getNextIndex(
            ScanNumberIndex currentScanIndex,
            int currentIndexVal,
            QVector<Eigen::MatrixX<int>> *connectedCentroidsMats,
            ScanNumberIndex *currentScanIndexReturn,
            int *currentIndexValReturn,
            ScanNumberIndex *nextScanIndexReturn,
            int *nextValReturn
    ) {

        int nextScanIndex = currentScanIndex;
        int nextVal = -1;

        if (nextScanIndex >= connectedCentroidsMats->size()) {
            *currentScanIndexReturn = currentScanIndex;
            *currentIndexValReturn = currentIndexVal;
            *nextScanIndexReturn = nextScanIndex;
            *nextValReturn = nextVal;
            return;
        }

        const Eigen::MatrixX<int> &currentCentroids = connectedCentroidsMats->at(nextScanIndex);

        for (int row = 0; row < currentCentroids.rows(); row++) {

            nextVal = currentCentroids.coeff(row, currentIndexVal);
            nextScanIndex += 1;

            if (nextVal < 0) {
                continue;
            }

            for (int rowIndexEraser = 0; rowIndexEraser < currentCentroids.rows(); ++rowIndexEraser) {
                (*connectedCentroidsMats)[currentScanIndex].coeffRef(rowIndexEraser, currentIndexVal) = -1;
            }

            break;
        }

        *currentScanIndexReturn = currentScanIndex;
        *currentIndexValReturn = currentIndexVal;
        *nextScanIndexReturn = nextScanIndex;
        *nextValReturn = nextVal;
    }


    Err groupConnectedCentroidsToHills(
            const QVector<Eigen::MatrixX<int>> &_connectedCentroidsMats,
            const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber,
            int minScanCount,
            QVector<FeatureFinderHill> *featureFinderHills
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(_connectedCentroidsMats); ree;
        featureFinderHills->clear();

        const QVector<ScanNumber> scanNumbers = scanPointsByScanNumber.keys().toVector();

        QVector<Eigen::MatrixX<int>> connectedCentroidsMats = _connectedCentroidsMats;

        const QList<ScanPoints> &scanPoints = scanPointsByScanNumber.values();

        for (ScanNumberIndex baseScanIndex = 0; baseScanIndex < connectedCentroidsMats.size(); baseScanIndex++) {

            const Eigen::MatrixX<int> &baseCentroids = connectedCentroidsMats.at(baseScanIndex);

            for (int baseCentroidsColInd = 0; baseCentroidsColInd< baseCentroids.cols(); ++baseCentroidsColInd) {

                for (int baseCentroidsRowInd = 0; baseCentroidsRowInd < baseCentroids.rows(); ++baseCentroidsRowInd) {

                    const int nextScanCentroidVal = baseCentroids.coeff(baseCentroidsRowInd, baseCentroidsColInd);
                    if (nextScanCentroidVal < 0) {
                        continue;
                    }

                    ScanPoint currentScanPoint = scanPoints.at(baseScanIndex).at(baseCentroidsColInd);

                    int currentIndexVal = baseCentroidsColInd;
                    ScanNumberIndex currentScanIndex = baseScanIndex;
                    int nextIndexVal = nextScanCentroidVal;
                    ScanNumberIndex nextScanIndex = currentScanIndex + 1 + baseCentroidsRowInd;

                    FeatureFinderHill featureFinderHill;
                    featureFinderHill.addPoint(
                            currentScanIndex,
                            scanNumbers.at(currentScanIndex),
                            currentScanPoint.x(),
                            currentScanPoint.y()
                    );

                    while (nextIndexVal != -1) {

                        getNextIndex(
                                nextScanIndex,
                                nextIndexVal,
                                &connectedCentroidsMats,
                                &currentScanIndex,
                                &currentIndexVal,
                                &nextScanIndex,
                                &nextIndexVal
                        );

                        if (nextIndexVal == -1) {
                            break;
                        }

                        currentScanPoint = scanPoints.at(currentScanIndex).at(currentIndexVal);
                        featureFinderHill.addPoint(
                                currentScanIndex,
                                scanNumbers.at(currentScanIndex),
                                currentScanPoint.x(),
                                currentScanPoint.y()
                        );
                    }

                    currentScanPoint = scanPoints.at(currentScanIndex).at(currentIndexVal);
                    featureFinderHill.addPoint(
                            currentScanIndex,
                            scanNumbers.at(currentScanIndex),
                            currentScanPoint.x(),
                            currentScanPoint.y()
                    );

                    if (featureFinderHill.scanCount() >= minScanCount) {
                        featureFinderHills->push_back(featureFinderHill);
                    }

                    break;
                }
            }
        }

        ERR_RETURN
    }


}//namespace
Err FeatureFinderHillBuilder::Private::buildHills(
        const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber,
        QVector<FeatureFinderHill> *featureFinderHills
) {

    QElapsedTimer et;
    et.start();

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPointsByScanNumber); ree;

    QVector<QVector<QVector<double>>> groupedMzVals;
    QVector<QVector<QVector<double>>> groupedIntensityVals;
    e = buildScanPointGroups(
            scanPointsByScanNumber.values().toVector(),
            &groupedMzVals,
            &groupedIntensityVals
    ); ree;
//    qDebug() << "Scans grouped in" << et.restart() << "mSec";

    QVector<Eigen::MatrixX<int>> connectedCentroidsMats;
    e = connectCentroidsInGroupedMzVals(
            groupedMzVals,
            m_featureFinderParameters.tolerancePPM,
            &connectedCentroidsMats
    ); ree;
//    qDebug() << connectedCentroidsMats.size() << "Centroid pairs connected" << et.restart() << "mSec";

    e = groupConnectedCentroidsToHills(
            connectedCentroidsMats,
            scanPointsByScanNumber,
            m_featureFinderParameters.minScanCount,
            featureFinderHills
    ); ree;
//    qDebug() << featureFinderHills->size() << "Centroids to hills" << et.restart() << "mSec";

    ERR_RETURN
}


namespace {


    struct RefineHillsInput {
        FeatureFinderHill featureFinderHill;
        FeatureFinderParameters params;
        int minScanCount = -1;
    };


    Err splitHill(
            const FeatureFinderHill &featureFinderHill,
            const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
            const QVector<double> &intensityVecSmoothed,
            int minHillScanCount,
            QVector<FeatureFinderHill> *refinedHills
            ) {

        ERR_INIT

        FeatureFinderHill refinedHill = featureFinderHill;
        refinedHill.updateIntensities(intensityVecSmoothed);

        for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {

            FeatureFinderHill splitHill = featureFinderHill;
            e = splitHill.refineHill(pii.first, pii.second); ree;

            if (splitHill.scanCount() < minHillScanCount) {
                continue;
            }

            refinedHills->push_back(splitHill);
        }

        ERR_RETURN
    }


    QPair<Err, QVector<FeatureFinderHill>> refineHillsLogic(const RefineHillsInput &rhi) {

        ERR_INIT

        PeakIntegratomaticParameters params;
        params.smoothCount = rhi.params.smoothCount;
        params.filterLength = rhi.params.filterLength;
        params.sigma = rhi.params.sigma;
        params.signalToNoiseRatio = rhi.params.signalToNoiseRatio;

        PeakIntegratomatic peakIntegratomatic;
        e = peakIntegratomatic.init(params); rree;

        QVector<PeakIntegrationIndexes> peakLimits;
        QVector<double> intensityVecSmoothed;
        e = peakIntegratomatic.findAllPeaksLimitsInXIC(
                rhi.featureFinderHill.intensities(),
                &peakLimits,
                &intensityVecSmoothed
        ); rree;

        const int minSplitVal = 2;
        if (peakLimits.size() < minSplitVal) {
            return {e, {rhi.featureFinderHill}};
        }

        QVector<FeatureFinderHill> newHills;
        e = splitHill(
                rhi.featureFinderHill,
                peakLimits,
                intensityVecSmoothed,
                rhi.minScanCount,
                &newHills
        ); rree;

        return {e, newHills};
    }


} //namespace
Err FeatureFinderHillBuilder::Private::refineHills(QVector<FeatureFinderHill> *featureFinderHills) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(*featureFinderHills); ree;

    QVector<RefineHillsInput> refineHillsInput;
    for (const FeatureFinderHill &ffh : (*featureFinderHills)) {
        RefineHillsInput rhi;
        rhi.params = m_featureFinderParameters;
        rhi.featureFinderHill = ffh;
        rhi.minScanCount = m_featureFinderParameters.minScanCount;
        refineHillsInput.push_back(rhi);
    }

    QVector<FeatureFinderHill> refinedHills;

    if (m_runParallel) {

        QFuture<QPair<Err, QVector<FeatureFinderHill>>> futures
                = QtConcurrent::mapped(refineHillsInput, refineHillsLogic);
        futures.waitForFinished();

        for (const QPair<Err, QVector<FeatureFinderHill>> &future: futures) {
            e = future.first; ree;
            refinedHills.append(future.second);
        }
    }

    else {

        for (const RefineHillsInput &rhi : refineHillsInput) {

            const QPair<Err, QVector<FeatureFinderHill>> refinedHillResult = refineHillsLogic(rhi);
            e = refinedHillResult.first; ree;

            refinedHills.append(refinedHillResult.second);
        }

    }

    *featureFinderHills = refinedHills;
    ERR_RETURN
}


void FeatureFinderHillBuilder::Private::setRunParallel(bool runParallel) {
    m_runParallel = runParallel;
}


///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////


FeatureFinderHillBuilder::FeatureFinderHillBuilder() : d_ptr(new Private()) {}

FeatureFinderHillBuilder::~FeatureFinderHillBuilder(){}

Err FeatureFinderHillBuilder::init(const FeatureFinderParameters &featureFinderParameters) {
    ERR_INIT
    e = d_ptr->init(featureFinderParameters); ree;
    ERR_RETURN
}


Err FeatureFinderHillBuilder::buildScanPointGroupsTest(
        const QVector<ScanPoints> &allScanPoints,
        QVector<QVector<QVector<double>>> *groupedMzVals,
        QVector<QVector<QVector<double>>> *groupedIntensityVals
) {

    ERR_INIT

    e = d_ptr->buildScanPointGroups(
            allScanPoints,
            groupedMzVals,
            groupedIntensityVals
    ); ree;

    ERR_RETURN
}


Err FeatureFinderHillBuilder::connectCentroidsInGroupedMzValsTest(
        const QVector<QVector<QVector<double>>> &groupedMzVals,
        double tolerancePPM,
        QVector<QVector<int>> *connectedCentroidsVecs
) {

    ERR_INIT

    QVector<Eigen::MatrixX<int>> connectedCentroidsMats;

    e = d_ptr->connectCentroidsInGroupedMzVals(
            groupedMzVals,
            tolerancePPM,
            &connectedCentroidsMats
    ); ree;

    for (const Eigen::MatrixX<int> &mat : connectedCentroidsMats) {

        QVector<int> vec;

        for (int row = 0; row < mat.rows(); row++) {

            const Eigen::VectorX<int> &v = mat.row(row);
            std::vector<int> vecReturn(v.data(), v.data() + v.size());
            vec.append(QVector<int>::fromStdVector(vecReturn));

        }

        connectedCentroidsVecs->push_back(vec);
    }

    ERR_RETURN
}


Err FeatureFinderHillBuilder::buildHills(
        const QMap<ScanNumber, ScanPoints> &scanPointsByScanNumber,
        QVector<FeatureFinderHill> *featureFinderHills
        ) {
    ERR_INIT
    e = d_ptr->buildHills(scanPointsByScanNumber, featureFinderHills); ree;
    ERR_RETURN
}


Err FeatureFinderHillBuilder::writeHillsToBatmassMzMrtFile(
        const QMap<ScanNumber, double> &scanNumberVsScanTime,
        const QVector<FeatureFinderHill> &featureFinderHills,
        const QString &destinationFilePath
        ) {

    ERR_INIT

    struct MzRtRow : public CSVReaderInputBase {
        double mzLo = -1.0;
        double mzHi = -1.0;
        double rtLo = -1.0;
        double rtHi = -1.0;

        QMap<QString, QVariant> map() override {
            return {
                {"mzLo", mzLo},
                {"mzHi", mzHi},
                {"rtLo", rtLo},
                {"rtHi", rtHi}
            };
        }
    };

    QVector<FeatureFinderHill> featureFinderHillsSorted = featureFinderHills;

    const auto sortLogic = [](const FeatureFinderHill &l, const FeatureFinderHill &r){

        if (MathUtils::tSame(l.mzMean(), r.mzMean())) {
            return l.minMaxScanNumber().first < r.minMaxScanNumber().first;
        }

        return l.mzMean() < r.mzMean();
    };
    std::sort(featureFinderHillsSorted.begin(), featureFinderHillsSorted.end(), sortLogic);

    QVector<MzRtRow> mzRtRows;
    for (const FeatureFinderHill &ffh : featureFinderHills) {

        const QPair<double, double> mzRange = ffh.mzMinMax();
        const QPair<int, int> scanNumberRange = ffh.minMaxScanNumber();

        MzRtRow mzRtRow;
        mzRtRow.mzLo = mzRange.first - 0.025;
        mzRtRow.mzHi = mzRange.second + 0.025;
        mzRtRow.rtLo = scanNumberVsScanTime.value(scanNumberRange.first);
        mzRtRow.rtHi = scanNumberVsScanTime.value(scanNumberRange.second);

        mzRtRows.push_back(mzRtRow);
    }

    const QVector<QSharedPointer<CSVReaderInputBase>> ptrs
            = CSVReaderInputBase::convertInputStructToSharedPointers(mzRtRows);

    CSVReader csvReader;
    e = csvReader.writeDataToCSV(
            destinationFilePath,
            ptrs
            ); ree;

    ERR_RETURN
}


void FeatureFinderHillBuilder::setRunParallel(bool runParallel) {
    d_ptr->setRunParallel(runParallel);
}

Err FeatureFinderHillBuilder::refineHills(QVector<FeatureFinderHill> *featureFinderHills) {
    ERR_INIT

    e = d_ptr->refineHills(featureFinderHills); ree;

    ERR_RETURN
}
