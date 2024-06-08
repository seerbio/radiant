//
// Created by anichols on 12/10/22.
//

#include "FeatureFinderHillBuilder.h"

#include "CSVReader.h"
#include "ErrorUtils.h"
#include "MathUtils.h"
#include "MsReaderBase.h"
#include "ParallelUtils.h"
//#include "PeakIntegratomatic.h"

#include <Eigen/Dense>

#include <QtConcurrent/QtConcurrent>

#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;


class Q_DECL_HIDDEN FeatureFinderHillBuilder::Private
{
public:

    using RTreeCoor = bg::model::point<double, 2, bg::cs::cartesian>;
    using RTreeBox = bg::model::box<RTreeCoor>;
    using RTreePoint = std::pair<RTreeBox, FeatureFinderHill*> ;
    using RTree = bgi::rtree<RTreePoint, bgi::dynamic_quadratic>;

    Private();
    ~Private();

    Err init(const FeatureFinderParameters &featureFinderParameters);

    Err buildScanPointGroups(
            const QVector<ScanPoints*> &allScanPoints,
            QVector<QVector<QVector<float>>> *groupedMzVals,
            QVector<QVector<QVector<float>>> *groupedIntensityVals
    );

    Err connectCentroidsInGroupedMzVals(
            const QVector<QVector<QVector<float>>> &groupedMzVals,
            double tolerancePPM,
            QVector<Eigen::MatrixX<int>> *connectedCentroids
    ) const;

    Err buildHills(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints);

    Err refineHills(bool useSmoothing);

    void setRunParallel(bool runParallel);

    Err initRTree();

    Err getRTreeLimits(
            double *mzMin,
            double *mzMax,
            FrameIndex *frameIndexMin,
            FrameIndex *frameIndexMax
    );

    Err getHills(
            double mzMin,
            double mzMax,
            FrameIndex frameIndexMin,
            FrameIndex frameIndexMax,
            QVector<FeatureFinderHill*> *featureFinderHills
    );

    Err featureFinderHills(QVector<FeatureFinderHill*> *featureFinderHills);

    Err getHills(
            double mzMin,
            double mzMax,
            QVector<FeatureFinderHill*> *featureFinderHills
    );

private:

    QVector<FeatureFinderHill> m_featureFinderHills;
    FeatureFinderParameters m_featureFinderParameters;
    bool m_runParallel;
    bool m_isInit;

    QMap<ScanNumberIndex, ScanNumber> m_scanNumberIndexVsScanNumber;

    RTree *m_rTree;
};

FeatureFinderHillBuilder::Private::Private()
: m_runParallel(false)
, m_isInit(false)
{}

FeatureFinderHillBuilder::Private::~Private(){
    delete m_rTree;
}

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
        QVector<ScanPoints*> allScanPointsTranch;
    };

    struct ScanGroupingParallelOutput {
        Err e = Error::eNoError;
        QVector<QVector<QVector<float>>> groupedMzVals;
        QVector<QVector<QVector<float>>> groupedIntensityVals;
    };

    QVector<ScanGroupingParallelInput> buildScanGroupingParallelInputs(
            const QVector<QVector<ScanPoints*>> &tranchedScanPoints,
            int skipScanCount
    ) {

        QVector<ScanGroupingParallelInput> scanGroupingParallelInputs;

        const int minAllowableGroupCount = 2;

        for (const QVector<ScanPoints*> &spVec : tranchedScanPoints) {

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

            QVector<QVector<float>> mzValsGroup;
            QVector<QVector<float>> intensityValsGroup;

            for (ScanNumberIndex nextScanIndex = scanIndex; nextScanIndex < stopNextScanIndex; ++nextScanIndex) {

                if (nextScanIndex >= sgpi.allScanPointsTranch.size()) {
                    break;
                }

                ScanPoints* scanPoints = sgpi.allScanPointsTranch.at(nextScanIndex);

                QVector<float> mzVals;
                QVector<float> intensityVals;
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
        const QVector<ScanPoints*> &allScanPoints,
        QVector<QVector<QVector<float>>> *groupedMzVals,
        QVector<QVector<QVector<float>>> *groupedIntensityVals
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(allScanPoints); ree;
    groupedMzVals->clear();
    groupedIntensityVals->clear();

    ScanGroupingParallelInput scanGroupingParallelInput;
    scanGroupingParallelInput.skipScanCount = m_featureFinderParameters.skipScanCount;
    scanGroupingParallelInput.allScanPointsTranch = allScanPoints;
    scanGroupingParallelInput.isLastTranche = true;

    const ScanGroupingParallelOutput output = scanGroupingParallelLogic(scanGroupingParallelInput);
    *groupedMzVals = output.groupedMzVals;
    *groupedIntensityVals = output.groupedIntensityVals;
    e = output.e; ree;

    ERR_RETURN
}

namespace {

    struct ConnectCentroidsInGroupedMzValsInput {
        QVector<QVector<float>> mzValsGroup;
        double tolerancePPM = -1.0;
    };


    QVector<ConnectCentroidsInGroupedMzValsInput> buildConnectCentroidsInGroupedMzValsInputVec(
            const QVector<QVector<QVector<float>>> &groupedMzVals,
            double tolerancePPM
    ) {

        QVector<ConnectCentroidsInGroupedMzValsInput> inputs;

        for (const QVector<QVector<float>> &mzValsGroup : groupedMzVals) {

            ConnectCentroidsInGroupedMzValsInput input;
            input.mzValsGroup = mzValsGroup;
            input.tolerancePPM = tolerancePPM;

            inputs.push_back(input);
        }

        return inputs;
    }


    Eigen::VectorX<int> findCentroidConnections(
            const QVector<float> &referenceScan,
            const QVector<float> &nextScan,
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

        const QVector<float> &referenceScan = input.mzValsGroup.front();

        const int returnMatRows = input.mzValsGroup.size() - 1;
        const int returnMatCols = referenceScan.size();
        Eigen::MatrixX<int> returnMat(returnMatRows, returnMatCols);

        for (ScanNumberIndex nextScanIndex = 1; nextScanIndex < input.mzValsGroup.size(); nextScanIndex++) {

            const QVector<float> &nextScan = input.mzValsGroup.at(nextScanIndex);
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
        const QVector<QVector<QVector<float>>> &groupedMzVals,
        double tolerancePPM,
        QVector<Eigen::MatrixX<int>> *connectedCentroids
) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(groupedMzVals); ree;

    QVector<ConnectCentroidsInGroupedMzValsInput> inputs = buildConnectCentroidsInGroupedMzValsInputVec(
            groupedMzVals,
            tolerancePPM
    );

    if (m_runParallel) {

        QFuture<Eigen::MatrixX<int>> futures
                = QtConcurrent::mapped(inputs, connectCentroidsInGroupedMzValsParallelLogic);

        futures.waitForFinished();

        for (const Eigen::MatrixX<int> &mat: futures) {
            connectedCentroids->push_back(mat);
        }
    }

    else {

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
                (*connectedCentroidsMats)[currentScanIndex].coeffRef(row, currentIndexVal) = -2;
                continue;
            }

            for (int rowIndexEraser = 0; rowIndexEraser < currentCentroids.rows(); ++rowIndexEraser) {
                (*connectedCentroidsMats)[currentScanIndex].coeffRef(rowIndexEraser, currentIndexVal) = -2;
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
            const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints,
            int minScanCount,
            QVector<FeatureFinderHill> *featureFinderHills
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(_connectedCentroidsMats); ree;
        featureFinderHills->clear();

        const QVector<ScanNumber> scanNumbers = scanNumberVsScanPoints.keys().toVector();

        //NOTE: This is not const because the next index function modifies it to zero out terminals
        QVector<Eigen::MatrixX<int>> connectedCentroidsMats = _connectedCentroidsMats;

        const QList<ScanPoints*> &scanPoints = scanNumberVsScanPoints.values();

        for (ScanNumberIndex baseScanIndex = 0; baseScanIndex < connectedCentroidsMats.size(); baseScanIndex++) {

            const Eigen::MatrixX<int> &baseCentroids = connectedCentroidsMats.at(baseScanIndex);

            for (int baseCentroidsColInd = 0; baseCentroidsColInd < baseCentroids.cols(); ++baseCentroidsColInd) {

                const Eigen::VectorX<int> colVec = baseCentroids.col(baseCentroidsColInd);
                const int colVecMaxForSinglePointAddition = colVec.maxCoeff();

                if (colVecMaxForSinglePointAddition == -1) {

                    ScanPoint currentScanPoint = scanPoints.at(baseScanIndex)->at(baseCentroidsColInd);
                    ScanNumberIndex currentScanIndex = baseScanIndex;

                    FeatureFinderHill featureFinderHill;
                    featureFinderHill.addPoint(
                            currentScanIndex,
                            scanNumbers.at(currentScanIndex),
                            currentScanPoint.x(),
                            currentScanPoint.y()
                    );

                    if (featureFinderHill.scanCount() >= minScanCount) {
                        featureFinderHills->push_back(featureFinderHill);
                    }

                    continue;
                }

                for (int baseCentroidsRowInd = 0; baseCentroidsRowInd < baseCentroids.rows(); ++baseCentroidsRowInd) {

                    const int nextScanCentroidVal = baseCentroids.coeff(baseCentroidsRowInd, baseCentroidsColInd);
                    if (nextScanCentroidVal < 0) {
                        continue;
                    }

                    ScanPoint currentScanPoint = scanPoints.at(baseScanIndex)->at(baseCentroidsColInd);

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

                    while (nextIndexVal >= 0) {

                        getNextIndex(
                                nextScanIndex,
                                nextIndexVal,
                                &connectedCentroidsMats,
                                &currentScanIndex,
                                &currentIndexVal,
                                &nextScanIndex,
                                &nextIndexVal
                        );

                        if (nextIndexVal <= -1) {
                            break;
                        }

                        currentScanPoint = scanPoints.at(currentScanIndex)->at(currentIndexVal);
                        featureFinderHill.addPoint(
                                currentScanIndex,
                                scanNumbers.at(currentScanIndex),
                                currentScanPoint.x(),
                                currentScanPoint.y()
                        );
                    }

                    currentScanPoint = scanPoints.at(currentScanIndex)->at(currentIndexVal);
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

    Err buildScanNumberIndexVsScanNumber(
            const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints,
            QMap<ScanNumberIndex, ScanNumber> *scanNumberIndexVsScanNumbers
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scanNumberVsScanPoints); ree;

        for (const ScanNumber &sn : scanNumberVsScanPoints.keys()) {
            scanNumberIndexVsScanNumbers->insert(scanNumberIndexVsScanNumbers->size(), sn);
        }

        ERR_RETURN
    }

}//namespace
Err FeatureFinderHillBuilder::Private::buildHills(
        const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints
) {

    QElapsedTimer et;
    et.start();

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanNumberVsScanPoints); ree;

    e = buildScanNumberIndexVsScanNumber(
            scanNumberVsScanPoints,
            &m_scanNumberIndexVsScanNumber
            ); ree;

    QVector<QVector<QVector<float>>> groupedMzVals;
    QVector<QVector<QVector<float>>> groupedIntensityVals;
    e = buildScanPointGroups(
            scanNumberVsScanPoints.values().toVector(),
            &groupedMzVals,
            &groupedIntensityVals
    ); ree;

    QVector<Eigen::MatrixX<int>> connectedCentroidsMats;
    e = connectCentroidsInGroupedMzVals(
            groupedMzVals,
            m_featureFinderParameters.tolerancePPM,
            &connectedCentroidsMats
    ); ree;

    e = groupConnectedCentroidsToHills(
            connectedCentroidsMats,
            scanNumberVsScanPoints,
            m_featureFinderParameters.minScanCount,
            &m_featureFinderHills
    ); ree;

    e = initRTree(); ree;

    ERR_RETURN
}

namespace {

    struct RefineHillsInput {
        FeatureFinderHill featureFinderHill;
        FeatureFinderParameters params;
        int minScanCount = -1;
    };


    Err splitHill(
            const QMap<ScanNumberIndex, ScanNumber> &scanNumberIndexVsScanNumber,
            const FeatureFinderHill &featureFinderHill,
            const QVector<double> &_intensityVecSmoothed,
            const QVector<PeakIntegrationIndexes> &peakIntegrationIndexes,
            int minHillScanCount,
            QVector<FeatureFinderHill> *refinedHills
            ) {

        ERR_INIT


        if (peakIntegrationIndexes.isEmpty()) {
            refinedHills->push_back(featureFinderHill);
            ERR_RETURN
        }

        for (const PeakIntegrationIndexes &pii : peakIntegrationIndexes) {

            const int distance = pii.second - pii.first + 1;

            const QVector<ScanNumberIndex> scanNumberIndexes
                = featureFinderHill.scanNumberIndexes().mid(pii.first, distance);
            const QVector<ScanNumber> scanNumbers
                = featureFinderHill.scanNumbers().mid(pii.first, distance);
            const QVector<float> intensities = featureFinderHill.intensities().mid(pii.first, distance);

            FeatureFinderHill splitHill;
            for (int i = 0; i < scanNumbers.size(); i++) {

                splitHill.addPoint(
                        scanNumberIndexes.at(i),
                        scanNumbers.at(i),
                        featureFinderHill.mzMean(),
                        intensities.at(i)
                        );
            }

            refinedHills->push_back(splitHill);
        }

        ERR_RETURN
    }

    QPair<Err, QVector<FeatureFinderHill>> refineHillsLogic(
            const QMap<ScanNumberIndex, ScanNumber> &scanNumberIndexVsScanNumber,
            const RefineHillsInput &rhi,
            bool useSmoothing
            ) {

        ERR_INIT

//        PeakIntegratomaticParameters params;
//        params.smoothCount = useSmoothing ? rhi.params.smoothCount : 0;
//        params.filterLength = rhi.params.filterLength;
//        params.sigma = rhi.params.sigma;
//        params.signalToNoiseRatio = rhi.params.signalToNoiseRatio;
//
//        PeakIntegratomatic peakIntegratomatic;
//        e = peakIntegratomatic.init(params); rree;
//
//        QVector<PeakIntegrationIndexes> peakLimits;
//        QVector<double> intensityVecSmoothed;
//        e = peakIntegratomatic.findAllPeaksLimitsInXIC(
//                rhi.featureFinderHill.intensities(),
//                &peakLimits,
//                &intensityVecSmoothed
//        ); rree;
//
//        if (peakLimits.isEmpty()) {
//
//            FeatureFinderHill ffh = rhi.featureFinderHill;
//
////            qDebug() << intensityVecSmoothed;
////            qDebug() << ffh.intensities();
////            qDebug() << ffh.scanNumberIndexes();
////            qDebug() << ffh.scanNumbers();
////            qDebug() << "**";
//
//            return {e, {rhi.featureFinderHill}};
//        }
//
        QVector<FeatureFinderHill> newHills;
//        e = splitHill(
//                scanNumberIndexVsScanNumber,
//                rhi.featureFinderHill,
//                intensityVecSmoothed,
//                peakLimits,
//                rhi.minScanCount,
//                &newHills
//        ); rree;

        return {e, newHills};
    }


} //namespace
Err FeatureFinderHillBuilder::Private::refineHills(bool useSmoothing) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_featureFinderHills); ree;

    const auto transformLogic = [&](const FeatureFinderHill &ffh){
        RefineHillsInput rhi;
        rhi.params = m_featureFinderParameters;
        rhi.featureFinderHill = ffh;
        rhi.minScanCount = m_featureFinderParameters.minScanCount;
        return rhi;
    };

    QVector<RefineHillsInput> refineHillsInput;
    refineHillsInput.reserve(m_featureFinderHills.size());
    std::transform(
            m_featureFinderHills.begin(),
            m_featureFinderHills.end(),
            std::back_inserter(refineHillsInput),
            transformLogic
            );

    QVector<FeatureFinderHill> refinedHills;

    if (m_runParallel) {

        const auto refineHillsLogicBinder = std::bind(
                refineHillsLogic,
                m_scanNumberIndexVsScanNumber,
                std::placeholders::_1,
                useSmoothing
        );

        QFuture<QPair<Err, QVector<FeatureFinderHill>>> futures
                = QtConcurrent::mapped(refineHillsInput, refineHillsLogicBinder);
        futures.waitForFinished();

        for (const QPair<Err, QVector<FeatureFinderHill>> &future: futures) {
            e = future.first; ree;
            refinedHills.append(future.second);
        }
    }

    else {

        for (const RefineHillsInput &rhi : refineHillsInput) {

            const QPair<Err, QVector<FeatureFinderHill>> refinedHillResult = refineHillsLogic(
                    m_scanNumberIndexVsScanNumber,
                    rhi,
                    useSmoothing
                    );

            e = refinedHillResult.first; ree;

            refinedHills.append(refinedHillResult.second);
        }

    }

    m_featureFinderHills.clear();
    m_featureFinderHills = refinedHills;

    e = initRTree(); ree;

    ERR_RETURN
}

void FeatureFinderHillBuilder::Private::setRunParallel(bool runParallel) {
    m_runParallel = runParallel;
}

Err FeatureFinderHillBuilder::Private::initRTree() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_featureFinderHills); ree;

    std::vector<RTreePoint> cloudLoader;
    cloudLoader.reserve(m_featureFinderHills.size());
    std::transform(
            m_featureFinderHills.begin(),
            m_featureFinderHills.end(),
            std::back_inserter(cloudLoader),
            [](FeatureFinderHill &ffh){

                const QPair<double, double> mzMinMax = ffh.mzMinMax();
                const QPair<ScanNumber, ScanNumber> intensityMinMax = ffh.scanNumberIndexMinMax();

                const RTreeBox box(
                        RTreeCoor(mzMinMax.first, static_cast<double>(intensityMinMax.first)),
                        RTreeCoor(mzMinMax.second, static_cast<double>(intensityMinMax.second))
                );
                const RTreePoint rTreePoint(box, &ffh);
                return rTreePoint;
            }
    );

    const int maxElements = 16;
    m_rTree = new RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    e = ErrorUtils::isFalse(m_rTree->empty()); ree;

    m_isInit = true;
    ERR_RETURN
}

Err FeatureFinderHillBuilder::Private::getRTreeLimits(
        double *mzMin,
        double *mzMax,
        FrameIndex *frameIndexMin,
        FrameIndex *frameIndexMax
) {

    ERR_INIT

    e = ErrorUtils::isFalse(m_rTree->empty()); ree

    *mzMin = m_rTree->bounds().min_corner().get<0>();
    *mzMax = m_rTree->bounds().max_corner().get<0>();
    *frameIndexMin = static_cast<FrameIndex>(m_rTree->bounds().min_corner().get<1>());
    *frameIndexMax = static_cast<FrameIndex>(m_rTree->bounds().max_corner().get<1>());

    ERR_RETURN
}

Err FeatureFinderHillBuilder::Private::getHills(
        double mzMin,
        double mzMax,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        QVector<FeatureFinderHill*> *featureFinderHills
) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;

    const RTreeBox queryBox(
            RTreeCoor(mzMin, static_cast<double>(frameIndexMin)),
            RTreeCoor(mzMax,static_cast<double>(frameIndexMax))
    );

    std::vector<RTreePoint> result;
    m_rTree->query(bgi::intersects(queryBox), std::back_inserter(result));

    for (const RTreePoint &rtp : result) {
        featureFinderHills->push_back(rtp.second);
    }

    ERR_RETURN
}

Err FeatureFinderHillBuilder::Private::getHills(
        double mzMin,
        double mzMax,
        QVector<FeatureFinderHill*> *featureFinderHills
        ) {

    ERR_INIT

    double mzMinRTree;
    double mzMaxRTree;
    FrameIndex frameIndexMinRTree;
    FrameIndex frameIndexMaxRTree;

    e = getRTreeLimits(
            &mzMinRTree,
            &mzMaxRTree,
            &frameIndexMinRTree,
            &frameIndexMaxRTree
            ); ree;

    e = getHills(
            mzMin,
            mzMax,
            frameIndexMinRTree,
            frameIndexMaxRTree,
            featureFinderHills
            ); ree;

    ERR_RETURN
}

Err FeatureFinderHillBuilder::Private::featureFinderHills(QVector<FeatureFinderHill *> *featureFinderHills) {
    ERR_INIT

    for (FeatureFinderHill &ffh : m_featureFinderHills) {
        featureFinderHills->push_back(&ffh);
    }

    ERR_RETURN
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

FeatureFinderHillBuilder::FeatureFinderHillBuilder() : d_ptr(QScopedPointer<Private>(new Private)) {}

FeatureFinderHillBuilder::~FeatureFinderHillBuilder(){}

Err FeatureFinderHillBuilder::init(const FeatureFinderParameters &featureFinderParameters) {
    ERR_INIT
    e = d_ptr->init(featureFinderParameters); ree;
    ERR_RETURN
}

Err FeatureFinderHillBuilder::buildScanPointGroupsTest(
        const QVector<ScanPoints*> &allScanPoints,
        QVector<QVector<QVector<float>>> *groupedMzVals,
        QVector<QVector<QVector<float>>> *groupedIntensityVals
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
        const QVector<QVector<QVector<float>>> &groupedMzVals,
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
            vec.append(QVector<int>(vecReturn.begin(), vecReturn.end()));

        }

        connectedCentroidsVecs->push_back(vec);
    }

    ERR_RETURN
}

Err FeatureFinderHillBuilder::buildHills(const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints) {
    ERR_INIT
    e = d_ptr->buildHills(scanNumberVsScanPoints); ree;
    ERR_RETURN
}

Err FeatureFinderHillBuilder::writeHillsToBatmassMzMrtFile(
        const QMap<ScanNumber, double> &scanNumberVsScanTime,
        const QVector<FeatureFinderHill> &featureFinderHills,
        const QString &destinationFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(destinationFilePath.contains(".mzrt.csv"));
    if (e != eNoError) {
        qDebug() << "to write hills and view in Batmass, file extension must be .mzrt.csv";
        rrr(eError)
    }

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
            return l.scanNumberMinMax().first < r.scanNumberMinMax().first;
        }

        return l.mzMean() < r.mzMean();
    };
    std::sort(featureFinderHillsSorted.begin(), featureFinderHillsSorted.end(), sortLogic);

    QVector<MzRtRow> mzRtRows;
    for (const FeatureFinderHill &ffh : featureFinderHills) {

        const QPair<double, double> mzRange = ffh.mzMinMax();
        const QPair<int, int> scanNumberRange = ffh.scanNumberMinMax();

        MzRtRow mzRtRow;
        mzRtRow.mzLo = mzRange.first - 0.025;
        mzRtRow.mzHi = mzRange.second + 0.025;
        mzRtRow.rtLo = scanNumberVsScanTime.value(scanNumberRange.first);
        mzRtRow.rtHi = scanNumberVsScanTime.value(scanNumberRange.second) + 0.1;

        mzRtRows.push_back(mzRtRow);
    }

    e = CSVReader::write(
            mzRtRows,
            destinationFilePath
            ); ree

    ERR_RETURN
}

void FeatureFinderHillBuilder::setRunParallel(bool runParallel) {
    d_ptr->setRunParallel(runParallel);
}

Err FeatureFinderHillBuilder::refineHills(bool useSmoothing) {
    ERR_INIT
    e = d_ptr->refineHills(useSmoothing); ree;
    ERR_RETURN
}

Err FeatureFinderHillBuilder::getHills(
        double mzMin,
        double mzMax,
        FrameIndex frameIndexMin,
        FrameIndex frameIndexMax,
        QVector<FeatureFinderHill*> *featureFinderHills
        ) {
    ERR_INIT
    e = d_ptr->getHills(
            mzMin,
            mzMax,
            frameIndexMin,
            frameIndexMax,
            featureFinderHills
            ); ree;
    ERR_RETURN
}

Err FeatureFinderHillBuilder::getHills(
        double mzMin,
        double mzMax,
        QVector<FeatureFinderHill*> *featureFinderHills
) {
    ERR_INIT
    e = d_ptr->getHills(
            mzMin,
            mzMax,
            featureFinderHills
    ); ree;
    ERR_RETURN
}

Err FeatureFinderHillBuilder::featureFinderHills(QVector<FeatureFinderHill *> *featureFinderHills) {
    ERR_INIT
    e = d_ptr->featureFinderHills(featureFinderHills); ree;
    ERR_RETURN
}
