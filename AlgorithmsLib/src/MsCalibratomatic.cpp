//
// Created by anichols on 4/16/23.
//

#include "MsCalibratomatic.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"
#include "TandemFragmentPredictotron.h"

#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/serialization/vector.hpp>

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
using rTreeCoor = bg::model::point<int, 1, bg::cs::cartesian>;
using rTreePeakBox = bg::model::box<rTreeCoor>;
using rTreeSearchBox = bg::model::box<rTreeCoor>;
using rTreePoint = std::pair<rTreePeakBox, QString> ;
using RTree = bgi::rtree<rTreePoint, bgi::dynamic_quadratic>;


MsCalibratomatic::MsCalibratomatic() : m_calPointK(10) {}

Err MsCalibratomatic::init(
        const QMap<QString, QString> &scoreVectorsVsScanFrameFilePaths,
        const PythiaParameters &pythiaParameters,
        int calPointK,
        FragLibraryTronDIA *fragLibraryTronDia
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    e = ErrorUtils::isNotEqual(calPointK, 0); ree;
    m_calPointK = calPointK;

    m_fragLibraryTronDia = fragLibraryTronDia;

    e = buildCalibrator(scoreVectorsVsScanFrameFilePaths); ree;

    qDebug() << "NearestNeighbors Treesize" << m_nnSearch.kdTreeSize();

    ERR_RETURN
}

Err MsCalibratomatic::buildCalibrator(
        const QMap<QString, QString> &scoreVectorsVsScanFrameFilePaths
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(scoreVectorsVsScanFrameFilePaths); ree;

    for (
        auto it = scoreVectorsVsScanFrameFilePaths.begin();
        it != scoreVectorsVsScanFrameFilePaths.end();
        it++
        ) {

        const QString &scoreVectorsFilePath = it.key();
        const QString &msFrameScansFilePath = it.value();

        qDebug() << "Processing files" << scoreVectorsFilePath;

        e = processLogicForFrameScores(
                scoreVectorsFilePath,
                msFrameScansFilePath
                ); ree;

    }

    e = loadCalibrationPointsToKDTree(); ree;

    ERR_RETURN
}

Err MsCalibratomatic::processLogicForFrameScores(
        const QString &scoreVectorsFilePath,
        const QString &msFrameScansFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_fragLibraryTronDia->isInit()); ree;

    m_scoreVectors.clear();
    m_frameIndexVsScanPoints.clear();
    m_peptideWithModsVsCharge.clear();
    m_topCandidatesInFrameIndex.clear();

    ParquetReader::read(
            scoreVectorsFilePath,
            &m_scoreVectors
    ); ree;

    e = buildPeptideSequenceWithModsVsCharge(); ree;

    QVector<MsFrameScanPointRows> msFrameScanPointRows;
    ParquetReader::read(
            msFrameScansFilePath,
            &msFrameScanPointRows
    ); ree;

    e = MsFrame::buildFrameIndexVsScanPoints(
            msFrameScanPointRows,
            &m_frameIndexVsScanPoints
            ); ree;

    const int topN = 1;
    e = getTopNCandidatesPerFrameIndex(
            topN,
            &m_topCandidatesInFrameIndex
            );ree;

    QVector<PeptideStringWithModsScoreResult> topScoringPSMs;
    e = getScoredPSMsUntilFirstDecoyIsFound(&topScoringPSMs); ree;

    e = buildCalibrationPoints(topScoringPSMs); ree;

    ERR_RETURN
}

Err MsCalibratomatic::buildPeptideSequenceWithModsVsCharge() {

    ERR_INIT
    m_peptideWithModsVsCharge.clear();
    e = ErrorUtils::isNotEmpty(m_scoreVectors); ree;

    for (const MsFrameScoreVectorReaderRow &row : m_scoreVectors) {
        e = ErrorUtils::isFalse(m_peptideWithModsVsCharge.contains(row.peptideStringWithMods)); ree;
        m_peptideWithModsVsCharge.insert(row.peptideStringWithMods, row.charge);
    }

    ERR_RETURN
}

namespace {

    RTree buildScorePeaksRTree(const QVector<MsFrameScoreVectorReaderRow> &scoreVectors) {

        std::vector<rTreePoint> cloudLoader;
        for (const MsFrameScoreVectorReaderRow &row : scoreVectors) {

            const FrameIndex frameIndexStart = row.scorePeakStart;
            const FrameIndex frameIndexEnd = row.scorePeakEnd;

            const rTreePeakBox peakBox(
                    rTreeCoor(frameIndexStart, 0),
                    rTreeCoor(frameIndexEnd, 0)
            );

            cloudLoader.emplace_back(peakBox, row.peptideStringWithMods);
        }

        const int maxElements = 16;
        return RTree(cloudLoader, bgi::dynamic_quadratic(maxElements));
    }

    Err getScorePeptidesByFrameIndex(
            FrameIndex frameIndex,
            const RTree &rtree,
            QStringList *peptidesFoundInFrameIndex
            ) {

        ERR_INIT

        peptidesFoundInFrameIndex->clear();

        const rTreeSearchBox queryBox(rTreeCoor(frameIndex, 0), rTreeCoor(frameIndex, 0));
        std::vector<rTreePoint> result;
        rtree.query(bgi::intersects(queryBox), std::back_inserter(result));

        for (const rTreePoint &rtp : result) {
            peptidesFoundInFrameIndex->push_back(rtp.second);
        }

        ERR_RETURN
    }

    Err buildPeptideWithModsVsScoreVec(
            const QVector<MsFrameScoreVectorReaderRow> &scoreVectors,
            QMap<PeptideStringWithMods, QVector<double>> *peptideWithModsVsScoreVec
            ) {

        ERR_INIT

        peptideWithModsVsScoreVec->clear();

        e = ErrorUtils::isNotEmpty(scoreVectors); ree;

        for (const MsFrameScoreVectorReaderRow &r : scoreVectors) {
            e = ErrorUtils::isFalse(peptideWithModsVsScoreVec->contains(r.peptideStringWithMods)); ree;
            peptideWithModsVsScoreVec->insert(r.peptideStringWithMods, r.scorePerFrameIndexOfTargetVec);
        }

        ERR_RETURN
    }

    const auto sortScoresAscLogic = [](
            const QPair<PeptideStringWithMods, double> &l,
            const QPair<PeptideStringWithMods, double> &r
    ){return l.second < r.second;};

    Err getFrameIndexPeptidesWithScoresDesc(
            FrameIndex frameIndex,
            const QStringList &peptidesFoundInFrameIndex,
            const QMap<PeptideStringWithMods, QVector<double>> &peptideWithModsVsScoreVec,
            QVector<QPair<PeptideStringWithMods, double>> *peptideScores
            ) {

        ERR_INIT

        QVector<double> scores;
        for (const QString &peptideStringWithMods : peptidesFoundInFrameIndex) {

            e = ErrorUtils::isTrue(peptideWithModsVsScoreVec.contains(peptideStringWithMods)); ree;
            const QVector<double> &scoresVec = peptideWithModsVsScoreVec.value(peptideStringWithMods);

            const QPair<PeptideStringWithMods, double> point = {peptideStringWithMods, scoresVec.at(frameIndex)};

            const double minScore = 0.001;
            if (point.second < minScore){
                continue;
            }
            scores.push_back(point.second);
            peptideScores->push_back(point);
        }

//        const double stDevs = 2.0;
//        const double thresholdScore = MathUtils::mean(scores) + (stDevs * MathUtils::stDev(scores));

        std::sort(peptideScores->rbegin(), peptideScores->rend(), sortScoresAscLogic);

        ERR_RETURN
    }

}//namespace
Err MsCalibratomatic::getTopNCandidatesPerFrameIndex(
        int topN,
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
) {

    ERR_INIT

    m_topCandidatesInFrameIndex.clear();

    const RTree scoresRTree = buildScorePeaksRTree(m_scoreVectors);

    QMap<PeptideStringWithMods, QVector<double>> peptideWithModsVsScoreVec;
    e = buildPeptideWithModsVsScoreVec(
            m_scoreVectors,
            &peptideWithModsVsScoreVec
    ); ree;

    const QList<FrameIndex> &frameIndexes = m_frameIndexVsScanPoints.keys();
    for (const FrameIndex frameIndex : frameIndexes) {

        QStringList peptidesFoundInFrameIndex;
        e = getScorePeptidesByFrameIndex(
                frameIndex,
                scoresRTree,
                &peptidesFoundInFrameIndex
        ); ree;

        QVector<QPair<PeptideStringWithMods, Score>> peptideWithScoresDesc;
        e = getFrameIndexPeptidesWithScoresDesc(
                frameIndex,
                peptidesFoundInFrameIndex,
                peptideWithModsVsScoreVec,
                &peptideWithScoresDesc
        ); ree;

        if (peptideWithScoresDesc.isEmpty()) {
            continue;
        }

        topN = std::min(topN, peptideWithScoresDesc.size());
        peptideWithScoresDesc.resize(topN);
        topCansInFrameIndex->insert(frameIndex, peptideWithScoresDesc);
    }

    ERR_RETURN
}

Err MsCalibratomatic::getScoredPSMsUntilFirstDecoyIsFound(QVector<PeptideStringWithModsScoreResult> *scoreNoDecoys) {

    ERR_INIT

    QVector<PeptideStringWithModsScoreResult> scores;
    for (auto it = m_topCandidatesInFrameIndex.begin(); it != m_topCandidatesInFrameIndex.end(); it++) {

        const FrameIndex frameIndex = it.key();
        const QVector<QPair<PeptideStringWithMods, Score>> &scorePair = it.value();

        for (const QPair<PeptideStringWithMods, Score> &pr : scorePair) {

            if (pr.first.isEmpty()) {
                continue;
            }

            PeptideStringWithModsScoreResult res;
            res.frameIndex = frameIndex;
            res.peptideStringWithMods = pr.first;
            res.score = pr.second;

            e = ErrorUtils::isTrue(m_peptideWithModsVsCharge.contains(res.peptideStringWithMods));ree;
            res.charge = m_peptideWithModsVsCharge.value(res.peptideStringWithMods);

            scores.push_back(res);
        }
    }

    using T = PeptideStringWithModsScoreResult;
    std::sort(scores.rbegin(), scores.rend(), [](const T &l, const T &r){return l.score < r.score;});

    for (const T &r : scores) {

        bool isDecoy;
        e = m_fragLibraryTronDia->peptideStringWithModsIsDecoy(
                r.peptideStringWithMods,
                &isDecoy
        ); ree;

        if (isDecoy) {
            break;
        }

        scoreNoDecoys->push_back(r);
    }

    ERR_RETURN
}

Err MsCalibratomatic::buildCalibrationPoints(const QVector<PeptideStringWithModsScoreResult> &scoresNoDecoys) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scoresNoDecoys); ree;

    for (const PeptideStringWithModsScoreResult &res : scoresNoDecoys) {

        const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
                res.peptideStringWithMods,
                res.charge
                );

        QVector<MS2Ion> theoTandemPrediction;
        e = m_fragLibraryTronDia->getMS2Ions(
                peptideSequenceChargeKey,
                &theoTandemPrediction
                ); ree
        e = ErrorUtils::isNotEmpty(theoTandemPrediction); ree;
        const ScanPoints theoTandemPredictionScanPoints
                = FragLibraryTronDIA::ms2IonsToScanPoints(theoTandemPrediction);

        const ScanPoints &scanPoints = m_frameIndexVsScanPoints.value(res.frameIndex);
        e = ErrorUtils::isNotEmpty(scanPoints); ree;

        const ExtractPoints extractPoints = MsUtils::extractPointsFromPoints(
                scanPoints,
                theoTandemPredictionScanPoints,
                m_params.ms2ExtractionWidthPPM
                );

        m_calibrationPoints[res.frameIndex].push_back(extractPoints);

    }

    ERR_RETURN
}

namespace {

    using DiffPPM = double;
    using Coors = QVector<double>;

    void buildNNInput(
            const QMap<FrameIndex, QVector<ExtractPoints>> &calibrationPoints,
            QVector<QPair<DiffPPM, Coors>> *valuesVsTreePoints,
            QVector<double> *ppmDiffVals
            ) {

        for (auto it = calibrationPoints.begin(); it != calibrationPoints.end(); it++) {

            const FrameIndex frameIndex = it.key();
            const auto frameIndexDouble = static_cast<double>(frameIndex);
            const QVector<ExtractPoints> &eps = it.value();

            for (const ExtractPoints &ep : eps) {

                for (const QPointF &p : ep.mzFoundVsSearched) {

                    const double mzFound = p.x();
                    if (mzFound < 0) {
                        continue;
                    }

                    const double mzTheo = p.y();
                    const double ppmDiff = 1e6 * ((mzFound - mzTheo) / mzTheo);

                    valuesVsTreePoints->push_back({ppmDiff, {frameIndexDouble, mzFound}});
                    ppmDiffVals->push_back(ppmDiff);
                }
            }
        }

    }

    void filterNNInput(
            const QVector<double> &ppmDiffVals,
            QVector<QPair<DiffPPM, Coors>> *valuesVsTreePoints,
            QVector<QPair<DiffPPM, Coors>> *valuesVsTreePointsRemoved
            ) {

        const double ppmMean = MathUtils::mean(ppmDiffVals);
        const double stDev = MathUtils::stDev(ppmDiffVals);

        const double stDevMultiplier = 3.0;
        const double ppmMin = ppmMean - (stDevMultiplier * stDev);
        const double ppmMax = ppmMean + (stDevMultiplier * stDev);

        const auto terminatorLogic = [ppmMin, ppmMax](const QPair<DiffPPM, Coors> &calPoint){
            return !(ppmMin <= calPoint.first && calPoint.first <= ppmMax);
        };

        const auto terminator = std::remove_if(
                valuesVsTreePoints->begin(),
                valuesVsTreePoints->end(),
                terminatorLogic
                );

        valuesVsTreePoints->erase(terminator, valuesVsTreePoints->end());

        const auto terminatorLogicRemoved = [ppmMin, ppmMax](const QPair<DiffPPM, Coors> &calPoint){
            return ppmMin <= calPoint.first && calPoint.first <= ppmMax;
        };

        const auto terminatorRemoved = std::remove_if(
                valuesVsTreePointsRemoved->begin(),
                valuesVsTreePointsRemoved->end(),
                terminatorLogicRemoved
        );

        valuesVsTreePointsRemoved->erase(terminatorRemoved, valuesVsTreePointsRemoved->end());

    }

}//namespace
Err MsCalibratomatic::loadCalibrationPointsToKDTree() {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_calibrationPoints); ree;

    QVector<QPair<DiffPPM, QVector<double>>> valuesVsTreePoints;
    QVector<QPair<DiffPPM, QVector<double>>> valuesVsTreePointsRemoved;
    QVector<double> ppmDiffVals;

    buildNNInput(
            m_calibrationPoints,
            &valuesVsTreePoints,
            &ppmDiffVals
            );

    valuesVsTreePointsRemoved = valuesVsTreePoints;

    qDebug() << "Calibration points count" << valuesVsTreePoints.size();
    filterNNInput(
            ppmDiffVals,
            &valuesVsTreePoints,
            &valuesVsTreePointsRemoved
            );
    qDebug() << "Calibration points filtered count" << valuesVsTreePoints.size();

    e = m_nnSearch.init(valuesVsTreePoints); ree;

    ERR_RETURN
}

namespace {

    QVector<QVector<double>> buildNearestNeighborsInput(const QMap<FrameIndex, ScanPoints> &indexVsScanPoints) {

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

                const double meanDiffPPM = MathUtils::mean(nnSearchResult.values);
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

Err MsCalibratomatic::init(
        const QString &calFilePath,
        const QString &matFilePath
        ) {

    ERR_INIT

    e = m_nnSearch.readNearestNeighbors(
            calFilePath,
            matFilePath
            ); ree;

    ERR_RETURN
}

Err MsCalibratomatic::writeCalibratomatic(
        const QString &msDataFilePath,
        QString *calibrationMatFilePath,
        QString *calibarationCalFilePath
        ) {

    ERR_INIT
    e = m_nnSearch.writeNearestNeighbors(
            msDataFilePath,
            calibrationMatFilePath,
            calibarationCalFilePath
            ); ree;
    ERR_RETURN
}
