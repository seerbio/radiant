//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretronProcessormatic.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"

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


Err MsFrameScoretronProcessormatic::processLogicForFrameScores(
        const QString &scoreVectorsFilePath,
        const QString &msFrameScansFilePath,
        int topNPSMs,
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
        ) {

    ERR_INIT

    QVector<MsFrameScoreVectorReaderRow> scoreVectors;
    ParquetReader::read(
            scoreVectorsFilePath,
            &scoreVectors
    ); ree;
    e = ErrorUtils::isNotEmpty(scoreVectors); ree;

    QVector<MsFrameScanPointRows> msFrameScanPointRows;
    ParquetReader::read(
            msFrameScansFilePath,
            &msFrameScanPointRows
    ); ree;
    e = ErrorUtils::isNotEmpty(msFrameScanPointRows); ree;

    QMap<FrameIndex, ScanPoints> frameIndexVsScanPoints;
    e = MsFrame::buildFrameIndexVsScanPoints(
            msFrameScanPointRows,
            &frameIndexVsScanPoints
    ); ree;

    e = getTopNCandidatesPerFrameIndex(
            scoreVectors,
            frameIndexVsScanPoints,
            topNPSMs,
            topCansInFrameIndex
    );ree;

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

        std::sort(peptideScores->rbegin(), peptideScores->rend(), sortScoresAscLogic);

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretronProcessormatic::getTopNCandidatesPerFrameIndex(
        const QVector<MsFrameScoreVectorReaderRow> &scoreVectors,
        const QMap<FrameIndex, ScanPoints> &frameIndexVsScanPoints,
        int topNPSMs,
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
        ) {

    ERR_INIT

    topCansInFrameIndex->clear();

    const RTree scoresRTree = buildScorePeaksRTree(scoreVectors);

    QMap<PeptideStringWithMods, QVector<double>> peptideWithModsVsScoreVec;
    e = buildPeptideWithModsVsScoreVec(
            scoreVectors,
            &peptideWithModsVsScoreVec
    ); ree;

    const QList<FrameIndex> &frameIndexes = frameIndexVsScanPoints.keys();
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

        const int topN = std::min(topNPSMs, peptideWithScoresDesc.size());

        peptideWithScoresDesc.resize(topN);
        topCansInFrameIndex->insert(frameIndex, peptideWithScoresDesc);

    }

    ERR_RETURN
}
