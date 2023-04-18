//
// Created by anichols on 4/16/23.
//

#include "MsCalibratomatic.h"

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


Err MsCalibratomatic::init(
        const PythiaParameters &pythiaParameters,
        FragLibraryTronDIA *fragLibraryTronDia
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

    m_fragLibraryTronDia = fragLibraryTronDia;

    ERR_RETURN
}

Err MsCalibratomatic::exec(
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

        e = processLogic(
                scoreVectorsFilePath,
                msFrameScansFilePath
                ); ree;

    }

    ERR_RETURN
}

Err MsCalibratomatic::processLogic(
        const QString &scoreVectorsFilePath,
        const QString &msFrameScansFilePath
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_fragLibraryTronDia->isInit()); ree;

    ParquetReader::read(
            scoreVectorsFilePath,
            &m_scoreVectors
    ); ree;

    QVector<MsFrameScanPointRows> msFrameScanPointRows;
    ParquetReader::read(
            msFrameScansFilePath,
            &msFrameScanPointRows
    ); ree;

    e = buildFrameIndexVsScanPoints(msFrameScanPointRows); ree;

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

Err MsCalibratomatic::buildFrameIndexVsScanPoints(const QVector<MsFrameScanPointRows> &msFrameScanPointRows) {

    ERR_INIT

    m_frameIndexVsScanPoints.clear();
    e = ErrorUtils::isNotEmpty(msFrameScanPointRows); ree;

    for (const MsFrameScanPointRows &row : msFrameScanPointRows) {

        ScanPoints scanPoints;
        e = ParallelUtils::zip(
                row.mzVals,
                row.intensityVals,
                &scanPoints
        ); ree

        m_frameIndexVsScanPoints.insert(row.frameIndex, scanPoints);
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



    }


    ERR_RETURN
}
