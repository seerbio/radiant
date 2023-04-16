//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretronProcessormatic.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"

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


Err MsFrameScoretronProcessormatic::init(const PythiaParameters &pythiaParameters) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    m_params = pythiaParameters;

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
//            qDebug() << "DELETE ME" << frameIndex << rtp.first.min_corner().get<0>()
//                     << rtp.first.max_corner().get<0>() << peptidesFoundInFrameIndex->back();
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
            ){
        return l.second < r.second;
    };

    void thresholdScoredPeptides(
            double threshold,
            QVector<QPair<PeptideStringWithMods, double>> *peptideScores
            ) {

        const auto terminatorLogic = [threshold](const QPair<PeptideStringWithMods, double> &s){
            return s.second < threshold;
        };

        const auto terminator = std::remove_if(
                peptideScores->begin(),
                peptideScores->end(),
                terminatorLogic
                );

        peptideScores->erase(terminator, peptideScores->end());
    }

    Err getFrameIndexPeptides(
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

        const double stDevs = 2.0;
        const double thresholdScore = MathUtils::mean(scores) + (stDevs * MathUtils::stDev(scores));

        std::sort(peptideScores->rbegin(), peptideScores->rend(), sortScoresAscLogic);

        ERR_RETURN
    }

    struct FrameIndexCandidate {
        double ogScore = -1.0;
        PeptideStringWithMods peptideStringWithMods;
        QVector<MS2Ion> tandemPredMs2Ions;
    };

    Err bindPredictionsToPeptides(
            const QVector<QPair<PeptideStringWithMods, double>> &peptidesWithScores,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &framePredictions,
            QVector<FrameIndexCandidate> *frameCadidates
            ) {

        ERR_INIT

        for (const QPair<PeptideStringWithMods, double> &pepScore : peptidesWithScores) {

            FrameIndexCandidate fic;
            fic.ogScore = pepScore.second;
            fic.peptideStringWithMods = pepScore.first;
            fic.tandemPredMs2Ions = framePredictions.value(fic.peptideStringWithMods);

//            e = ErrorUtils::isNotEmpty(fic.tandemPredMs2Ions); ree;
            frameCadidates->push_back(fic);
        }

        ERR_RETURN
    }


    Err rescoreFrameIndex(
            const ScanPoints &scanPoints,
            const QVector<FrameIndexCandidate> &frameIndexCandidates,
            double ppmTol
            ) {

        ERR_INIT

        const auto transformLogic = [](const MS2Ion &ion) {return QPointF(ion.mz, ion.intensity);};

        e = ErrorUtils::isNotEmpty(scanPoints); ree;

        for (const FrameIndexCandidate &fic : frameIndexCandidates) {

            QVector<QPointF> mzToExtract;
            std::transform(
                    fic.tandemPredMs2Ions.begin(),
                    fic.tandemPredMs2Ions.end(),
                    std::back_inserter(mzToExtract),
                    transformLogic
                    );

            const ExtractPoints extractedPoints = MsUtils::extractPointsFromPoints(
                    scanPoints,
                    mzToExtract,
                    ppmTol
            );

            qDebug() << extractedPoints.mzFoundVsSearched;
            qDebug() << "";
            qDebug() << extractedPoints.intensityFoundVsSearched;
            qDebug() << "*******";

        }

        ERR_RETURN
    }

    Err processLogic(
            const MsFrameScanPointRows &row,
            const RTree &scoresRTree,
            const QMap<PeptideStringWithMods, QVector<double>> &peptideWithModsVsScoreVec,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &framePredictions,
            double ppmTol,
            QVector<QPair<PeptideStringWithMods, double>> *_peptideScores
            ) {

        ERR_INIT

//        peptideScores->clear();

        const FrameIndex frameIndex = row.frameIndex;
        ScanPoints scanPoints;
        e = ParallelUtils::zip(
                row.mzVals,
                row.intensityVals,
                &scanPoints
                ); ree;

        QStringList peptidesFoundInFrameIndex;
        e = getScorePeptidesByFrameIndex(
                frameIndex,
                scoresRTree,
                &peptidesFoundInFrameIndex
                ); ree;

        QVector<QPair<PeptideStringWithMods, double>> peptideScores;
        e = getFrameIndexPeptides(
                frameIndex,
                peptidesFoundInFrameIndex,
                peptideWithModsVsScoreVec,
                &peptideScores
                ); ree;

        QVector<FrameIndexCandidate> frameCadidates;
        e = bindPredictionsToPeptides(
                peptideScores,
                framePredictions,
                &frameCadidates
                ); ree;

        e = rescoreFrameIndex(
                scanPoints,
                frameCadidates,
                ppmTol
                ); ree;

//        const int minThresholdCount = 10;
//        if (peptideScores->size() > minThresholdCount) {
//
//            thresholdScoredPeptides(
//                    thresholdScore,
//                    peptideScores
//            );
//
//        }

        ERR_RETURN
    }

}//namespace
Err MsFrameScoretronProcessormatic::rescoreMsFrame(
        const QVector<MsFrameScanPointRows> &msFrameScanPointRows,
        const QVector<MsFrameScoreVectorReaderRow> &scoreVectors,
        const QMap<PeptideStringWithMods, QVector<MS2Ion>> &framePredictions
        ){

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msFrameScanPointRows); ree;
    e = ErrorUtils::isNotEmpty(scoreVectors); ree;
    e = ErrorUtils::isNotEmpty(framePredictions); ree;

    RTree scoresRTree = buildScorePeaksRTree(scoreVectors);

    QMap<PeptideStringWithMods, QVector<double>> peptideWithModsVsScoreVec;
    e = buildPeptideWithModsVsScoreVec(
            scoreVectors,
            &peptideWithModsVsScoreVec
            ); ree;

    for (const MsFrameScanPointRows &row : msFrameScanPointRows) {

        if (row.frameIndex != 128) {
            continue;
        }

        QVector<QPair<PeptideStringWithMods, double>> peptideScores;
        e = processLogic(
                row,
                scoresRTree,
                peptideWithModsVsScoreVec,
                framePredictions,
                m_params.ms2ExtractionWidthPPM,
                &peptideScores
                ); ree;

        if (peptideScores.isEmpty()) {
            continue;
        }

        qDebug() << peptideScores;
    }

    ERR_RETURN
}

