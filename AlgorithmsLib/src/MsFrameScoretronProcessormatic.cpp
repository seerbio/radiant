//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretronProcessormatic.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "ParallelUtils.h"
#include "ParquetReader.h"
#include "TandemSpectraDeconvolvotron.h"

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


Err MsFrameScoretronProcessormatic::init(
        const QMap<PeptideStringWithMods, QVector<MS2Ion>> &fragPreds,
        const MsFrame &msFrame,
        const PythiaParameters &params,
        const QString &scoreVectorsFilePath
) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isTrue(msFrame.isValid()); ree;
    e = ErrorUtils::isNotEmpty(fragPreds); ree;
    e = ErrorUtils::fileExists(scoreVectorsFilePath); ree;

    m_params = params;
    m_fragPreds = fragPreds;
    m_msFrame = msFrame;
    m_scoreVectorsFilePath = scoreVectorsFilePath;

    ERR_RETURN
}

Err MsFrameScoretronProcessormatic::processLogicForFrameScores(
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex,
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> *topCansInFrameIndexVsDiscScore
        ) {

    ERR_INIT

    topCansInFrameIndex->clear();

    e = ErrorUtils::isTrue(m_msFrame.isValid()); ree;

    QVector<MsFrameScoreVectorReaderRow> scoreVectors;
    ParquetReader::read(
            m_scoreVectorsFilePath,
            &scoreVectors
    ); ree;
    e = ErrorUtils::isNotEmpty(scoreVectors); ree;

    e = getTopNCandidatesPerFrameIndex(
            scoreVectors,
            m_msFrame.frameIndexVsScanPoints(),
            topCansInFrameIndex
    );ree;


    e = calculateDiscriminateScoreForFrameIndexes(
            *topCansInFrameIndex,
            topCansInFrameIndexVsDiscScore
    ); ree;

    ERR_RETURN
}

Err MsFrameScoretronProcessormatic::calculateDiscriminateScoreForFrameIndexes(
        const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> &topCansInFrameIndex,
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore >>> *topCansInFrameIndexVsDiscScore
) {

    ERR_INIT

    for (auto it = topCansInFrameIndex.begin(); it != topCansInFrameIndex.end(); it++) {

        const FrameIndex frameIndex = it.key();
        const QVector<QPair<PeptideStringWithMods, Score>> &peptideStringWithModsScore = it.value();

        e = calculateDiscriminateScoreForFrame(
                peptideStringWithModsScore,
                frameIndex,
                topCansInFrameIndexVsDiscScore
        ); ree;

//#define DEBUG_DISC
#ifdef DEBUG_DISC
        if (frameIndex == 128) {
            qDebug() << frameIndex;
            const auto &discScores = topCansInFrameIndexVsDiscScore->value(frameIndex);
            for (const QPair<PeptideStringWithMods, DiscScore > &pr : discScores) {
                qDebug() << pr << m_fragPredsIsDecoy.value(pr.first) << pr.second;
            }
            qDebug() << "************";
        }
#endif

    }

    ERR_RETURN
}

namespace {

    Err getFrameIndexMs2Ions(
            double mzMax,
            const QVector<QPair<PeptideStringWithMods, Score>> &peptideStringWithModsScore,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &fragPreds,
            QMap<PeptideStringWithMods, QVector<MS2Ion>> *scanPreds
    ) {

        ERR_INIT

        for (const QPair<PeptideStringWithMods, Score> &pepScore : peptideStringWithModsScore) {

            const PeptideStringWithMods &peptideStringWithMods = pepScore.first;

            QVector<MS2Ion> ms2Ions = fragPreds.value(peptideStringWithMods);

            const double mzMin = 0.0;
            FragLibReader::filterMs2IonsByMz(
                    mzMin,
                    mzMax,
                    &ms2Ions
            );

            e = ErrorUtils::isTrue(fragPreds.contains(peptideStringWithMods)); ree;
            scanPreds->insert(peptideStringWithMods, ms2Ions);
        }

        ERR_RETURN
    }

    void sortDiscScoresDesc(
            QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> *frameIndexVsPeptideStringWithModsDiscScore
    ) {

        for (
                auto it = frameIndexVsPeptideStringWithModsDiscScore->begin();
                it != frameIndexVsPeptideStringWithModsDiscScore->end();
                it++
                ) {

            QVector<QPair<PeptideStringWithMods, Score>> &scores = it.value();
            using Sort = QPair<PeptideStringWithMods, Score>;
            std::sort(scores.rbegin(), scores.rend(), [](const Sort &l, const Sort &r){
                return l.second < r.second;
            });

        }

    }

}//namespace
Err MsFrameScoretronProcessormatic::calculateDiscriminateScoreForFrame(
        const QVector<QPair<PeptideStringWithMods, Score>> &peptideStringWithModsScore,
        FrameIndex frameIndex,
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> *frameIndexVsPeptideStringWithModsDiscScore
) {

    ERR_INIT

    QMap<PeptideStringWithMods, QVector<MS2Ion>> scanPreds;
    e = getFrameIndexMs2Ions(
            m_params.mzMaxDataStructure,
            peptideStringWithModsScore,
            m_fragPreds,
            &scanPreds
    ); ree;

    const QMap<FrameIndex, ScanPoints> &frame = m_msFrame.frameIndexVsScanPoints();
    ScanPoints scanPoints = frame.value(frameIndex);

    QMap<PeptideStringWithMods, DiscScore> pepSeqVsWeight;

    TandemSpectraDeconvolvotron deconvolvotron;
    const int precision = 2;
    const int maxIters = 50;
    const double stopTolerance = 1e-8;
    e = deconvolvotron.init(
            precision,
            m_params.mzMaxDataStructure,
            maxIters,
            stopTolerance
    ); ree;

    e = deconvolvotron.deconvolveTandemSpectra(
            scanPoints,
            scanPreds,
            &pepSeqVsWeight
    ); ree;

    for (auto itt = pepSeqVsWeight.begin(); itt != pepSeqVsWeight.end(); itt++) {

        const PeptideStringWithMods &peptideStringWithMods = itt.key();
        const double discriminateScore = itt.value();

        (*frameIndexVsPeptideStringWithModsDiscScore)[frameIndex].push_back({peptideStringWithMods, discriminateScore});
    }

    sortDiscScoresDesc(frameIndexVsPeptideStringWithModsDiscScore);

    ERR_RETURN
}

namespace {

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

    Err buildFrameIndexVsApexPeptideStringWithModsVec(
            const QVector<MsFrameScoreVectorReaderRow> &scoreVectors,
            QMap<FrameIndex, QVector<PeptideStringWithMods>> *frameIndexVsApexPeptideStringWithModsVec
            ) {

        ERR_INIT

        frameIndexVsApexPeptideStringWithModsVec->clear();

        e = ErrorUtils::isNotEmpty(scoreVectors); ree;

        for (const MsFrameScoreVectorReaderRow &row : scoreVectors){
            (*frameIndexVsApexPeptideStringWithModsVec)[row.frameIndexMaxScore].push_back(row.peptideStringWithMods);
        }

        ERR_RETURN
    }

    Err getFrameIndexPeptidesWithScoresDesc(
            FrameIndex frameIndex,
            const QVector<PeptideStringWithMods> &peptidesFoundInFrameIndex,
            const QMap<PeptideStringWithMods, QVector<double>> &peptideWithModsVsScoreVec,
            QVector<QPair<PeptideStringWithMods, double>> *peptideScores
    ) {

        ERR_INIT

        const auto sortScoresAscLogic = [](
                const QPair<PeptideStringWithMods, double> &l,
                const QPair<PeptideStringWithMods, double> &r
        ){return l.second < r.second;};

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
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> *topCansInFrameIndex
        ) {

    ERR_INIT

    topCansInFrameIndex->clear();

    QMap<FrameIndex, QVector<PeptideStringWithMods>> frameIndexVsApexPeptideStringWithModsVec;
    e = buildFrameIndexVsApexPeptideStringWithModsVec(
            scoreVectors,
            &frameIndexVsApexPeptideStringWithModsVec
            ); ree;

    QMap<PeptideStringWithMods, QVector<double>> peptideWithModsVsScoreVec;
    e = buildPeptideWithModsVsScoreVec(
            scoreVectors,
            &peptideWithModsVsScoreVec
            ); ree;

    const QList<FrameIndex> &frameIndexes = frameIndexVsScanPoints.keys();
    for (const FrameIndex frameIndex : frameIndexes) {

        if (!frameIndexVsApexPeptideStringWithModsVec.contains(frameIndex)){
            continue;
        }

        const QVector<PeptideStringWithMods> peptidesFoundInFrameIndex
                = frameIndexVsApexPeptideStringWithModsVec.value(frameIndex);

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

        //TODO double check this return here drewholio
//        qDebug() << peptideWithScoresDesc.size();
        const int topN = std::min(m_params.returnPSMTopN, peptideWithScoresDesc.size());
        peptideWithScoresDesc.resize(topN);

        topCansInFrameIndex->insert(frameIndex, peptideWithScoresDesc);

    }

    ERR_RETURN
}




