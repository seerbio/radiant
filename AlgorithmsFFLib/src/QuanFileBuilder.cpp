//
// Created by andrewnichols on 8/15/24.
//

#include "QuanFileBuilder.h"

#include "CandidateScores.h"
#include "EigenUtils.h"
#include "MsFrame.h"
#include "ParquetReader.h"
#include "TurboXIC.h"
#include "TurboXIC2D.h"
#include "QuanReader.h"

#include "Eigen/Dense"

#include <QElapsedTimer>
#include <QtConcurrent/QtConcurrent>

namespace {

    constexpr int topN = 12;

    QMap<MzTargetKey, QVector<CandidateScores*>> buildMzTargetKeyVsCandidateScoresPntrs(
        const QVector<CandidateScores*> &candidateScores
        ) {

        QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs;

        for (CandidateScores *cs : candidateScores) {
            mzTargetKeyVsCandidateScoresPntrs[cs->targetKey].push_back(cs);
        }

        return mzTargetKeyVsCandidateScoresPntrs;
    }

    struct ParallelInput {
        MzTargetKey targetKey;
        QVector<CandidateScores*> candidateScores;
    };

    Err buildParallelInputs(
        const QMap<MzTargetKey, QVector<CandidateScores*>> &mzTargetKeyVsCandidateScoresPntrs,
        QVector<ParallelInput> *parallelInputs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzTargetKeyVsCandidateScoresPntrs); ree;

        parallelInputs->clear();

        for (const MzTargetKey &mzk : mzTargetKeyVsCandidateScoresPntrs.keys()) {

            ParallelInput pi;
            pi.targetKey = mzk;
            pi.candidateScores = mzTargetKeyVsCandidateScoresPntrs.value(pi.targetKey);

            parallelInputs->push_back(pi);
        }

        ERR_RETURN
    }

    QPair<FrameIndex, FrameIndex> findFrameIndexMinMax(const QMap<Index, XICPoints> &indexVsXICPoints) {

        QPair<FrameIndex, FrameIndex> frameIndexMinMax = {std::numeric_limits<FrameIndex>::max(), 0};

        for (const XICPoints &xicPoints : indexVsXICPoints) {

            if (xicPoints.empty()) {
                continue;
            }

            const auto minMaxLogic = [](const XICPoint &l, const XICPoint &r) {return l.scanNumber < r.scanNumber;};

            const auto frameIndexMinMaxIter = std::minmax_element(
                xicPoints.begin(),
                xicPoints.end(),
                minMaxLogic
                );

            frameIndexMinMax.first = std::min(frameIndexMinMax.first, frameIndexMinMaxIter.first->scanNumber);
            frameIndexMinMax.second = std::max(frameIndexMinMax.second, frameIndexMinMaxIter.second->scanNumber);
        }

        return frameIndexMinMax;
    }

    Err buildXICMatrix(
        const QMap<Index, XICPoints> &indexVsXICPoints,
        Eigen::MatrixX<float> *xicMatrix
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(indexVsXICPoints); ree;

        const QPair<FrameIndex, FrameIndex> frameIndexMinMax = findFrameIndexMinMax(indexVsXICPoints);

        e = ErrorUtils::isAboveThreshold(
            frameIndexMinMax.second,
            frameIndexMinMax.first,
            ErrorUtilsParam::IncludeThreshold
            );
        if (e != eNoError) {
            qDebug() << indexVsXICPoints;
        }; ree;

        const int cols = indexVsXICPoints.size() + 1;
        const int rows = frameIndexMinMax.second - frameIndexMinMax.first + 1;

        xicMatrix->resize(rows, cols);
        xicMatrix->setZero();

        QVector<FrameIndex> frameIndexes(rows);
        std::iota(frameIndexes.begin(), frameIndexes.end(), frameIndexMinMax.first);

        const Eigen::VectorX<float> frameIndexesVec = EigenUtils::convertQVectorToEigenVector(frameIndexes).cast<float>();

        xicMatrix->col(0) = frameIndexesVec;

        for (auto it = indexVsXICPoints.begin(); it != indexVsXICPoints.end(); ++it) {

            const int col = it.key() + 1;
            const XICPoints &xicPoints = it.value();

            for (const XICPoint &xicPoint : xicPoints) {
                const int row = xicPoint.scanNumber - frameIndexMinMax.first;
                xicMatrix->coeffRef(row, col) = xicPoint.intensity;
            }
        }

        ERR_RETURN
    }

    Err buildQuanReaderRow(
        const CandidateScores *cs,
        QuanReaderRow *qr
        ) {

        ERR_INIT

        const QVector<MS2Ion> ms2Ions = cs->isDecoy
                                      ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
                                      : cs->targetDecoyCandidatePair->ms2IonsTarget();

        constexpr int maxIonsSize = 12;
        QVector<float> mzSearchedVals(maxIonsSize, -1.0f);

        for (int i = 0; i < std::min(maxIonsSize, ms2Ions.size()); i++) {
            mzSearchedVals[i] = ms2Ions.at(i).mz;
        }

        qr->peptideStringWithMods = cs->targetDecoyCandidatePair->peptideStringWithMods();
        qr->charge = cs->targetDecoyCandidatePair->charge();
        qr->targetKey = cs->targetKey;
        qr->mzSearched1 = mzSearchedVals.at(0);
        qr->mzSearched2 = mzSearchedVals.at(1);
        qr->mzSearched3 = mzSearchedVals.at(2);
        qr->mzSearched4 = mzSearchedVals.at(3);
        qr->mzSearched5 = mzSearchedVals.at(4);
        qr->mzSearched6 = mzSearchedVals.at(5);
        qr->mzSearched7 = mzSearchedVals.at(6);
        qr->mzSearched8 = mzSearchedVals.at(7);
        qr->mzSearched9 = mzSearchedVals.at(8);
        qr->mzSearched10 = mzSearchedVals.at(9);
        qr->mzSearched11 = mzSearchedVals.at(10);
        qr->mzSearched12 = mzSearchedVals.at(11);

        qr->classifierScore = static_cast<float>(cs->classifierScore);
        qr->discScore = static_cast<float>(cs->discriminantScore);
        qr->qValue = static_cast<float>(cs->qValue);
        qr->mzInterferences = cs->mzInterferences;
        qr->isDecoy = cs->isDecoy;

        qr->scanTime = cs->scanTime;
        qr->scanTimeStart = cs->scanTimeStart;
        qr->scanTimeEnd = cs->scanTimeEnd;
        qr->cosineSimSum100 = cs->featuresArray[Features::CosineSimSum100Top12];
        qr->intensityValMz1 = cs->featuresArray[Features::IntensityFoundMax1];
        qr->intensityValMz2 = cs->featuresArray[Features::IntensityFoundMax2];
        qr->intensityValMz3 = cs->featuresArray[Features::IntensityFoundMax3];
        qr->intensityValMz4 = cs->featuresArray[Features::IntensityFoundMax4];
        qr->intensityValMz5 = cs->featuresArray[Features::IntensityFoundMax5];
        qr->intensityValMz6 = cs->featuresArray[Features::IntensityFoundMax6];
        qr->intensityValMz7 = cs->featuresArray[Features::IntensityFoundMax7];
        qr->intensityValMz8 = cs->featuresArray[Features::IntensityFoundMax8];
        qr->intensityValMz9 = cs->featuresArray[Features::IntensityFoundMax9];
        qr->intensityValMz10 = cs->featuresArray[Features::IntensityFoundMax10];
        qr->intensityValMz11 = cs->featuresArray[Features::IntensityFoundMax11];
        qr->intensityValMz12 = cs->featuresArray[Features::IntensityFoundMax12];

        qr->scanTimeAlt = cs->scanTimeAlt;
        qr->scanTimeStartAlt = cs->scanTimeStartAlt;
        qr->scanTimeEndAlt = cs->scanTimeEndAlt;
        qr->cosineSimSum100Alt = cs->cosineSimSum100Alt;
        qr->scanTimePredicted = cs->scanTimePredicted;

        for (int i = 0; i < cs->intensityValsAlt.size(); i++) {

            switch (i) {
                case 0:
                    qr->intensityValMzAlt1 = cs->intensityValsAlt.at(i);
                    break;
                case 1:
                    qr->intensityValMzAlt2 = cs->intensityValsAlt.at(i);
                    break;
                case 2:
                    qr->intensityValMzAlt3 = cs->intensityValsAlt.at(i);
                    break;
                case 3:
                    qr->intensityValMzAlt4 = cs->intensityValsAlt.at(i);
                    break;
                case 4:
                    qr->intensityValMzAlt5 = cs->intensityValsAlt.at(i);
                    break;
                case 5:
                    qr->intensityValMzAlt6 = cs->intensityValsAlt.at(i);
                    break;
                case 6:
                    qr->intensityValMzAlt7 = cs->intensityValsAlt.at(i);
                    break;
                case 7:
                    qr->intensityValMzAlt8 = cs->intensityValsAlt.at(i);
                    break;
                case 8:
                    qr->intensityValMzAlt9 = cs->intensityValsAlt.at(i);
                    break;
                case 9:
                    qr->intensityValMzAlt10 = cs->intensityValsAlt.at(i);
                    break;
                case 10:
                    qr->intensityValMzAlt11 = cs->intensityValsAlt.at(i);
                    break;
                case 11:
                    qr->intensityValMzAlt12 = cs->intensityValsAlt.at(i);
                    break;
                default:
                    rrr(eValueError);
            }

        }

        ERR_RETURN
    }

    QPair<Err, QVector<QuanReaderRow>> parallelLogic(const ParallelInput &pi) {

        ERR_INIT

        QVector<QuanReaderRow> quanReaderRows;

        for (const CandidateScores *cs : pi.candidateScores) {

            QuanReaderRow quanReaderRow;
            e = buildQuanReaderRow(
                cs,
                &quanReaderRow
            ); rree;

            quanReaderRows.push_back(quanReaderRow);
        }

        return {e, quanReaderRows};
    }

}//namespace
Err QuanFileBuilder::buildQuanFile(
    const QVector<CandidateScores*> &candidateScores,
    const QString& outputFilePath
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;
    e = ErrorUtils::isNotEmpty(outputFilePath); ree;

    QElapsedTimer et;
    et.start();

    const QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs
                                        = buildMzTargetKeyVsCandidateScoresPntrs(candidateScores);

    QVector<ParallelInput> parallelInputs;
    e = buildParallelInputs(
        mzTargetKeyVsCandidateScoresPntrs,
        &parallelInputs
        ); ree;

    QVector<QuanReaderRow> quanReaderRows;

#define RUN_PARALLEL_BUILDER
#ifdef RUN_PARALLEL_BUILDER
    QFuture<QPair<Err, QVector<QuanReaderRow>>> futures = QtConcurrent::mapped(
        parallelInputs,
        parallelLogic
        );
    futures.waitForFinished();

    for (const QPair<Err, QVector<QuanReaderRow>> &pr : futures) {
        e = pr.first; ree;
        quanReaderRows.append(pr.second);
    }
#else
    for (const ParallelInput &pi : parallelInputs) {
        QPair<Err, QVector<QuanReaderRow>> errVsQuanReaderRows = parallelLogic(pi);
        e = errVsQuanReaderRows.first;
        quanReaderRows.append(errVsQuanReaderRows.second);
    }
#endif

    e = ParquetReader::write(quanReaderRows, outputFilePath); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "Quan file built in"
            << et.elapsed()
            << "mSec";

    ERR_RETURN
}
