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
        MsFrame *msFrame = nullptr;
        QVector<CandidateScores*> candidateScores;
        float ppmTol = -1.0;
    };

    Err buildParallelInputs(
        const QMap<MzTargetKey, QVector<CandidateScores*>> &mzTargetKeyVsCandidateScoresPntrs,
        const QMap<MzTargetKey, MsFrame*> &mzTargetKeyVsMsFrame,
        float ppmTol,
        QVector<ParallelInput> *parallelInputs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(mzTargetKeyVsCandidateScoresPntrs); ree;
        e = ErrorUtils::isNotEmpty(mzTargetKeyVsMsFrame); ree;

        parallelInputs->clear();

        for (const MzTargetKey &mzk : mzTargetKeyVsCandidateScoresPntrs.keys()) {

            e = ErrorUtils::contains(mzk, mzTargetKeyVsMsFrame); ree;

            ParallelInput pi;
            pi.targetKey = mzk;
            pi.candidateScores = mzTargetKeyVsCandidateScoresPntrs.value(pi.targetKey);
            pi.msFrame = mzTargetKeyVsMsFrame.value(pi.targetKey);
            pi.ppmTol = ppmTol;

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
        const QMap<Index, XICPoints> &indexVsXICPoints,
        MsFrame *msFrame,
        const CandidateScores *cs,
        QuanReaderRow *qr
        ) {

        ERR_INIT

        qr->peptideStringWithMods = cs->targetDecoyCandidatePair->peptideStringWithMods();
        qr->charge = cs->targetDecoyCandidatePair->charge();
        qr->targetKey = cs->targetKey;
        qr->mzSearched1 = cs->featuresArray[CandidateScores::Features::MzSearched1];
        qr->mzSearched2 = cs->featuresArray[CandidateScores::Features::MzSearched2];
        qr->mzSearched3 = cs->featuresArray[CandidateScores::Features::MzSearched3];
        qr->mzSearched4 = cs->featuresArray[CandidateScores::Features::MzSearched4];
        qr->mzSearched5 = cs->featuresArray[CandidateScores::Features::MzSearched5];
        qr->mzSearched6 = cs->featuresArray[CandidateScores::Features::MzSearched6];
        qr->mzSearched7 = cs->featuresArray[CandidateScores::Features::MzSearched7];
        qr->mzSearched8 = cs->featuresArray[CandidateScores::Features::MzSearched8];
        qr->mzSearched9 = cs->featuresArray[CandidateScores::Features::MzSearched9];
        qr->mzSearched10 = cs->featuresArray[CandidateScores::Features::MzSearched10];
        qr->mzSearched11 = cs->featuresArray[CandidateScores::Features::MzSearched11];
        qr->mzSearched12 = cs->featuresArray[CandidateScores::Features::MzSearched12];
        qr->classifierScore = static_cast<float>(cs->classifierScore);
        qr->discScore = static_cast<float>(cs->discriminantScore);
        qr->qValue = static_cast<float>(cs->qValue);
        qr->mzInterferences = cs->mzInterferences;
        qr->isDecoy = cs->isDecoy;

        if (qr->isDecoy) {
            ERR_RETURN
        }

        Eigen::MatrixX<float> xicMatrix;
        e = buildXICMatrix(indexVsXICPoints, &xicMatrix); ree;

        for (int col = 0; col < xicMatrix.cols(); col++) {

            const Eigen::VectorX<float> &matCol = xicMatrix.col(col);
            const QVector<float> matColVec = EigenUtils::convertEigenVectorToQVector(matCol);

            switch (col) {
            case 0:
                std::transform(
                    matColVec.begin(),
                    matColVec.end(),
                    std::back_inserter(qr->scanTimes),
                    [msFrame](float f){return msFrame->scanTimeFromFrameIndex(static_cast<int>(f));}
                    );
                break;
            case 1:
                qr->intensityValsMz1 = matColVec;
                break;
            case 2:
                qr->intensityValsMz2 = matColVec;
                break;
            case 3:
                qr->intensityValsMz3 = matColVec;
                break;
            case 4:
                qr->intensityValsMz4 = matColVec;
                break;
            case 5:
                qr->intensityValsMz5 = matColVec;
                break;
            case 6:
                qr->intensityValsMz6 = matColVec;
                break;
            case 7:
                qr->intensityValsMz7 = matColVec;
                break;
            case 8:
                qr->intensityValsMz8 = matColVec;
                break;
            case 9:
                qr->intensityValsMz9 = matColVec;
                break;
            case 10:
                qr->intensityValsMz10 = matColVec;
                break;
            case 11:
                qr->intensityValsMz11 = matColVec;
                break;
            case 12:
                qr->intensityValsMz12 = matColVec;
                break;
            default:
                rrr(eValueError)
            }
        }

        ERR_RETURN
    }

    QPair<Err, QVector<QuanReaderRow>> parallelLogic(const ParallelInput &pi) {

        ERR_INIT
        TurboXIC2D turboXIC;
        e = turboXIC.init(pi.msFrame->frameIndexVsScanPoints());

        QVector<QuanReaderRow> quanReaderRows;

        for (const CandidateScores *cs : pi.candidateScores) {

            e = ErrorUtils::isAboveThreshold(
                cs->frameIndexEnd,
                cs->frameIndexStart,
                ErrorUtilsParam::ExcludeThreshold
                ); rree;

            if (cs->isDecoy) {
                QMap<Index, XICPoints> indexVsXICPointsDecoy;
                QuanReaderRow quanReaderRowDecoy;
                e = buildQuanReaderRow(
                    indexVsXICPointsDecoy,
                    pi.msFrame,
                    cs,
                    &quanReaderRowDecoy
                ); rree;
                quanReaderRows.push_back(quanReaderRowDecoy);
                continue;
            }

            const int peakLength = cs->frameIndexEnd - cs->frameIndexStart;
            const int bufferLength = std::max(1, static_cast<int>(std::round(peakLength / 2.0)));

            const FrameIndex frameIndexStart = cs->frameIndexStart - bufferLength;
            const FrameIndex frameIndexEnd = cs->frameIndexEnd + bufferLength;

            QVector<MS2Ion> ms2Ions = cs->isDecoy
                            ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
                            : cs->targetDecoyCandidatePair->ms2IonsTarget();

            ms2Ions.resize(std::min(topN, ms2Ions.size()));

            QMap<Index, XICPoints> indexVsXICPoints;

            for (const MS2Ion &ms2Ion : ms2Ions) {

                const float mzTol = MathUtils::calculatePPM(ms2Ion.mz, pi.ppmTol);
                const float mzStart = ms2Ion.mz - mzTol;
                const float mzEnd = ms2Ion.mz + mzTol;

                const XICPoints xicPoints = turboXIC.extractPointsXIC(
                    mzStart,
                    mzEnd,
                    static_cast<float>(frameIndexStart),
                    static_cast<float>(frameIndexEnd)
                    );

                indexVsXICPoints.insert(indexVsXICPoints.size(), xicPoints);
            }

            const int xicPointsCount = std::accumulate(
                indexVsXICPoints.begin(),
                indexVsXICPoints.end(),
                0,
                [](int sum, const XICPoints &x){return sum + x.size();}
                );

            if(xicPointsCount < 1) {
                continue;
            }

            QuanReaderRow quanReaderRow;
            e = buildQuanReaderRow(
                indexVsXICPoints,
                pi.msFrame,
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
    const QMap<MzTargetKey, MsFrame*> &mzTargetKeyVsMsFrame,
    const QString& outputFilePath,
    float ppmTolerance
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;
    e = ErrorUtils::isNotEmpty(outputFilePath); ree;
    e = ErrorUtils::isNotEmpty(mzTargetKeyVsMsFrame); ree;

    QElapsedTimer et;
    et.start();

    const QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs
                                        = buildMzTargetKeyVsCandidateScoresPntrs(candidateScores);

    QVector<ParallelInput> parallelInputs;
    e = buildParallelInputs(
        mzTargetKeyVsCandidateScoresPntrs,
        mzTargetKeyVsMsFrame,
        ppmTolerance,
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

    qDebug() << "Quan file built in" << et.elapsed() << "mSec";

    ERR_RETURN
}
