//
// Created by anichols on 4/14/23.
//

#include "MsFrameScoretron.h"

#include "BiophysicalCalcs.h"
#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MsFrameScoretronProcessormatic.h"
#include "MsFrameScoreVectorReader.h"
#include "MsReaderParquet.h"
#include "ParallelUtils.h"
#include "PeakIntegratomatic.h"
#include "TandemSpectraDeconvolvotron.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>


namespace {

    Err peptideStringWithModsFromPeptideSequenceChargeKey(
            const PeptideSequenceChargeKey &peptideSequenceChargeKey,
            PeptideStringWithMods *peptideStringWithMods
            ){

        ERR_INIT

        const int expectedSplitSize = 2;

        const QStringList peptideSequenceChargeKeySplit = peptideSequenceChargeKey.split(
                S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP,
                QString::SkipEmptyParts
                );

        e = ErrorUtils::isEqual(
                peptideSequenceChargeKeySplit.size(),
                expectedSplitSize); ree;

        *peptideStringWithMods = peptideSequenceChargeKeySplit.front();

        ERR_RETURN
    }

    void filterMs2IonsByMz(
            double mzStart,
            double mzEnd,
            QVector<MS2Ion> *ms2Ions
            ) {

        const auto terminatorLogic = [mzStart, mzEnd](const MS2Ion &ion){
            return !(mzStart <= ion.mz && ion.mz <= mzEnd);
        };

        const auto terminator = std::remove_if(
                ms2Ions->begin(),
                ms2Ions->end(),
                terminatorLogic
                );

        ms2Ions->erase(terminator, ms2Ions->end());
    }

    void getTopNMostIntenseMs2Ions(
            int topNMs2Ions,
            QVector<MS2Ion> *ms2Ions
            ) {

        const auto sortIntensityAsc = [](const MS2Ion &l, const MS2Ion &r){
            return l.intensity < r.intensity;
        };

        std::sort(ms2Ions->rbegin(), ms2Ions->rend(), sortIntensityAsc);

        topNMs2Ions = std::min(topNMs2Ions, ms2Ions->size());

        ms2Ions->resize(topNMs2Ions);

        const auto sortMzAsc = [](const MS2Ion &l, const MS2Ion &r) {
            return l.mz < r.mz;
        };

        std::sort(ms2Ions->begin(), ms2Ions->end(), sortMzAsc);
    }

    struct ScoringMatrices {
        Eigen::MatrixX<double> scoringMatrixMz;
        Eigen::MatrixX<double> scoringMatrixIntensity;
        Eigen::MatrixX<double> theoMatrixMatrixIntensity;
    };

    Eigen::MatrixX<double> buildTheoMat(
            const QVector<MS2Ion> &ms2Ions,
            int frameScanCount,
            int topNMs2Ions
    ) {

        Eigen::MatrixX<double> theoMat(frameScanCount, topNMs2Ions);
        theoMat.setZero();

        Eigen::VectorX<double> theoIntensities(topNMs2Ions);
        for (int i = 0; i < topNMs2Ions; i++) {
            theoIntensities(i) = ms2Ions.at(i).intensity;
        }
        for (int i = 0; i < frameScanCount; i++) {
            theoMat.row(i) = theoIntensities;
        }

        return theoMat;
    }

    ScoringMatrices buildTargetScoringMatrices(
            const PythiaParameters &params,
            const QVector<MS2Ion> &ms2IonsTandemPred,
            int frameScanCount,
            int topNMs2Ions,
            double ppmTolerance,
            TurboXIC *turboXic
    ) {

        Eigen::MatrixX<double> scoringMatIntensity(frameScanCount, topNMs2Ions);
        scoringMatIntensity.setZero();

        Eigen::MatrixX<double> scoringMatMz(frameScanCount, topNMs2Ions);
        scoringMatMz.setOnes();

        QVector<MS2Ion> ms2Ions = ms2IonsTandemPred;

        filterMs2IonsByMz(
                params.mzMinDataStructure,
                params.mzMaxDataStructure,
                &ms2Ions
                );

        getTopNMostIntenseMs2Ions(
                params.topNMs2Ions,
                &ms2Ions
                );

        const Eigen::MatrixX<double> theoMat = buildTheoMat(
                ms2Ions,
                frameScanCount,
                topNMs2Ions
        );

        for (int colIdx = 0; colIdx < ms2Ions.size(); colIdx++) {

            const double mz = ms2Ions.at(colIdx).mz;
            const double massTol = MathUtils::calculatePPM(
                    mz,
                    ppmTolerance
            );

            const double mzMin = mz - massTol;
            const double mzMax = mz + massTol;
            const FrameIndex frameIndexMin = 0;
            const FrameIndex frameIndexMax = frameScanCount;

            const XICPoints xic = turboXic->extractPoints(
                    mzMin,
                    mzMax,
                    frameIndexMin,
                    frameIndexMax
            );

            if (xic.isEmpty()) {
                continue;
            }

            for (auto it = xic.begin(); it != xic.end(); it++) {
                const FrameIndex frameIndex = it.key();
                const Intensity intensity = it.value();

                scoringMatIntensity(frameIndex, colIdx) = intensity;
                scoringMatMz(frameIndex, colIdx) = mz;
            }
        }

        ScoringMatrices sm;
        sm.scoringMatrixIntensity = scoringMatIntensity;
        sm.theoMatrixMatrixIntensity = theoMat;
        sm.scoringMatrixMz = scoringMatMz;

        return sm;
    }

    Err buildPerFrameIndexScoreVectors(
            const ScoringMatrices &scoringMatrices,
            FrameIndexScoreResultOfTarget *frameIndexScoreResultOfTarget
    ) {

        ERR_INIT

        const double numberOfSearchedFrags = static_cast<double>(scoringMatrices.scoringMatrixIntensity.cols());

        const Eigen::VectorX<double> cosineSimPerFrameIndexOfTarget = EigenUtils::rowWiseCosineSimilarOfMatrices(
                scoringMatrices.scoringMatrixIntensity,
                scoringMatrices.theoMatrixMatrixIntensity
        );

        const Eigen::VectorX<double> mzProductPerFrameIndexOfTarget = scoringMatrices.scoringMatrixMz.rowwise().prod();

        const Eigen::VectorX<double> intensityPerFrameIndexOfTarget = scoringMatrices.scoringMatrixIntensity.rowwise().sum();

        const Eigen::VectorX<int> foundIonsPerFrameIndexOfTarget
                = (scoringMatrices.scoringMatrixIntensity.array() > 0.0).rowwise().count().cast<int>();

        const Eigen::VectorX<double> fractionFoundIonsPerFrameIndexOfTarget
                = foundIonsPerFrameIndexOfTarget.cast<double>() / numberOfSearchedFrags;

        const Eigen::VectorX<double> foundIonsPerFrameIndexOfTargetFactorial
                = foundIonsPerFrameIndexOfTarget.unaryExpr([](int x) { return MathUtils::factorial(x); }).cast<double>();

        const Eigen::VectorX<double> scorePerFrameIndexOfTarget = (foundIonsPerFrameIndexOfTargetFactorial.array()
                                                                  * fractionFoundIonsPerFrameIndexOfTarget.array()
                                                                  * mzProductPerFrameIndexOfTarget.array().sqrt()
                                                                  * cosineSimPerFrameIndexOfTarget.array()).sqrt().sqrt();

        frameIndexScoreResultOfTarget->cosineSimPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(cosineSimPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget->intensityPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(intensityPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget->foundIonsPerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(foundIonsPerFrameIndexOfTarget);

        frameIndexScoreResultOfTarget->scorePerFrameIndexOfTargetVec
                = EigenUtils::convertEigenVectorToQVector(scorePerFrameIndexOfTarget);

        const double stopThresholdFraction = 0.1;
        const int filterLength = 7;
        const int smoothCount = 1;
        e = PeakIntegratomatic::simpleIntegrator(
                frameIndexScoreResultOfTarget->scorePerFrameIndexOfTargetVec,
                stopThresholdFraction,
                filterLength,
                smoothCount,
                &frameIndexScoreResultOfTarget->bestScorePeakLimits
        ); ree;

        ERR_RETURN
    }

    Err scoreFrameTargets(
            const MsFrame &frame,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPreds,
            const PythiaParameters &params,
            QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
            ) {

        ERR_INIT

        qDebug() << "Building matrix" << frame.uniqueMsInfoScanKey();
        qDebug() << "Frame size" << frame.scanCount();

        const QMap<FrameIndex, ScanPoints> scanPoints = frame.frameIndexVsScanPoints();

        TurboXIC turboXic;
        e = turboXic.init(scanPoints); ree;

        for (auto it = tandemPreds.begin(); it != tandemPreds.end(); it++) {

            const PeptideStringWithMods &peptideStringWithMods = it.key();
            const QVector<MS2Ion> &ms2IonsTandemPred = it.value();

            const ScoringMatrices targetScoringMatrices = buildTargetScoringMatrices(
                    params,
                    ms2IonsTandemPred,
                    frame.scanCount(),
                    params.topNMs2Ions,
                    params.ms2ExtractionWidthPPM,
                    &turboXic
            );

            FrameIndexScoreResultOfTarget frameIndexScoreResultOfTarget;
            e = buildPerFrameIndexScoreVectors(
                    targetScoringMatrices,
                    &frameIndexScoreResultOfTarget
            ); ree;

            frameIndexScoreResultOfTarget.charge = BiophysicalCalcs::calculateChargeFromSequence(
                    peptideStringWithMods,
                    params.aminoAcids,
                    frame.meanPrecursorRange()
                    );

            pepStrWModsVsFrameIndexScoreResultOfTargets->insert(peptideStringWithMods, frameIndexScoreResultOfTarget);
        }

        ERR_RETURN
    }



}//namespace
QPair<Err, QVector<PSMsReaderRow>> MsFrameScoretron::scoreCandidates(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop,
        bool applySmooth2D
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(params.isValid()); rree;
    m_params = params;

    e = buildFragIonLibForTargetMz(
            params,
            fragLibFilePath,
            mzTargetStartStop,
            &m_fragPreds,
            &m_fragPredsIsDecoy
            ); rree;
    e = ErrorUtils::isNotEmpty(m_fragPreds); rree;

    e = buildMsFrame(
            msDataFilePath,
            uniqueMsInfoScanKey,
            params,
            mzTargetStartStop,
            applySmooth2D,
            &m_msFrame
            ); rree;

//#define WRITE_FRAME
#ifdef WRITE_FRAME
//    for (int i = 0; i < 412; i++) {
//        qDebug() << i << m_msFrame.scanNumberFromFrameIndex(i);
//    }
    const QString framesFilePath = msDataFilePath + ".frameScans";
    e = m_msFrame.writeFrameScans(framesFilePath);
    qDebug() << e << "Drewholio" << framesFilePath;
#endif

    qDebug() << uniqueMsInfoScanKey << m_msFrame.scanCount() << "Candidates:" << m_fragPreds.size();

    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> pepStrWModsVsFrameIndexScoreResultOfTargets;
    e = processFrameLogic(
            {m_msFrame, m_fragPreds},
            params,
            &pepStrWModsVsFrameIndexScoreResultOfTargets
    ); rree

    filterByFoundMzCount(
            params.minFoundMzPeaks,
            &pepStrWModsVsFrameIndexScoreResultOfTargets
            );

    const QString outputFilePathFrameScores = msDataFilePath + "." + m_msFrame.uniqueMsInfoScanKey() + ".frameScores";
    e = writeFrameTargetScoreVectors(
            pepStrWModsVsFrameIndexScoreResultOfTargets,
            outputFilePathFrameScores
    ); rree;


    /////////

    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> topCansInFrameIndexVsScore;
    e = MsFrameScoretronProcessormatic::processLogicForFrameScores(
            outputFilePathFrameScores,
            m_msFrame,
            params.returnPSMTopN,
            &topCansInFrameIndexVsScore
            ); rree;

    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> topCansInFrameIndexVsDiscScore;
    e = calculateDiscriminateScoreForFrameIndexes(
            topCansInFrameIndexVsScore,
            &topCansInFrameIndexVsDiscScore
            ); rree;

    QVector<PSMsReaderRow> psmsReaderRows;
    e = buildPSMsReaderRows(
            pepStrWModsVsFrameIndexScoreResultOfTargets,
            topCansInFrameIndexVsScore,
            topCansInFrameIndexVsDiscScore,
            &psmsReaderRows
            ); rree;

    return {e, psmsReaderRows};
}

Err MsFrameScoretron::buildFragIonLibForTargetMz(
        const PythiaParameters &params,
        const QString &fragLibUri,
        const QPair<double, double> &mzTargetStartStop,
        QMap<PeptideStringWithMods, QVector<MS2Ion>> *outputIons,
        QMap<PeptideStringWithMods, bool> *outputDecoy
) {

    ERR_INIT

    const int monoOffset = 0;
    for (Charge charge = params.chargeStateMin; charge <= params.chargeStateMax; ++charge) {

        const double massStart = BiophysicalCalcs::calculateMassFromThomson(
                mzTargetStartStop.first,
                charge,
                monoOffset
        );

        const double massEnd = BiophysicalCalcs::calculateMassFromThomson(
                mzTargetStartStop.second,
                charge,
                monoOffset
        );

        FragLibReader fragLibReader;
        e = fragLibReader.init(fragLibUri); ree;

        QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> peptideSequenceChargeKeyVsMS2Ions;
        QMap<PeptideSequenceChargeKey, bool> peptideSequenceChargeKeyVsIsDecoy;
        e = fragLibReader.getMS2Ions(
                massStart,
                massEnd,
                &peptideSequenceChargeKeyVsMS2Ions,
                &peptideSequenceChargeKeyVsIsDecoy
        ); ree;

        for (auto it = peptideSequenceChargeKeyVsMS2Ions.begin(); it != peptideSequenceChargeKeyVsMS2Ions.end(); it++) {

            const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
            QVector<MS2Ion> &ms2Ions = it.value();

            PeptideStringWithMods peptideStringWithMods;
            e = peptideStringWithModsFromPeptideSequenceChargeKey(
                    peptideSequenceChargeKey,
                    &peptideStringWithMods
            ); ree;

            e = ErrorUtils::isTrue(peptideSequenceChargeKeyVsIsDecoy.contains(peptideSequenceChargeKey)); ree;

            outputIons->insert(peptideStringWithMods, ms2Ions);
            outputDecoy->insert(peptideStringWithMods, peptideSequenceChargeKeyVsIsDecoy.value(peptideSequenceChargeKey));
        }

    }

    ERR_RETURN
}

Err MsFrameScoretron::buildMsFrame(
        const QString &msDataFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const PythiaParameters &params,
        QPair<double, double> mzTargetStartStop,
        bool applySmooth2D,
        MsFrame *msFrame
) {

    ERR_INIT

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(
            msDataFilePath,
            MsParquetReaderNamespace::PERCURSOR_TARGET_MZ,
            {mzTargetStartStop.first, mzTargetStartStop.second}
    ); ree;

    const QMap<ScanNumber, ScanPoints> targetScanPoints = msReaderParquet.getScanPoints();
    e = ErrorUtils::isNotEmpty(targetScanPoints); ree;

    e = msFrame->init(
            params,
            uniqueMsInfoScanKey,
            targetScanPoints,
            mzTargetStartStop
    ); ree;

    if (applySmooth2D) {
//        e = msFrame->denoiseFrame();
//        e = msFrame->deisotopeFrame();
//        e = msFrame->smoothFrame();
//        e = msFrame->gaussianSmooth2D(); ree;
    }

    ERR_RETURN
}

Err MsFrameScoretron::processFrameLogic(
        const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk,
        const PythiaParameters &params,
        QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
) {

    ERR_INIT

    const MsFrame &frame = chunk.first;
    const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPreds = chunk.second;

    const QPair<double, double> &mzTargetStartEnd = frame.precursorMzTargetStartEnd();
    qDebug() << "Processing window" << mzTargetStartEnd.first << mzTargetStartEnd.second;

    e = scoreFrameTargets(
            frame,
            tandemPreds,
            params,
            pepStrWModsVsFrameIndexScoreResultOfTargets
    ); ree;

    ERR_RETURN
}

void MsFrameScoretron::filterByFoundMzCount(
        int minFoundMzCount,
        QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> *pepStrWModsVsFrameIndexScoreResultOfTargets
) {

    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> pepStrWModsVsFrameIndexScoreResultOfTargetsFiltered;
    for (
            auto it = pepStrWModsVsFrameIndexScoreResultOfTargets->begin();
            it != pepStrWModsVsFrameIndexScoreResultOfTargets->end();
            it++
            ) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const FrameIndexScoreResultOfTarget &frameIndexScoreResultOfTarget = it.value();

        const QVector<int> &founds = frameIndexScoreResultOfTarget.foundIonsPerFrameIndexOfTargetVec;
        const int maxFoundCount = *std::max_element(founds.begin(), founds.end());

        if (maxFoundCount < minFoundMzCount) {
            continue;
        }

        pepStrWModsVsFrameIndexScoreResultOfTargetsFiltered.insert(peptideStringWithMods, frameIndexScoreResultOfTarget);
    }

    *pepStrWModsVsFrameIndexScoreResultOfTargets = pepStrWModsVsFrameIndexScoreResultOfTargetsFiltered;
}

Err MsFrameScoretron::writeFrameTargetScoreVectors(
        const QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> &pepStrWModsVsFrameIndexScoreResultOfTargets,
        const QString &outputFilePath
) {

    ERR_INIT

    int counter = 0;
    QVector<MsFrameScoreVectorReaderRow> msFrameScoreVectorReaderRows;
    for (auto it = pepStrWModsVsFrameIndexScoreResultOfTargets.begin(); it != pepStrWModsVsFrameIndexScoreResultOfTargets.end(); it++) {

        const PeptideStringWithMods &peptideStringWithMods = it.key();
        const FrameIndexScoreResultOfTarget &frameIndexScoreResultOfTarget = it.value();

        MsFrameScoreVectorReaderRow row;
        row.peptideStringWithMods = peptideStringWithMods;
        row.cosineSimPerFrameIndexOfTargetVec = frameIndexScoreResultOfTarget.cosineSimPerFrameIndexOfTargetVec;
        row.foundIonsPerFrameIndexOfTargetVec = frameIndexScoreResultOfTarget.foundIonsPerFrameIndexOfTargetVec;
        row.scorePerFrameIndexOfTargetVec = frameIndexScoreResultOfTarget.scorePerFrameIndexOfTargetVec;
        row.intensityPerFrameIndexOfTargetVec = frameIndexScoreResultOfTarget.intensityPerFrameIndexOfTargetVec;
        row.scorePeakStart = frameIndexScoreResultOfTarget.bestScorePeakLimits.first;
        row.scorePeakEnd = frameIndexScoreResultOfTarget.bestScorePeakLimits.second;
        row.charge = frameIndexScoreResultOfTarget.charge;
        msFrameScoreVectorReaderRows.push_back(row);
        counter++;
    }

    e = ParquetReader::write(
            msFrameScoreVectorReaderRows,
            outputFilePath
    ); ree;

    qDebug() << "Found Targets:" << counter;
    qDebug() << "Frame written to:" << outputFilePath;

    ERR_RETURN
}

Err MsFrameScoretron::calculateDiscriminateScoreForFrameIndexes(
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
            const PythiaParameters &params,
            const QVector<QPair<PeptideStringWithMods, Score>> &peptideStringWithModsScore,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &fragPreds,
            QMap<PeptideStringWithMods, QVector<MS2Ion>> *scanPreds
            ) {

        ERR_INIT

        for (const QPair<PeptideStringWithMods, Score> &pepScore : peptideStringWithModsScore) {

            const PeptideStringWithMods &peptideStringWithMods = pepScore.first;

            QVector<MS2Ion> ms2Ions = fragPreds.value(peptideStringWithMods);

            const double mzMin = 0.0;
            filterMs2IonsByMz(
                    mzMin,
                    params.mzMaxDataStructure,
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
Err MsFrameScoretron::calculateDiscriminateScoreForFrame(
        const QVector<QPair<PeptideStringWithMods, Score>> &peptideStringWithModsScore,
        FrameIndex frameIndex,
        QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> *frameIndexVsPeptideStringWithModsDiscScore
        ) {

    ERR_INIT

    QMap<PeptideStringWithMods, QVector<MS2Ion>> scanPreds;
    e = getFrameIndexMs2Ions(
            m_params,
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

    void deleteUnfoundPoints(QVector<QPointF> *points) {

        const auto terminatorLogic = [](const QPointF &p){
            return p.x() < 0;
        };

        const auto terminator = std::remove_if(
                points->begin(),
                points->end(),
                terminatorLogic
                );

        points->erase(terminator, points->end());
    }

    QVector<QPointF> buildExtractPoints(
            const ScanPoints &scanPoints,
            double mzMax,
            double extractPPM,
            bool extractIntensity,
            QVector<MS2Ion> *ms2Ions
            ) {

        const double mzMin = 0.0;

        filterMs2IonsByMz(
                mzMin,
                mzMax,
                ms2Ions
        );

        QVector<QPointF> ms2IonsQPoints;
        std::transform(
                ms2Ions->begin(),
                ms2Ions->end(),
                std::back_inserter(ms2IonsQPoints),
                [](const MS2Ion &ion){return QPointF(ion.mz, ion.intensity);}
        );

        ExtractPoints extractPoints = MsUtils::extractPointsFromPoints(
                scanPoints,
                ms2IonsQPoints,
                extractPPM
                );

        if (extractIntensity) {
            deleteUnfoundPoints(&extractPoints.intensityFoundVsSearched);
            return extractPoints.intensityFoundVsSearched;
        }

        deleteUnfoundPoints(&extractPoints.mzFoundVsSearched);
        return extractPoints.mzFoundVsSearched;
    }

}//namespace
Err MsFrameScoretron::buildPSMsReaderRows(
        const QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> &pepStrWModsVsFrameIndexScoreResultOfTargets,
        const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> &topCansInFrameIndex,
        const QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, DiscScore>>> &topCansInFrameIndexVsDiscScore,
        QVector<PSMsReaderRow> *psmsReaderRows
) {

    ERR_INIT

    const QList<FrameIndex> &scoreKeys = topCansInFrameIndex.keys();
    const QList<FrameIndex> &discScoreKeys = topCansInFrameIndexVsDiscScore.keys();

    e = ErrorUtils::isNotEmpty(scoreKeys); ree;
    e = ErrorUtils::isTrue(scoreKeys == discScoreKeys); ree;
    
    for (FrameIndex frameIndex : scoreKeys) {

        QMap<PeptideStringWithMods, PSMsReaderRow> psmReaderRowsForFrame;

        const QVector<QPair<PeptideStringWithMods, Score>> &frameCandidateScores
                = topCansInFrameIndex.value(frameIndex);

        int counter = 1;
        for (const QPair<PeptideStringWithMods, Score> &scorePr : frameCandidateScores) {

            const PeptideStringWithMods &peptideStringWithMods = scorePr.first;
            const FrameIndexScoreResultOfTarget &frameIndexScoreResultOfTarget
                    = pepStrWModsVsFrameIndexScoreResultOfTargets.value(peptideStringWithMods);

            PSMsReaderRow &psmsReaderRow = psmReaderRowsForFrame[peptideStringWithMods];
            psmsReaderRow.score = scorePr.second;
            psmsReaderRow.frameRankScore = counter++;
            psmsReaderRow.frameIndex = frameIndex;
            psmsReaderRow.scanNumber = m_msFrame.scanNumberFromFrameIndex(frameIndex);
            psmsReaderRow.charge = frameIndexScoreResultOfTarget.charge;
            psmsReaderRow.peptideStringWithMods = peptideStringWithMods;
            psmsReaderRow.isDecoy = m_fragPredsIsDecoy.value(peptideStringWithMods);
            psmsReaderRow.uniqueMsInfoScanKey = m_msFrame.uniqueMsInfoScanKey();

            const ScanPoints &scanPoints = m_msFrame.getScanPointsByScanNumber(psmsReaderRow.scanNumber);

            QVector<MS2Ion> ms2Ions = m_fragPreds.value(peptideStringWithMods);
            const QVector<QPointF> mzFoundVsTheo = buildExtractPoints(
                    scanPoints,
                    m_params.mzMaxDataStructure,
                    m_params.ms2ExtractionWidthPPM,
                    false,
                    &ms2Ions
                    );
            const QPair<QVector<double>, QVector<double>> extractedMz = ParallelUtils::unZip(mzFoundVsTheo);
            psmsReaderRow.mzFound = extractedMz.first;
            psmsReaderRow.mzSearched = extractedMz.second;


            const QVector<QPointF> intensityFoundVsTheo = buildExtractPoints(
                    scanPoints,
                    m_params.mzMaxDataStructure,
                    m_params.ms2ExtractionWidthPPM,
                    true,
                    &ms2Ions
                    );
            const QPair<QVector<double>, QVector<double>> extractedIntensity = ParallelUtils::unZip(intensityFoundVsTheo);
            psmsReaderRow.intensityFound = extractedIntensity.first;
            psmsReaderRow.intensitySearched = extractedIntensity.second;
        }

        const QVector<QPair<PeptideStringWithMods, DiscScore>> &frameCandidateDiscScores
                = topCansInFrameIndexVsDiscScore.value(frameIndex);

        counter = 1;
        for (const QPair<PeptideStringWithMods, DiscScore> &discScorePr : frameCandidateDiscScores) {

            const PeptideStringWithMods &peptideStringWithMods = discScorePr.first;

            PSMsReaderRow &psmsReaderRow = psmReaderRowsForFrame[peptideStringWithMods];
            psmsReaderRow.discScore = discScorePr.second;
            psmsReaderRow.frameRankDiscScore = counter++;
        }

        psmsReaderRows->append(psmReaderRowsForFrame.values().toVector());
    }

    ERR_RETURN
}