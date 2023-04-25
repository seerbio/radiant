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
#include "PeakIntegratomatic.h"
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

    Err buildFragIonLibForTarget(
            const PythiaParameters &params,
            const QString &fragLibUri,
            const QPair<double, double> &mzTargetStartStop,
            QMap<PeptideStringWithMods, QVector<MS2Ion>> *output
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
            e = fragLibReader.getMS2Ions(
                    massStart,
                    massEnd,
                    &peptideSequenceChargeKeyVsMS2Ions
                    ); ree;

            for (auto it = peptideSequenceChargeKeyVsMS2Ions.begin(); it != peptideSequenceChargeKeyVsMS2Ions.end(); it++) {

                const PeptideSequenceChargeKey &peptideSequenceChargeKey = it.key();
                QVector<MS2Ion> &ms2Ions = it.value();

                PeptideStringWithMods peptideStringWithMods;
                e = peptideStringWithModsFromPeptideSequenceChargeKey(
                        peptideSequenceChargeKey,
                        &peptideStringWithMods
                        ); ree;

                filterMs2IonsByMz(
                        params.mzMinDataStructure,
                        params.mzMaxDataStructure,
                        &ms2Ions
                        );

                getTopNMostIntenseMs2Ions(
                        params.topNMs2Ions,
                        &ms2Ions
                        );

                output->insert(peptideStringWithMods, ms2Ions);
            }

        }

        ERR_RETURN
    }

    Err buildMsFrame(
            const QString &msDataFilePath,
            const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
            const PythiaParameters &params,
            const QPair<double, double> mzTargetStartStop,
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

        const bool denoise = true;
        const bool deisotope = true;
        const bool smooth = false;
        e = msFrame->preprocessMsFrame(
                denoise,
                deisotope,
                smooth
        ); ree;

        ERR_RETURN
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
            const QVector<MS2Ion> &ms2IonsTandemPred,
            int frameScanCount,
            int topNMs2Ions,
            double featureFinderTolerancePPM,
            TurboXIC *turboXic
    ) {

        Eigen::MatrixX<double> scoringMatIntensity(frameScanCount, topNMs2Ions);
        scoringMatIntensity.setZero();

        Eigen::MatrixX<double> scoringMatMz(frameScanCount, topNMs2Ions);
        scoringMatMz.setOnes();

        const Eigen::MatrixX<double> theoMat = buildTheoMat(
                ms2IonsTandemPred,
                frameScanCount,
                topNMs2Ions
        );

        for (int colIdx = 0; colIdx < ms2IonsTandemPred.size(); colIdx++) {

            const double mz = ms2IonsTandemPred.at(colIdx).mz;
            const double massTol = MathUtils::calculatePPM(
                    mz,
                    featureFinderTolerancePPM
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

    struct FrameIndexScoreResultOfTarget {
        QVector<double> cosineSimPerFrameIndexOfTargetVec;
        QVector<double> intensityPerFrameIndexOfTargetVec;
        QVector<int> foundIonsPerFrameIndexOfTargetVec;
        QVector<double> scorePerFrameIndexOfTargetVec;
        PeakIntegrationIndexes bestScorePeakLimits = {-1, -1};
        int charge = -1;
    };

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

    Err writeFrameTargetScoreVectors(
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

    Err scoreFrameTargets(
            const MsFrame &frame,
            const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPreds,
            const PythiaParameters &params,
            const QString &msDataFilePath,
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
                    ms2IonsTandemPred,
                    frame.scanCount(),
                    params.topNMs2Ions,
                    params.featureFinderTolerancePPM,
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

    Err processFrameLogic(
            const QPair<MsFrame, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &chunk,
            const PythiaParameters &params,
            const QString &msDataFilePath,
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
                msDataFilePath,
                pepStrWModsVsFrameIndexScoreResultOfTargets
        ); ree;

        ERR_RETURN
    }

    void filterByFoundMzCount(
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

}//namespace
QPair<Err, QPair<UniqueMsInfoScanKey, QString>> MsFrameScoretron::scoreCandidates(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
        ) {

    ERR_INIT

    QMap<PeptideStringWithMods, QVector<MS2Ion>> fragPreds;
    e = buildFragIonLibForTarget(
            params,
            fragLibFilePath,
            mzTargetStartStop,
            &fragPreds
            ); rree;
    e = ErrorUtils::isNotEmpty(fragPreds); rree;

    MsFrame msFrame;
    e = buildMsFrame(
            msDataFilePath,
            uniqueMsInfoScanKey,
            params,
            mzTargetStartStop,
            &msFrame
            ); rree;

    qDebug() << uniqueMsInfoScanKey << msFrame.scanCount() << fragPreds.size();

    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> pepStrWModsVsFrameIndexScoreResultOfTargets;
    e = processFrameLogic(
            {msFrame, fragPreds},
            params,
            msDataFilePath,
            &pepStrWModsVsFrameIndexScoreResultOfTargets
    ); rree

    filterByFoundMzCount(
            params.minFoundMzPeaks,
            &pepStrWModsVsFrameIndexScoreResultOfTargets
            );

    const QString outputFilePath = msDataFilePath + "." + msFrame.uniqueMsInfoScanKey() + ".frameScores";
//    e = writeFrameTargetScoreVectors(
//            pepStrWModsVsFrameIndexScoreResultOfTargets,
//            outputFilePath
//    ); rree;

    QMap<FrameIndex, QVector<QPair<PeptideStringWithMods, Score>>> topCansInFrameIndex;
    e = MsFrameScoretronProcessormatic::processLogicForFrameScores(
            outputFilePath,
            msFrame,
            params.returnPSMTopN,
            &topCansInFrameIndex
            ); rree;

    return {e, {uniqueMsInfoScanKey, outputFilePath}};
}
