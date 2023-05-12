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
            theoIntensities(i) = ms2Ions.at(i).y();
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

        FragLibReader::filterMs2IonsByMz(
                params.mzMinDataStructure,
                params.mzMaxDataStructure,
                &ms2Ions
                );

        FragLibReader::getTopNMostIntenseMs2Ions(
                params.topNMs2Ions,
                &ms2Ions
                );

        const Eigen::MatrixX<double> theoMat = buildTheoMat(
                ms2Ions,
                frameScanCount,
                topNMs2Ions
        );

        for (int colIdx = 0; colIdx < ms2Ions.size(); colIdx++) {

            const double mz = ms2Ions.at(colIdx).x();
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

                if (scoringMatIntensity.coeff(frameIndex, colIdx) > intensity) {
                    continue;
                }

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

        const Eigen::VectorX<double> scorePerFrameIndexOfTarget = (
                foundIonsPerFrameIndexOfTargetFactorial.array()
                * fractionFoundIonsPerFrameIndexOfTarget.array()
                * mzProductPerFrameIndexOfTarget.array().sqrt()
                * cosineSimPerFrameIndexOfTarget.array()
                ).sqrt().sqrt();

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
Err MsFrameScoretron::buildFrameScoreVectors(QString *frameScoreVectorsFilePath) {

    ERR_INIT

//#define WRITE_FRAME
#ifdef WRITE_FRAME
//    for (int i = 0; i < 412; i++) {
//        qDebug() << i << m_msFrame.scanNumberFromFrameIndex(i);
//    }
    const QString framesFilePath = msDataFilePath + ".frameScans";
    e = m_msFrame.writeFrameScans(framesFilePath);
    qDebug() << e << "Drewholio" << framesFilePath;
#endif

    QMap<PeptideStringWithMods, FrameIndexScoreResultOfTarget> pepStrWModsVsFrameIndexScoreResultOfTargets;
    e = processFrameLogic(
            {m_msFrame, m_fragPreds},
            m_params,
            &pepStrWModsVsFrameIndexScoreResultOfTargets
    ); ree

    filterByFoundMzCount(
            m_params.minFoundMzPeaks,
            &pepStrWModsVsFrameIndexScoreResultOfTargets
            );

    *frameScoreVectorsFilePath = m_msDataFilePath + "." + m_msFrame.uniqueMsInfoScanKey() + ".frameScores";
    e = writeFrameTargetScoreVectors(
            pepStrWModsVsFrameIndexScoreResultOfTargets,
            *frameScoreVectorsFilePath
    ); ree;


    ERR_RETURN
}

Err MsFrameScoretron::buildFragIonLibForTargetMz(const QString &fragLibUri) {

    ERR_INIT

    const int monoOffset = 0;
    for (Charge charge = m_params.chargeStateMin; charge <= m_params.chargeStateMax; ++charge) {

        const double massStart = BiophysicalCalcs::calculateMassFromThomson(
                m_mzTargetStartStop.first,
                charge,
                monoOffset
        );

        const double massEnd = BiophysicalCalcs::calculateMassFromThomson(
                m_mzTargetStartStop.second,
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

            m_fragPreds.insert(peptideStringWithMods, ms2Ions);
            m_fragPredsIsDecoy.insert(peptideStringWithMods, peptideSequenceChargeKeyVsIsDecoy.value(peptideSequenceChargeKey));
        }

    }

    ERR_RETURN
}

Err MsFrameScoretron::buildMsFrame() {

    ERR_INIT

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(
            m_msDataFilePath,
            MsParquetReaderNamespace::PERCURSOR_TARGET_MZ,
            {m_mzTargetStartStop.first, m_mzTargetStartStop.second}
    ); ree;

    const QMap<ScanNumber, ScanPoints> targetScanPoints = msReaderParquet.getScanPoints();
    e = ErrorUtils::isNotEmpty(targetScanPoints); ree;

    e = m_msFrame.init(
            m_params,
            m_uniqueMsInfoScanKey,
            targetScanPoints,
            m_mzTargetStartStop
    ); ree;

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
        row.frameIndexMaxScore = MathUtils::findMaxIndexInVector(frameIndexScoreResultOfTarget.scorePerFrameIndexOfTargetVec);
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

Err MsFrameScoretron::init(
        const PythiaParameters &params,
        const QString &msDataFilePath,
        const QString &fragLibFilePath,
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        const QPair<double, double> &mzTargetStartStop
        ) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    e = ErrorUtils::fileExists(fragLibFilePath); ree;
    e = ErrorUtils::isTrue(params.isValid()); ree;
    e = ErrorUtils::isNotEmpty(uniqueMsInfoScanKey); ree;
    e = ErrorUtils::isAboveThreshold(
            mzTargetStartStop.second,
            mzTargetStartStop.first,
            ErrorUtilsParam::ExcludeThreshold
            ); ree;

    m_params = params;
    m_msDataFilePath = msDataFilePath;
    m_uniqueMsInfoScanKey = uniqueMsInfoScanKey;
    m_mzTargetStartStop = mzTargetStartStop;

    e = buildFragIonLibForTargetMz(fragLibFilePath); ree;
    e = ErrorUtils::isNotEmpty(m_fragPreds); ree;

    e = buildMsFrame(); ree;
    e = ErrorUtils::isAboveThreshold(
            m_msFrame.scanCount(),
            0,
            ErrorUtilsParam::ExcludeThreshold
            );ree;

    qDebug() << "TargetKey" << uniqueMsInfoScanKey
             << "Scan Count" << m_msFrame.scanCount()
             << "Candidate Count:" << m_fragPreds.size();

    ERR_RETURN
}
