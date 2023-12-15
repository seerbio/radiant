//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAFFWorkflow.h"

#include "ClassifierWeightsManager.h"
#include "EigenUtils.h"
#include "FastaFileToPeptidesListWorkFlow.h"
#include "FastaReader.h"
#include "FDRCLassifierNeuralNet.h"
#include "FragLibReader.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
#include "ProteinDigestomatic.h"

#include <QElapsedTimer>


PythiaDIAFFWorkflow::PythiaDIAFFWorkflow() : m_minTopNMs2Ions(6) {}

Err PythiaDIAFFWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &fastaUri
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::fileExists(fragLibUri); ree;

    if (!fastaUri.isEmpty()) {
        e = ErrorUtils::fileExists(fastaUri); ree;
        m_fastaUri = fastaUri;
    }

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;

    const double massMin
        = pythiaParameters.peptideLengthMin * Molecule(MolecularFormulas::alanineFormula).monoisotopicMass();

    const double massMax
            = pythiaParameters.peptideLengthMax * Molecule(MolecularFormulas::tryptophanFormula).monoisotopicMass();

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            m_fragLibUri,
            massMin,
            massMax,
            &fragLibReaderRows
            ); ree;

    e = m_targetDecoyCandidatePairManager.init(
            m_pythiaParameters,
            &fragLibReaderRows
            ); ree;

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath); ree;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    e = msReaderPointerAcc.ptr->collateMS2MzTargetFrames(
            &diaTargetFrames
            ); ree;

    const int msLevel = 1;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
    e = msReaderPointerAcc.ptr->getScanPoints(msLevel, &scanNumberVsScanTimeMS1); ree;

    e = m_targetDecoyCandidatePairScoretron.init(
            m_pythiaParameters,
            scanNumberVsScanTimeMS1,
            &msReaderPointerAcc,
            &diaTargetFrames
            ); ree;

    if (!m_pythiaParameters.turboMode) {

        e = buildCalibration(
                &msReaderPointerAcc,
                &diaTargetFrames,
                &scanNumberVsScanTimeMS1
        ); ree;

        e = optimizeParameters(&msReaderPointerAcc); ree;
    }

    e = mainAnalysis(&msReaderPointerAcc); ree;

    QVector<CandidateScores*> candidateScoreClassifierPntrs;
    e = applyNeuralNetClassifier(
            msReaderPointerAcc.ptr->scanTimeMinMax(),
            &candidateScoreClassifierPntrs
            ); ree;

    qDebug() << "Updating" << candidateScoreClassifierPntrs.size() << "PSMs";
    e = updateProteinGroupAnnotation(
            m_fastaUri,
            &candidateScoreClassifierPntrs
            ); ree;

    QVector<CandidateScores> candidateScoreClassifier;
    std::transform(
            candidateScoreClassifierPntrs.begin(),
            candidateScoreClassifierPntrs.end(),
            std::back_inserter(candidateScoreClassifier),
            [](const CandidateScores *cs){return *cs;}
            );

    const QString resultsFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
    e = ParquetReader::write(candidateScoreClassifier, resultsFilePath); ree;

//    limit what scores are collected depending on extendedScores or NeuralNetScores is set

    ERR_RETURN
}

namespace {

    void filterDuplicateCandidateScoresByDiscriminantScore(QVector<CandidateScores*> *candidateScores) {

        QMap<QString, CandidateScores*> keyVsCandidatesFoundBest;
        for (CandidateScores *cs : *candidateScores) {

            const QString key = cs->peptideStringWithMods + QString::number(cs->charge);
            if (keyVsCandidatesFoundBest.contains(key)) {
                continue;
            }

            keyVsCandidatesFoundBest.insert(key, cs);
        }

        *candidateScores = keyVsCandidatesFoundBest.values().toVector();
    }

    Err buildMsCalibrationReaderRows(
            const QVector<CandidateScores*> &_candidateScores,
            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(_candidateScores); ree;

        qDebug() << _candidateScores.size() << "Found for recalibartion";

        QVector<CandidateScores*> candidateScoresFiltered = _candidateScores;
        filterDuplicateCandidateScoresByDiscriminantScore(&candidateScoresFiltered);
        qDebug() << candidateScoresFiltered.size() << "Found for recalibartion filtered";

        const auto msCalibrationReaderRowsInsertLogic = [](CandidateScores *cs){

            MsCalibarationReaderRow row;
            row.peptideStringWithMods = cs->peptideStringWithMods;
            row.iRTPredicted = static_cast<float>(cs->iRTPredicted);
            row.scanTime = cs->scanTime;
            row.mzSearchedVec = cs->mzSearchedVec;
            row.mzFoundMeanVec = cs->mzFoundMeanVec;
            row.mzFoundStDevVec = cs->mzFoundStDevVec;

            return row;
        };

        std::transform(
                candidateScoresFiltered.begin(),
                candidateScoresFiltered.end(),
                std::back_inserter(*msCalibrationReaderRows),
                msCalibrationReaderRowsInsertLogic
        );

//#define WRITE_CALIBRATION_ROWS_TS
#ifdef WRITE_CALIBRATION_ROWS_TS
        qDebug() << "ACHTUNG, ACHTUNG, ACHTUNG!!!! WRITING CAL FILE IN:"; einfo;
        const QString filename = QStringLiteral("/home/anichols/temp/cal2.prq");
        e = ParquetReader::write(*msCalibrationReaderRows, filename); ree;
#endif

        ERR_RETURN
    }


}//namespace
Err PythiaDIAFFWorkflow::buildCalibration(
        MsReaderPointerAcc *msReaderPointerAcc,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    const double calibrationTrainingFraction = 1.0;
    const bool useExtendedScores = false;
    const bool useNeuralNetworkScores = false;
    const int minTrainingCountTranche = 50;

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos(),
            calibrationTrainingFraction,
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

    const int targetCount = std::accumulate(
            mzTargetKeyVsTargetDecoyCandidatePointers.begin(),
            mzTargetKeyVsTargetDecoyCandidatePointers.end(),
            0,
            [](int sum, const QVector<TargetDecoyCandidatePair*> &tv){return sum + tv.size();}
            );

    const double sizePerTranche = 10000.0;
    const int trancheSize = std::max(static_cast<int>(targetCount / sizePerTranche), 1);
    qDebug() << "Tranch size for calibration:" << trancheSize << "target count:" << targetCount << "sizePerTranche:" << static_cast<int>(sizePerTranche);

    QVector<QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>>> mzTargetKeyVsTargetDecoyCandidatePointersTranched;
    e = ParallelUtils::trancheMapValueVectorsByKeyForParallelization(
            mzTargetKeyVsTargetDecoyCandidatePointers,
            trancheSize,
            &mzTargetKeyVsTargetDecoyCandidatePointersTranched
            ); ree;

    const int topNMS2IonsCalibration = 6;

    for (int i = 0; i < mzTargetKeyVsTargetDecoyCandidatePointers.size(); i++) {

        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> localTranches;

        for (int j = 0; j <= i; j++) {
            const QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> &tranche = mzTargetKeyVsTargetDecoyCandidatePointersTranched.at(j);
            for (auto it = tranche.begin(); it != tranche.end(); it++) {
                const MzTargetKey &mzTargetKey = it.key();
                const QVector<TargetDecoyCandidatePair*> &tdcps = it.value();
                localTranches[mzTargetKey].append(tdcps);
            }
        }

        m_candidateScores.clear();

        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2IonsCalibration,
                m_msCalibratomatic,
                &localTranches,
                &m_candidateScores
                ); ree;

        const QPair<double, double> scanTimeMinMax = msReaderPointerAcc->ptr->scanTimeMinMax();
        e = setDiscriminantScoreForCandidates(
                scanTimeMinMax,
                useExtendedScores,
                useNeuralNetworkScores,
                topNMS2IonsCalibration
        ); ree;

        QVector<CandidateScores*> candidateScoresPntrs;
        e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;
        e = setQValueForCandidates(QValueScoreType::DiscriminantScore, &candidateScoresPntrs); ree;

        QMap<QString, int> fdrVsCount;
        e = FDRCLassifierNeuralNet::outputFDRResults(
                m_candidateScores,
                true,
                &fdrVsCount
        ); ree;

        const double fdrThreshold = 0.1;
        int targetCountBelowFDRThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                m_candidateScores,
                fdrThreshold,
                &targetCountBelowFDRThreshold
        ); ree;

        std::sort(candidateScoresPntrs.rbegin(), candidateScoresPntrs.rend(), [](CandidateScores *l, CandidateScores *r){
            return l->discriminateScore < r->discriminateScore;
        });

        const int idealTrainingCountAtGivenFDR = 1000;
        int minTrainingCount = std::max(minTrainingCountTranche, targetCountBelowFDRThreshold);
        minTrainingCount = std::min(minTrainingCount, idealTrainingCountAtGivenFDR);
        qDebug() << "Training RT count 10% FDR:" << minTrainingCount;

        candidateScoresPntrs.resize(minTrainingCount);
        QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
        e = buildMsCalibrationReaderRows(candidateScoresPntrs, &msCalibrationReaderRows); ree;

        const int numberOfStDevs = 3;

        e = m_msCalibratomatic.initRtOnly(msCalibrationReaderRows); ree;
        qDebug() << "scanTimeWindowStDev x 3:" << m_msCalibratomatic.scanTimeStDev(numberOfStDevs);

        const QString onePercenFDRKey = "1";
        if (fdrVsCount.value(onePercenFDRKey) >= 200) {

            msCalibrationReaderRows.resize(fdrVsCount.value(onePercenFDRKey));
            e = m_msCalibratomatic.initMzOnly(msCalibrationReaderRows); ree;

            e = recalibrateMzVals(
                    diaTargetFrames,
                    scanNumberVsScanTimeMS1,
                    &m_targetDecoyCandidatePairScoretron,
                    msReaderPointerAcc
            ); ree;
        }

        if (minTrainingCount >= idealTrainingCountAtGivenFDR) {
            qDebug() << "scanTimeWindowStDev x 3:" << m_msCalibratomatic.scanTimeStDev(numberOfStDevs);
            break;
        }
    }

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
        const QVector<MsScanInfo> &msScanInfos,
        double selectionFraction,
        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> *mzTargetKeyVsTargetDecoyCandidatePointers
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(msScanInfos); ree;
    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    mzTargetKeyVsTargetDecoyCandidatePointers->clear();

    for (const MsScanInfo &msScanInfo : msScanInfos) {

        if (msScanInfo.msLevel < 2) {
            continue;
        }

        QVector<TargetDecoyCandidatePair*> targetDecoyPointers;
        e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
                msScanInfo.precursorTargetMz - (msScanInfo.isoWindowLower + m_pythiaParameters.precursorExtractionWindowThomsons),
                msScanInfo.precursorTargetMz + (msScanInfo.isoWindowLower + m_pythiaParameters.precursorExtractionWindowThomsons),
                selectionFraction,
                &targetDecoyPointers
                ); ree;

        mzTargetKeyVsTargetDecoyCandidatePointers->insert(msScanInfo.targetKey(), targetDecoyPointers);
        if (m_pythiaParameters.verbosity > 1) {
            qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyPointers.size() << "targetDecoyPairs found";
        }

    }

    ERR_RETURN
}

namespace {

    struct TargetDecoyCandidateScores {
        CandidateScores *candidateScoresTarget;
        CandidateScores *candidateScoresDecoy;
        ScoresTargets *scoresTarget;
        ScoresDecoys  *scoresDecoy;
    };

    QString buildCandidateKey(const CandidateScores &candidateScores) {
        const QString decoyToString = candidateScores.isDecoy ? "_1" : "_0";
        return candidateScores.peptideStringWithMods + QString::number(candidateScores.charge) + candidateScores.targetKey + decoyToString;
    }

    struct BuildClassiferParallelInput {
        QVector<ScoresTargets*> scoresTargets;
        QVector<ScoresDecoys*> scoresDecoys;
    };

    struct BuildClassifierParallelOutput {
        QVector<QVector<double>> A;
        QVector<double> b;
    };

    Err buildParallelInput(
            const QVector<ScoresTargets*> &scoresTargets,
            const QVector<ScoresDecoys*> &scoresDecoys,
            QVector<BuildClassiferParallelInput> *inputs
            ) {
        ERR_INIT

        QVector<QVector<ScoresTargets*>> scoresTargetsTranched;
        QVector<QVector<ScoresDecoys*>> scoresDecoysTranched;

        e = ParallelUtils::trancheVectorForParallelization(
                scoresTargets,
                ParallelUtils::numberOfAvailableSystemProcessors(),
                &scoresTargetsTranched
                ); ree;

        e = ParallelUtils::trancheVectorForParallelization(
                scoresDecoys,
                ParallelUtils::numberOfAvailableSystemProcessors(),
                &scoresDecoysTranched
                ); ree;

        e = ErrorUtils::isEqual(scoresTargetsTranched.size(), scoresDecoysTranched.size()); ree;
        for (int i = 0; i < scoresDecoysTranched.size(); i++) {
            BuildClassiferParallelInput input;
            input.scoresTargets = scoresTargetsTranched.at(i);
            input.scoresDecoys = scoresDecoysTranched.at(i);
            inputs->push_back(input);
        }

        ERR_RETURN
    }

    QPair<Err, BuildClassifierParallelOutput> buildClassifierDataParallel(const BuildClassiferParallelInput &input) {

        ERR_INIT

        BuildClassifierParallelOutput output;

        e = ClassifierWeightsManager::buildDataClassifier1(
                input.scoresTargets,
                input.scoresDecoys,
                &output.A,
                &output.b
        ); rree;

        return {e, output};
    }

    Err processParallelResults(
            const QPair<Err, BuildClassifierParallelOutput> &result,
            int inputsSize,
            QVector<QVector<double>> *A,
            QVector<double> *b
            ) {

        ERR_INIT
        e = result.first; ree;

        const QVector<QVector<double>> &ALocal = result.second.A;
        const QVector<double> bLocal = result.second.b;

        if (A->isEmpty() || b->isEmpty()) {
            A->resize(ALocal.size());
            for(int row = 0; row < A->size(); row++) {
                QVector<double> v(ALocal.front().size(), 0.0);
                (*A)[row] = v;
            }
            b->resize(bLocal.size());
        }

        for (int row = 0; row < A->size(); row++) {
            for (int col = 0; col < A->front().size(); col++) {
                (*A)[row][col] += ALocal[row][col] / inputsSize;
            }
        }

        for (int row = 0; row < bLocal.size(); row++) {
            (*b)[row] += bLocal[row] / inputsSize;
        }

        ERR_RETURN
    }

    QPair<Err, QVector<double>> scoreVectorParallel(
            const CandidateScores &candidateScores,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            int theoMzIonsSize,
            const QPair<double, double> &scanTimeMinMax
            ) {

        ERR_INIT

        QVector<double> vec;
        e = CandidateScores::buildScoreVector(
                candidateScores,
                useExtendedScores,
                useNeuralNetworkScores,
                theoMzIonsSize,
                scanTimeMinMax,
                &vec
                ); rree;

        return {e, vec};
    }


}//namespace
Err PythiaDIAFFWorkflow::setDiscriminantScoreForCandidates(
        const QPair<double, double> &scanTimeMinMax,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int theoMzIonsSize
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(m_candidateScores.isEmpty()); ree;

    QElapsedTimer et;
    et.start();

    const auto scoreBuildBinder = std::bind(
            scoreVectorParallel,
            std::placeholders::_1,
            useExtendedScores,
            useNeuralNetworkScores,
            theoMzIonsSize,
            scanTimeMinMax
    );

    //TODO clean this up. Abstract sub processes into methods.

    QFuture<QPair<Err, QVector<double>>> futuresScoreBuilder = QtConcurrent::mapped(
            m_candidateScores,
            scoreBuildBinder
            );
    futuresScoreBuilder.waitForFinished();

    QVector<QVector<double>> scoreVectors;
    for (const QPair<Err, QVector<double>> &res : futuresScoreBuilder) {
        e = res.first; ree;
        scoreVectors.push_back(res.second);
    }

    e = ErrorUtils::isEqual(m_candidateScores.size(), scoreVectors.size()); ree;

    QMap<QString, TargetDecoyCandidateScores> keyVsTargetDecoyCandidateScores;

    for(int i = 0; i < scoreVectors.size(); i++) {

        CandidateScores *cs = &m_candidateScores[i];

        const QString key = cs->peptideStringWithMods + QString::number(cs->charge) + cs->targetKey;
        QVector<double> *sv = &scoreVectors[i];

        if (cs->isDecoy) {
            keyVsTargetDecoyCandidateScores[key].scoresDecoy = sv;
            keyVsTargetDecoyCandidateScores[key].candidateScoresDecoy = cs;
            continue;
        }

        keyVsTargetDecoyCandidateScores[key].candidateScoresTarget = cs;
        keyVsTargetDecoyCandidateScores[key].scoresTarget = sv;
    }
    qDebug() << "build key vs scores" << et.restart() << "mSec";

    QVector<ScoresTargets*> scoresTargets;
    QVector<ScoresDecoys*> scoresDecoys;
    QVector<CandidateScores*> candidateScoresTargetsPtrs;
    QVector<CandidateScores*> candidateScoresDecoysPtrs;

    for (auto it = keyVsTargetDecoyCandidateScores.begin(); it != keyVsTargetDecoyCandidateScores.end(); it++) {

        const QString &key = it.key();
        const TargetDecoyCandidateScores &tdcs = it.value();

        e = ErrorUtils::isEqual(tdcs.scoresTarget->size(), tdcs.scoresDecoy->size());
        if (e != eNoError) {
            qDebug() << "target decoys not paired for key" << key;
            rrr(eValueError);
        }

        scoresTargets.push_back(tdcs.scoresTarget);
        scoresDecoys.push_back(tdcs.scoresDecoy);
        candidateScoresTargetsPtrs.push_back(tdcs.candidateScoresTarget);
        candidateScoresDecoysPtrs.push_back(tdcs.candidateScoresDecoy);
    }
    qDebug() << "build separate" << et.restart() << "mSec";

    QVector<BuildClassiferParallelInput> inputs;
    e = buildParallelInput(
            scoresTargets,
            scoresDecoys,
            &inputs
    ); ree;

    QFuture<QPair<Err, BuildClassifierParallelOutput>> futures = QtConcurrent::mapped(
            inputs,
            buildClassifierDataParallel
            );
    futures.waitForFinished();

    QVector<QVector<double>> A;
    QVector<double> b;
    for (const QPair<Err, BuildClassifierParallelOutput> &result : futures) {
        e = processParallelResults(
                result,
                inputs.size(),
                &A,
                &b
                ); ree;
    }
    qDebug() << "build classifier" << et.restart() << "mSec";

    QVector<double> weights;
    e = ClassifierWeightsManager::fitWeights(A, b, &weights); ree;
    qDebug() << "fit weights" << et.restart() << "mSec";

    qDebug() << "Weights:" << weights;
    qDebug() << "b:" << b;

    QVector<double> discScoreTargets;
    e = ClassifierWeightsManager::applyWeights(scoresTargets, weights, &discScoreTargets); ree;

    QVector<double> discScoreDecoys;
    e = ClassifierWeightsManager::applyWeights(scoresDecoys, weights, &discScoreDecoys); ree;

    qDebug() << "apply weights scores" << et.restart() << "mSec";

    e = ErrorUtils::isEqual(scoresTargets.size(), scoresDecoys.size()); ree;
    e = ErrorUtils::isEqual(scoresTargets.size(), candidateScoresTargetsPtrs.size()); ree;
    e = ErrorUtils::isEqual(scoresDecoys.size(), candidateScoresDecoysPtrs.size()); ree;
    e = ErrorUtils::isEqual(discScoreTargets.size(), scoresTargets.size());
    e = ErrorUtils::isEqual(discScoreDecoys.size(), scoresDecoys.size());

    for (int i = 0; i < scoresDecoys.size(); i++) {

        CandidateScores* csTarget = candidateScoresTargetsPtrs.at(i);
        csTarget->discriminateScore = discScoreTargets.at(i);

        CandidateScores* csDecoy = candidateScoresDecoysPtrs.at(i);
        csDecoy->discriminateScore = discScoreDecoys.at(i);

    }
    qDebug() << "Setting scores" << et.restart() << "mSec";

    ERR_RETURN

}

Err PythiaDIAFFWorkflow::recalibrateMzVals(
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1,
        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
        MsReaderPointerAcc *msReaderPointerAcc
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msCalibratomatic.isInit()); ree;
    e = ErrorUtils::isFalse(diaTargetFrames->isEmpty()); ree;
    e = ErrorUtils::isFalse(scanNumberVsScanTimeMS1->isEmpty()); ree;

    qDebug() << "Recalibrating mz vals";

    QElapsedTimer et;
    et.start();

    const QList<MzTargetKey> &diaTargetFrameKeys = diaTargetFrames->keys();

    for (const MzTargetKey &k: diaTargetFrameKeys) {
        e = m_msCalibratomatic.recalibrateScanPoints(diaTargetFrames->value(k)); ree;
    }

    qDebug() << "Recalibrating MS1 mz vals frame";
    e = m_msCalibratomatic.recalibrateScanPoints(scanNumberVsScanTimeMS1); ree;

    e = targetDecoyCandidatePairScoretron->init(
            m_pythiaParameters,
            *scanNumberVsScanTimeMS1,
            msReaderPointerAcc,
            diaTargetFrames
            ); ree;

    qDebug() << "Points recalibrated in" << et.elapsed() << "mSec";

    ERR_RETURN
}

namespace {

    Err buildsetQValueForCandidateScoresInputs(
            const QValueScoreType &qValueScoreType,
            QVector<CandidateScores*> *candidateScores,
            QHash<PeptideSequenceChargeKey, double> *identifierVsTargets,
            QHash<PeptideSequenceChargeKey, double> *identifierVsDecoys

    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
        identifierVsTargets->clear();
        identifierVsDecoys->clear();

        for (const CandidateScores *cs : *candidateScores) {


            const PeptideSequenceChargeKey peptideSequenceChargeKey = buildCandidateKey(*cs);

            double classifierScore = -1.0;
            if (qValueScoreType == QValueScoreType::NNClassifierScore) {
                classifierScore = cs->classifierScore;
            }
            else {
                classifierScore = cs->discriminateScore;
            }


            if (cs->isDecoy) {
                identifierVsDecoys->insert(peptideSequenceChargeKey, classifierScore);
                continue;
            }

            identifierVsTargets->insert(peptideSequenceChargeKey, classifierScore);

        }

        ERR_RETURN
    }

    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            const QHash<PeptideSequenceChargeKey, double> &identifierVsQValue,
            const QHash<PeptideSequenceChargeKey, double> &identifierVsDecoyRatio,
            QVector<CandidateScores*> *candidateScores
    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

        for (CandidateScores *cs : *candidateScores) {

            if (cs->isDecoy) {
                continue;
            }

            const PeptideSequenceChargeKey peptideSequenceChargeKey = buildCandidateKey(*cs);

            e = ErrorUtils::isTrue(identifierVsQValue.contains(peptideSequenceChargeKey)); ree;
            e = ErrorUtils::isTrue(identifierVsDecoyRatio.contains(peptideSequenceChargeKey)); ree;

            cs->qValue = identifierVsQValue.value(peptideSequenceChargeKey);
            cs->decoyRatio = identifierVsDecoyRatio.value(peptideSequenceChargeKey);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::setQValueForCandidates(
        const QValueScoreType &qValueScoreType,
        QVector<CandidateScores*> *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    QElapsedTimer et;
    et.start();

    QHash<PeptideSequenceChargeKey, double> identifierVsTargets;
    QHash<PeptideSequenceChargeKey, double> identifierVsDecoys;
    e = buildsetQValueForCandidateScoresInputs(
            qValueScoreType,
            candidateScores,
            &identifierVsTargets,
            &identifierVsDecoys
    ); ree;

    QHash<PeptideSequenceChargeKey, double> identifierVsQValue;
    QHash<PeptideSequenceChargeKey, double> identifierVsDecoyRatio;
    e = MathUtils::calculateQValue(
            identifierVsTargets,
            identifierVsDecoys,
            &identifierVsQValue,
            &identifierVsDecoyRatio
    ); ree;

    e = setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            identifierVsQValue,
            identifierVsDecoyRatio,
            candidateScores
    ); ree;

    ERR_RETURN
}

namespace {

    struct DOEResult {
        double ppm = -1.0;
        double scanTimeStDev = -1.0;
//        double cosineSimAnchor = -1.0;
        int fdrCount = -1;
    };

    Err buildDOE(
            const PythiaParameters &pythiaParameters,
            double mzPPMStDev,
            double scanTimeStDev,
            QVector<PythiaParameters> *pythiaParametersExperiments
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(mzPPMStDev > 0.0); ree;
        e = ErrorUtils::isTrue(scanTimeStDev > 0.0); ree;
        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

//        const QVector<QVector<double>> experiments = {
//                {1.5, 1.0,  2.0},
//                {3.5, 1.0,  2.0},
//                {1.5,  3.0,  2.0},
//                {3.5,  3.0,  2.0},
//                {1.5,  2.0, 1.0},
//                {3.5,  2.0, 1.0},
//                {1.5,  2.0,  3.0},
//                {3.5,  2.0,  3.0},
//                {2.5, 1.0, 1.0},
//                {2.5,  3.0, 1.0},
//                {2.5, 1.0,  3.0},
//                {2.5,  3.0,  3.0},
//                {2.5,  2.0,  2.0}
//        };

//        const QVector<QVector<double>> experiments = {
//                {2.5, 1.0},
//                {3.5, 1.0},
//                {2.5,  3.0},
//                {3.5,  3.0},
//                {2.5,  2.0},
//                {3.5,  2.0}
//        };

        const QVector<QVector<double>> experiments = {
                {2.0, 2.0},
                {2.5,  2.0},
                {3.5,  2.0},
                {4.5, 2.0},
                {5.5, 2.0}
        };

        for (const QVector<double> &exp : experiments) {

            PythiaParameters params = pythiaParameters;
            params.ms2ExtractionWidthPPM = mzPPMStDev * exp.at(0);
            params.scanTimeWindowMinutes = scanTimeStDev * exp.at(1);

//            switch (static_cast<int>(exp.at(2))) {
//                case 1:
//                    params.cosineSimToAnchorThreshold = 0.9;
//                    break;
//                case 2:
//                    params.cosineSimToAnchorThreshold = 0.935;
//                    break;
//                case 3:
//                    params.cosineSimToAnchorThreshold = 0.97;
//                    break;
//                default:
//                    params.cosineSimToAnchorThreshold = pythiaParameters.cosineSimToAnchorThreshold;
//            }

            pythiaParametersExperiments->push_back(params);
        }

        for (const PythiaParameters &pp : *pythiaParametersExperiments) {
            qDebug() << "ppmTol" << pp.ms2ExtractionWidthPPM
                     << "scanTimeWinMin" << pp.scanTimeWindowMinutes
                     << "cosineSimSumAnchThreshold" << pp.cosineSimToAnchorThreshold;
        }

        ERR_RETURN
    }


    Err findMostFrequentValue(
            const QVector<QPair<QString, int>> &countsVector,
            double *value
            ) {
        ERR_INIT

        e = ErrorUtils::isNotEmpty(countsVector); ree;

        if (countsVector.size() == 1) {
            e = ErrorUtils::toDouble(countsVector.front().first, value); ree;
            ERR_RETURN
        }

        QVector<QPair<QString, int>> cv = countsVector;

        std::sort(cv.rbegin(), cv.rend(), [](const QPair<QString, int> &l, const QPair<QString, int> &r){
            return l.second < r.second;
        });

        if (cv.at(0).second == cv.at(1).second) {
            double v1;
            e = ErrorUtils::toDouble(cv.at(0).first, &v1); ree;

            double v2;
            e = ErrorUtils::toDouble(cv.at(1).first, &v2); ree;

            * value = (v1 + v2) / 2;
            ERR_RETURN
        }

        e = ErrorUtils::toDouble(cv.front().first, value); ree;

        ERR_RETURN
    }

    Err getTopFrequencyParameters(
            QVector<DOEResult> *results,
            double *ppmSetting
//            double *scanTimeWidthSetting
//            double *cosineSimSetting
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(*results); ree;

        std::sort(
                results->rbegin(),
                results->rend(),
                [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
                );

        for (const DOEResult &r : *results) {
//            qDebug() << r.ppm << r.scanTimeStDev << r.cosineSimAnchor << r.fdrCount;
            qDebug() << r.ppm << r.scanTimeStDev << r.fdrCount;
        }

        *ppmSetting = results->front().ppm;
//        *scanTimeWidthSetting = results->front().scanTimeStDev;

//        const int topNResults = 3;
//        results->resize(topNResults);
//
//        QMap<QString, int> ppmCounts;
//        QMap<QString, int> scanTimeWidthCounts;
////        QMap<QString, int> cosineSimCounts;
//        for (const DOEResult &r : *results) {
//
//            const QString ppmString = QString::number(r.ppm);
//            const QString scanTimeWidthString = QString::number(r.scanTimeStDev);
////            const QString cosineSimString = QString::number(r.cosineSimAnchor);
//
//            ppmCounts[ppmString]++;
//            scanTimeWidthCounts[scanTimeWidthString]++;
////            cosineSimCounts[cosineSimString]++;
//        }
//
//        const QVector<QPair<QString, int>> ppmCountsVector = ParallelUtils::convertMapToVectorPairs(ppmCounts);
//        const QVector<QPair<QString, int>> scanTimeWidthCountsVector = ParallelUtils::convertMapToVectorPairs(scanTimeWidthCounts);
////        const QVector<QPair<QString, int>> cosineSimCountsVector = ParallelUtils::convertMapToVectorPairs(cosineSimCounts);
//
//        e = findMostFrequentValue(ppmCountsVector, ppmSetting); ree;
//        e = findMostFrequentValue(scanTimeWidthCountsVector, scanTimeWidthSetting); ree;
////        e = findMostFrequentValue(cosineSimCountsVector, cosineSimSetting); ree;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::optimizeParameters(MsReaderPointerAcc *msReaderPointerAcc) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

    const int topNMS2IonsOptimization = 6;
    qDebug() << "Using top:" << topNMS2IonsOptimization << "fragments for optimization";

    const double desiredTestSize = std::min(20000.0, static_cast<double>(m_targetDecoyCandidatePairManager.targetsCount()));
    const double selectionFraction = desiredTestSize / m_targetDecoyCandidatePairManager.targetsCount();
    qDebug() << "Optimization selection fraction" << selectionFraction;

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos(),
            selectionFraction,
            &mzTargetKeyVsTargetDecoyCandidatePointers
    ); ree;

    QVector<PythiaParameters> pythiaParametersExperiments;
    e = buildDOE(
            m_pythiaParameters,
            m_msCalibratomatic.mzStDev(),
            m_msCalibratomatic.scanTimeStDev(),
            &pythiaParametersExperiments
            ); ree;

    const bool useExtendedScores = false;
    const bool useNeuralNetworkScores = false;

    QVector<DOEResult> results;
    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {

        qDebug() << "ppmTol" << pythiaParams.ms2ExtractionWidthPPM
                 << "scanTimeWinMin" << pythiaParams.scanTimeWindowMinutes
                 << "cosineSimSumAnchThreshold" << pythiaParams.cosineSimToAnchorThreshold;

        e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(pythiaParams); ree;
        qDebug() << "STarting opt";

        m_candidateScores.clear();
        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2IonsOptimization,
                m_msCalibratomatic,
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                &m_candidateScores
        ); ree;

        const QPair<double, double> scanTimeMinMax = msReaderPointerAcc->ptr->scanTimeMinMax();
        e = setDiscriminantScoreForCandidates(
                scanTimeMinMax,
                useExtendedScores,
                useNeuralNetworkScores,
                topNMS2IonsOptimization
        ); ree;

        QVector<CandidateScores*> candidateScoresPntrs;
        e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;
        e = setQValueForCandidates(QValueScoreType::DiscriminantScore, &candidateScoresPntrs); ree;

        QMap<QString, int> fdrVsCount;
        e = FDRCLassifierNeuralNet::outputFDRResults(
                m_candidateScores,
                true,
                &fdrVsCount
        ); ree;
        qDebug() << "Ending opt";

        const double fdrThreshold = 0.1;
        int targetCountAboveFDRQValueThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                m_candidateScores,
                fdrThreshold,
                &targetCountAboveFDRQValueThreshold
                ); ree;

        DOEResult res;
        res.ppm = pythiaParams.ms2ExtractionWidthPPM;
//        res.scanTimeStDev = pythiaParams.scanTimeWindowMinutes;
//        res.cosineSimAnchor = pythiaParams.cosineSimToAnchorThreshold;
        res.fdrCount = targetCountAboveFDRQValueThreshold;
        results.push_back(res);
    }

    e = getTopFrequencyParameters(
            &results,
            &m_pythiaParameters.ms2ExtractionWidthPPM
//            &m_pythiaParameters.scanTimeWindowMinutes
//            &m_pythiaParameters.cosineSimToAnchorThreshold
            ); ree;

    e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(m_pythiaParameters); ree;

    qDebug() << "Optimal ppm setting:" << m_pythiaParameters.ms2ExtractionWidthPPM;
    qDebug() << "Optimal scanTimeWindow setting:" << m_pythiaParameters.scanTimeWindowMinutes;
    qDebug() << "Optimal cosineSimSum setting:" << m_pythiaParameters.cosineSimToAnchorThreshold;

//    //TODO replace this with response surface derived DOE parameters when you figure out how to do it.
////    const DOEResult bestParametersFDR = *std::max_element(
////            results.rbegin(),
////            results.rend(),
////            [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
////    );
////
////    m_pythiaParameters.ms2ExtractionWidthPPM = bestParametersFDR.ppm;
////    m_pythiaParameters.scanTimeWindowMinutes = bestParametersFDR.scanTimeStDev;
////    m_pythiaParameters.cosineSimToAnchorThreshold = bestParametersFDR.cosineSimAnchor;
//    ///////////////////////////////////////////////////////////////////////////////////////////////////

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::mainAnalysis(MsReaderPointerAcc *msReaderPointerAcc) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

    m_candidateScores.clear();

    const bool useExtendedScores = !m_pythiaParameters.turboMode;
    const bool useNeuralNetworkScores = !m_pythiaParameters.turboMode;

    const int topNMs2IonsMainAnalysis = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions))
            );

    qDebug() << "Using top:" << topNMs2IonsMainAnalysis << "fragments for main analysis";

    const double selectionFractionBypass = -1.0;

    QElapsedTimer et;
    et.start();

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos(),
            selectionFractionBypass,
            &mzTargetKeyVsTargetDecoyCandidatePointers
    ); ree;
    qDebug() << "Targets fetched" << et.restart() << "mSec";

    e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
            topNMs2IonsMainAnalysis,
            m_msCalibratomatic,
            &mzTargetKeyVsTargetDecoyCandidatePointers,
            &m_candidateScores
    ); ree;
    qDebug() << "Targets scored" << et.restart() << "mSec";

    const QPair<double, double> scanTimeMinMax = msReaderPointerAcc->ptr->scanTimeMinMax();
    e = setDiscriminantScoreForCandidates(
            scanTimeMinMax,
            useExtendedScores,
            useNeuralNetworkScores,
            topNMs2IonsMainAnalysis
    ); ree;
    qDebug() << "Targets discriminated" << et.restart() << "mSec";

    QVector<CandidateScores*> candidateScoresPntrs;
    e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;
    e = setQValueForCandidates(QValueScoreType::DiscriminantScore, &candidateScoresPntrs); ree;
    qDebug() << "Targets qval'd" << et.restart() << "mSec";

    QMap<QString, int> fdrVsCount;
    e = FDRCLassifierNeuralNet::outputFDRResults(
            m_candidateScores,
            true,
            &fdrVsCount
    ); ree;
    qDebug() << "Targets resulted" << et.restart() << "mSec";

    const double fdrThreshold = 0.1;
    int targetCountBelowFDRThreshold;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
            m_candidateScores,
            fdrThreshold,
            &targetCountBelowFDRThreshold
    ); ree;
    qDebug() << "Targets counted" << et.restart() << "mSec";



    std::sort(candidateScoresPntrs.rbegin(), candidateScoresPntrs.rend(), [](CandidateScores *l, CandidateScores *r){
        return l->discriminateScore < r->discriminateScore;
    });
    qDebug() << "Targets sorted" << et.restart() << "mSec";

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::buildCandidateScoresPtrs(QVector<CandidateScores*> *candidateScoresPntrs) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_candidateScores); ree;

    std::transform(
            m_candidateScores.begin(),
            m_candidateScores.end(),
            std::back_inserter(*candidateScoresPntrs),
            [](CandidateScores &cs){return &cs;}
    );

    ERR_RETURN
}

namespace {

    Err filterScoredCandidatesTo50PercentFDR(
            QVector<CandidateScores*> *candidateScoresTargetsAndDecoys,
            QVector<CandidateScores*> *candidateScoresTargetsAndDecoys50PercentFDRFiltered
            ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoresTargetsAndDecoys->isEmpty()); ree;

        std::sort(
                candidateScoresTargetsAndDecoys->rbegin(),
                candidateScoresTargetsAndDecoys->rend(),
                [](const CandidateScores *l, const CandidateScores *r){return l->discriminateScore < r->discriminateScore;}
        );

        const double fdrThreshold = 0.5;
        for (CandidateScores *csp : *candidateScoresTargetsAndDecoys) {

            candidateScoresTargetsAndDecoys50PercentFDRFiltered->push_back(csp);

            if (csp->qValue >= fdrThreshold && !csp->isDecoy) {
                break;
            }
        }

        std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
        std::shuffle(
                candidateScoresTargetsAndDecoys50PercentFDRFiltered->begin(),
                candidateScoresTargetsAndDecoys50PercentFDRFiltered->end(),
                rng
                );

        ERR_RETURN
    }

    Err minMaxScaleScores(
            const QVector<KarnnNNTarget> &karnnNNTargets,
            QVector<KarnnNNTarget> *karnnNNTargetsNorm
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(karnnNNTargets); ree;

        QVector<QVector<double>> vecs;
        std::transform(
                karnnNNTargets.begin(),
                karnnNNTargets.end(),
                std::back_inserter(vecs),
                [](const KarnnNNTarget &kt){return kt.scoreVec;}
                );

        Eigen::MatrixX<double> mat = EigenUtils::convertQVectorsToEigenMatrix(vecs);
        EigenUtils::minMaxScaleMatrix(&mat);

        const QVector<QVector<double>> vecsNorm = EigenUtils::convertEigenMatrixToQVectors(mat);

        e = ErrorUtils::isEqual(vecsNorm.size(), karnnNNTargets.size()); ree;

        for (int i = 0; i < vecsNorm.size(); i++) {
            KarnnNNTarget ktNew = karnnNNTargets.at(i);
            ktNew.scoreVec = vecsNorm.at(i);
            karnnNNTargetsNorm->push_back(ktNew);
        }

        ERR_RETURN
    }

    Err buildKarnnNNTargetsNormalized(
            const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            const QPair<double, double> scanTimeMinMax,
            int topNMs2Ions,
            QVector<KarnnNNTarget> *karnnNNTargetsNorm
    ){

        ERR_INIT

        e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoys50PercentFDRFiltered); ree;
        e = ErrorUtils::isTrue(scanTimeMinMax.second > scanTimeMinMax.first); ree;

        QVector<KarnnNNTarget> karnnNNTargets;
        for (int i = 0; i < candidateScoresTargetsAndDecoys50PercentFDRFiltered.size(); i++) {
            const CandidateScores *cs = candidateScoresTargetsAndDecoys50PercentFDRFiltered.at(i);
            KarnnNNTarget karnnNnTarget;
            karnnNnTarget.seq = cs->peptideStringWithMods;
            karnnNnTarget.isDecoy = cs->isDecoy;
            karnnNnTarget.index = i;

            e = CandidateScores::buildScoreVector(
                    *cs,
                    true,
                    true,
                    topNMs2Ions,
                    scanTimeMinMax,
                    &karnnNnTarget.scoreVec
            ); ree;
            karnnNNTargets.push_back(karnnNnTarget);
        }

        e = minMaxScaleScores(karnnNNTargets, karnnNNTargetsNorm); ree;

        ERR_RETURN
    }

    Err predictNNScores(
            const QVector<KarnnNNTarget> &karnnNNTargetsNorm,
            QVector<float> *predictions
            ) {

        ERR_INIT

        const int batchSize = std::min(50, std::max(1, static_cast<int>(karnnNNTargetsNorm.size() / 100.0)));
        const int maxIters = 1;
        qDebug() << "Batch Size:" << batchSize << "MaxIters:" << maxIters;

        QVector<QVector<float>> xData;
        QVector<float> yData;
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            xData.push_back(QVector<float>(kt.scoreVec.begin(), kt.scoreVec.end()));
            yData.push_back(kt.isDecoy ? 1.0 : 0.0);
        }

        const int baggingSize = 12;
        const float learningRate = 0.003;
        const int epochs = 1;
        FDRCLassifierNeuralNet fdrcLassifierNeuralNet;
        e = fdrcLassifierNeuralNet.init(
                epochs,
                baggingSize,
                batchSize,
                learningRate
        ); ree;

        e = fdrcLassifierNeuralNet.exec(xData, yData, predictions); ree;

        ERR_RETURN
    }

    Err processPredictions(
            const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            const QVector<float> &predictions,
            bool reportDecoys,
            QVector<KarnnNNTarget> *karnnNNTargetsNorm,
            QVector<CandidateScores*> *candidateScoreClassifier
            ) {

        ERR_INIT

        for (int i = 0; i < predictions.size(); i++) {
            (*karnnNNTargetsNorm)[i].nnScore = predictions.at(i);
        }

        std::sort(
                karnnNNTargetsNorm->begin(),
                karnnNNTargetsNorm->end(),
                [](const KarnnNNTarget &l, const KarnnNNTarget &r){return l.nnScore < r.nnScore;}
        );

        int counter = 0;
        int falsePositives = 0;
        const double fdrCutoff = 0.008;
        for (const KarnnNNTarget &rp : *karnnNNTargetsNorm) {

            CandidateScores *candidateScoresNew = candidateScoresTargetsAndDecoys50PercentFDRFiltered.at(rp.index);
            candidateScoresNew->classifierScore = rp.nnScore;
            candidateScoreClassifier->push_back(candidateScoresNew);

            ++counter;

            if (rp.nnScore > 0.5 || (falsePositives / static_cast<double>(counter)) > fdrCutoff) {
                if (!reportDecoys) {
                    break;
                }
                else {
                    counter--;
                    continue;
                }
            }

            if (rp.isDecoy){
                falsePositives++;
            }
        }

        qDebug() << "False Pos" << falsePositives << "Total" << counter << "FDR 0.5 nnScore cuttoff" << falsePositives / (counter + 0.0);

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::applyNeuralNetClassifier(
        const QPair<double, double> &scanTimeMinMax,
        QVector<CandidateScores*> *candidateScoreClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_candidateScores); ree;

    QVector<CandidateScores*> candidateScoresTargetsAndDecoys;
    e = buildCandidateScoresPtrs(&candidateScoresTargetsAndDecoys); ree;

    QVector<CandidateScores*> candidateScoresTargetsAndDecoys50PercentFDRFiltered;
    e = filterScoredCandidatesTo50PercentFDR(
            &candidateScoresTargetsAndDecoys,
            &candidateScoresTargetsAndDecoys50PercentFDRFiltered
            ); ree;

//#define BYPASS_NEURAL_NET
#ifdef BYPASS_NEURAL_NET
    std::sort(
            candidateScoresTargetsAndDecoysShuffled.rbegin(),
            candidateScoresTargetsAndDecoysShuffled.rend(),
            [](const CandidateScores &l, const CandidateScores &r){return l.discriminateScore < r.discriminateScore;}
            );

    int counter = 0;
    int decoyCounter = 0;
    for (const CandidateScores &cs : candidateScoresTargetsAndDecoysShuffled) {

        if (cs.isDecoy) {
            decoyCounter++;
        }

        const double qVal = static_cast<double>(decoyCounter) / ++counter;
        if (qVal > 0.01) {
            break;
        }

        candidateScoreClassifier->push_back(cs);
    }
#else

    const int totalCount = candidateScoresTargetsAndDecoys50PercentFDRFiltered.size();
    const int decoyCount = static_cast<int>(std::count_if(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered.begin(),
            candidateScoresTargetsAndDecoys50PercentFDRFiltered.end(),
            [](const CandidateScores* cs){return cs->isDecoy;}
            ));

    qDebug() << "target vs decoy count" << totalCount - decoyCount << decoyCount
             << "total" << totalCount;

    QVector<KarnnNNTarget> karnnNNTargetsNorm;
    e = buildKarnnNNTargetsNormalized(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            scanTimeMinMax,
            m_pythiaParameters.topNMs2Ions,
            &karnnNNTargetsNorm
            ); ree;

    QVector<float> predictions;
    e = predictNNScores(karnnNNTargetsNorm, &predictions); ree;


    e = processPredictions(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            predictions,
            m_pythiaParameters.reportDecoys,
            &karnnNNTargetsNorm,
            candidateScoreClassifier
            ); ree;

    e = setQValueForCandidates(QValueScoreType::NNClassifierScore, candidateScoreClassifier); ree
#endif

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::updateProteinGroupAnnotation(
        const QString &fastaFilePath,
        QVector<CandidateScores*> *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    FastaReader fastaReader;
    e = fastaReader.parseFastaFile(fastaFilePath); ree;

    QVector<PeptideSequence> peptideSequences;
    QMap<PeptideStringWithMods, QVector<FastaEntry>> peptideStringWithModsVsFastaEntries;
    e = FastaFileToPeptidesListWorkFlow::digestFastaEntries(
            m_pythiaParameters,
            fastaReader.fastaEntries(),
            &peptideSequences,
            &peptideStringWithModsVsFastaEntries
            ); ree;

    QMap<PeptideStringWithMods, QVector<FastaEntry>> peptideStringWithModsVsFastaEntriesLeucinesReplaced;
    for (auto it = peptideStringWithModsVsFastaEntries.begin(); it != peptideStringWithModsVsFastaEntries.end(); it++) {
        QString peptideSeqReplacedLeucines = it.key();
        peptideSeqReplacedLeucines = peptideSeqReplacedLeucines.replace('L', 'J').replace('I', 'J');
        peptideStringWithModsVsFastaEntriesLeucinesReplaced.insert(peptideSeqReplacedLeucines, it.value());
    }

    for (int i = 0; i < candidateScores->size(); i++) {

        CandidateScores *cs = (*candidateScores)[i];

        QString peptideSeqReplacedLeucines = cs->peptideStringWithMods;
        peptideSeqReplacedLeucines = peptideSeqReplacedLeucines.replace('L', 'J').replace('I', 'J');

        const QVector<FastaEntry> &fastaEntries = peptideStringWithModsVsFastaEntriesLeucinesReplaced.value(peptideSeqReplacedLeucines);

        QStringList fastaDescriptions;
        std::transform(
                fastaEntries.begin(),
                fastaEntries.end(),
                std::back_inserter(fastaDescriptions),
                [](const FastaEntry &fe){return fe.fastaDescription;}
                );

        cs->proteinGroup = fastaDescriptions.join(';');

    }

    ERR_RETURN
}
