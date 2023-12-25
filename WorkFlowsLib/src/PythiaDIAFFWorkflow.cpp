//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAFFWorkflow.h"

#include "CandidateScores.h"
#include "ClassifierWeightsManager.h"
#include "EigenUtils.h"
#include "FastaFileToPeptidesListWorkFlow.h"
#include "FastaReader.h"
#include "FDRCLassifierNeuralNet.h"
#include "FragLibReader.h"
#include "MsFrame.h"
#include "MsReaderMzML.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
#include "ProteinDigestomatic.h"
#include "TurboXIC.h"

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
                [](const CandidateScores *l, const CandidateScores *r){
                    return l->discriminantScore < r->discriminantScore;
                }
        );

        const double fdrThreshold = 0.5;
        int counter = 0;
        for (CandidateScores *csp : *candidateScoresTargetsAndDecoys) {

            counter++;

            if (csp->qValue >= fdrThreshold && !csp->isDecoy) {
                break;
            }
        }

        *candidateScoresTargetsAndDecoys50PercentFDRFiltered = *candidateScoresTargetsAndDecoys;
        candidateScoresTargetsAndDecoys50PercentFDRFiltered->resize(static_cast<int>(counter));

        std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);

        const int shuffleCount = 3;
        for (int i = 0; i < shuffleCount; i++) {
            std::shuffle(
                    candidateScoresTargetsAndDecoys50PercentFDRFiltered->begin(),
                    candidateScoresTargetsAndDecoys50PercentFDRFiltered->end(),
                    rng
            );
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath); ree;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    msReaderPointerAcc.ptr->collateMS2MzTargetFrames(&diaTargetFrames);

    const int msLevel = 1;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanPointsMS1;
    e = msReaderPointerAcc.ptr->getScanPoints(msLevel, &scanNumberVsScanPointsMS1); ree;

    QMap<ScanNumber, ScanPoints*> ms1FramePtrs;
    for (auto it = scanNumberVsScanPointsMS1.begin(); it != scanNumberVsScanPointsMS1.end(); it++) {
        ms1FramePtrs.insert(it.key(), &it.value());
    }

    MsFrame msFrameMS1;
    e = msFrameMS1.init(ms1FramePtrs, msReaderPointerAcc.ptr->getScanNumberVsScanTime()); ree;

    TurboXIC turboXICMS1;
    e = turboXICMS1.init(msFrameMS1.frameIndexVsScanPoints()); ree;

    e = m_targetDecoyCandidatePairScoretron.init(
            m_pythiaParameters,
            scanNumberVsScanPointsMS1,
            &msReaderPointerAcc,
            &diaTargetFrames,
            &turboXICMS1
            ); ree;

    if (!m_pythiaParameters.turboMode) {

        QVector<CandidateScores*> candidateScoresTrainings;
        e = buildCalibration(
                &msReaderPointerAcc,
                &diaTargetFrames,
                &scanNumberVsScanPointsMS1,
                &candidateScoresTrainings
                ); ree;


        e = msReaderPointerAcc.ptr->getScanPoints(msLevel, &scanNumberVsScanPointsMS1); ree;
        ms1FramePtrs.clear();
        for (auto it = scanNumberVsScanPointsMS1.begin(); it != scanNumberVsScanPointsMS1.end(); it++) {
            ms1FramePtrs.insert(it.key(), &it.value());
        }
        e = msFrameMS1.init(ms1FramePtrs, msReaderPointerAcc.ptr->getScanNumberVsScanTime()); ree;
        e = turboXICMS1.init(msFrameMS1.frameIndexVsScanPoints()); ree;

        e = m_targetDecoyCandidatePairScoretron.init(
                m_pythiaParameters,
                scanNumberVsScanPointsMS1,
                &msReaderPointerAcc,
                &diaTargetFrames,
                &turboXICMS1
                ); ree;

        e = optimizeParameters(
                candidateScoresTrainings,
                &msReaderPointerAcc
                ); ree;
    }

    int targetCountBelowFDRThreshold;
    e = mainAnalysis(
            &msReaderPointerAcc,
            &targetCountBelowFDRThreshold
            ); ree;

    QVector<CandidateScores*> candidateScoresTargetsAndDecoys;
    e = buildCandidateScoresPtrs(&candidateScoresTargetsAndDecoys); ree;

    QVector<CandidateScores*> candidateScoresTargetsAndDecoys50PercentFDRFiltered;
    e = filterScoredCandidatesTo50PercentFDR(
            &candidateScoresTargetsAndDecoys,
            &candidateScoresTargetsAndDecoys50PercentFDRFiltered
    ); ree;

    qDebug() << "Analyzing" << candidateScoresTargetsAndDecoys50PercentFDRFiltered.size() << "for filtering";
    const int ionsSharedToReject = 4;
    e = removeInterferingCandidates(ionsSharedToReject, &candidateScoresTargetsAndDecoys50PercentFDRFiltered); ree;
    qDebug() << candidateScoresTargetsAndDecoys50PercentFDRFiltered.size() << "after filtering";

    QVector<CandidateScores*> candidateScoreClassifierPntrs;
    const int seedFirstTry = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST;
    e = applyNeuralNetClassifier(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            seedFirstTry,
            &candidateScoreClassifierPntrs
            ); ree;

    if (candidateScoreClassifierPntrs.size() <= targetCountBelowFDRThreshold) {

        std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
        std::shuffle(
                candidateScoresTargetsAndDecoys50PercentFDRFiltered.begin(),
                candidateScoresTargetsAndDecoys50PercentFDRFiltered.end(),
                rng
        );

        const int seedSecondTry = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST + 111;
        candidateScoreClassifierPntrs.clear();
        e = applyNeuralNetClassifier(
                candidateScoresTargetsAndDecoys50PercentFDRFiltered,
                seedSecondTry,
                &candidateScoreClassifierPntrs
                ); ree;

        if (candidateScoreClassifierPntrs.size() <= targetCountBelowFDRThreshold){
            QVector<CandidateScores*> candidateScoresPntrs;
            buildCandidateScoresPtrs(&candidateScoresPntrs);

            std::sort(candidateScoresPntrs.rbegin(), candidateScoresPntrs.rend(), [](CandidateScores *l, CandidateScores *r){
                return l->discriminantScore < r->discriminantScore;
            });

            candidateScoresPntrs.resize(targetCountBelowFDRThreshold);
            candidateScoreClassifierPntrs = candidateScoresPntrs;
        }

    }

    qDebug() << "Updating" << candidateScoreClassifierPntrs.size() << "PSMs";
    e = updateProteinGroupAnnotation(
            m_fastaUri,
            &candidateScoreClassifierPntrs
            ); ree;

    int counter = 0;
    int decoys = 0;
    int entrap = 0;
    for (CandidateScores *cs : candidateScoreClassifierPntrs) {
        counter++;

        if (cs->proteinGroup.contains("_ARATH") && !cs->proteinGroup.contains("_HUMAN")) {
            entrap++;
        }

        if (cs->isDecoy) {
            decoys++;
        }
        if (decoys/(double)counter >= 0.01) {
            break;
        }

    }
    qDebug() << "Counter:" << counter << "Decoys:" <<  decoys << "Entrap:" << entrap << "Entrap%" << entrap / (double)counter;


//    QVector<CandidateScores> candidateScoreClassifier;
//    std::transform(
//            candidateScoreClassifierPntrs.begin(),
//            candidateScoreClassifierPntrs.end(),
//            std::back_inserter(candidateScoreClassifier),
//            [](const CandidateScores *cs){return *cs;}
//            );

//    const QString resultsFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
//    e = ParquetReader::write(candidateScoreClassifier, resultsFilePath); ree;

//    limit what scores are collected depending on extendedScores or NeuralNetScores is set

    ERR_RETURN
}

namespace {

    void filterDuplicateCandidateScoresByDiscriminantScore(QVector<CandidateScores*> *candidateScores) {

        QMap<QString, CandidateScores*> keyVsCandidatesFoundBest;
        for (CandidateScores *cs : *candidateScores) {

            const QString key
                = cs->targetDecoyCandidatePair->peptideStringWithMods() + QString::number(cs->targetDecoyCandidatePair->charge());
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

        const int top6 = 6;
        const auto msCalibrationReaderRowsInsertLogic = [](CandidateScores *cs){

            MsCalibarationReaderRow row;
            row.peptideStringWithMods = cs->targetDecoyCandidatePair->peptideStringWithMods();
            row.iRTPredicted = static_cast<float>(cs->targetDecoyCandidatePair->iRt());
            row.scanTime = cs->scanTime;
            row.mzSearchedVec = cs->featuresArray.mid(CandidateScores::Features::MzSearched1, top6);
            row.mzFoundMeanVec = cs->featuresArray.mid(CandidateScores::Features::MzFoundMean1, top6);
            row.mzFoundStDevVec = cs->featuresArray.mid(CandidateScores::Features::MzFoundStDev1, top6);

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
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1,
        QVector<CandidateScores*> *candidateScoresForTrainings
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

    const auto sizePerTranche = static_cast<double>(m_pythiaParameters.trancheSizeMax);
    const int trancheSize = std::max(static_cast<int>(targetCount / sizePerTranche), 1);
    qDebug() << "Tranch size for calibration:" << trancheSize << "target count:" << targetCount << "sizePerTranche:" << static_cast<int>(sizePerTranche);

    QVector<QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>>> mzTargetKeyVsTargetDecoyCandidatePointersTranched;
    e = ParallelUtils::trancheMapValueVectorsByKeyForParallelization(
            mzTargetKeyVsTargetDecoyCandidatePointers,
            trancheSize,
            &mzTargetKeyVsTargetDecoyCandidatePointersTranched
            ); ree;

    const int topNMS2IonsCalibration = 6;

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> topScoresOfTranches;

    for (int i = 0; i < mzTargetKeyVsTargetDecoyCandidatePointers.size(); i++) {

        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> topScoresOfTranchesAndTrancheCombined = mzTargetKeyVsTargetDecoyCandidatePointersTranched.at(i);;
        for (auto it = topScoresOfTranches.begin(); it != topScoresOfTranches.end(); it++) {
            const MzTargetKey &mzTargetKey = it.key();
            const QVector<TargetDecoyCandidatePair*> &tdcps = it.value();
            topScoresOfTranchesAndTrancheCombined[mzTargetKey].append(tdcps);
        }

        m_candidateScores.clear();

        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2IonsCalibration,
                msReaderPointerAcc->ptr->scanTimeMinMax(),
                m_msCalibratomatic,
                &topScoresOfTranchesAndTrancheCombined,
                &m_candidateScores
                ); ree

        std::sort(
                m_candidateScores.rbegin(),
                m_candidateScores.rend(),
                [](const CandidateScores &l, const CandidateScores &r){
                    return l.featuresArray[CandidateScores::Features::CosineSimSum100]
                         < r.featuresArray[CandidateScores::Features::CosineSimSum100];
                }
                );

        const QPair<double, double> scanTimeMinMax = msReaderPointerAcc->ptr->scanTimeMinMax();
        e = setDiscriminantScoreForCandidates(
                useExtendedScores,
                useNeuralNetworkScores
        ); ree;

        QVector<CandidateScores*> candidateScoresPntrs;
        e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;
        e = setQValueForCandidates(QValueScoreType::DiscriminantScore, &candidateScoresPntrs); ree;

        std::sort(candidateScoresPntrs.rbegin(), candidateScoresPntrs.rend(), [](CandidateScores *l, CandidateScores *r){
            return l->discriminantScore < r->discriminantScore;
        });

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

        const int idealTrainingCountAtGivenFDR = 1000;
        int minTrainingCount = std::max(minTrainingCountTranche, targetCountBelowFDRThreshold);
        minTrainingCount = std::min(minTrainingCount, idealTrainingCountAtGivenFDR);
        qDebug() << "Training RT count 10% FDR:" << minTrainingCount;

        candidateScoresPntrs.resize(minTrainingCount);

        QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
        e = buildMsCalibrationReaderRows(candidateScoresPntrs, &msCalibrationReaderRows); ree;
        topScoresOfTranches.clear();
        for (CandidateScores* p : candidateScoresPntrs) {
            topScoresOfTranches[p->targetKey].push_back(p->targetDecoyCandidatePair);
        }

        const int numberOfStDevs = 3;
        e = m_msCalibratomatic.initRtOnly(msCalibrationReaderRows); ree;
        qDebug() << "scanTimeWindowStDev x" << numberOfStDevs <<":" << m_msCalibratomatic.scanTimeStDev(numberOfStDevs);

        if (minTrainingCount >= idealTrainingCountAtGivenFDR) {

            const int minMzTrainingSize = 200;
            const QString twoPercenFDRKey = "2";
            if (fdrVsCount.value(twoPercenFDRKey) >= minMzTrainingSize) {

                msCalibrationReaderRows.resize(fdrVsCount.value(twoPercenFDRKey));
                e = m_msCalibratomatic.initMzOnly(msCalibrationReaderRows); ree;
                e = recalibrateMzVals(
                        diaTargetFrames,
                        scanNumberVsScanTimeMS1,
                        &m_targetDecoyCandidatePairScoretron,
                        msReaderPointerAcc
                        ); ree;

            }

            QVector<CandidateScores*> candidateScoresPntrsOpt;
            e = buildCandidateScoresPtrs(&candidateScoresPntrsOpt); ree;
            std::sort(candidateScoresPntrsOpt.rbegin(), candidateScoresPntrsOpt.rend(), [](CandidateScores *l, CandidateScores *r){
                return l->discriminantScore < r->discriminantScore;
            });

            candidateScoresPntrsOpt.resize(idealTrainingCountAtGivenFDR * 2);

            *candidateScoresForTrainings = candidateScoresPntrsOpt;
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
                msScanInfo.precursorTargetMz + (msScanInfo.isoWindowUpper + m_pythiaParameters.precursorExtractionWindowThomsons),
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
        return candidateScores.targetDecoyCandidatePair->peptideStringWithMods()
                + QString::number(candidateScores.targetDecoyCandidatePair->charge()) + candidateScores.targetKey + decoyToString;
    }

    struct BuildClassiferParallelInput {
        QVector<ScoresTargets*> scoresTargets;
        QVector<ScoresDecoys*> scoresDecoys;
    };

    struct BuildClassifierParallelOutput {
        QVector<QVector<float>> A;
        QVector<float> b;
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
            QVector<QVector<float>> *A,
            QVector<float> *b
            ) {

        ERR_INIT
        e = result.first; ree;

        const QVector<QVector<float>> &ALocal = result.second.A;
        const QVector<float> bLocal = result.second.b;

        if (A->isEmpty() || b->isEmpty()) {
            A->resize(ALocal.size());
            for(int row = 0; row < A->size(); row++) {
                QVector<float> v(ALocal.front().size(), 0.0);
                (*A)[row] = v;
            }
            b->resize(bLocal.size());
        }

        for (int row = 0; row < A->size(); row++) {
            for (int col = 0; col < A->front().size(); col++) {
                (*A)[row][col] += ALocal[row][col] / static_cast<float>(inputsSize);
            }
        }

        for (int row = 0; row < bLocal.size(); row++) {
            (*b)[row] += bLocal[row] / static_cast<float>(inputsSize);
        }

        ERR_RETURN
    }

    QPair<Err, QVector<float>> scoreVectorParallel(
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            CandidateScores* candidateScores
            ) {

        ERR_INIT


        if (useNeuralNetworkScores) {
            QVector<float> &vec = candidateScores->featuresArray;
            return {e, vec};
        }
        else if (useExtendedScores) {
            QVector<float> vec = candidateScores->featuresArray.mid(0, CandidateScores::Features::CosineSimToAnchor1);
            return {e, vec};
        }
        else {
            QVector<float> vec = candidateScores->featuresArray.mid(0, 5);
            return {e, vec};
        }

        return {eError, {}};
    }

}//namespace
Err PythiaDIAFFWorkflow::setDiscriminantScoreForCandidates(
        bool useExtendedScores,
        bool useNeuralNetworkScores
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(m_candidateScores.isEmpty()); ree;

    QElapsedTimer et;
    et.start();

    const auto scoreBuildBinder = std::bind(
            scoreVectorParallel,
            useExtendedScores,
            useNeuralNetworkScores,
            std::placeholders::_1
    );

    //TODO clean this up. Abstract sub processes into methods.

    QVector<CandidateScores*> candidateScoresPntrs;
    e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;

    QFuture<QPair<Err, QVector<float>>> futuresScoreBuilder = QtConcurrent::mapped(
            candidateScoresPntrs,
            scoreBuildBinder
            );
    futuresScoreBuilder.waitForFinished();

    QVector<QVector<float>> scoreVectors;
    for (const QPair<Err, QVector<float>> &res : futuresScoreBuilder) {
        e = res.first; ree;
        scoreVectors.push_back(res.second);
    }

    e = ErrorUtils::isEqual(m_candidateScores.size(), scoreVectors.size()); ree;

    QMap<QString, TargetDecoyCandidateScores> keyVsTargetDecoyCandidateScores;

    for(int i = 0; i < scoreVectors.size(); i++) {

        CandidateScores *cs = &m_candidateScores[i];

        const QString key = cs->targetDecoyCandidatePair->peptideStringWithMods()
                                + QString::number(cs->targetDecoyCandidatePair->charge()) + cs->targetKey;
        QVector<float> *sv = &scoreVectors[i];

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

    QVector<QVector<float>> A;
    QVector<float> b;
    for (const QPair<Err, BuildClassifierParallelOutput> &result : futures) {
        e = processParallelResults(
                result,
                inputs.size(),
                &A,
                &b
                ); ree;
    }
    qDebug() << "build classifier" << et.restart() << "mSec";

    QVector<float> weights;
    e = ClassifierWeightsManager::fitWeights(A, b, &weights); ree;
    qDebug() << "fit weights" << et.restart() << "mSec";

    qDebug() << "Weights:" << weights;
    qDebug() << "b:" << b;

    QVector<float> discScoreTargets;
    e = ClassifierWeightsManager::applyWeights(scoresTargets, weights, &discScoreTargets); ree;

    QVector<float> discScoreDecoys;
    e = ClassifierWeightsManager::applyWeights(scoresDecoys, weights, &discScoreDecoys); ree;

    qDebug() << "apply weights scores" << et.restart() << "mSec";

    e = ErrorUtils::isEqual(scoresTargets.size(), scoresDecoys.size()); ree;
    e = ErrorUtils::isEqual(scoresTargets.size(), candidateScoresTargetsPtrs.size()); ree;
    e = ErrorUtils::isEqual(scoresDecoys.size(), candidateScoresDecoysPtrs.size()); ree;
    e = ErrorUtils::isEqual(discScoreTargets.size(), scoresTargets.size());
    e = ErrorUtils::isEqual(discScoreDecoys.size(), scoresDecoys.size());

    for (int i = 0; i < scoresDecoys.size(); i++) {

        CandidateScores* csTarget = candidateScoresTargetsPtrs.at(i);
        csTarget->discriminantScore = discScoreTargets.at(i);

        CandidateScores* csDecoy = candidateScoresDecoysPtrs.at(i);
        csDecoy->discriminantScore = discScoreDecoys.at(i);
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
                classifierScore = cs->discriminantScore;
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
                {1.5, 2.0},
                {2.0, 2.0},
                {2.5,  2.0},
                {3.5,  2.0},
                {4.5, 2.0}
        };

        for (const QVector<double> &exp : experiments) {

            PythiaParameters params = pythiaParameters;
            params.ms2ExtractionWidthPPM = mzPPMStDev * exp.at(0);

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
            qDebug() << "ppmTol" << pp.ms2ExtractionWidthPPM;
        }

        ERR_RETURN
    }

    Err getTopFrequencyParameters(
            QVector<DOEResult> *results,
            double *ppmSetting
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(*results); ree;

        std::sort(
                results->rbegin(),
                results->rend(),
                [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
                );

        for (const DOEResult &r : *results) {
            qDebug() << r.ppm << r.scanTimeStDev << r.fdrCount;
        }

        *ppmSetting = results->front().ppm;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::optimizeParameters(
        const QVector<CandidateScores*> &candidateScoresTrainings,
        MsReaderPointerAcc *msReaderPointerAcc
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

    const int topNMS2IonsOptimization = 6;
    qDebug() << "Using top:" << topNMS2IonsOptimization << "fragments for optimization";


    qDebug() << "Optimization selection fraction" << candidateScoresTrainings.size();

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    for (CandidateScores *cs : candidateScoresTrainings) {
        mzTargetKeyVsTargetDecoyCandidatePointers[cs->targetKey].push_back(cs->targetDecoyCandidatePair);
//        qDebug() << cs->targetDecoyCandidatePair->peptideStringWithMods() << cs->targetKey;
    }

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

        qDebug() << "ppmTol" << pythiaParams.ms2ExtractionWidthPPM;

        e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(pythiaParams); ree;
        qDebug() << "STarting opt";

        m_candidateScores.clear();

        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2IonsOptimization,
                msReaderPointerAcc->ptr->scanTimeMinMax(),
                m_msCalibratomatic,
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                &m_candidateScores
        ); ree;

        const QPair<double, double> scanTimeMinMax = msReaderPointerAcc->ptr->scanTimeMinMax();
        e = setDiscriminantScoreForCandidates(
                useExtendedScores,
                useNeuralNetworkScores
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
        res.fdrCount = targetCountAboveFDRQValueThreshold;
        results.push_back(res);
    }

    e = getTopFrequencyParameters(
            &results,
            &m_pythiaParameters.ms2ExtractionWidthPPM
            ); ree;

    e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(m_pythiaParameters); ree;

    qDebug() << "Optimal ppm setting:" << m_pythiaParameters.ms2ExtractionWidthPPM;
    qDebug() << "Optimal scanTimeWindow setting:" << m_msCalibratomatic.scanTimeStDev(m_pythiaParameters.scanTimeWindowStDevs);


    ERR_RETURN
}

Err PythiaDIAFFWorkflow::mainAnalysis(
        MsReaderPointerAcc *msReaderPointerAcc,
        int *targetCountBelowFDRThresholdOnePercent
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

    m_candidateScores.clear();

    const bool useExtendedScores = false;
    const bool useNeuralNetworkScores = false;

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
            msReaderPointerAcc->ptr->scanTimeMinMax(),
            m_msCalibratomatic,
            &mzTargetKeyVsTargetDecoyCandidatePointers,
            &m_candidateScores
    ); ree;
    qDebug() << "Targets scored" << et.restart() << "mSec";

    e = setDiscriminantScoreForCandidates(
            useExtendedScores,
            useNeuralNetworkScores
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

    const double fdrThresholdOnePercent = 0.01;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
            m_candidateScores,
            fdrThresholdOnePercent,
            targetCountBelowFDRThresholdOnePercent
    ); ree;
    qDebug() << "Targets counted" << et.restart() << "mSec";

    std::sort(candidateScoresPntrs.rbegin(), candidateScoresPntrs.rend(), [](CandidateScores *l, CandidateScores *r){
        return l->discriminantScore < r->discriminantScore;
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

Err PythiaDIAFFWorkflow::removeInterferingCandidates(
        int ionsSharedToReject,
        QVector<CandidateScores*> *candidates
        ) {

    ERR_INIT

    QMap<ScanNumber, QVector<CandidateScores*>> scanNumberVsCandidateScoresPtr;
    for (int i = 0; i < candidates->size(); i++) {
        CandidateScores *cs = candidates->at(i);
        scanNumberVsCandidateScoresPtr[cs->scanNumber].push_back(cs);
    }

    QVector<CandidateScores*> filteredCandidates;
    for (auto it = scanNumberVsCandidateScoresPtr.begin(); it != scanNumberVsCandidateScoresPtr.end(); it++) {

        const ScanNumber scanNumber = it.key();
        QVector<CandidateScores*> &scanNumberCandidates = it.value();
        const int foundCandidatesCountInScan = scanNumberCandidates.size();

        std::sort(
                scanNumberCandidates.rbegin(),
                scanNumberCandidates.rend(),
                [](CandidateScores *l, CandidateScores *r){return l->discriminantScore < r->discriminantScore;}
                );

        if (foundCandidatesCountInScan == 1) {
            filteredCandidates.push_back(scanNumberCandidates.front());
            continue;
        }

        const int cols = static_cast<int>(std::round(m_pythiaParameters.mzMaxDataStructure - m_pythiaParameters.mzMinDataStructure));
        const int rows = foundCandidatesCountInScan;
        Eigen::MatrixX<float> mat(rows, cols);
        mat.setZero();

        int top6MzValsFound = 6;

        for (int row = 0; row < rows; row++) {

            const CandidateScores* cs = scanNumberCandidates.at(row);

            const QVector<float> top6MzValsFoundArr = cs->featuresArray.mid(CandidateScores::Features::MzFoundMean1, top6MzValsFound);
            for (float mz : top6MzValsFoundArr) {
                const int col = static_cast<int>(std::round(mz));

                if (col < 0 || col >= cols) {
                    continue;
                }
                mat.coeffRef(row, col) = 1.0f;
            }
        }

        QSet<int> excludeIndexes;
        for (int row = 0; row < rows - 1; row++) {

            const Eigen::VectorX<float> &vecBase = mat.row(row);

            int rowNext = row + 1;
            const Eigen::MatrixX<float> &matNext = mat.block(rowNext, 0, rows - rowNext, cols);

            const Eigen::VectorX<float> dotProds =  matNext * vecBase;

            if (static_cast<int>(dotProds.maxCoeff()) < ionsSharedToReject) {
                continue;
            }

            for (int i = 0; i < dotProds.size(); i++) {
                const float sharedIons = dotProds.coeff(i);
                if (static_cast<int>(sharedIons) >= ionsSharedToReject) {
                    excludeIndexes.insert(row + i + 1);
                }
            }
        }

        for (int i = 0; i < scanNumberCandidates.size(); i++) {
            CandidateScores *cs = scanNumberCandidates.at(i);
            if (excludeIndexes.contains(i)) {
                continue;
            }
            filteredCandidates.push_back(cs);
        }
    }

    candidates->clear();
    *candidates = filteredCandidates;

    ERR_RETURN
}

namespace {

    Err minMaxScaleScores(
            const QVector<KarnnNNTarget> &karnnNNTargets,
            QVector<KarnnNNTarget> *karnnNNTargetsNorm
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(karnnNNTargets); ree;

        QVector<QVector<float>> vecs;
        std::transform(
                karnnNNTargets.begin(),
                karnnNNTargets.end(),
                std::back_inserter(vecs),
                [](const KarnnNNTarget &kt){return kt.scoreVec;}
                );

        Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(vecs);
        EigenUtils::minMaxScaleMatrix(&mat);

        const QVector<QVector<float>> vecsNorm = EigenUtils::convertEigenMatrixToQVectors(mat);

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
            QVector<KarnnNNTarget> *karnnNNTargetsNorm
    ){

        ERR_INIT

        e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoys50PercentFDRFiltered); ree;

        QVector<KarnnNNTarget> karnnNNTargets;
        for (int i = 0; i < candidateScoresTargetsAndDecoys50PercentFDRFiltered.size(); i++) {
            const CandidateScores *cs = candidateScoresTargetsAndDecoys50PercentFDRFiltered.at(i);
            KarnnNNTarget karnnNnTarget;
            karnnNnTarget.seq = cs->targetDecoyCandidatePair->peptideStringWithMods();
            karnnNnTarget.isDecoy = cs->isDecoy;
            karnnNnTarget.index = i;
            karnnNnTarget.scoreVec = cs->featuresArray;
            karnnNNTargets.push_back(karnnNnTarget);
        }

        e = minMaxScaleScores(karnnNNTargets, karnnNNTargetsNorm); ree;

        ERR_RETURN
    }

    Err predictNNScores(
            const QVector<KarnnNNTarget> &karnnNNTargetsNorm,
            int seed,
            QVector<float> *predictions
            ) {

        ERR_INIT

        const int batchSize = std::min(50, std::max(1, static_cast<int>(karnnNNTargetsNorm.size() / 100.0)));
        const int maxIters = 1;
        qDebug() << "Batch Size:" << batchSize << "MaxIters:" << maxIters;

        QVector<QVector<float>> xData;
        QVector<float> yData;
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVec);
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

        e = fdrcLassifierNeuralNet.exec(
                xData,
                yData,
                seed,
                predictions
                ); ree;

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
        const double fdrCutoff = 0.0075; //TODO make this settable
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
        const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys50PercentFDRFiltered,
        int seed,
        QVector<CandidateScores*> *candidateScoreClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_candidateScores); ree;

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
            &karnnNNTargetsNorm
            ); ree;

    QVector<float> predictions;
    e = predictNNScores(
            karnnNNTargetsNorm,
            seed,
            &predictions
            ); ree;

    e = processPredictions(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            predictions,
            m_pythiaParameters.reportDecoys,
            &karnnNNTargetsNorm,
            candidateScoreClassifier
            ); ree;

    e = setQValueForCandidates(QValueScoreType::NNClassifierScore, candidateScoreClassifier); ree

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

        QString peptideSeqReplacedLeucines = cs->targetDecoyCandidatePair->peptideStringWithMods();
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
