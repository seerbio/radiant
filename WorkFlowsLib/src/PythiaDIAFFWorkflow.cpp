//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAFFWorkflow.h"

#include "CandidateScores.h"
#include "ClassifierWeightsManager.h"
#include "Deconvolvotron.h"
#include "DiscriminantScoretron.h"
#include "EigenUtils.h"
#include "FastaReader.h"
#include "FDRCLassifierNeuralNet.h"
#include "FragLibReader.h"
#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.h"
#include "MsReaderPointerAcc.h"
#include "PythiaDIAFFWorkflowAlgos/OptimizeMassAccuracyPPMSettertron.h"
#include "ParallelUtils.h"
#include "PeptideStringWithMods.h"
#include "QuanFileBuilder.h"
#include "QuanTransitionRefinertron.h"
#include "QValueSettertron.h"
#include "SequenceSubstringSearchomatic.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

#include <QElapsedTimer>

#include "PythiaDIAFFWorkflowAlgos/PythiaDIAFFWorkflowSharedMethods.h"


PythiaDIAFFWorkflow::PythiaDIAFFWorkflow()
: m_minTopNMs2Ions(6)
, m_minTrainingCountTranche(50)
{}

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
    m_pythiaParameters.print();

    QVector<FragLibReaderRow> fragLibReaderRows;
    e = FragLibReader::getFragLibReaderRows(
            m_fragLibUri,
            &fragLibReaderRows
            ); ree;

    e = m_targetDecoyCandidatePairManager.init(
            m_pythiaParameters,
            &fragLibReaderRows
            ); ree;

    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(&m_targetDecoyPairPntrs); ree

    ERR_RETURN
}

namespace {

    Err filterScoredCandidatesForNeuralNet(
            int minMs2FragCount,
            QVector<CandidateScores*> *candidateScoresTargetsAndDecoys
            ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoresTargetsAndDecoys->isEmpty()); ree;

        const auto terminatorLogic = [minMs2FragCount](CandidateScores *cs) {
            return cs->featuresArray[Features::CosineSimSum100] < static_cast<float>(minMs2FragCount);
        };
        const auto terminator = std::remove_if(candidateScoresTargetsAndDecoys->begin(), candidateScoresTargetsAndDecoys->end(), terminatorLogic);
        candidateScoresTargetsAndDecoys->erase(terminator, candidateScoresTargetsAndDecoys->end());

        std::sort(
                candidateScoresTargetsAndDecoys->rbegin(),
                candidateScoresTargetsAndDecoys->rend(),
                [](const CandidateScores *l, const CandidateScores *r){
                    return l->discriminantScore < r->discriminantScore;
                });

        int counter = 0;
        for (const CandidateScores *csp : *candidateScoresTargetsAndDecoys) {

            counter++;
            if (constexpr double fdrTrainingThreshold = 0.85; csp->qValue >= fdrTrainingThreshold && !csp->isDecoy) {
                break;
            }
        }

        candidateScoresTargetsAndDecoys->resize(counter);

        std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);

        const int shuffleCount = 3;
        for (int i = 0; i < shuffleCount; i++) {
            std::shuffle(
                    candidateScoresTargetsAndDecoys->begin(),
                    candidateScoresTargetsAndDecoys->end(),
                    rng
            );
        }

        ERR_RETURN
    }

    Err populateAltIdTargetKeysLogic(const QVector<CandidateScores*> &candidateScoreses) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(candidateScoreses); ree;

        for (CandidateScores *csOG : candidateScoreses) {

            for (const CandidateScores *csAlt : candidateScoreses) {

                if (csOG == csAlt) {
                    continue;
                }

                if (csOG->targetKey == csAlt->targetKey) {
                    csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                            - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                        / csOG->featuresArray[Features::CosineSimSum100];
                    csOG->featuresArray[Features::DiscriminantScore] = csOG->discriminantScore;
                    continue;
                }

                switch (csAlt->targetDecoyCandidatePair->charge()) {
                case 1:
                    if (!csAlt->isDecoy) {
                        csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                            = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                                - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                            / csOG->featuresArray[Features::CosineSimSum100];
                        csOG->featuresArray[Features::AltTargetKeyIdTimeDeltaCharge1_1] = std::abs((csOG->scanTime - csAlt->scanTime) / csOG->scanTime);
                        break;
                    }
                    csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                            - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                        / csOG->featuresArray[Features::CosineSimSum100];
                    break;

                case 2:
                    if (!csAlt->isDecoy) {
                        csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                            = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                                - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                            / csOG->featuresArray[Features::CosineSimSum100];
                        csOG->featuresArray[Features::AltTargetKeyIdTimeDeltaCharge2_1] = std::abs((csOG->scanTime - csAlt->scanTime) / csOG->scanTime);
                        break;
                    }
                    csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                            - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                        / csOG->featuresArray[Features::CosineSimSum100];
                    break;

                case 3:
                    if (!csAlt->isDecoy) {
                        csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                            = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                                - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                            / csOG->featuresArray[Features::CosineSimSum100];
                        csOG->featuresArray[Features::AltTargetKeyIdTimeDeltaCharge3_1] = std::abs((csOG->scanTime - csAlt->scanTime) / csOG->scanTime);
                        break;
                    }
                    csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                            - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                        / csOG->featuresArray[Features::CosineSimSum100];
                    break;

                case 4:
                    if (!csAlt->isDecoy) {
                        csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                            = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                                - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                            / csOG->featuresArray[Features::CosineSimSum100];
                        csOG->featuresArray[Features::AltTargetKeyIdTimeDeltaCharge4_1] = std::abs((csOG->scanTime - csAlt->scanTime) / csOG->scanTime);
                        break;
                    }
                    csOG->featuresArray[Features::AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[Features::CosineSimSum100] * csOG->featuresArray[Features::CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[Features::CosineSim100MS1])
                                            - (csAlt->featuresArray[Features::CosineSimSum100]  * csAlt->featuresArray[Features::CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[Features::CosineSim100MS1]))
                                        / csOG->featuresArray[Features::CosineSimSum100];
                    break;

                    default:
                        rrr(eValueError);
                }
            }
        }

        ERR_RETURN
    }

    Err populateAltIdTargetKeys(
        int threadCount,
        QVector<CandidateScores*> *candidateScoresPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoresPntrs->isEmpty()); ree;

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PopulateAltIdTargetKeys";

        QMap<PeptideStringWithMods, QVector<CandidateScores*>> pepStrWModsVsCandScoresEntries;
        for (CandidateScores *cs : *candidateScoresPntrs) {
            pepStrWModsVsCandScoresEntries[cs->targetDecoyCandidatePair->peptideStringWithMods()].push_back(cs);
        }

        for (const QVector<CandidateScores*> &csPntrs : pepStrWModsVsCandScoresEntries) {
            e = populateAltIdTargetKeysLogic(csPntrs); ree;
        }

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PopulateAltIdTargetKeys Finished";

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {
    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    MsReaderPointerAcc msReaderPointerAcc;

    // m_pythiaParameters.useLazyLoading = true;
    // m_pythiaParameters.ms2ExtractionWidthPPMOverride = 7.5;
    msReaderPointerAcc.setUseLazyLoading(m_pythiaParameters.useLazyLoading);

    e = msReaderPointerAcc.openFile(msDataFilePath); ree;
    msReaderPointerAcc.ptr->printSize();

    e = m_targetDecoyCandidatePairScoretron.init(
            m_pythiaParameters,
            &msReaderPointerAcc
            ); ree;

    {
        MsCalibratomaticSettertron msCalibratomaticSettertron;
        e = msCalibratomaticSettertron.init(
            &m_pythiaParameters,
            &msReaderPointerAcc,
            &m_targetDecoyCandidatePairManager,
            &m_targetDecoyCandidatePairScoretron
            ); ree;
        e = msCalibratomaticSettertron.buildCalibration(&m_msCalibratomatic); ree;
    }

    if (m_pythiaParameters.ms2ExtractionWidthPPMOverride > 0.0) {
        m_pythiaParameters.ms2ExtractionWidthPPM = m_pythiaParameters.ms2ExtractionWidthPPMOverride;
        m_pythiaParameters.ms1ExtractionWidthPPM = m_pythiaParameters.ms1ExtractionWidthPPMOverride;

        if (m_pythiaParameters.ms1ExtractionWidthPPMOverride < 0.0) {
            m_pythiaParameters.ms1ExtractionWidthPPM = m_pythiaParameters.ms2ExtractionWidthPPM;
        }
        e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(m_pythiaParameters); ree;

        qDebug()
        << qPrintable(S_GLOBAL_TIMER.elapsed())
        << "Ms2 Accuracy overridden to" << m_pythiaParameters.ms2ExtractionWidthPPM;
        qDebug()
        << qPrintable(S_GLOBAL_TIMER.elapsed())
        << "Ms1 Accuracy overridden to" << m_pythiaParameters.ms1ExtractionWidthPPM;
    }
    else {
         OptimizeMassAccuracyPPMSettertron msOptimizeMassAccuracyPPMSettertron;
         e = msOptimizeMassAccuracyPPMSettertron.initExec(
                &msReaderPointerAcc,
                &m_msCalibratomatic,
                &m_pythiaParameters,
                &m_targetDecoyCandidatePairScoretron,
                &m_targetDecoyPairPntrs
                ); ree;
    }

    int targetCountBelowFDRThreshold;
    e = mainAnalysis(
        &msReaderPointerAcc,
        &targetCountBelowFDRThreshold
     ); ree;

    QVector<CandidateScores*> candidateScoresTargetsAndDecoys;
    e = PythiaDIAFFWorkflowSharedMethods::buildCandidateScoresPtrs(
        m_candidateScores,
        &candidateScoresTargetsAndDecoys
        ); ree;

    e = populateAltIdTargetKeys(
        m_pythiaParameters.threadCount,
        &candidateScoresTargetsAndDecoys
        ); ree;

    QVector<CandidateScores*> candidateScoreClassifierPntrs;
    e = applyNeuralNetClassifier(
            candidateScoresTargetsAndDecoys,
            S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST,
            &candidateScoreClassifierPntrs
            ); ree;

    int targetCountBelowFDRThresholdOnePercent;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
        candidateScoreClassifierPntrs,
        m_pythiaParameters.percentFDR / 100.0,
        &targetCountBelowFDRThresholdOnePercent
        ); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "Pre Neural Net PSMs Count" << targetCountBelowFDRThreshold
            << "| Post Neural Net Count PSMs" << targetCountBelowFDRThresholdOnePercent;

    const bool candidateScoresSortedHiLo = std::is_sorted(
        candidateScoreClassifierPntrs.begin(),
        candidateScoreClassifierPntrs.end(),
        [](const CandidateScores *l, const CandidateScores *r){return l->classifierScore < r->classifierScore;}
        );
    e = ErrorUtils::isTrue(candidateScoresSortedHiLo); ree;

    if (!m_pythiaParameters.reportDecoys) {

        const auto terminatorLogicFDRFilter
            = [&](CandidateScores *cs){return cs->qValue > m_pythiaParameters.percentFDR / 100.0;};

        const auto terminator = std::remove_if(
                candidateScoreClassifierPntrs.begin(),
                candidateScoreClassifierPntrs.end(),
                terminatorLogicFDRFilter
                );

        candidateScoreClassifierPntrs.erase(terminator, candidateScoreClassifierPntrs.end());
    }
    else {

        int counter = 0;
        int decoyCounter = 0;
        for (const CandidateScores *candidateScores : candidateScoreClassifierPntrs) {

            counter++;
            if (candidateScores->isDecoy) {
                decoyCounter++;
            }

            constexpr double fdrThreshold = 0.5;
            if (decoyCounter / static_cast<double>(counter) >= fdrThreshold) {
                break;
            }
        }

        candidateScoreClassifierPntrs.resize(counter);
    }
    m_pythiaParameters.writePythiaDIA = true;
    if (m_pythiaParameters.writePythiaDIA) {

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Annotating" << candidateScoreClassifierPntrs.size() << "PSMs";
        e = updateProteinGroupAnnotation(
                m_fastaUri,
                &candidateScoreClassifierPntrs
                ); ree;
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Annotation finished";

        QVector<CandidateScoresReaderRow> candidateScoreReaderRows;
        std::transform(
                candidateScoreClassifierPntrs.begin(),
                candidateScoreClassifierPntrs.end(),
                std::back_inserter(candidateScoreReaderRows),
                [](const CandidateScores *cs){return CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);}
                );

        const QString resultsFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
        e = ParquetReader::write(candidateScoreReaderRows, resultsFilePath); ree;

//TODO delete this after dev is done.
#define REPORT_ENTRAP
#ifdef REPORT_ENTRAP
        if (m_pythiaParameters.verbosity > -1) {
            int counter = 0;
            int decoys = 0;
            int entrap = 0;

            for (CandidateScores *cs : candidateScoreClassifierPntrs) {
                counter++;

                if (cs->proteinGroup.contains("_ARATH") && !cs->proteinGroup.contains("_HUMAN") && !cs->isDecoy) {
                    entrap++;
                }

                if (cs->isDecoy) {
                    decoys++;
                }
                if (decoys/static_cast<double>(counter) >= 0.01) {
                    break;
                }

            }
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                    << "Alt:" << targetCountBelowFDRThresholdOnePercent
                    << "| Counter:" << counter
                    << "| Decoys:" <<  decoys
                    << "| Entrap:" << entrap
                    << "| Entrap%" << entrap / (double)counter;
        }
#endif
    }
    // constexpr int frameIndexBuffer = 1;
    // QuanTransitionRefinertron quanTransitionRefinertron(m_pythiaParameters.ms2ExtractionWidthPPM, frameIndexBuffer);
    // const QString quanFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_QUAN_FILE_EXTENSION;
    // e = QuanFileBuilder::buildQuanFile(
    //     candidateScoreClassifierPntrs,
    //     quanFilePath
    //     ); ree;

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::mainAnalysis(
        const MsReaderPointerAcc *msReaderPointerAcc,
        int *targetCountBelowFDRThresholdOnePercent
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

    m_candidateScores.clear();

    constexpr bool useExtendedScores = true;
    constexpr bool useNeuralNetworkScores = false;
    constexpr int topNMs2IonsMainAnalysis = 12;
    constexpr bool useTopNIntegrationsParameter = false;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "Using top:"
            << topNMs2IonsMainAnalysis
            << "fragments for main analysis";

    QElapsedTimer et;
    et.start();

    const QVector<MsScanInfo> uniqueMsScanInfos = msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = PythiaDIAFFWorkflowSharedMethods::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            m_targetDecoyPairPntrs,
            m_pythiaParameters,
            uniqueMsScanInfos,
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

    qDebug()
    << qPrintable(S_GLOBAL_TIMER.elapsed())
    << "TargetKeys size:" << mzTargetKeyVsTargetDecoyCandidatePointers.size();

    QMap<MzTargetKey, TurboXIC*> nullToBuildTurboXICInParallelLoop;

    constexpr int splitter = 2;
    const int threadCount = uniqueMsScanInfos.size() < m_pythiaParameters.threadCount
                          ? std::min(uniqueMsScanInfos.size() * splitter, m_pythiaParameters.threadCount)
                          : m_pythiaParameters.threadCount;

    m_weights = DiscriminantScoretron::defaultWeights(useExtendedScores, useNeuralNetworkScores);

    constexpr float minPeakCount = 3.9;
    m_candidateScores.clear();
    e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
            topNMs2IonsMainAnalysis,
            m_msCalibratomatic,
            minPeakCount,
            threadCount,
            useExtendedScores,
            useNeuralNetworkScores,
            useTopNIntegrationsParameter,
            nullToBuildTurboXICInParallelLoop,
            m_weights,
            &mzTargetKeyVsTargetDecoyCandidatePointers,
            &m_candidateScores
            ); ree
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Targets scored" << et.restart() << "mSec";

    QVector<CandidateScores*> candidateScoresVecBatchPntrs;
    QMap<int, int> fdrVsCounts;
    QVector<float> weights;
    e = PythiaDIAFFWorkflowSharedMethods::processBatch(
        m_candidateScores,
        m_pythiaParameters,
        useExtendedScores,
        useNeuralNetworkScores,
        &candidateScoresVecBatchPntrs,
        &fdrVsCounts,
        &weights
        ); ree;

    QString fdrString;
    e = FDRCLassifierNeuralNet::outPutFDRCounts(fdrVsCounts, &fdrString); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << qPrintable(fdrString);
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Targets resulted" << et.restart() << "mSec";

    constexpr double fdrThresholdOnePercent = 0.01;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
            m_candidateScores,
            fdrThresholdOnePercent,
            targetCountBelowFDRThresholdOnePercent
    ); ree;

    PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(&candidateScoresVecBatchPntrs);
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Targets sorted" << et.restart() << "mSec";

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
                [](const KarnnNNTarget &kt){return kt.scoreVecNormalized;}
                );

        Eigen::MatrixX<float> mat = EigenUtils::convertQVectorsToEigenMatrix(vecs);
        EigenUtils::minMaxScaleMatrix(&mat);

        const QVector<QVector<float>> vecsNorm = EigenUtils::convertEigenMatrixToQVectors(mat);

        e = ErrorUtils::isEqual(vecsNorm.size(), karnnNNTargets.size()); ree;

        for (int i = 0; i < vecsNorm.size(); i++) {
            KarnnNNTarget ktNew = karnnNNTargets.at(i);
            ktNew.scoreVecNormalized = vecsNorm.at(i);
            karnnNNTargetsNorm->push_back(ktNew);
        }

        ERR_RETURN
    }

    Err buildKarnnNNTargetsNormalized(
            const QVector<CandidateScores*> &candidateScoresTargetsAndDecoysFDRFiltered,
            QVector<KarnnNNTarget> *karnnNNTargetsNorm
    ){

        ERR_INIT

        e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoysFDRFiltered); ree;

        QVector<KarnnNNTarget> karnnNNTargets;
        karnnNNTargets.reserve(candidateScoresTargetsAndDecoysFDRFiltered.size());
        for (int i = 0; i < candidateScoresTargetsAndDecoysFDRFiltered.size(); i++) {
            CandidateScores *cs = candidateScoresTargetsAndDecoysFDRFiltered.at(i);
            KarnnNNTarget karnnNnTarget;
            karnnNnTarget.candidateScores = cs;
            karnnNnTarget.scoreVecNormalized = DiscriminantScoretron::scoreVectorLogic(true, true, cs);

            karnnNNTargets.push_back(karnnNnTarget);
        }

        e = minMaxScaleScores(karnnNNTargets, karnnNNTargetsNorm); ree;

        ERR_RETURN
    }

    Err trainNeuralNetwork(
            const QVector<KarnnNNTarget> &karnnNNTargetsNorm,
            int seed,
            int threadCount,
            int verbosity,
            FDRCLassifierNeuralNet *fdrcLassifierNeuralNet
            ) {

        ERR_INIT

        constexpr int batchSizeMin = 500;
        const int batchSize = std::min(batchSizeMin, std::max(1, static_cast<int>(karnnNNTargetsNorm.size() / 100.0)));

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Batch Size:" << batchSize;

        QVector<QVector<float>> xData;
        QVector<float> yData;
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVecNormalized);
            yData.push_back(kt.candidateScores->isDecoy ? 1.0 : 0.0);
        }

        constexpr int baggingSize = 4;
        constexpr float learningRate = 0.003;
        constexpr int epochs = 3; //TODO make this settable

        e = fdrcLassifierNeuralNet->init(
                epochs,
                baggingSize,
                batchSize,
                learningRate,
                threadCount
        ); ree;

        e = fdrcLassifierNeuralNet->trainClassifier(
                xData,
                yData,
                seed,
                verbosity
                ); ree;

        ERR_RETURN
    }

    Err predictClassifierScores(
        const QVector<KarnnNNTarget> &karnnNNTargetsNorm,
        FDRCLassifierNeuralNet *fdrcLassifierNeuralNet,
        QVector<float> *predictions
        ) {

        ERR_INIT

        QVector<QVector<float>> xData;
        QVector<float> yData;
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVecNormalized);
            yData.push_back(kt.candidateScores->isDecoy ? 1.0 : 0.0);
        }

        e = fdrcLassifierNeuralNet->predictBaggedClassifiers(
            xData,
            predictions
            ); ree;

        ERR_RETURN
    }

    Err processPredictions(
            const QVector<float> &predictions,
            QVector<KarnnNNTarget> *karnnNNTargetsNorm
            ) {

        ERR_INIT

        e = ErrorUtils::isFalse(karnnNNTargetsNorm->isEmpty()); ree;

        for (int i = 0; i < predictions.size(); i++) {
            (*karnnNNTargetsNorm)[i].candidateScores->classifierScore = predictions.at(i);
        }

        ERR_RETURN
    }

    Err serialFilterByValue(
        double threshold,
        QVector<CandidateScores*> *candidateScoreses
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoreses->isEmpty()); ree;

        std::sort(
            candidateScoreses->begin(),
            candidateScoreses->end(),
            [](const CandidateScores *l, const CandidateScores *r){
                return l->classifierScore < r->classifierScore;
            });

        int counter = 0;
        for (const CandidateScores *csp : *candidateScoreses) {

            counter++;
            if (csp->qValue >= threshold && !csp->isDecoy) {
                break;
            }
        }

        candidateScoreses->resize(counter);

        ERR_RETURN
    }

    Err subsetKarnnNNTargetsForTraining(
        const QVector<KarnnNNTarget> &karnnNNTargetsNorm,
        QVector<KarnnNNTarget> *karnnNNTargetsNormTrain
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(karnnNNTargetsNorm); ree;

        *karnnNNTargetsNormTrain = karnnNNTargetsNorm;
        std::sort(
            karnnNNTargetsNormTrain->rbegin(),
            karnnNNTargetsNormTrain->rend(),
            [](const KarnnNNTarget &l, const KarnnNNTarget &r) {
                return l.candidateScores->discriminantScore < r.candidateScores->discriminantScore;
            });

        int counter = 0;
        for (const KarnnNNTarget &knt : *karnnNNTargetsNormTrain) {
            counter++;
            if (constexpr double fdrTrainingThreshold = 0.65; knt.candidateScores->qValue >= fdrTrainingThreshold && !knt.candidateScores->isDecoy) {
                break;
            }
        }
        karnnNNTargetsNormTrain->resize(counter);

        std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
        constexpr int shuffleCount = 3;
        for (int i = 0; i < shuffleCount; i++) {
            std::shuffle(
                    karnnNNTargetsNormTrain->begin(),
                    karnnNNTargetsNormTrain->end(),
                    rng
            );
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::applyNeuralNetClassifier(
        const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys,
        int seed,
        QVector<CandidateScores*> *candidateScoreClassifier
        ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_candidateScores); ree;

    QVector<CandidateScores*> candidateScoresTargetsAndDecoysNeuralNet = candidateScoresTargetsAndDecoys;
    e = filterScoredCandidatesForNeuralNet(
        m_pythiaParameters.minMs2FragCount,
        &candidateScoresTargetsAndDecoysNeuralNet
        ); ree;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Analyzing" << candidateScoresTargetsAndDecoysNeuralNet.size() << "for filtering";

// #define WRITENN
#ifdef WRITENN
    QFile file("/home/anichols/Desktop/scores.csv");
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {

    }

    QTextStream out(&file);

    for (const CandidateScores &cs : m_candidateScores) {
        for (float value : cs.featuresArray) {
            out << value << ",";
        }
        out << "\n";
    }
    file.close();
#endif

    candidateScoreClassifier->clear();

    const int totalCount = candidateScoresTargetsAndDecoysNeuralNet.size();
    const int decoyCount = static_cast<int>(std::count_if(
            candidateScoresTargetsAndDecoysNeuralNet.begin(),
            candidateScoresTargetsAndDecoysNeuralNet.end(),
            [](const CandidateScores* cs){return cs->isDecoy;}
            ));

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "target vs decoy count for neural net training"
             << totalCount - decoyCount << ":" << decoyCount
             << "total" << totalCount;

    QVector<KarnnNNTarget> karnnNNTargetsNorm;
    e = buildKarnnNNTargetsNormalized(
            candidateScoresTargetsAndDecoysNeuralNet,
            &karnnNNTargetsNorm
            ); ree;

    QVector<KarnnNNTarget> karnnNNTargetsNormTrain;
    e = subsetKarnnNNTargetsForTraining(
        karnnNNTargetsNorm,
        &karnnNNTargetsNormTrain
        ); ree;

    FDRCLassifierNeuralNet fdrClassifierNeuralNet;
    e = trainNeuralNetwork(
            karnnNNTargetsNormTrain,
            seed,
            m_pythiaParameters.threadCount,
            m_pythiaParameters.verbosity,
            &fdrClassifierNeuralNet
            ); ree;

    QVector<float> predictions;
    e = predictClassifierScores(
        karnnNNTargetsNorm,
        &fdrClassifierNeuralNet,
        &predictions
        ); ree;

    e = processPredictions(
            predictions,
            &karnnNNTargetsNorm
            ); ree;

    *candidateScoreClassifier = candidateScoresTargetsAndDecoysNeuralNet;

    std::sort(
        candidateScoreClassifier->begin(),
        candidateScoreClassifier->end(),
        [](CandidateScores *l, CandidateScores *r){return l->classifierScore < r->classifierScore;}
        );

    e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::NNClassifierScore,
            candidateScoreClassifier
            ); ree

    constexpr double fdrQValThreshold = 0.5;
    e = serialFilterByValue(
        fdrQValThreshold,
        candidateScoreClassifier
        ); ree;

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::updateProteinGroupAnnotation(
        const QString &fastaFilePath,
        QVector<CandidateScores*> *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    SequenceSubstringSearchomatic sequenceSubstringSearchomatic;

    FastaReader fastaReader;
    e = fastaReader.parseFastaFile(fastaFilePath); ree;
    e = sequenceSubstringSearchomatic.init(fastaReader.fastaEntries()); ree;

    QStringList searchSequences;
    std::transform(
        candidateScores->begin(),
        candidateScores->end(),
        std::back_inserter(searchSequences),
        [](CandidateScores *cs){return cs->targetDecoyCandidatePair->peptideString();}
        );

    QVector<QString> proteinGroupsAll;
    e = sequenceSubstringSearchomatic.findProteinGroups(searchSequences, &proteinGroupsAll);
    e = ErrorUtils::isEqual(
        proteinGroupsAll.size(),
        candidateScores->size()
        ); ree;

    for (int i = 0; i < proteinGroupsAll.size(); i++) {

        if (candidateScores->at(i)->isDecoy) {
            QStringList proteinGroupSplit = proteinGroupsAll.at(i).split(";");
            for (QString &pg : proteinGroupSplit) {
                pg = pg.replace('>', ">decoy_");
            }

            const QString updatedProteinGroup = proteinGroupSplit.join(";");
            (*candidateScores)[i]->proteinGroup = updatedProteinGroup;
            continue;
        }

        (*candidateScores)[i]->proteinGroup = proteinGroupsAll.at(i);
    }

    ERR_RETURN
}

// namespace {
//
//     struct SpectrumCentricParallelInput {
//         MzTargetKey mzTargetKey;
//         QVector<CandidateScores*> candidateScoresPntrs;
//         QMap<ScanNumber, ScanPoints*> diaTargetFrame;
//         PythiaParameters pythiaParameters;
//         MsCalibratomatic msCalibratomatic;
//         QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
//     };
//
//     Err buildSpectrumCentricParallelInput(
//         const QMap<MzTargetKey, QVector<CandidateScores*>> &mzTargetKeyVsCandidateScoresPointers,
//         const QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> &diaTargetFrames,
//         const MsCalibratomatic &msCalibratomatic,
//         const PythiaParameters &pythiaParameters,
//         const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
//         QVector<SpectrumCentricParallelInput> *spectrumCentricParallelInputs
//         ) {
//
//         ERR_INIT
//         e = ErrorUtils::isNotEmpty(mzTargetKeyVsCandidateScoresPointers); ree;
//         e = ErrorUtils::isNotEmpty(diaTargetFrames); ree;
//
//         for(auto it = mzTargetKeyVsCandidateScoresPointers.begin(); it != mzTargetKeyVsCandidateScoresPointers.end(); ++it) {
//             SpectrumCentricParallelInput spectrumCentricParallelInput;
//             spectrumCentricParallelInput.mzTargetKey = it.key();
//             spectrumCentricParallelInput.candidateScoresPntrs = it.value();
//             spectrumCentricParallelInput.msCalibratomatic = msCalibratomatic;
//             spectrumCentricParallelInput.pythiaParameters = pythiaParameters;
//             spectrumCentricParallelInput.scanNumberVsScanTime = scanNumberVsScanTime;
//
//             e = ErrorUtils::contains(spectrumCentricParallelInput.mzTargetKey, diaTargetFrames); ree;
//             spectrumCentricParallelInput.diaTargetFrame = diaTargetFrames.value(spectrumCentricParallelInput.mzTargetKey);
//             spectrumCentricParallelInputs->push_back(spectrumCentricParallelInput);
//         }
//
//         ERR_RETURN
//     }
//
//     QPair<Err, QVector<QPair<CandidateScores*, DeconvolvotronResult>>> spectrumCentricParallelLogic(const SpectrumCentricParallelInput &scpi) {
//
//         ERR_INIT
//
//         e = ErrorUtils::isNotEmpty(scpi.candidateScoresPntrs); rree;
//         e = ErrorUtils::isNotEmpty(scpi.diaTargetFrame); rree;
//         e = ErrorUtils::isNotEmpty(scpi.mzTargetKey); rree;
//
//         QElapsedTimer et;
//         et.start();
//
//         SpectrumCentricMzTargetFrameSearch spectrumCentricMzTargetFrameSearch;
//         e = spectrumCentricMzTargetFrameSearch.init(
//             scpi.pythiaParameters,
//             scpi.msCalibratomatic,
//             scpi.diaTargetFrame,
//             scpi.candidateScoresPntrs,
//             scpi.scanNumberVsScanTime
//             ); rree;
//
//         QVector<QPair<CandidateScores*, DeconvolvotronResult>> candidateScoresPntrsVsScore;
//         e = spectrumCentricMzTargetFrameSearch.assignIdsToScans(&candidateScoresPntrsVsScore); rree;
//
//         qDebug() << "Processed MzTargetKey" << scpi.mzTargetKey << et.elapsed() << "mSec";
//
//         return {e, candidateScoresPntrsVsScore};
//     }
//
// }//namespace
// Err PythiaDIAFFWorkflow::spectrumCentricSearch(
//     const QVector<CandidateScores*> &candidateScoresPntrs,
//     const MsCalibratomatic &msCalibratomatic,
//     const MsReaderPointerAcc *msReaderPointerAcc
//     ) const {
//
//     ERR_INIT
//     e = ErrorUtils::isNotEmpty(m_targetDecoyPairPntrs); ree;
//
//     QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs;
//     for (CandidateScores *cs : candidateScoresPntrs) {
//         mzTargetKeyVsCandidateScoresPntrs[cs->targetKey].push_back(cs);
//     }
//
//     QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
//     e = msReaderPointerAcc->ptr->collateMS2MzTargetFrames(&diaTargetFrames); ree;
//
//     QVector<SpectrumCentricParallelInput> spectrumCentricParallelInputs;
//     e = buildSpectrumCentricParallelInput(
//         mzTargetKeyVsCandidateScoresPntrs,
//         diaTargetFrames,
//         msCalibratomatic,
//         m_pythiaParameters,
//         msReaderPointerAcc->ptr->getScanNumberVsScanTime(),
//         &spectrumCentricParallelInputs
//         ); ree;
//
// #define PARALLEL_DCON
// #ifdef PARALLEL_DCON
//     QFuture<QPair<Err, QVector<QPair<CandidateScores*, DeconvolvotronResult>>>> futures = QtConcurrent::mapped(
//         spectrumCentricParallelInputs,
//         spectrumCentricParallelLogic
//         ); ree;
//     futures.waitForFinished();
//
//     for (const QPair<Err, QVector<QPair<CandidateScores*, DeconvolvotronResult>>> &res : futures) {
//         e = res.first; ree;
//         const QVector<QPair<CandidateScores*, DeconvolvotronResult>> &prs = res.second;
//     }
// #else
//     for (const SpectrumCentricParallelInput &inp : spectrumCentricParallelInputs) {
//         QPair<Err, QVector<QPair<CandidateScores*, DeconvolvotronResult>>> res = spectrumCentricParallelLogic(inp); ree;
//         e = res.first; ree;
//         qDebug() << res.second.size() << "SDLFKDJSL";
//     }
// #endif
//
//     ERR_RETURN
// }
