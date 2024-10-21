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

    const double massMin
        = pythiaParameters.peptideLengthMin * Molecule(MolecularFormulas::glycineFormula).monoisotopicMass();

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

    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(&m_targetDecoyPairPntrs); ree

    ERR_RETURN
}

namespace {

    Err filterScoredCandidatesForNeuralNet(
            int minMs2FragCount,
            QVector<CandidateScores*> *candidateScoresTargetsAndDecoys,
            QVector<CandidateScores*> *candidateScoresTargetsAndDecoys50PercentFDRFiltered
            ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoresTargetsAndDecoys->isEmpty()); ree;

        const auto terminatorLogic = [minMs2FragCount](CandidateScores *cs) {
            return cs->featuresArray[CandidateScores::Features::CosineSimSum100] < minMs2FragCount;
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
            if (constexpr double fdrThreshold = 0.5; csp->qValue >= fdrThreshold && !csp->isDecoy) {
                break;
            }
        }

        *candidateScoresTargetsAndDecoys50PercentFDRFiltered = *candidateScoresTargetsAndDecoys;
        candidateScoresTargetsAndDecoys50PercentFDRFiltered->resize(counter);

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

    Err populateAltIdTargetKeys(QVector<CandidateScores*> *candidateScoresPntrs) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoresPntrs->isEmpty()); ree;

        struct CandidateScoresEntry {
            int charge = -1;
            float cosineSimSum100 = -1.0;
            float scanTime = -1.0;
            MzTargetKey mzTargetKey;
            bool isDecoy = false;
        };

        QMap<PeptideStringWithMods, QVector<CandidateScoresEntry>> pepStrWModsVsCandScoresEntries;

        for (CandidateScores *cs : *candidateScoresPntrs) {

            CandidateScoresEntry cse;
            cse.mzTargetKey = cs->targetKey;
            cse.charge = cs->targetDecoyCandidatePair->charge();
            cse.scanTime = cs->scanTime;
            cse.cosineSimSum100 = cs->featuresArray[CandidateScores::Features::CosineSimSum100];
            cse.isDecoy = cs->isDecoy;

            pepStrWModsVsCandScoresEntries[cs->targetDecoyCandidatePair->peptideStringWithMods()].push_back(cse);
        }

        for (QVector<CandidateScoresEntry> &cse : pepStrWModsVsCandScoresEntries) {
            std::sort(cse.begin(), cse.end(), [](const CandidateScoresEntry &l, const CandidateScoresEntry &r){
                if (l.charge == r.charge) {

                    if (l.mzTargetKey == r.mzTargetKey) {
                        return l.cosineSimSum100 > r.cosineSimSum100;
                    }

                    return l.mzTargetKey < r.mzTargetKey;
                }
                return l.charge< r.charge;
            });
        }

        for (CandidateScores *cs : *candidateScoresPntrs) {

            const QVector<CandidateScoresEntry> &csEntries
                    = pepStrWModsVsCandScoresEntries.value(cs->targetDecoyCandidatePair->peptideStringWithMods());

            for(const CandidateScoresEntry &cse : csEntries) {

                constexpr float scoreThreshold = 0.01;

                if (cse.mzTargetKey == cs->targetKey) {

                    switch (cs->targetDecoyCandidatePair->charge()) {
                        case 1:
                            if(!MathUtils::tSame(cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG], cs->featuresArray[CandidateScores::Features::CosineSimSum100])
                                && cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG] > scoreThreshold) {
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_1] = cse.cosineSimSum100;
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_1] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                                break;
                            }

                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG] = cs->featuresArray[CandidateScores::Features::CosineSimSum100];
                            break;

                        case 2:
                            if(!MathUtils::tSame(cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG], cs->featuresArray[CandidateScores::Features::CosineSimSum100])
                               && cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG] > scoreThreshold) {
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_1] = cse.cosineSimSum100;
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_1] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                                break;
                            }

                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG] = cs->featuresArray[CandidateScores::Features::CosineSimSum100];
                            break;
                        case 3:
                            if(!MathUtils::tSame(cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG], cs->featuresArray[CandidateScores::Features::CosineSimSum100])
                               && cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG] > scoreThreshold) {
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_1] = cse.cosineSimSum100;
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_1] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                                break;
                            }

                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG] = cs->featuresArray[CandidateScores::Features::CosineSimSum100];
                            break;
                        case 4:
                            if(!MathUtils::tSame(cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG], cs->featuresArray[CandidateScores::Features::CosineSimSum100])
                               && cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG] > scoreThreshold) {
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_1] = cse.cosineSimSum100;
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_1] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                                break;
                            }

                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG] = cs->featuresArray[CandidateScores::Features::CosineSimSum100];
                            break;
                        default:
                            rrr(eValueError);
                    }

                    continue;
                }

                switch (cse.charge) {

                    case 1:
                        if (cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_2] > scoreThreshold) {
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_3] = cse.cosineSimSum100;
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_3] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                            break;
                        }

                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_2] = cse.cosineSimSum100;
                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_2] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                        break;
                    case 2:
                        if (cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_2] > scoreThreshold) {
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_3] = cse.cosineSimSum100;
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_3] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                            break;
                        }

                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_2] = cse.cosineSimSum100;
                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_2] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                        break;
                    case 3:
                        if (cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2] > scoreThreshold) {
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_3] = cse.cosineSimSum100;
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_3] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                            break;
                        }

                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2] = cse.cosineSimSum100;
                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_2] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                        break;
                    case 4:
                        if (cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_2] > scoreThreshold) {
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_3] = cse.cosineSimSum100;
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_3] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                            break;
                        }

                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_2] = cse.cosineSimSum100;
                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_2] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                        break;

                    default:
                        rrr(eValueError);
                }

            }
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {
    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    MsReaderPointerAcc msReaderPointerAcc;

    e = msReaderPointerAcc.openFile(msDataFilePath); ree;
    msReaderPointerAcc.ptr->printSize();

    e = m_targetDecoyCandidatePairScoretron.init(
            m_pythiaParameters,
            &msReaderPointerAcc
            ); ree;

    MsCalibratomaticSettertron msCalibratomaticSettertron;
    e = msCalibratomaticSettertron.init(
        &m_pythiaParameters,
        &msReaderPointerAcc,
        &m_targetDecoyCandidatePairManager,
        &m_targetDecoyCandidatePairScoretron
        ); ree;
    e = msCalibratomaticSettertron.buildCalibration(&m_msCalibratomatic); ree;

    if (m_msCalibratomatic.isInitRT()) {

        OptimizeMassAccuracyPPMSettertron msOptimizeMassAccuracyPPMSettertron;
        e = msOptimizeMassAccuracyPPMSettertron.initExec(
            &msReaderPointerAcc,
            &m_msCalibratomatic,
            &m_pythiaParameters,
            &m_targetDecoyCandidatePairScoretron,
            &m_targetDecoyPairPntrs
            ); ree;

        m_weights = msOptimizeMassAccuracyPPMSettertron.weights();
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

    QVector<CandidateScores*> candidateScoresTargetsAndDecoysNeuralNet;
    e = filterScoredCandidatesForNeuralNet(
        m_pythiaParameters.minMs2FragCount,
        &candidateScoresTargetsAndDecoys,
        &candidateScoresTargetsAndDecoysNeuralNet
        ); ree;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Analyzing" << candidateScoresTargetsAndDecoysNeuralNet.size() << "for filtering";

    e = populateAltIdTargetKeys(&candidateScoresTargetsAndDecoysNeuralNet); ree;

    QVector<CandidateScores*> candidateScoreClassifierPntrs;
    e = applyNeuralNetClassifier(
            candidateScoresTargetsAndDecoysNeuralNet,
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

    QElapsedTimer etTrans;
    etTrans.start();
    QVector<CandidateScores*> candidateScoreClassifierPntrsFDRFiltered = candidateScoreClassifierPntrs;
    const auto terminator = std::remove_if(
        candidateScoreClassifierPntrsFDRFiltered.begin(),
        candidateScoreClassifierPntrsFDRFiltered.end(),
        [&](const CandidateScores *cs){return cs->isDecoy || cs->qValue >= m_pythiaParameters.percentFDR / 100.0f;}
        );
    candidateScoreClassifierPntrsFDRFiltered.erase(terminator, candidateScoreClassifierPntrsFDRFiltered.end());

    constexpr int frameIndexBuffer = 1;
    QuanTransitionRefinertron quanTransitionRefinertron(m_pythiaParameters.ms2ExtractionWidthPPM, frameIndexBuffer);
    e = quanTransitionRefinertron.refineTransitions(candidateScoreClassifierPntrsFDRFiltered); ree;

    const QString quanFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_QUAN_FILE_EXTENSION;
    e = QuanFileBuilder::buildQuanFile(
        candidateScoreClassifierPntrs,
        quanFilePath
        ); ree;

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::mainAnalysis(
        MsReaderPointerAcc *msReaderPointerAcc,
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

    QMap<MzTargetKey, TurboXIC*> nullToBuildTurboXICInParallelLoop;

    constexpr int splitter = 2;
    const int threadCount = uniqueMsScanInfos.size() < m_pythiaParameters.threadCount
                          ? std::min(uniqueMsScanInfos.size() * splitter, m_pythiaParameters.threadCount)
                          : m_pythiaParameters.threadCount;

    m_weights = m_weights.isEmpty()
              ? DiscriminantScoretron::defaultWeights(useExtendedScores, useNeuralNetworkScores)
              : m_weights;

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
            karnnNnTarget.seq = cs->targetDecoyCandidatePair->peptideStringWithMods();
            karnnNnTarget.isDecoy = cs->isDecoy;
            karnnNnTarget.index = i;
            karnnNnTarget.scoreVec = DiscriminantScoretron::scoreVectorLogic(true, true, cs);

            karnnNNTargets.push_back(karnnNnTarget);
        }

        e = minMaxScaleScores(karnnNNTargets, karnnNNTargetsNorm); ree;

        ERR_RETURN
    }

    Err predictNNScores(
            const QVector<KarnnNNTarget> &karnnNNTargetsNorm,
            int seed,
            int threadCount,
            int verbosity,
            QVector<float> *predictions
            ) {

        ERR_INIT

        constexpr int batchSizeMin = 500;
        const int batchSize = std::min(batchSizeMin, std::max(1, static_cast<int>(karnnNNTargetsNorm.size() / 100.0)));

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Batch Size:" << batchSize;

        QVector<QVector<float>> xData;
        QVector<float> yData;
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVec);
            yData.push_back(kt.isDecoy ? 1.0 : 0.0);
        }

        constexpr int baggingSize = 6;
        constexpr float learningRate = 0.003;
        constexpr int epochs = 3; //TODO make this settable
        FDRCLassifierNeuralNet fdrcLassifierNeuralNet;
        e = fdrcLassifierNeuralNet.init(
                epochs,
                baggingSize,
                batchSize,
                learningRate,
                threadCount
        ); ree;

        e = fdrcLassifierNeuralNet.exec(
                xData,
                yData,
                seed,
                verbosity,
                predictions
                ); ree;

        ERR_RETURN
    }

    Err processPredictions(
            const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            const QVector<float> &predictions,
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

        std::transform(
                karnnNNTargetsNorm->begin(),
                karnnNNTargetsNorm->end(),
                std::back_inserter(*candidateScoreClassifier),
                [candidateScoresTargetsAndDecoys50PercentFDRFiltered](const KarnnNNTarget &kt){
                    CandidateScores *candidateScoresNew = candidateScoresTargetsAndDecoys50PercentFDRFiltered.at(kt.index);
                    candidateScoresNew->classifierScore = kt.nnScore;
                    return candidateScoresNew;
                }
                );

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::applyNeuralNetClassifier(
        const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys50PercentFDRFiltered,
        int seed,
        QVector<CandidateScores*> *candidateScoreClassifier
        ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_candidateScores); ree;

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

    const int totalCount = candidateScoresTargetsAndDecoys50PercentFDRFiltered.size();
    const int decoyCount = static_cast<int>(std::count_if(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered.begin(),
            candidateScoresTargetsAndDecoys50PercentFDRFiltered.end(),
            [](const CandidateScores* cs){return cs->isDecoy;}
            ));

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "target vs decoy count for neural net training"
             << totalCount - decoyCount << ":" << decoyCount
             << "total" << totalCount;

//#define PRINT_AVERAGES
#ifdef PRINT_AVERAGES
    Eigen::VectorX<float> vec(CandidateScores::Features::FeaturesSize);
    vec.setZero();
    for (const CandidateScores *cs : candidateScoresTargetsAndDecoys50PercentFDRFiltered) {
        vec += EigenUtils::convertQVectorToEigenVector(cs->featuresArray);
    }
    vec /= candidateScoresTargetsAndDecoys50PercentFDRFiltered.size();
    for(int i = 0; i < vec.size(); i++) {
        qDebug() << i << vec.coeff(i);
    }
    einfo;
#endif

    QVector<KarnnNNTarget> karnnNNTargetsNorm;
    e = buildKarnnNNTargetsNormalized(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            &karnnNNTargetsNorm
            ); ree;

    QVector<float> predictions;
    e = predictNNScores(
            karnnNNTargetsNorm,
            seed,
            m_pythiaParameters.threadCount,
            m_pythiaParameters.verbosity,
            &predictions
            ); ree;

    e = processPredictions(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered,
            predictions,
            &karnnNNTargetsNorm,
            candidateScoreClassifier
            ); ree;

    e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::NNClassifierScore,
            candidateScoreClassifier
            ); ree

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
    e = ErrorUtils::isEqual(proteinGroupsAll.size(), candidateScores->size()); ree;

    for (int i = 0; i < proteinGroupsAll.size(); i++) {
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
