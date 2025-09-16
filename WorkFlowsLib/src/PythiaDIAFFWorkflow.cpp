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
#include "IonMobilitron.h"
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
, m_calibratomaticFeatures(DiscriminantScoretron::featuresCalibration())
, m_ppmOptimizationFeatures(DiscriminantScoretron::featuresOptimization())
, m_neuralNetFeatures(DiscriminantScoretron::featuresNeuralNetwork())
{}

PythiaDIAFFWorkflow::~PythiaDIAFFWorkflow() {

    for (FeatureFinderHillBuilder *h : m_scanNumberVsFeatureFinderHillBuildersPntrsTIMS) {
        delete h;
    }

}

Err PythiaDIAFFWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &fastaUri,
        const QString &outputFolderPath
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::fileExists(fragLibUri); ree;

    if (!fastaUri.isEmpty()) {
        e = ErrorUtils::fileExists(fastaUri); ree;
        m_fastaUri = fastaUri;
    }

    ParallelUtils::printSystemDetails();

    m_pythiaParameters = pythiaParameters;

/***** DEV OVERRIDES *****/
// #define DEV_OVERRIDES
#ifdef DEV_OVERRIDES
    // m_pythiaParameters.useLazyLoading = true;
    // m_pythiaParameters.ms2ExtractionWidthPPMOverride = 100;
    // m_pythiaParameters.peakCenter = 4;
    // m_pythiaParameters.writePythiaDIA = false;
    m_pythiaParameters.reannotate = true;
		// m_pythiaParameters.baggingSize = 12;
		// m_pythiaParameters.epochs = 12;
		// m_pythiaParameters.nodesFraction = 0.5;
    // m_pythiaParameters.baggingSize = 4;
    // m_pythiaParameters.threadCount = 8;
    // m_pythiaParameters.shortReport = true;
    // m_pythiaParameters.calibrationTrainingVolume = 100000;
    // m_pythiaParameters.ionsSharedToReject = 4;
    // m_pythiaParameters.filterLengthIntegration = 7;
    // m_pythiaParameters.maxAnchorColumnIndex = 6;
    // m_pythiaParameters.minMs2FragCount = 2;
    // m_pythiaParameters.stopThresholdFractionMS2 = 0.75;
    // m_pythiaParameters.calibrationTrainingVolume = 1000;

    qDebug() << "ACTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
#endif
/**************************/

    m_fragLibUri = fragLibUri;
    m_outputFolderPath = outputFolderPath;
    m_pythiaParameters.print();

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Reading library";

    e = FragLibReader::getFragLibReaderRows(
            m_fragLibUri,
            &m_fragLibReaderRows
            ); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Finished reading library";

    e = m_targetDecoyCandidatePairManager.init(
            m_pythiaParameters,
            &m_fragLibReaderRows
            ); ree;

    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(&m_targetDecoyPairPntrs); ree //TODO HERHHERHEHREH

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
            return cs->featuresArray[CosineSimSum100] < static_cast<float>(minMs2FragCount);
        };
        const auto terminator = std::remove_if(
        	candidateScoresTargetsAndDecoys->begin(),
        	candidateScoresTargetsAndDecoys->end(),
        	terminatorLogic
        	);
        candidateScoresTargetsAndDecoys->erase(terminator, candidateScoresTargetsAndDecoys->end());

        PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(candidateScoresTargetsAndDecoys);

    	e = QValueSettertron::setQValueForCandidates(
			QValueSettertron::QValueScoreType::DiscriminantScore,
			candidateScoresTargetsAndDecoys
			); ree;

    	QMap<int, int> fdrVsCount;
    	e = FDRCLassifierNeuralNet::outputFDRResults(
			*candidateScoresTargetsAndDecoys,
			true,
			&fdrVsCount
			); ree;

    	constexpr int fdrTrainingThresholdInt = 50;

        candidateScoresTargetsAndDecoys->resize(fdrVsCount.value(fdrTrainingThresholdInt));

        std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);

        constexpr int shuffleCount = 3;
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
                    csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                            - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                        / csOG->featuresArray[CosineSimSum100];
                    csOG->featuresArray[DiscriminantScore] = csOG->discriminantScore;
                    continue;
                }

                switch (csAlt->targetDecoyCandidatePair->charge()) {
                case 1:
                    if (!csAlt->isDecoy) {
                        csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                            = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                                - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                            / csOG->featuresArray[CosineSimSum100];
                        csOG->featuresArray[AltTargetKeyIdTimeDeltaCharge1_1] = std::abs((csOG->scanTime - csAlt->scanTime) / csOG->scanTime);
                        break;
                    }
                    csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                            - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                        / csOG->featuresArray[CosineSimSum100];
                    break;

                case 2:
                    if (!csAlt->isDecoy) {
                        csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                            = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                                - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                            / csOG->featuresArray[CosineSimSum100];
                        csOG->featuresArray[AltTargetKeyIdTimeDeltaCharge2_1] = std::abs((csOG->scanTime - csAlt->scanTime) / csOG->scanTime);
                        break;
                    }
                    csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                            - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                        / csOG->featuresArray[CosineSimSum100];
                    break;

                case 3:
                    if (!csAlt->isDecoy) {
                        csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                            = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                                - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                            / csOG->featuresArray[CosineSimSum100];
                        csOG->featuresArray[AltTargetKeyIdTimeDeltaCharge3_1] = std::abs((csOG->scanTime - csAlt->scanTime) / csOG->scanTime);
                        break;
                    }
                    csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                            - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                        / csOG->featuresArray[CosineSimSum100];
                    break;

                case 4:
                    if (!csAlt->isDecoy) {
                        csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                            = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                                - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                            / csOG->featuresArray[CosineSimSum100];
                        csOG->featuresArray[AltTargetKeyIdTimeDeltaCharge4_1] = std::abs((csOG->scanTime - csAlt->scanTime) / csOG->scanTime);
                        break;
                    }
                    csOG->featuresArray[AltTargetKeyIdDiscScoreChargeOG_alt]
                                        = ((csOG->featuresArray[CosineSimSum100] * csOG->featuresArray[CosineSimSpectrumOverTimeCubed] * csOG->featuresArray[CosineSim100MS1])
                                            - (csAlt->featuresArray[CosineSimSum100]  * csAlt->featuresArray[CosineSimSpectrumOverTimeCubed]  * csAlt->featuresArray[CosineSim100MS1]))
                                        / csOG->featuresArray[CosineSimSum100];
                    break;

                    default:
                        rrr(eValueError);
                }
            }
        }

        ERR_RETURN
    }

    Err populateAltIdTargetKeys(
        QVector<CandidateScores*> *candidateScoresPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoresPntrs->isEmpty()); ree;

        PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(candidateScoresPntrs);

        QHash<PeptideStringWithMods, QVector<CandidateScores*>> pepStrWModsVsCandScoresEntries;
        for (CandidateScores *cs : *candidateScoresPntrs) {
            pepStrWModsVsCandScoresEntries[cs->targetDecoyCandidatePair->peptideStringWithMods()].push_back(cs);
        }

        for (const QVector<CandidateScores*> &csPntrs : pepStrWModsVsCandScoresEntries) {
            e = populateAltIdTargetKeysLogic(csPntrs); ree;
        }

        ERR_RETURN
    }

    Err writePythiaDIA(
        const QVector<CandidateScores*> &candidateScoresPntrs,
        bool writeShortReport,
        const QString &outputFolderPath,
        MsReaderPointerAcc *msReaderPointerAcc
        ) {
        ERR_INIT

        if (writeShortReport) {

            QVector<CandidateScoresReaderRowTrunc> candidateScoreReaderRows;
            std::transform(
                    candidateScoresPntrs.begin(),
                    candidateScoresPntrs.end(),
                    std::back_inserter(candidateScoreReaderRows),
                    [](const CandidateScores *cs){return CandidateScoresReaderRowTrunc::buildCandidateScoresReaderRow(cs);}
                    );

            QString resultsFilePath = msReaderPointerAcc->ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;

            if (!outputFolderPath.isEmpty()) {
                const QFileInfo fileInfo(resultsFilePath);
                const QString msDataFileName = fileInfo.fileName();
                resultsFilePath = outputFolderPath + msDataFileName;
            }
            e = ParquetReader::write(candidateScoreReaderRows, resultsFilePath); ree;

            ERR_RETURN
        }

        QVector<CandidateScoresReaderRow> candidateScoreReaderRows;
        std::transform(
                candidateScoresPntrs.begin(),
                candidateScoresPntrs.end(),
                std::back_inserter(candidateScoreReaderRows),
                [](const CandidateScores *cs){return CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);}
                );

        QString resultsFilePath = msReaderPointerAcc->ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;

        if (!outputFolderPath.isEmpty()) {
            const QFileInfo fileInfo(resultsFilePath);
            const QString msDataFileName = fileInfo.fileName();
            resultsFilePath = outputFolderPath + msDataFileName;
        }
        e = ParquetReader::write(candidateScoreReaderRows, resultsFilePath); ree;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    MsReaderPointerAcc msReaderPointerAcc;

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
            m_calibratomaticFeatures,
            &m_pythiaParameters,
            &msReaderPointerAcc,
            &m_targetDecoyCandidatePairManager,
            &m_targetDecoyCandidatePairScoretron,
            false
            ); ree;
        e = msCalibratomaticSettertron.buildCalibration(&m_msCalibratomatic);

        if (e != eNoError) {

            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Retrying calibration w/ exclusion";

            MsCalibratomaticSettertron msCalibratomaticSettertronBackup;
            e = msCalibratomaticSettertronBackup.init(
                m_calibratomaticFeatures,
                &m_pythiaParameters,
                &msReaderPointerAcc,
                &m_targetDecoyCandidatePairManager,
                &m_targetDecoyCandidatePairScoretron,
                true
                ); ree;
            e = msCalibratomaticSettertronBackup.buildCalibration(&m_msCalibratomatic);
            if (e != eNoError) {
                qWarning() << "Could not find enough significant PSMs for calibration";
                eee_absorb;
            }
        }
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
            m_ppmOptimizationFeatures,
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
        m_candidateScorePairs,
        &candidateScoresTargetsAndDecoys
        ); ree;

    e = populateAltIdTargetKeys(&candidateScoresTargetsAndDecoys); ree;

    // if (msReaderPointerAcc.ptr->isTIMS()) {
    //
    //     e = IonMobilitron::assignIonMobilityValues(
    //         m_pythiaParameters,
    //         candidateScoresTargetsAndDecoys,
    //         &m_scanNumberVsFeatureFinderHillBuildersPntrsTIMS
    //         ); ree;
    //
    //     e = predictIonMobilityIndexes(candidateScoresTargetsAndDecoys); ree;
    // }

// #define WRITE_DISC_RESULTS
#ifdef WRITE_DISC_RESULTS
    const QString resultsFilePath = "disc_score_results_" + QString::number(m_pythiaParameters.threadCount) +".prq";
    QVector<CandidateScoresReaderRow> candidateScoreReaderRows;
    std::transform(
            candidateScoresTargetsAndDecoys.begin(),
            candidateScoresTargetsAndDecoys.end(),
            std::back_inserter(candidateScoreReaderRows),
            [](const CandidateScores *cs){return CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);}
            );
    e = ParquetReader::write(candidateScoreReaderRows, resultsFilePath); ree;
#endif

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
        [](const CandidateScores *l, const CandidateScores *r) {
            if (MathUtils::tSame(l->classifierScore, r->classifierScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                return l->discriminantScore > r->discriminantScore;
            }
            return l->classifierScore < r->classifierScore;
        });
    e = ErrorUtils::isTrue(candidateScoresSortedHiLo); ree;

    filterDecoysOrNot(&candidateScoreClassifierPntrs);

    if (m_pythiaParameters.reannotate) {
        e = updateProteinGroupAnnotation(
                m_fastaUri,
                targetCountBelowFDRThresholdOnePercent,
                &candidateScoreClassifierPntrs
                ); ree;
    }

    if (m_pythiaParameters.writePythiaDIA) {
        e = writePythiaDIA(
            candidateScoreClassifierPntrs,
            m_pythiaParameters.shortReport,
            m_outputFolderPath,
            &msReaderPointerAcc
            ); ree;
    }

    QString quanFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_QUAN_FILE_EXTENSION;
    if (!m_outputFolderPath.isEmpty()) {
        const QFileInfo fileInfo(quanFilePath);
        const QString msDataFileName = fileInfo.fileName();
        quanFilePath = m_outputFolderPath + msDataFileName;
    }
    e = QuanFileBuilder::buildQuanFile(
        candidateScoreClassifierPntrs,
        quanFilePath
        ); ree;

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::setCalibratomaticFeatures(const QVector<Features> &features) {
    ERR_INIT
    e = ErrorUtils::isNotEmpty(features); ree;
    m_calibratomaticFeatures = features;
    ERR_RETURN
}

Err PythiaDIAFFWorkflow::setPPMOptimizationFeatures(const QVector<Features> &features) {
    ERR_INIT
    e = ErrorUtils::isNotEmpty(features); ree;
    m_ppmOptimizationFeatures = features;
    ERR_RETURN
}

Err PythiaDIAFFWorkflow::setNeuralNetFeatures(const QVector<Features> &features) {
    ERR_INIT
    e = ErrorUtils::isNotEmpty(features); ree;
    m_neuralNetFeatures = features;
    ERR_RETURN
}

ResultsSummary PythiaDIAFFWorkflow::resultsSummary() const {
    return m_resultsSummary;
}

Err PythiaDIAFFWorkflow::mainAnalysis(
        const MsReaderPointerAcc *msReaderPointerAcc,
        int *targetCountBelowFDRThresholdOnePercent
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

    m_candidateScorePairs.clear();

    constexpr int topNMs2IonsMainAnalysis = 12;
    constexpr bool useTopNIntegrationsParameter = false;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "Using top:"
            << topNMs2IonsMainAnalysis
            << "fragments for main analysis";

    QElapsedTimer et;
    et.start();

    const QVector<MsScanInfo> uniqueMsScanInfos = msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    qDebug()
    << qPrintable(S_GLOBAL_TIMER.elapsed())
    << "TargetKeys size:" << uniqueMsScanInfos.size();

    constexpr int splitter = 2;
    const int threadCount = uniqueMsScanInfos.size() < m_pythiaParameters.threadCount
                          ? std::min(uniqueMsScanInfos.size() * splitter, m_pythiaParameters.threadCount)
                          : m_pythiaParameters.threadCount;

    m_weights = DiscriminantScoretron::defaultWeights(m_ppmOptimizationFeatures);

    constexpr float minPeakCount = 3.9;
    m_candidateScorePairs.clear();
    e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
            m_ppmOptimizationFeatures,
            topNMs2IonsMainAnalysis,
            m_msCalibratomatic,
            minPeakCount,
            threadCount,
            useTopNIntegrationsParameter,
            uniqueMsScanInfos,
            m_weights,
            &m_targetDecoyPairPntrs,
            &m_candidateScorePairs
            ); ree
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Targets scored" << et.restart() << "mSec";

    QVector<CandidateScores*> candidateScoresVecBatchPntrs;
    QMap<int, int> fdrVsCounts;
    QVector<float> weights;
    e = PythiaDIAFFWorkflowSharedMethods::processBatch(
        m_ppmOptimizationFeatures,
        m_candidateScorePairs,
        m_pythiaParameters,
        &candidateScoresVecBatchPntrs,
        &fdrVsCounts,
        &weights
        ); ree;

    QString fdrString;
    e = FDRCLassifierNeuralNet::outPutFDRCounts(fdrVsCounts, &fdrString); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Targets resulted" << et.restart() << "mSec";
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << qPrintable(fdrString);

    constexpr double fdrThresholdOnePercent = 0.01;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
            candidateScoresVecBatchPntrs,
            fdrThresholdOnePercent,
            targetCountBelowFDRThresholdOnePercent
    ); ree;

    PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(&candidateScoresVecBatchPntrs);

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
        const QVector<Features> &neuralNetFeatures,
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

// #define WRITE_KARNNN_NORM_TO_FILE
#ifdef WRITE_KARNNN_NORM_TO_FILE
            karnnNnTarget.scoreVecNormalized = = CandidateScores::selectFeaturesArrayFeatures(
                    cs->featuresArray,
                    DiscriminantScoretron::featuresNeuralNetwork()
                    );
#else
        karnnNnTarget.scoreVecNormalized = CandidateScores::selectFeaturesArrayFeatures(
            cs->featuresArray,
            neuralNetFeatures
            );
#endif

            karnnNNTargets.push_back(karnnNnTarget);
        }

        e = minMaxScaleScores(karnnNNTargets, karnnNNTargetsNorm); ree;

        ERR_RETURN
    }

	QPair<Err, FDRCLassifierNeuralNet> trainNeuralNetworkLogic(
		const QPair<QVector<QVector<float>>, QVector<float>> &trainingVecs,
		const PythiaParameters &pythiaParameters,
		int batchSize
		) {
	    ERR_INIT
    	FDRCLassifierNeuralNet fdrcLassifierNeuralNet;

    	constexpr int baggingOverride = 1;
    	e = fdrcLassifierNeuralNet.init(
				pythiaParameters.epochs,
				baggingOverride,
				batchSize,
				pythiaParameters.learningRate,
				pythiaParameters.nodesFraction,
				pythiaParameters.threadCount
		); rree;

    	e = fdrcLassifierNeuralNet.trainClassifier(
				trainingVecs.first,
				trainingVecs.second,
				S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST,
				pythiaParameters.verbosity
				); rree;

    	return {e, fdrcLassifierNeuralNet};
    }

    Err trainNeuralNetwork(
            const PythiaParameters &pythiaParameters,
            const QVector<QVector<KarnnNNTarget>> &karnnNNTargetsNormTranched,
            int seed,
            QVector<FDRCLassifierNeuralNet> *fdrcLassifierNeuralNets,
            QVector<QVector<KarnnNNTarget>> *inferenceKarnnVecs
            ) {

        ERR_INIT

    	e = ErrorUtils::isTrue(karnnNNTargetsNormTranched.size() > 1); ree;

    	const int trainingSize = std::accumulate(
    		karnnNNTargetsNormTranched.begin(),
    		karnnNNTargetsNormTranched.end(),
    		0,
    		[](int total, const QVector<KarnnNNTarget> &karnnNNTargets){return total + karnnNNTargets.size();}
    		);

        constexpr int batchSizeMin = 500;
        const int batchSize = std::min(batchSizeMin, std::max(1, static_cast<int>(trainingSize / 100.0)));
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Batch Size:" << batchSize;

    	QVector<QPair<QVector<QVector<float>>, QVector<float>>> trainingVecs(karnnNNTargetsNormTranched.size());
    	inferenceKarnnVecs->resize(karnnNNTargetsNormTranched.size());

    	for (int i = 0; i < karnnNNTargetsNormTranched.size(); i++) {
    		for (int j = 0; j < karnnNNTargetsNormTranched.size(); j++) {

    			const QVector<KarnnNNTarget> &karnnNNTargetsNorm = karnnNNTargetsNormTranched.at(j);

    			QVector<QVector<float>> xData;
    			QVector<float> yData;

    			if (i == j) {
    				(*inferenceKarnnVecs)[i] = karnnNNTargetsNorm;
    				continue;
    			}

    			for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
    				xData.push_back(kt.scoreVecNormalized);
    				yData.push_back(kt.candidateScores->isDecoy ? 1.0 : 0.0);
    			}

    			trainingVecs[i].first.append(xData);
    			trainingVecs[i].second.append(yData);
    		}
    	}

// #define WRITENN_NORM
#ifdef WRITENN_NORM
        QFile file("nn_train_data.csv");
        if (!file.open(QIODevice::WriteOnly)) {
            rrr(eFileError)
        }
        QDataStream out(&file);
        out << xData;
        out << yData;
        file.close();
#endif

		if (pythiaParameters.parallelNeuralNets) {
			const auto loadLogicBinder = std::bind(
				trainNeuralNetworkLogic,
				std::placeholders::_1,
				pythiaParameters,
				batchSize
				);

			QFuture<QPair<Err, FDRCLassifierNeuralNet>> future = QtConcurrent::mapped(
				trainingVecs,
				loadLogicBinder
				);
			future.waitForFinished();

			for (const QPair<Err, FDRCLassifierNeuralNet> &result : future) {
				e = result.first; ree;
				fdrcLassifierNeuralNets->push_back(result.second);
			}
		}
		else {
			for (int i = 0; i < karnnNNTargetsNormTranched.size(); i++) {
				const QPair<Err, FDRCLassifierNeuralNet> result = trainNeuralNetworkLogic(
					trainingVecs[i],
					pythiaParameters,
					batchSize
					);

				e = result.first; ree;
				fdrcLassifierNeuralNets->push_back(result.second);
			}
		}

        ERR_RETURN
    }

    QPair<Err, QVector<float>> predictClassiferScoresLogic(
        const QVector<KarnnNNTarget> &karnnNNTargetsNorm,
        FDRCLassifierNeuralNet *fdrcLassifierNeuralNet
        ) {

        ERR_INIT

        QVector<float> predictions;

        QVector<QVector<float>> xData;
        xData.reserve(karnnNNTargetsNorm.size());
        QVector<float> yData;
        yData.reserve(karnnNNTargetsNorm.size());
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVecNormalized);
            yData.push_back(kt.candidateScores->isDecoy ? 1.0 : 0.0);
        }

        e = fdrcLassifierNeuralNet->predictBaggedClassifiers(
            xData,
            &predictions
            ); rree;

        // Normalize predictions to a standard range, using a linear transformation
        // that fixes a 1% FDR cutoff and the median of decoy scores (in log space).

        // 1. Compute 1% FDR cutoff

        // Create temporary candidate scores for FDR calculation using QValueSettertron
        QVector<CandidateScores*> tempCandidateScores;
        tempCandidateScores.reserve(karnnNNTargetsNorm.size());

        for (int i = 0; i < karnnNNTargetsNorm.size(); ++i) {
            // Create a temporary copy of the candidate with updated classifier score
            CandidateScores *tempCs = new CandidateScores(*karnnNNTargetsNorm[i].candidateScores);
            tempCs->classifierScore = predictions[i];
            tempCandidateScores.push_back(tempCs);
        }

        // Use QValueSettertron to calculate q-values based on classifier scores
        e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::NNClassifierScore,
            &tempCandidateScores
        ); rree;

        // Find the highest score where q-value <= 1% FDR
        constexpr double targetFDR = 0.01;
        float fdrCutoff = 0.0f;

        for (const CandidateScores *cs : tempCandidateScores) {
            if (!cs->isDecoy && cs->qValue <= targetFDR) {
                fdrCutoff = std::max(fdrCutoff, static_cast<float>(cs->classifierScore));
            }
        }

        // Clean up temporary candidates
        for (CandidateScores *cs : tempCandidateScores) {
            delete cs;
        }
        tempCandidateScores.clear();

        // 2. Compute median of decoys

        const int decoyCount = std::accumulate(yData.begin(), yData.end(), 0,
            [](int sum, float val) { return sum + static_cast<int>(val); });
        QVector<float> decoyScores;
        decoyScores.reserve(decoyCount);
        for (int i = 0; i < predictions.size(); ++i) {
            if (yData[i] == 1.0f) {
                decoyScores.push_back(predictions[i]);
            }
        }
        const double decoyMedian = static_cast<float>(MathUtils::median(decoyScores));

        // 3. Compute slope and intercept of linear transformation (in log space)
        const double logCutoff = std::log(fdrCutoff);
        const double m = 1.0f / (std::log(decoyMedian) - logCutoff);
        const double b = -m * logCutoff;

        // 4. Apply linear transformation to predictions:
        //
        //    For each score:
        // a. Convert to log space
        // b. Apply linear transformation
        // c. Convert back from log space
        //
        // Future enhancement: use Eigen to vectorize this operation
        // and employ fused-multiply-add (FMA) intrinsic.

        QVector<float> normPredictions;
        normPredictions.reserve(predictions.size());

        for (int i = 0; i < predictions.size(); ++i) {
            double logScore = std::log(predictions[i]);
            double transformedLogScore = m * logScore + b;
            normPredictions[i] = std::exp(transformedLogScore);
        }

        return {e, normPredictions};
    }

    Err predictClassifierScores(
		const QVector<QVector<KarnnNNTarget>> &inferenceKarnnVecs,
        QVector<FDRCLassifierNeuralNet> &fdrcLassifierNeuralNets,
        QVector<QVector<float>> *predictions
        ) {

        ERR_INIT

    	predictions->resize(inferenceKarnnVecs.size());
    	for (int i = 0; i < inferenceKarnnVecs.size(); i++) {
    		const QPair<Err, QVector<float>> result = predictClassiferScoresLogic(
    			inferenceKarnnVecs[i],
    			&fdrcLassifierNeuralNets[i]
    			);
    		e = result.first; ree;
    		(*predictions)[i] = result.second;
    	}

        // const auto applyWeightsBinder = std::bind(
        //     predictClassiferScoresLogic,
        //     std::placeholders::_1,
        //     fdrcLassifierNeuralNet
        //     );

        // QVector<QVector<KarnnNNTarget>> karnnNNTargetsNormVecsTranched;
        // e = ParallelUtils::trancheVectorForParallelizationInOrder(
        //     karnnNNTargetsNorm,
        //     ParallelUtils::numberOfAvailableSystemProcessors(),
        //     0,
        //     &karnnNNTargetsNormVecsTranched
        //     ); ree;
        //
        // QFuture<QPair<Err, QVector<float>>> future = QtConcurrent::mapped(
        //     karnnNNTargetsNormVecsTranched,
        //     applyWeightsBinder
        //     );
        // future.waitForFinished();
        //
        // for (const QPair<Err, QVector<float>> &pr : future) {
        //     e = pr.first; ree;
        //     predictions->append(pr.second);
        // }

        ERR_RETURN
    }

    Err processPredictions(
        const QVector<float> &predictions,
        QVector<KarnnNNTarget> *karnnNNTargetsNorm,
        int foldIndex
    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(karnnNNTargetsNorm->isEmpty()); ree;

        for (int i = 0; i < predictions.size(); i++) {
            (*karnnNNTargetsNorm)[i].candidateScores->classifierScore = predictions.at(i);
            (*karnnNNTargetsNorm)[i].candidateScores->classifierFold = foldIndex;
        }

        ERR_RETURN
    }

    Err serialFilterByValue(
        double threshold,
        QVector<CandidateScores*> *candidateScoreses
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoreses->isEmpty()); ree;

        PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersClassifierScoreAsc(candidateScoreses);

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

	void logNeuralNetStats(const QVector<KarnnNNTarget> &karnnNNTargetsNorm) {
    	const int totalCount = karnnNNTargetsNorm.size();
    	const int decoyCount = static_cast<int>(std::count_if(
				karnnNNTargetsNorm.begin(),
				karnnNNTargetsNorm.end(),
				[](const KarnnNNTarget &kt){return kt.candidateScores->isDecoy;}
				));

    	qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
				 << "target vs decoy count for neural net training"
				 << totalCount - decoyCount << ":" << decoyCount
				 << "total" << totalCount;

// #define WRITE_KARNNN_NORM_TO_FILE
#ifdef WRITE_KARNNN_NORM_TO_FILE
    	e = updateProteinGroupAnnotation(
			m_fastaUri,
			0,
			&candidateScoresTargetsAndDecoysNeuralNet
			); ree;

    	const QString filenameNN = "kareNN_" + QString::number(m_pythiaParameters.threadCount) + ".dat";
    	QFile file(filenameNN);
    	if (!file.open(QIODevice::WriteOnly)) {
    		rrr(eFileError);
    	}

    	QDataStream out(&file);
    	out << karnnNNTargetsNormTrain.size();

    	for (const KarnnNNTarget &kt : karnnNNTargetsNormTrain) {
    		float label =kt.candidateScores->isDecoy ? 1.0 : 0.0;
    		out << kt.scoreVecNormalized;
    		out << label;
    		out << kt.candidateScores->proteinGroup;
    	}
    	file.close();
#endif
    }

}//namespace
Err PythiaDIAFFWorkflow::applyNeuralNetClassifier(
        const QVector<CandidateScores*> &candidateScoresTargetsAndDecoys,
        int seed,
        QVector<CandidateScores*> *candidateScoreClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_candidateScorePairs); ree;
    e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoys); ree;

    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> targetDecoyCandidateScorePairsPntrs;
    for (QPair<CandidateScoresTarget, CandidateScoresDecoy> &pr : m_candidateScorePairs) {
        targetDecoyCandidateScorePairsPntrs.push_back({&pr.first, &pr.second});
    }

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

// #define WRITE_CANDIDATES_TO_FILE
#ifdef WRITE_CANDIDATES_TO_FILE
    QVector<CandidateScoresReaderRow> candidateScoresReaderRows;
    for (const CandidateScores *cs : candidateScoresTargetsAndDecoysNeuralNet) {
        CandidateScoresReaderRow csrr = CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);
        candidateScoresReaderRows.push_back(csrr);
    }
    const QString filename = "candidates_" + QString::number(m_pythiaParameters.threadCount) + ".prq";
    e = ParquetReader::write(candidateScoresReaderRows, filename); ree;
#endif

    candidateScoreClassifier->clear();

    QVector<KarnnNNTarget> karnnNNTargetsNorm;
    e = buildKarnnNNTargetsNormalized(
        m_neuralNetFeatures,
        candidateScoresTargetsAndDecoysNeuralNet,
        &karnnNNTargetsNorm
        ); ree;
	logNeuralNetStats(karnnNNTargetsNorm);

	QVector<QVector<KarnnNNTarget>> karnnNNTargetsNormTranched;
	e = ParallelUtils::trancheVectorForParallelization(
		karnnNNTargetsNorm,
		m_pythiaParameters.baggingSize,
		&karnnNNTargetsNormTranched
		); ree;

    QVector<FDRCLassifierNeuralNet> fdrClassifierNeuralNets;
	QVector<QVector<KarnnNNTarget>> inferenceKarnnVecs;
    e = trainNeuralNetwork(
            m_pythiaParameters,
            karnnNNTargetsNormTranched,
            seed,
            &fdrClassifierNeuralNets,
            &inferenceKarnnVecs
            ); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Inference start";
    QVector<QVector<float>> predictions;
    e = predictClassifierScores(
        inferenceKarnnVecs,
        fdrClassifierNeuralNets,
        &predictions
        ); ree;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Inference end";

	for (int i = 0; i < predictions.size(); i++) {
	    processPredictions(predictions[i], &inferenceKarnnVecs[i], i); ree;
	}

    *candidateScoreClassifier = candidateScoresTargetsAndDecoysNeuralNet;

// #define WRITE_RESULTS_TO_FILE
#ifdef WRITE_RESULTS_TO_FILE
    QVector<CandidateScoresReaderRow> candidateScoresReaderRowsResults;
    for (const CandidateScores *cs : candidateScoresTargetsAndDecoysNeuralNet) {
        CandidateScoresReaderRow csrr = CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);
        candidateScoresReaderRowsResults.push_back(csrr);
    }
    const QString filenameResults = "results_" + QString::number(m_pythiaParameters.threadCount) + ".prq";
    e = ParquetReader::write(candidateScoresReaderRowsResults, filenameResults); ree;
#endif

    e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::NNClassifierScore,
            &targetDecoyCandidateScorePairsPntrs
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
        int targetCountBelowFDRThresholdOnePercent,
        QVector<CandidateScores*> *candidateScores
        ) {

    ERR_INIT

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Annotating" << candidateScores->size() << "PSMs";

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    FastaReader fastaReader;
    e = fastaReader.parseFastaFile(fastaFilePath); ree;

    const QVector<FastaEntry> fastaEntries = fastaReader.fastaEntries().values().toVector();

    QVector<QVector<FastaEntry>> fastaEntriesTranched;
    e = ParallelUtils::trancheVectorForParallelization(
        fastaEntries,
        m_pythiaParameters.threadCount,
        &fastaEntriesTranched
        ); ree;

    SequenceSubstringSearchomatic sequenceSubstringSearchomatic;
    e = sequenceSubstringSearchomatic.init(fastaEntries); ree;

    QStringList searchSequences;
    std::transform(
        candidateScores->begin(),
        candidateScores->end(),
        std::back_inserter(searchSequences),
        [](CandidateScores *cs){return cs->targetDecoyCandidatePair->peptideString();}
        );

    QVector<QString> proteinGroupsAll;
    e = sequenceSubstringSearchomatic.findProteinGroups(searchSequences, &proteinGroupsAll); ree;
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

//TODO delete this after dev is done.
#define REPORT_ENTRAP
#ifdef REPORT_ENTRAP

    QVector<CandidateScores*> candidateScoresCopy = *candidateScores;

    int counter = 0;
    int decoys = 0;
    int entrap = 0;
    int fmrCorrectCount = 0;

    PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(&candidateScoresCopy);

    for (CandidateScores *cs : candidateScoresCopy) {
        counter++;

        if (decoys/static_cast<double>(counter) >= 0.01) {
            break;
        }

        if (cs->isDecoy) {
            decoys++;
            continue;
        }

        if (!cs->proteinGroup.contains("_HUMAN")
            && !cs->isDecoy
            && cs->proteinGroup.contains("_ARATH")
            && !cs->proteinGroup.isEmpty()
            ) {
            entrap++;
            }

        if (entrap / static_cast<double>(counter) < 0.01) {
            fmrCorrectCount++;
        }
    }
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "By Discriminant Score:"
            << "Alt:" << targetCountBelowFDRThresholdOnePercent
            << "| Counter:" << counter
            << "| Decoys:" <<  decoys
            << "| Entrap:" << entrap
            << "| Counter FMR Corrected:" << fmrCorrectCount
            << "| Entrap%:" << entrap / static_cast<double>(counter);

    m_resultsSummary.psmCountLDA = counter;
    m_resultsSummary.decoyCountLDA = decoys;
    m_resultsSummary.entrapCountLDA = entrap;
    m_resultsSummary.psmCountCorrectedFMRLDA = fmrCorrectCount;
    m_resultsSummary.eFDRLDA = entrap / static_cast<double>(counter);

    PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersClassifierScoreAsc(&candidateScoresCopy);

    counter = 0;
    decoys = 0;
    entrap = 0;
    fmrCorrectCount = 0;

    for (CandidateScores *cs : candidateScoresCopy) {
        counter++;

        if (decoys/static_cast<double>(counter) >= 0.01) {
            break;
        }

        if (cs->isDecoy) {
            decoys++;
            continue;
        }

        if (!cs->proteinGroup.contains("_HUMAN")
            && !cs->isDecoy
            && cs->proteinGroup.contains("_ARATH")
            && !cs->proteinGroup.isEmpty()
            ) {
            entrap++;
            }

        if (entrap / static_cast<double>(counter) < 0.01) {
            fmrCorrectCount++;
        }
    }
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "By Classifier Score:"
            << "Alt:" << targetCountBelowFDRThresholdOnePercent
            << "| Counter:" << counter
            << "| Decoys:" <<  decoys
            << "| Entrap:" << entrap
            << "| Counter FMR Corrected:" << fmrCorrectCount
            << "| Entrap%:" << entrap / static_cast<double>(counter);

    m_resultsSummary.altPSMCountNeuralNet = targetCountBelowFDRThresholdOnePercent;
    m_resultsSummary.psmCountNeuralNet = counter;
    m_resultsSummary.decoyCountNeuralNet = decoys;
    m_resultsSummary.entrapCountNeuralNet = entrap;
    m_resultsSummary.psmCountCorrectedFMRNeuralNet = fmrCorrectCount;
    m_resultsSummary.eFDRNeuralNet = entrap / static_cast<double>(counter);
#endif

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Annotation finished";

    ERR_RETURN
}

void PythiaDIAFFWorkflow::filterDecoysOrNot(QVector<CandidateScores*> *candidateScoreClassifierPntrs) const {

    if (!m_pythiaParameters.reportDecoys) {

        const auto terminatorLogicFDRFilter
            = [&](CandidateScores *cs){return cs->qValue > m_pythiaParameters.percentFDR / 100.0;};

        const auto terminator = std::remove_if(
                candidateScoreClassifierPntrs->begin(),
                candidateScoreClassifierPntrs->end(),
                terminatorLogicFDRFilter
                );

        candidateScoreClassifierPntrs->erase(terminator, candidateScoreClassifierPntrs->end());
    }
    else {

        int counter = 0;
        int decoyCounter = 0;
        for (const CandidateScores *candidateScores : *candidateScoreClassifierPntrs) {

            counter++;
            if (candidateScores->isDecoy) {
                decoyCounter++;
            }

            if (
                constexpr double fdrThreshold = 0.5;
                decoyCounter / static_cast<double>(counter) >= fdrThreshold
                ) {
                break;
            }
        }

        candidateScoreClassifierPntrs->resize(counter);
    }
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
