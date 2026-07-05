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
#include "IdLevelQValueAnnotator.h"
#include "IonMobilitron.h"
#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.h"
#include "MsReaderPointerAcc.h"
#include "PythiaDIAFFWorkflowAlgos/OptimizeMassAccuracyPPMSettertron.h"
#include "ParallelUtils.h"
#include "PeptideStringWithMods.h"
#include "QuanTransitionRefinertron.h"
#include "QValueSettertron.h"
#include "SequenceSubstringSearchomatic.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

#include <QElapsedTimer>
#include <QFile>
#include <QFileInfo>
#include <QSet>
#include <QTextStream>

#include <algorithm>

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
    // m_pythiaParameters.ms2ExtractionWidthPPMOverride = 10;
    // m_pythiaParameters.peakCenter = 4;
    // m_pythiaParameters.writePythiaDIA = false;
    m_pythiaParameters.reannotate = true;
	m_pythiaParameters.baggingSize = 10;
	m_pythiaParameters.epochs = 12;
	// m_pythiaParameters.parallelNeuralNets = true;
		// m_pythiaParameters.nodesFraction = 0.5;
    // m_pythiaParameters.baggingSize = 4;
    // m_pythiaParameters.threadCount = 8;
    // m_pythiaParameters.shortReport = true;
    // m_pythiaParameters.calibrationTrainingVolume = 100000;
    // m_pythiaParameters.ionsSharedToReject = 4;
    // m_pythiaParameters.filterLengthIntegration = 7;
    // m_pythiaParameters.maxAnchorColumnIndex = 6;
    // m_pythiaParameters.minMs2FragCount = 2;
    // m_pythiaParameters.stopThresholdFractionMS2 = 0.65;
    // m_pythiaParameters.calibrationTrainingVolume = 1000;

    qDebug() << "ACHTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACHTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACHTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACHTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACHTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
    qDebug() << "ACHTUNG!!! TURN OVERRIDES OFF IN PRODUCTION!!!!";
#endif
/**************************/

    m_fragLibUri = fragLibUri;
    m_outputFolderPath = outputFolderPath;
    m_pythiaParameters.print();

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Reading library";

    e = FragLibReader::getFragLibReaderRows(
            m_fragLibUri,
            m_pythiaParameters.useAlternativeDecoys,
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

    constexpr int TIMS_NEURAL_NET_AUTO_INFERENCE_CANDIDATE_LIMIT = 200000;

    Err filterScoredCandidatesForNeuralNet(
            int minMs2FragCount,
            int neuralNetCandidateLimit,
            bool stratifyByTargetKey,
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
    	const int minTrainingVol = std::min(neuralNetCandidateLimit, candidateScoresTargetsAndDecoys->size());
        const int targetTrainingVol = std::max(fdrVsCount.value(fdrTrainingThresholdInt), minTrainingVol);

        if (stratifyByTargetKey && targetTrainingVol < candidateScoresTargetsAndDecoys->size()) {
            QMap<MzTargetKey, QVector<CandidateScores*>> targetKeyVsCandidateScores;
            for (CandidateScores *candidateScores : *candidateScoresTargetsAndDecoys) {
                if (candidateScores == nullptr) {
                    continue;
                }
                targetKeyVsCandidateScores[candidateScores->targetKey].push_back(candidateScores);
            }

            QVector<CandidateScores*> candidateScoresStratified;
            candidateScoresStratified.reserve(targetTrainingVol);
            bool appended = true;
            int rank = 0;
            while (candidateScoresStratified.size() < targetTrainingVol && appended) {
                appended = false;
                for (auto it = targetKeyVsCandidateScores.begin();
                     it != targetKeyVsCandidateScores.end() && candidateScoresStratified.size() < targetTrainingVol;
                     ++it) {
                    if (rank >= it.value().size()) {
                        continue;
                    }

                    candidateScoresStratified.push_back(it.value().at(rank));
                    appended = true;
                }
                ++rank;
            }

            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMS neural-net candidate pool stratified by target key"
                     << "initial_rows" << candidateScoresTargetsAndDecoys->size()
                     << "selected_rows" << candidateScoresStratified.size()
                     << "target_keys" << targetKeyVsCandidateScores.size()
                     << "target_rows" << targetTrainingVol;
            *candidateScoresTargetsAndDecoys = candidateScoresStratified;
        }
        else {
            candidateScoresTargetsAndDecoys->resize(targetTrainingVol);
        }

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

            QString resultsFilePath = msReaderPointerAcc->ptr->filePath() + S_GLOBAL_SETTINGS.DOT_RADIANT_DIA_FILE_EXTENSION;

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

        QString resultsFilePath = msReaderPointerAcc->ptr->filePath() + S_GLOBAL_SETTINGS.DOT_RADIANT_DIA_FILE_EXTENSION;

        if (!outputFolderPath.isEmpty()) {
            const QFileInfo fileInfo(resultsFilePath);
            const QString msDataFileName = fileInfo.fileName();
            resultsFilePath = outputFolderPath + msDataFileName;
        }
        e = ParquetReader::write(candidateScoreReaderRows, resultsFilePath); ree;

        ERR_RETURN
    }

    bool isCandidateScoresWritable(const CandidateScores *candidateScores) {
        if (candidateScores == nullptr || candidateScores->targetDecoyCandidatePair == nullptr) {
            return false;
        }

        constexpr int requiredOutputIonCount = 12;
        return candidateScores->integrations.size() >= requiredOutputIonCount
               && candidateScores->ionLabels.size() >= requiredOutputIonCount;
    }

    int removeNonWritableCandidateScores(QVector<CandidateScores*> *candidateScoresPntrs) {
        if (candidateScoresPntrs == nullptr || candidateScoresPntrs->isEmpty()) {
            return 0;
        }

        const int originalSize = candidateScoresPntrs->size();
        const auto terminator = std::remove_if(
            candidateScoresPntrs->begin(),
            candidateScoresPntrs->end(),
            [](const CandidateScores *candidateScores) {
                return !isCandidateScoresWritable(candidateScores);
            }
            );
        candidateScoresPntrs->erase(terminator, candidateScoresPntrs->end());
        return originalSize - candidateScoresPntrs->size();
    }

    QString cleanTsvField(QString value) {
        value.replace('\t', ' ');
        value.replace('\n', ' ');
        value.replace('\r', ' ');
        return value;
    }

    Err writeCandidateScoresDebug(
        const QVector<CandidateScores*> &candidateScoresPntrs,
        const QString &suffix,
        const QString &outputFolderPath,
        const MsReaderPointerAcc *msReaderPointerAcc
        ) {

        ERR_INIT

        if (candidateScoresPntrs.isEmpty()
            || msReaderPointerAcc == nullptr
            || msReaderPointerAcc->ptr.isNull()) {
            ERR_RETURN
        }

        QString resultsFilePath = msReaderPointerAcc->ptr->filePath()
                                  + S_GLOBAL_SETTINGS.DOT_RADIANT_DIA_FILE_EXTENSION
                                  + suffix;

        if (!outputFolderPath.isEmpty()) {
            const QFileInfo fileInfo(resultsFilePath);
            resultsFilePath = outputFolderPath + fileInfo.fileName();
        }

        QFile file(resultsFilePath);
        if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
            rrr(eFileError);
        }

        QTextStream out(&file);
        out << "PairPointer\tPairKey\tCandidateRowIsDecoy\tLibraryEntryIsDecoy\t"
            << "PeptideStringWithMods\tCharge\tTargetKey\tIsDecoy\tClassifierScore\tDiscriminantScore\tQValue\t"
            << "ClassifierFold\t"
            << "CosineSimSum100\tCosineSimSum45\tCosineSimSpectrumOverTimeCubed\tCosineSim100MS1\t"
            << "IonMobilityLibrary\tIonMobilityFound\tIonMobilityIndex\tIonMobilityIndexStart\tIonMobilityIndexEnd\t"
            << "IonMobilityDeltaAbs\tMs2IonMobilityWeightedDeltaAbs\tMs2IonMobilityApexDeltaAbsMean\t"
            << "Ms2IonMobilityMatchedIonFraction\tMs2IonMobilityFwhmMean\t"
            << "Ms2IonMobilityRtCosineMean\tMs2IonMobilityRtCosineStDev\t"
            << "Ms2IonMobilityRtApexAgreementFraction\tScanTime\tScanTimeDeltaAbs\t"
            << "TotalIntensityLog\tProteinGroup";
        for (int i = 1; i <= 12; ++i) {
            out << "\tMzSearched" << i;
        }
        out << '\n';

        int writtenCount = 0;
        for (const CandidateScores *cs : candidateScoresPntrs) {
            if (cs == nullptr
                || cs->targetDecoyCandidatePair == nullptr) {
                continue;
            }

            const QString pairPointer = QString::number(
                reinterpret_cast<quintptr>(cs->targetDecoyCandidatePair),
                16
                );
            const QString pairKey = cs->targetDecoyCandidatePair->peptideStringWithMods()
                                    + QStringLiteral("|z")
                                    + QString::number(cs->targetDecoyCandidatePair->charge())
                                    + QStringLiteral("|")
                                    + cs->targetKey
                                    + QStringLiteral("|")
                                    + pairPointer;

            out << pairPointer << '\t'
                << cleanTsvField(pairKey) << '\t'
                << cs->isDecoy << '\t'
                << cs->targetDecoyCandidatePair->isDecoy() << '\t'
                << cleanTsvField(cs->targetDecoyCandidatePair->peptideStringWithMods()) << '\t'
                << cs->targetDecoyCandidatePair->charge() << '\t'
                << cleanTsvField(cs->targetKey) << '\t'
                << (cs->isDecoy || cs->targetDecoyCandidatePair->isDecoy()) << '\t'
                << cs->classifierScore << '\t'
                << cs->discriminantScore << '\t'
                << cs->qValue << '\t'
                << cs->classifierFold << '\t'
                << cs->featuresArray[CosineSimSum100] << '\t'
                << cs->featuresArray[CosineSimSum45] << '\t'
                << cs->featuresArray[CosineSimSpectrumOverTimeCubed] << '\t'
                << cs->featuresArray[CosineSim100MS1] << '\t'
                << cs->targetDecoyCandidatePair->iIM() << '\t'
                << cs->imDriftTime << '\t'
                << cs->ionMobilityIndex << '\t'
                << cs->ionMobilityIndexStart << '\t'
                << cs->ionMobilityIndexEnd << '\t'
                << cs->featuresArray[IonMobilityDeltaAbs] << '\t'
                << cs->featuresArray[Ms2IonMobilityWeightedDeltaAbs] << '\t'
                << cs->featuresArray[Ms2IonMobilityApexDeltaAbsMean] << '\t'
                << cs->featuresArray[Ms2IonMobilityMatchedIonFraction] << '\t'
                << cs->featuresArray[Ms2IonMobilityFwhmMean] << '\t'
                << cs->featuresArray[Ms2IonMobilityRtCosineMean] << '\t'
                << cs->featuresArray[Ms2IonMobilityRtCosineStDev] << '\t'
                << cs->featuresArray[Ms2IonMobilityRtApexAgreementFraction] << '\t'
                << cs->scanTime << '\t'
                << cs->featuresArray[ScanTimeDeltaAbs] << '\t'
                << cs->featuresArray[TotalIntensityLog] << '\t'
                << cleanTsvField(cs->proteinGroup);
            const QVector<MS2Ion> &ms2Ions = cs->isDecoy
                ? cs->targetDecoyCandidatePair->ms2IonsDecoy()
                : cs->targetDecoyCandidatePair->ms2IonsTarget();
            for (int i = 0; i < 12; ++i) {
                out << '\t' << (i < ms2Ions.size() ? ms2Ions.at(i).mz : -1.0);
            }
            out << '\n';
            ++writtenCount;
        }

        file.close();

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "Wrote candidate debug table"
                 << resultsFilePath
                 << writtenCount
                 << "rows from"
                 << candidateScoresPntrs.size()
                 << "target/decoy rows";

        ERR_RETURN
    }

    bool isIonMobilityOnlyFeature(const Features feature) {
        switch (feature) {
        case Ms1IntensityFoundApex100IM:
        case IonMobilityDelta:
        case IonMobilityDeltaAbs:
        case IonMobilityPdAbs:
        case Ms2IonMobilityWeightedDelta:
        case Ms2IonMobilityWeightedDeltaAbs:
        case Ms2IonMobilityApexDeltaAbsMean:
        case Ms2IonMobilityApexDeltaAbsStDev:
        case Ms2IonMobilityMatchedIonFraction:
        case Ms2IonMobilityFwhmMean:
        case Ms2IonMobilityFwhmStDev:
        case Ms2IonMobilityRtCosineMean:
        case Ms2IonMobilityRtCosineStDev:
        case Ms2IonMobilityRtApexAgreementFraction:
            return true;
        default:
            return false;
        }
    }

    void removeIonMobilityOnlyFeatures(QVector<Features> *features) {
        const auto terminator = std::remove_if(
            features->begin(),
            features->end(),
            isIonMobilityOnlyFeature
            );
        features->erase(terminator, features->end());
    }

    void appendFeatureIfMissing(QVector<Features> *features, const Features feature) {
        if (!features->contains(feature)) {
            features->push_back(feature);
        }
    }

    void configureWorkflowFeaturesForReader(
        bool isTIMS,
        QVector<Features> *calibratomaticFeatures,
        QVector<Features> *ppmOptimizationFeatures,
        QVector<Features> *neuralNetFeatures
        ) {

        *calibratomaticFeatures = DiscriminantScoretron::featuresCalibration();
        *ppmOptimizationFeatures = DiscriminantScoretron::featuresOptimization();
        *neuralNetFeatures = DiscriminantScoretron::featuresNeuralNetwork();

        if (isTIMS) {
            appendFeatureIfMissing(ppmOptimizationFeatures, Ms2IonMobilityRtCosineMean);
            appendFeatureIfMissing(ppmOptimizationFeatures, Ms2IonMobilityRtCosineStDev);
            appendFeatureIfMissing(ppmOptimizationFeatures, Ms2IonMobilityRtApexAgreementFraction);

            appendFeatureIfMissing(neuralNetFeatures, Ms2IonMobilityRtCosineMean);
            appendFeatureIfMissing(neuralNetFeatures, Ms2IonMobilityRtCosineStDev);
            appendFeatureIfMissing(neuralNetFeatures, Ms2IonMobilityRtApexAgreementFraction);
            return;
        }

        removeIonMobilityOnlyFeatures(calibratomaticFeatures);
        removeIonMobilityOnlyFeatures(ppmOptimizationFeatures);
        removeIonMobilityOnlyFeatures(neuralNetFeatures);
    }

}//namespace
Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;
    MsReaderPointerAcc msReaderPointerAcc;

    msReaderPointerAcc.setUseLazyLoading(m_pythiaParameters.useLazyLoading);

    const bool hasAnalysisScanTimeRange
        = m_pythiaParameters.analysisScanTimeMin >= 0.0
          && m_pythiaParameters.analysisScanTimeMax > m_pythiaParameters.analysisScanTimeMin;
    const QFileInfo msDataFileInfo(msDataFilePath);
    const bool useBrukerEarlyScanTimeFilter
        = hasAnalysisScanTimeRange
          && msDataFileInfo.isDir()
          && msDataFileInfo.suffix().compare(
              S_GLOBAL_SETTINGS.BRUKER_FILE_EXTENSION,
              Qt::CaseInsensitive
              ) == 0;

    if (useBrukerEarlyScanTimeFilter) {
        e = msReaderPointerAcc.openFile(
            msDataFilePath,
            QStringLiteral("scanTime"),
            {
                m_pythiaParameters.analysisScanTimeMin,
                m_pythiaParameters.analysisScanTimeMax
            }
            ); ree;
    }
    else {
        e = msReaderPointerAcc.openFile(msDataFilePath); ree;
    }

    if (hasAnalysisScanTimeRange && !useBrukerEarlyScanTimeFilter) {
        e = msReaderPointerAcc.ptr->restrictScanTimeRange(
            static_cast<ScanTime>(m_pythiaParameters.analysisScanTimeMin),
            static_cast<ScanTime>(m_pythiaParameters.analysisScanTimeMax)
            ); ree;
    }
    msReaderPointerAcc.ptr->printSize();

    configureWorkflowFeaturesForReader(
        msReaderPointerAcc.ptr->isTIMS(),
        &m_calibratomaticFeatures,
        &m_ppmOptimizationFeatures,
        &m_neuralNetFeatures
        );

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
    bool usedDiscriminantFallback = false;
    e = applyNeuralNetClassifier(
            candidateScoresTargetsAndDecoys,
            &msReaderPointerAcc,
            S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST,
            &candidateScoreClassifierPntrs,
            &usedDiscriminantFallback
            ); ree;

    const int removedNonWritableCandidateScores = removeNonWritableCandidateScores(&candidateScoreClassifierPntrs);
    if (removedNonWritableCandidateScores > 0) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "Removed non-writable candidate score rows"
                 << removedNonWritableCandidateScores;
    }

    if (candidateScoreClassifierPntrs.isEmpty()) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "No writable candidate score rows after neural-net filtering; skipping result write";
        ERR_RETURN
    }

    int targetCountBelowFDRThresholdOnePercent;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
        candidateScoreClassifierPntrs,
        m_pythiaParameters.percentFDR / 100.0,
        &targetCountBelowFDRThresholdOnePercent
        ); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Pre Neural Net PSMs Count" << targetCountBelowFDRThreshold
             << "| Post Neural Net Count PSMs" << targetCountBelowFDRThresholdOnePercent;

    const bool candidateScoresSortedHiLo = usedDiscriminantFallback
        ? std::is_sorted(
            candidateScoreClassifierPntrs.begin(),
            candidateScoreClassifierPntrs.end(),
            [](const CandidateScores *l, const CandidateScores *r) {
                if (MathUtils::tSame(l->discriminantScore, r->discriminantScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                    if (MathUtils::tSame(l->featuresArray[CosineSimSum100], r->featuresArray[CosineSimSum100], S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                        if (MathUtils::tSame(l->featuresArray[CosineSimSpectrumOverTime], r->featuresArray[CosineSimSpectrumOverTime], S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                            return l->isDecoy > r->isDecoy;
                        }
                        return l->featuresArray[CosineSimSpectrumOverTime] > r->featuresArray[CosineSimSpectrumOverTime];
                    }
                    return l->featuresArray[CosineSimSum100] > r->featuresArray[CosineSimSum100];
                }
                return l->discriminantScore > r->discriminantScore;
            })
        : std::is_sorted(
            candidateScoreClassifierPntrs.begin(),
            candidateScoreClassifierPntrs.end(),
            [](const CandidateScores *l, const CandidateScores *r) {
                if (MathUtils::tSame(l->classifierScore, r->classifierScore, S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL)) {
                    return l->discriminantScore > r->discriminantScore;
                }
                return l->classifierScore < r->classifierScore;
            });
    e = ErrorUtils::isTrue(candidateScoresSortedHiLo); ree;

    if (m_pythiaParameters.reannotate) {
        e = updateProteinGroupAnnotation(
                m_fastaUri,
                targetCountBelowFDRThresholdOnePercent,
                &candidateScoreClassifierPntrs
                ); ree;
    }

    const bool useLocalRtIdLevelQValues = !msReaderPointerAcc.ptr.isNull() && msReaderPointerAcc.ptr->isTIMS();
    e = IdLevelQValueAnnotator::annotate(
        &candidateScoreClassifierPntrs,
        !usedDiscriminantFallback,
        useLocalRtIdLevelQValues ? m_pythiaParameters.timsLocalFdrRtBinSeconds : 0.0
        ); ree;

    filterDecoysOrNot(&candidateScoreClassifierPntrs);

    if (m_pythiaParameters.writeRadiantDIA) {
        e = writePythiaDIA(
            candidateScoreClassifierPntrs,
            m_pythiaParameters.shortReport,
            m_outputFolderPath,
            &msReaderPointerAcc
            ); ree;
    }

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

Err PythiaDIAFFWorkflow::rescoreTimsFilteredCandidatesForNeuralNet(
    const MsReaderPointerAcc *msReaderPointerAcc,
    QVector<CandidateScores*> *candidateScoresTargetsAndDecoysNeuralNet,
    QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> *rescoredCandidateScorePairs,
    QVector<Features> *neuralNetFeatures
    ) {

    ERR_INIT

    if (msReaderPointerAcc == nullptr
        || msReaderPointerAcc->ptr.isNull()
        || !msReaderPointerAcc->ptr->isTIMS()
        || candidateScoresTargetsAndDecoysNeuralNet == nullptr
        || candidateScoresTargetsAndDecoysNeuralNet->isEmpty()
        || rescoredCandidateScorePairs == nullptr
        || neuralNetFeatures == nullptr) {
        ERR_RETURN
    }

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> targetKeyVsTargetDecoyCandidatePairs;
    QSet<TargetDecoyCandidatePair*> seenTargetDecoyCandidatePairs;

    const int maxTimsSecondStageCandidateRows = m_pythiaParameters.timsSecondStageCandidateRowLimit;
    const int maxTimsSecondStageUniquePrecursors = m_pythiaParameters.timsSecondStageUniquePrecursorLimit;
    int inspectedCandidateRows = 0;
    for (const CandidateScores *candidateScores : *candidateScoresTargetsAndDecoysNeuralNet) {
        if (inspectedCandidateRows >= maxTimsSecondStageCandidateRows
            || seenTargetDecoyCandidatePairs.size() >= maxTimsSecondStageUniquePrecursors) {
            break;
        }

        ++inspectedCandidateRows;
        if (candidateScores == nullptr || candidateScores->targetDecoyCandidatePair == nullptr) {
            continue;
        }

        TargetDecoyCandidatePair *targetDecoyCandidatePair = candidateScores->targetDecoyCandidatePair;
        if (seenTargetDecoyCandidatePairs.contains(targetDecoyCandidatePair)) {
            continue;
        }

        seenTargetDecoyCandidatePairs.insert(targetDecoyCandidatePair);
        targetKeyVsTargetDecoyCandidatePairs[candidateScores->targetKey].push_back(targetDecoyCandidatePair);
    }

    if (targetKeyVsTargetDecoyCandidatePairs.isEmpty()) {
        ERR_RETURN
    }

    QVector<Features> rescoringFeatures = m_ppmOptimizationFeatures;
    appendFeatureIfMissing(&rescoringFeatures, Ms2IonMobilityRtCosineMean);
    appendFeatureIfMissing(&rescoringFeatures, Ms2IonMobilityRtCosineStDev);
    appendFeatureIfMissing(&rescoringFeatures, Ms2IonMobilityRtApexAgreementFraction);

    appendFeatureIfMissing(neuralNetFeatures, Ms2IonMobilityRtCosineMean);
    appendFeatureIfMissing(neuralNetFeatures, Ms2IonMobilityRtCosineStDev);
    appendFeatureIfMissing(neuralNetFeatures, Ms2IonMobilityRtApexAgreementFraction);

    const QVector<float> rescoringWeights = DiscriminantScoretron::defaultWeights(rescoringFeatures);
    constexpr int topNMs2IonsTimsSecondStage = 8;
    constexpr bool useTopNIntegrationsParameter = false;
    constexpr float minPeakCountTims = 2.9f;
    const int threadCount = std::max(1, std::min(m_pythiaParameters.threadCount, 8));

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "TIMS 4D second-stage rescoring"
             << "unique_precursors" << seenTargetDecoyCandidatePairs.size()
             << "selected_candidate_rows" << inspectedCandidateRows
             << "candidate_rows" << candidateScoresTargetsAndDecoysNeuralNet->size()
             << "candidate_row_limit" << maxTimsSecondStageCandidateRows
             << "unique_precursor_limit" << maxTimsSecondStageUniquePrecursors
             << "target_keys" << targetKeyVsTargetDecoyCandidatePairs.size();

    QElapsedTimer timer;
    timer.start();
    rescoredCandidateScorePairs->clear();
    m_targetDecoyCandidatePairScoretron.setUseAdaptiveTimsMobilityCentering(true);
    e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
        rescoringFeatures,
        topNMs2IonsTimsSecondStage,
        m_msCalibratomatic,
        minPeakCountTims,
        threadCount,
        useTopNIntegrationsParameter,
        QMap<MzTargetKey, TurboXIC*>(),
        rescoringWeights,
        &targetKeyVsTargetDecoyCandidatePairs,
        rescoredCandidateScorePairs
        );
    m_targetDecoyCandidatePairScoretron.setUseAdaptiveTimsMobilityCentering(false);
    ree;

    if (rescoredCandidateScorePairs->isEmpty()) {
        ERR_RETURN
    }

    QHash<TargetDecoyCandidatePair*, QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> rescoredScoresByCandidate;
    rescoredScoresByCandidate.reserve(rescoredCandidateScorePairs->size());
    for (QPair<CandidateScoresTarget, CandidateScoresDecoy> &rescoredPair : *rescoredCandidateScorePairs) {
        TargetDecoyCandidatePair *targetDecoyCandidatePair = rescoredPair.first.targetDecoyCandidatePair;
        if (targetDecoyCandidatePair == nullptr) {
            targetDecoyCandidatePair = rescoredPair.second.targetDecoyCandidatePair;
        }
        if (targetDecoyCandidatePair == nullptr) {
            continue;
        }
        rescoredScoresByCandidate.insert(targetDecoyCandidatePair, {&rescoredPair.first, &rescoredPair.second});
    }

    int replacedCandidateRows = 0;
    for (CandidateScores *&candidateScores : *candidateScoresTargetsAndDecoysNeuralNet) {
        if (candidateScores == nullptr || candidateScores->targetDecoyCandidatePair == nullptr) {
            continue;
        }

        const auto rescoredIt = rescoredScoresByCandidate.constFind(candidateScores->targetDecoyCandidatePair);
        if (rescoredIt == rescoredScoresByCandidate.constEnd()) {
            continue;
        }

        CandidateScores *replacementCandidateScores
            = candidateScores->isDecoy ? rescoredIt.value().second : rescoredIt.value().first;
        if (!isCandidateScoresWritable(replacementCandidateScores)) {
            continue;
        }

        candidateScores = replacementCandidateScores;
        ++replacedCandidateRows;
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "TIMS 4D second-stage rescored"
             << "pairs" << rescoredCandidateScorePairs->size()
             << "replaced_candidate_rows" << replacedCandidateRows
             << "candidate_rows" << candidateScoresTargetsAndDecoysNeuralNet->size()
             << "msec" << timer.elapsed();

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::mainAnalysis(
        const MsReaderPointerAcc *msReaderPointerAcc,
        int *targetCountBelowFDRThresholdOnePercent
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

	    m_candidateScorePairs.clear();
        m_timsSecondStageCandidateScorePairs.clear();

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

    const float minPeakCount = msReaderPointerAcc->ptr->isTIMS() ? 2.9f : 3.9f;
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
        &weights,
        msReaderPointerAcc->ptr->isTIMS()
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
    bool is_nn_decoy(const CandidateScores *candidate_scores) {
        return candidate_scores->isDecoy || candidate_scores->targetDecoyCandidatePair->isDecoy();
    }

    bool is_nn_decoy(const KarnnNNTarget &kt) {
        return is_nn_decoy(kt.candidateScores);
    }

    float get_nn_decoy_label(const KarnnNNTarget &kt) {
        // Note: for decoy-origin candidates, label both elements of the pair as decoys
        return is_nn_decoy(kt) ? 1.0f : 0.0f;
    }

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

    struct NeuralNetworkTrainingFoldInput {
        QPair<QVector<QVector<float>>, QVector<float>> trainingVecs;
        int seed = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST;
    };

	QPair<Err, FDRCLassifierNeuralNet> trainNeuralNetworkLogic(
        const NeuralNetworkTrainingFoldInput &trainingFoldInput,
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
				pythiaParameters.focalLossGamma,
				pythiaParameters.threadCount
		); rree;

    	e = fdrcLassifierNeuralNet.trainClassifier(
                trainingFoldInput.trainingVecs.first,
                trainingFoldInput.trainingVecs.second,
                trainingFoldInput.seed,
				pythiaParameters.verbosity
				); rree;

    	return {e, fdrcLassifierNeuralNet};
    }

    Err trainNeuralNetwork(
            const PythiaParameters &pythiaParameters,
            const QVector<QVector<KarnnNNTarget>> &karnnNNTargetsNormTranched,
            int seed,
            bool useFoldSpecificSeeds,
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

        QVector<NeuralNetworkTrainingFoldInput> trainingFoldInputs(karnnNNTargetsNormTranched.size());
    	inferenceKarnnVecs->resize(karnnNNTargetsNormTranched.size());

    	for (int i = 0; i < karnnNNTargetsNormTranched.size(); i++) {
            trainingFoldInputs[i].seed = useFoldSpecificSeeds
                                             ? seed + i + 1
                                             : S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST;
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
    				yData.push_back(get_nn_decoy_label(kt));
    			}

                trainingFoldInputs[i].trainingVecs.first.append(xData);
                trainingFoldInputs[i].trainingVecs.second.append(yData);
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
                trainingFoldInputs,
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
                    trainingFoldInputs[i],
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
        FDRCLassifierNeuralNet *fdrcLassifierNeuralNet,
        bool normalizePredictions,
        bool keepRawPredictionsWhenNoFdrCutoff
        ) {

        ERR_INIT

        QVector<float> predictions;

        QVector<QVector<float>> xData;
        xData.reserve(karnnNNTargetsNorm.size());
        QVector<float> yData;
        yData.reserve(karnnNNTargetsNorm.size());
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVecNormalized);
            yData.push_back(get_nn_decoy_label(kt));
        }

        e = fdrcLassifierNeuralNet->predictBaggedClassifiers(
            xData,
            &predictions
            ); rree;

        if (!normalizePredictions) {
            return {e, predictions};
        }

        // Normalize predictions to a standard range, using a linear transformation
        // that fixes a 1% FDR cutoff and the median of decoy scores (in log space).

        // 1. Compute 1% FDR cutoff

        // Create temporary candidate scores for FDR calculation using QValueSettertron
        QVector<CandidateScores> tempCandidateScores;
        tempCandidateScores.reserve(karnnNNTargetsNorm.size());

        for (int i = 0; i < karnnNNTargetsNorm.size(); ++i) {
            // Create a temporary copy of the candidate with updated classifier score
            CandidateScores tempCs = *karnnNNTargetsNorm[i].candidateScores;
            tempCs.classifierScore = predictions[i];
            tempCandidateScores.push_back(tempCs);
        }

        QVector<CandidateScores*> tempCandidateScoresPntrs;
        std::transform(
            tempCandidateScores.begin(),
            tempCandidateScores.end(),
            std::back_inserter(tempCandidateScoresPntrs),
            [](CandidateScores &cs){return &cs;}
        );

        e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::NNClassifierScore,
            &tempCandidateScoresPntrs
        ); rree;

        constexpr double targetFDR = 0.01;
        float fdrCutoff = 0.0f;

        for (const CandidateScores *cs : tempCandidateScoresPntrs) {
            if (!is_nn_decoy(cs) && cs->qValue <= targetFDR) {
                fdrCutoff = std::max(fdrCutoff, static_cast<float>(cs->classifierScore));
            }
        }

        if (fdrCutoff <= 0.0f) {
            if (keepRawPredictionsWhenNoFdrCutoff) {
                qWarning() << "No targets found at 1% FDR threshold, keeping raw neural-net predictions";
                return {e, predictions};
            }

            qWarning() << "No targets found at 1% FDR threshold, using minimum value";
            fdrCutoff = std::numeric_limits<float>::min();
        }

        const int decoyCount = std::count_if(
            tempCandidateScoresPntrs.begin(),
            tempCandidateScoresPntrs.end(),
            [](CandidateScores *cs){ return is_nn_decoy(cs); }
            );

        e = ErrorUtils::isTrue(decoyCount > 0); rree;

        {
            // Sanity check
            const int checkDecoyCount = std::count_if(
                yData.begin(),
                yData.end(),
                [](float y){ return y == 1.0f; }
                );

            e = ErrorUtils::isEqual(decoyCount, checkDecoyCount); rree;
        }

        QVector<float> decoyScores;
        decoyScores.reserve(decoyCount);
        for (int i = 0; i < predictions.size(); ++i) {
            if (yData[i] == 1.0f) {
                decoyScores.push_back(predictions[i]);
            }
        }

        e = ErrorUtils::isEqual(decoyScores.size(), decoyCount); rree;

        const double decoyMedian = MathUtils::median(decoyScores);
        e = ErrorUtils::isTrue(decoyMedian > 0.0); rree;

        const double logCutoff = std::log(fdrCutoff);
        const double logDecoyMedian = std::log(decoyMedian);
        const double denominator = logDecoyMedian - logCutoff;

        e = ErrorUtils::isFalse(MathUtils::tZero(denominator)); rree;

        const double m = 1.0 / denominator;
        const double b = -logCutoff / denominator;

        // 4. Apply linear transformation to predictions:
        //
        //    For each score:
        // a. Convert to log space
        // b. Apply linear transformation
        // c. Convert back from log space
        //
        // Future enhancement: use Eigen to vectorize this operation
        // and employ fused-multiply-add (FMA) intrinsic.

        QVector<float> normPredictions(predictions.size(), 0.0f);

        for (int i = 0; i < predictions.size(); ++i) {
            // Validate prediction value before log transformation
            float prediction = predictions[i];
            if (prediction <= 0.0f) {
                continue;
            }

            double logScore = std::log(prediction);
            double transformedLogScore = m * logScore + b;
            normPredictions[i] = std::exp(transformedLogScore);
        }

        return {e, normPredictions};
    }

    Err predictClassifierScores(
		const QVector<QVector<KarnnNNTarget>> &inferenceKarnnVecs,
        QVector<FDRCLassifierNeuralNet> &fdrcLassifierNeuralNets,
        bool normalizePredictions,
        bool keepRawPredictionsWhenNoFdrCutoff,
        QVector<QVector<float>> *predictions
        ) {

	    ERR_INIT

    	predictions->resize(inferenceKarnnVecs.size());
    	for (int i = 0; i < inferenceKarnnVecs.size(); i++) {
            if (inferenceKarnnVecs[i].isEmpty()) {
                (*predictions)[i].clear();
                continue;
            }

    		const QPair<Err, QVector<float>> result = predictClassiferScoresLogic(
    			inferenceKarnnVecs[i],
    			&fdrcLassifierNeuralNets[i],
                normalizePredictions,
                keepRawPredictionsWhenNoFdrCutoff
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

    Err predictClassifierScoresWithAllNets(
        const QVector<KarnnNNTarget> &inferenceKarnnVec,
        QVector<FDRCLassifierNeuralNet> &fdrcLassifierNeuralNets,
        bool normalizePredictions,
        bool keepRawPredictionsWhenNoFdrCutoff,
        QVector<float> *predictions
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(inferenceKarnnVec.isEmpty()); ree;
        e = ErrorUtils::isFalse(fdrcLassifierNeuralNets.isEmpty()); ree;

        QVector<double> predictionSums(inferenceKarnnVec.size(), 0.0);
        int modelCount = 0;
        for (FDRCLassifierNeuralNet &fdrcLassifierNeuralNet : fdrcLassifierNeuralNets) {
            const QPair<Err, QVector<float>> result = predictClassiferScoresLogic(
                inferenceKarnnVec,
                &fdrcLassifierNeuralNet,
                normalizePredictions,
                keepRawPredictionsWhenNoFdrCutoff
                );
            if (result.first != eNoError || result.second.size() != inferenceKarnnVec.size()) {
                continue;
            }

            for (int i = 0; i < result.second.size(); ++i) {
                predictionSums[i] += result.second.at(i);
            }
            ++modelCount;
        }

        predictions->clear();
        predictions->reserve(inferenceKarnnVec.size());
        if (modelCount <= 0) {
            predictions->fill(1.0f, inferenceKarnnVec.size());
            ERR_RETURN
        }

        for (double predictionSum : predictionSums) {
            predictions->push_back(static_cast<float>(predictionSum / modelCount));
        }

        ERR_RETURN
    }

    Err processPredictions(
        const QVector<float> &predictions,
        QVector<KarnnNNTarget> *karnnNNTargetsNorm,
        int foldIndex
    ) {

        ERR_INIT

        if (karnnNNTargetsNorm->isEmpty() && predictions.isEmpty()) {
            ERR_RETURN
        }
        e = ErrorUtils::isFalse(karnnNNTargetsNorm->isEmpty()); ree;
        e = ErrorUtils::isEqual(predictions.size(), karnnNNTargetsNorm->size()); ree;

        for (int i = 0; i < predictions.size(); i++) {
            (*karnnNNTargetsNorm)[i].candidateScores->classifierScore = predictions.at(i);
            (*karnnNNTargetsNorm)[i].candidateScores->classifierFold = foldIndex;
        }

        ERR_RETURN
    }

    Err serialFilterByValue(
        double threshold,
        QVector<CandidateScores*> *candidateScoreses,
        const bool useClassifierScore
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoreses->isEmpty()); ree;

        if (useClassifierScore) {
            PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersClassifierScoreAsc(candidateScoreses);
        }
        else {
            PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(candidateScoreses);
        }

        int counter = 0;
        for (const CandidateScores *csp : *candidateScoreses) {

            counter++;
            if (csp->qValue >= threshold && !csp->isDecoy && !csp->targetDecoyCandidatePair->isDecoy()) {
                break;
            }
        }

        candidateScoreses->resize(counter);

        ERR_RETURN
    }

    Err filterByQValueThreshold(
        double threshold,
        QVector<CandidateScores*> *candidateScoreses,
        const bool useClassifierScore
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScoreses->isEmpty()); ree;

        int writeIndex = 0;
        QSet<TargetDecoyCandidatePair*> acceptedCandidatePairs;
        acceptedCandidatePairs.reserve(candidateScoreses->size());
        for (CandidateScores *csp : *candidateScoreses) {
            if (csp == nullptr
                || csp->targetDecoyCandidatePair == nullptr
                || csp->qValue > threshold) {
                continue;
            }
            acceptedCandidatePairs.insert(csp->targetDecoyCandidatePair);
        }

        for (CandidateScores *csp : *candidateScoreses) {
            if (csp == nullptr
                || csp->targetDecoyCandidatePair == nullptr
                || !acceptedCandidatePairs.contains(csp->targetDecoyCandidatePair)) {
                continue;
            }
            (*candidateScoreses)[writeIndex] = csp;
            ++writeIndex;
        }
        candidateScoreses->resize(writeIndex);

        if (useClassifierScore) {
            PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersClassifierScoreAsc(candidateScoreses);
        }
        else {
            PythiaDIAFFWorkflowSharedMethods::sortCandidatePointersDiscScoreDesc(candidateScoreses);
        }

        ERR_RETURN
    }

    Err applyTimsHighEvidenceGlobalQValueFilter(
        QVector<CandidateScores*> *candidateScores,
        QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> *targetDecoyCandidateScorePairs,
        const QValueSettertron::QValueScoreType &qValueScoreType,
        const PythiaParameters &pythiaParameters
        ) {

        ERR_INIT

        if (candidateScores == nullptr
            || candidateScores->isEmpty()
            || targetDecoyCandidateScorePairs == nullptr
            || targetDecoyCandidateScorePairs->isEmpty()) {
            ERR_RETURN
        }

        const float minCosineSimSum100 = pythiaParameters.timsHighEvidenceMinCosineSimSum100;
        const float minCosineSimSpectrumOverTimeCubed = pythiaParameters.timsHighEvidenceMinCosineSimSpectrumOverTimeCubed;
        const float maxScanTimeDeltaAbs = pythiaParameters.timsHighEvidenceMaxScanTimeDeltaAbs;

        auto filteredPairsForThresholds = [targetDecoyCandidateScorePairs](
            const float minCosineSimSum100Threshold,
            const float minCosineSimSpectrumOverTimeCubedThreshold,
            const float maxScanTimeDeltaAbsThreshold
            ) {
            QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> selectedPairs;
            selectedPairs.reserve(targetDecoyCandidateScorePairs->size());
            for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pair : *targetDecoyCandidateScorePairs) {
                CandidateScoresTarget *targetScores = pair.first;
                CandidateScoresDecoy *decoyScores = pair.second;
                if (targetScores == nullptr || decoyScores == nullptr) {
                    continue;
                }

                if (targetScores->featuresArray[CosineSimSum100] < minCosineSimSum100Threshold
                    || targetScores->featuresArray[CosineSimSpectrumOverTimeCubed] < minCosineSimSpectrumOverTimeCubedThreshold
                    || targetScores->featuresArray[ScanTimeDeltaAbs] > maxScanTimeDeltaAbsThreshold) {
                    continue;
                }

                selectedPairs.push_back(pair);
            }

            return selectedPairs;
        };

        if (pythiaParameters.timsHighEvidenceFilterSweep) {
            struct SweepThresholds {
                float minCosineSimSum100;
                float minCosineSimSpectrumOverTimeCubed;
                float maxScanTimeDeltaAbs;
            };

            QVector<SweepThresholds> sweepThresholds = {
                {3.6f, 0.15f, 120.0f},
                {3.8f, 0.20f, 120.0f},
                {4.0f, 0.20f, 90.0f},
                {4.0f, 0.25f, 90.0f},
                {4.2f, 0.30f, 70.0f},
                {4.4f, 0.35f, 70.0f},
                {minCosineSimSum100, minCosineSimSpectrumOverTimeCubed, maxScanTimeDeltaAbs},
            };

            for (const SweepThresholds &thresholds : sweepThresholds) {
                QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> sweepPairs = filteredPairsForThresholds(
                    thresholds.minCosineSimSum100,
                    thresholds.minCosineSimSpectrumOverTimeCubed,
                    thresholds.maxScanTimeDeltaAbs
                );

                if (sweepPairs.isEmpty()) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMS high-evidence threshold sweep"
                             << "pairs" << 0
                             << "q01_targets" << 0
                             << "q01_decoys" << 0
                             << "min_cos_sum100" << thresholds.minCosineSimSum100
                             << "min_cos_time_cubed" << thresholds.minCosineSimSpectrumOverTimeCubed
                             << "max_scan_time_delta_abs" << thresholds.maxScanTimeDeltaAbs;
                    continue;
                }

                e = QValueSettertron::setQValueForCandidates(
                    qValueScoreType,
                    &sweepPairs,
                    true,
                    true,
                    pythiaParameters.timsLocalFdrRtBinSeconds
                    ); ree;

                const int q01Targets = static_cast<int>(std::count_if(
                    sweepPairs.begin(),
                    sweepPairs.end(),
                    [](const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pair) {
                        return pair.first != nullptr
                               && !is_nn_decoy(pair.first)
                               && pair.first->qValue <= 0.01;
                    }));

                int q01Decoys = 0;
                for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pair : sweepPairs) {
                    if (pair.first != nullptr
                        && is_nn_decoy(pair.first)
                        && pair.first->qValue <= 0.01) {
                        ++q01Decoys;
                    }
                    if (pair.second != nullptr
                        && is_nn_decoy(pair.second)
                        && pair.second->qValue <= 0.01) {
                        ++q01Decoys;
                    }
                }

                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMS high-evidence threshold sweep"
                         << "pairs" << sweepPairs.size()
                         << "q01_targets" << q01Targets
                         << "q01_decoys" << q01Decoys
                         << "min_cos_sum100" << thresholds.minCosineSimSum100
                         << "min_cos_time_cubed" << thresholds.minCosineSimSpectrumOverTimeCubed
                         << "max_scan_time_delta_abs" << thresholds.maxScanTimeDeltaAbs;
            }
        }

        QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> filteredPairs = filteredPairsForThresholds(
            minCosineSimSum100,
            minCosineSimSpectrumOverTimeCubed,
            maxScanTimeDeltaAbs
            );

        if (filteredPairs.isEmpty()) {
            ERR_RETURN
        }

        e = QValueSettertron::setQValueForCandidates(
            qValueScoreType,
            &filteredPairs,
            true,
            true,
            pythiaParameters.timsLocalFdrRtBinSeconds
            ); ree;

        const int originalRows = candidateScores->size();
        candidateScores->clear();
        candidateScores->reserve(filteredPairs.size() * 2);
        for (const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pair : filteredPairs) {
            candidateScores->push_back(pair.first);
            candidateScores->push_back(pair.second);
        }

        const int targetRows = static_cast<int>(std::count_if(
            candidateScores->begin(),
            candidateScores->end(),
            [](const CandidateScores *candidateScore) {
                return candidateScore != nullptr && !is_nn_decoy(candidateScore);
            }));
        const int decoyRows = candidateScores->size() - targetRows;

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMS high-evidence global q-value filter"
                 << "rows" << originalRows << "->" << candidateScores->size()
                 << "pairs" << filteredPairs.size()
                 << "targets" << targetRows
                 << "decoys" << decoyRows
                 << "min_cos_sum100" << minCosineSimSum100
                 << "min_cos_time_cubed" << minCosineSimSpectrumOverTimeCubed
                 << "max_scan_time_delta_abs" << maxScanTimeDeltaAbs
                 << "local_rt_bin_sec" << pythiaParameters.timsLocalFdrRtBinSeconds;

        ERR_RETURN
    }

		void logNeuralNetStats(const QVector<KarnnNNTarget> &karnnNNTargetsNorm) {
    	const int totalCount = karnnNNTargetsNorm.size();
    	const int decoyCount = static_cast<int>(std::count_if(
				karnnNNTargetsNorm.begin(),
				karnnNNTargetsNorm.end(),
				[](const KarnnNNTarget &kt) { return is_nn_decoy(kt); }
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
    		float label = is_nn_decoy(kt) ? 1.0 : 0.0;
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
        const MsReaderPointerAcc *msReaderPointerAcc,
        int seed,
        QVector<CandidateScores*> *candidateScoreClassifier,
        bool *usedDiscriminantFallback
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_candidateScorePairs); ree;
    e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoys); ree;
    if (usedDiscriminantFallback != nullptr) {
        *usedDiscriminantFallback = false;
    }

    const bool isTimsRun = msReaderPointerAcc != nullptr
                           && !msReaderPointerAcc->ptr.isNull()
                           && msReaderPointerAcc->ptr->isTIMS();
    QVector<CandidateScores*> candidateScoresTargetsAndDecoysNeuralNet = candidateScoresTargetsAndDecoys;
    QVector<CandidateScores*> candidateScoresForDiscriminantFallback = candidateScoresTargetsAndDecoys;
    m_timsSecondStageCandidateScorePairs.clear();
	QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> targetDecoyCandidateScorePairsPntrs;
    QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> discriminantFallbackTargetDecoyPairPntrs;

    auto completeCandidateRowsToTargetDecoyPairs =
        [this](QVector<CandidateScores*> *candidateScoreRows) {
            if (candidateScoreRows == nullptr || candidateScoreRows->isEmpty()) {
                return;
            }

            QHash<TargetDecoyCandidatePair*, QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> pairByCandidate;
            pairByCandidate.reserve(m_candidateScorePairs.size());
            for (QPair<CandidateScoresTarget, CandidateScoresDecoy> &pr : m_candidateScorePairs) {
                TargetDecoyCandidatePair *candidatePair = pr.first.targetDecoyCandidatePair;
                if (candidatePair == nullptr) {
                    candidatePair = pr.second.targetDecoyCandidatePair;
                }
                if (candidatePair == nullptr) {
                    continue;
                }

                pairByCandidate.insert(candidatePair, {&pr.first, &pr.second});
            }

            QSet<CandidateScores*> existingRows;
            QSet<TargetDecoyCandidatePair*> selectedCandidatePairs;
            existingRows.reserve(candidateScoreRows->size());
            selectedCandidatePairs.reserve(candidateScoreRows->size() / 2);
            for (CandidateScores *candidateScores : *candidateScoreRows) {
                if (candidateScores == nullptr || candidateScores->targetDecoyCandidatePair == nullptr) {
                    continue;
                }

                existingRows.insert(candidateScores);
                selectedCandidatePairs.insert(candidateScores->targetDecoyCandidatePair);
            }

            const int initialRowCount = candidateScoreRows->size();
            int appendedComplementRows = 0;
            for (TargetDecoyCandidatePair *candidatePair : selectedCandidatePairs) {
                const auto pairIt = pairByCandidate.constFind(candidatePair);
                if (pairIt == pairByCandidate.constEnd()) {
                    continue;
                }

                CandidateScoresTarget *targetScores = pairIt.value().first;
                CandidateScoresDecoy *decoyScores = pairIt.value().second;
                if (targetScores != nullptr && !existingRows.contains(targetScores)) {
                    candidateScoreRows->push_back(targetScores);
                    existingRows.insert(targetScores);
                    ++appendedComplementRows;
                }
                if (decoyScores != nullptr && !existingRows.contains(decoyScores)) {
                    candidateScoreRows->push_back(decoyScores);
                    existingRows.insert(decoyScores);
                    ++appendedComplementRows;
                }
            }

            if (appendedComplementRows > 0) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMS neural-net candidate pool completed target-decoy pairs"
                         << "initial_rows" << initialRowCount
                         << "appended_rows" << appendedComplementRows
                         << "final_rows" << candidateScoreRows->size();
            }
        };

    auto rebuildTargetDecoyCandidateScorePairPointers =
        [&targetDecoyCandidateScorePairsPntrs](QVector<QPair<CandidateScoresTarget, CandidateScoresDecoy>> &candidateScorePairs) {
            targetDecoyCandidateScorePairsPntrs.clear();
            targetDecoyCandidateScorePairsPntrs.reserve(candidateScorePairs.size());
            for (QPair<CandidateScoresTarget, CandidateScoresDecoy> &pr : candidateScorePairs) {
                targetDecoyCandidateScorePairsPntrs.push_back({&pr.first, &pr.second});
            }
        };
    auto rebuildTargetDecoyCandidateScorePairPointersFromCandidateRows =
        [&targetDecoyCandidateScorePairsPntrs](const QVector<CandidateScores*> &candidateScoreRows) {
            QHash<TargetDecoyCandidatePair*, QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> pairByCandidate;
            pairByCandidate.reserve(candidateScoreRows.size() / 2);
            for (CandidateScores *candidateScores : candidateScoreRows) {
                if (candidateScores == nullptr || candidateScores->targetDecoyCandidatePair == nullptr) {
                    continue;
                }

                QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pair
                    = pairByCandidate[candidateScores->targetDecoyCandidatePair];
                if (candidateScores->isDecoy) {
                    pair.second = candidateScores;
                }
                else {
                    pair.first = candidateScores;
                }
            }

            targetDecoyCandidateScorePairsPntrs.clear();
            targetDecoyCandidateScorePairsPntrs.reserve(pairByCandidate.size());
            for (auto it = pairByCandidate.begin(); it != pairByCandidate.end(); ++it) {
                if (it.value().first != nullptr && it.value().second != nullptr) {
                    targetDecoyCandidateScorePairsPntrs.push_back(it.value());
                }
            }
        };

    rebuildTargetDecoyCandidateScorePairPointers(m_candidateScorePairs);

    if (m_pythiaParameters.writeFullCandidateDebug) {
        e = writeCandidateScoresDebug(
            candidateScoresTargetsAndDecoysNeuralNet,
            QStringLiteral(".scored_candidates.tsv"),
            m_outputFolderPath,
            msReaderPointerAcc
            ); ree;
    }

    const int neuralNetInferenceCandidateLimit =
        isTimsRun
            ? std::max(
                  m_pythiaParameters.neuralNetCandidateLimit,
                  m_pythiaParameters.timsNeuralNetInferenceCandidateLimit > 0
                      ? m_pythiaParameters.timsNeuralNetInferenceCandidateLimit
                      : TIMS_NEURAL_NET_AUTO_INFERENCE_CANDIDATE_LIMIT
                  )
            : m_pythiaParameters.neuralNetCandidateLimit;

    QVector<CandidateScores*> candidateScoresForNeuralNetTraining;
    if (isTimsRun && neuralNetInferenceCandidateLimit > m_pythiaParameters.neuralNetCandidateLimit) {
        candidateScoresForNeuralNetTraining = candidateScoresTargetsAndDecoys;
        e = filterScoredCandidatesForNeuralNet(
            m_pythiaParameters.minMs2FragCount,
            m_pythiaParameters.neuralNetCandidateLimit,
            true,
            &candidateScoresForNeuralNetTraining
            ); ree;
        completeCandidateRowsToTargetDecoyPairs(&candidateScoresForNeuralNetTraining);

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMS neural-net split pool"
                 << "training_rows" << candidateScoresForNeuralNetTraining.size()
                 << "inference_limit" << neuralNetInferenceCandidateLimit;
    }

    e = filterScoredCandidatesForNeuralNet(
        m_pythiaParameters.minMs2FragCount,
        neuralNetInferenceCandidateLimit,
        isTimsRun,
		        &candidateScoresTargetsAndDecoysNeuralNet
		        ); ree;

    if (isTimsRun) {
        completeCandidateRowsToTargetDecoyPairs(&candidateScoresTargetsAndDecoysNeuralNet);
        candidateScoresForDiscriminantFallback = candidateScoresTargetsAndDecoysNeuralNet;
        rebuildTargetDecoyCandidateScorePairPointersFromCandidateRows(candidateScoresForDiscriminantFallback);
        discriminantFallbackTargetDecoyPairPntrs = targetDecoyCandidateScorePairsPntrs;
    }

	    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Analyzing" << candidateScoresTargetsAndDecoysNeuralNet.size() << "for filtering";

    QVector<Features> neuralNetFeatures = m_neuralNetFeatures;
    constexpr bool enableTimsSecondStageRescoring = true;
    if (msReaderPointerAcc != nullptr
        && !msReaderPointerAcc->ptr.isNull()
        && msReaderPointerAcc->ptr->isTIMS()
        && enableTimsSecondStageRescoring) {
        e = rescoreTimsFilteredCandidatesForNeuralNet(
            msReaderPointerAcc,
            &candidateScoresTargetsAndDecoysNeuralNet,
            &m_timsSecondStageCandidateScorePairs,
            &neuralNetFeatures
            ); ree;

        if (!m_timsSecondStageCandidateScorePairs.isEmpty()) {
            rebuildTargetDecoyCandidateScorePairPointersFromCandidateRows(candidateScoresTargetsAndDecoysNeuralNet);
        }
    }
    else if (msReaderPointerAcc != nullptr
             && !msReaderPointerAcc->ptr.isNull()
             && msReaderPointerAcc->ptr->isTIMS()) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMS 4D second-stage rescoring skipped; using main IMS-aware candidate scores";
    }

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
        neuralNetFeatures,
        candidateScoresTargetsAndDecoysNeuralNet,
        &karnnNNTargetsNorm
        ); ree;

    QVector<KarnnNNTarget> karnnNNTargetsNormTrain = karnnNNTargetsNorm;
    QVector<KarnnNNTarget> karnnNNTargetsNormInferOnly;
    if (!candidateScoresForNeuralNetTraining.isEmpty()) {
        QSet<TargetDecoyCandidatePair*> trainingCandidatePairs;
        trainingCandidatePairs.reserve(candidateScoresForNeuralNetTraining.size() / 2);
        for (CandidateScores *candidateScores : candidateScoresForNeuralNetTraining) {
            if (candidateScores == nullptr || candidateScores->targetDecoyCandidatePair == nullptr) {
                continue;
            }
            trainingCandidatePairs.insert(candidateScores->targetDecoyCandidatePair);
        }

        karnnNNTargetsNormTrain.clear();
        karnnNNTargetsNormTrain.reserve(candidateScoresForNeuralNetTraining.size());
        karnnNNTargetsNormInferOnly.reserve(karnnNNTargetsNorm.size());
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            if (kt.candidateScores != nullptr
                && kt.candidateScores->targetDecoyCandidatePair != nullptr
                && trainingCandidatePairs.contains(kt.candidateScores->targetDecoyCandidatePair)) {
                karnnNNTargetsNormTrain.push_back(kt);
            }
            else {
                karnnNNTargetsNormInferOnly.push_back(kt);
            }
        }

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMS neural-net split pool normalized"
                 << "training_rows" << karnnNNTargetsNormTrain.size()
                 << "infer_only_rows" << karnnNNTargetsNormInferOnly.size()
                 << "total_inference_rows" << karnnNNTargetsNorm.size();
    }

	logNeuralNetStats(karnnNNTargetsNormTrain);

		QVector<QVector<KarnnNNTarget>> karnnNNTargetsNormTranched;
	e = ParallelUtils::trancheVectorForParallelization(
			karnnNNTargetsNormTrain,
			m_pythiaParameters.baggingSize,
			&karnnNNTargetsNormTranched
			); ree;

    QVector<FDRCLassifierNeuralNet> fdrClassifierNeuralNets;
	QVector<QVector<KarnnNNTarget>> inferenceKarnnVecs;
    e = trainNeuralNetwork(
            m_pythiaParameters,
            karnnNNTargetsNormTranched,
            seed,
            isTimsRun,
            &fdrClassifierNeuralNets,
            &inferenceKarnnVecs
            ); ree;

    const bool useMonotonePairQValues = msReaderPointerAcc != nullptr
                                        && !msReaderPointerAcc->ptr.isNull()
                                        && msReaderPointerAcc->ptr->isTIMS();
    const bool useTargetKeyStratifiedQValues = useMonotonePairQValues;
    const double localRtBinSecondsForQValues = useMonotonePairQValues
                                                   ? m_pythiaParameters.timsLocalFdrRtBinSeconds
                                                   : 0.0;
    const bool normalizeNeuralNetPredictions = isTimsRun
                                                   ? false
                                                   : m_pythiaParameters.normalizeNeuralNetPredictions;
    if (isTimsRun && m_pythiaParameters.normalizeNeuralNetPredictions) {
        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "TIMS neural-net prediction normalization disabled";
    }

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Inference start";
    QVector<QVector<float>> predictions;
    e = predictClassifierScores(
        inferenceKarnnVecs,
        fdrClassifierNeuralNets,
        normalizeNeuralNetPredictions,
        useMonotonePairQValues,
        &predictions
        ); ree;
	    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Inference end";

		for (int i = 0; i < predictions.size(); i++) {
		    processPredictions(predictions[i], &inferenceKarnnVecs[i], i); ree;
		}

        if (!karnnNNTargetsNormInferOnly.isEmpty()) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMS inference-only candidate prediction start"
                     << karnnNNTargetsNormInferOnly.size();
            QVector<float> inferOnlyPredictions;
            e = predictClassifierScoresWithAllNets(
                karnnNNTargetsNormInferOnly,
                fdrClassifierNeuralNets,
                normalizeNeuralNetPredictions,
                useMonotonePairQValues,
                &inferOnlyPredictions
                ); ree;
            e = processPredictions(
                inferOnlyPredictions,
                &karnnNNTargetsNormInferOnly,
                m_pythiaParameters.baggingSize
                ); ree;
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMS inference-only candidate prediction end";
        }

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
            &targetDecoyCandidateScorePairsPntrs,
            useMonotonePairQValues,
            useTargetKeyStratifiedQValues,
            localRtBinSecondsForQValues
            ); ree

    const double targetFdrThreshold = m_pythiaParameters.percentFDR / 100.0;
    int neuralNetTargetCount = 0;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
        candidateScoresTargetsAndDecoysNeuralNet,
        targetFdrThreshold,
        &neuralNetTargetCount
        ); ree;

    const bool canUseDiscriminantFallback = useMonotonePairQValues;
	    if (canUseDiscriminantFallback) {
	        QVector<CandidateScores*> discriminantCandidates = candidateScoresForDiscriminantFallback;
            QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> discriminantPairPntrs =
                discriminantFallbackTargetDecoyPairPntrs.isEmpty()
                    ? targetDecoyCandidateScorePairsPntrs
                    : discriminantFallbackTargetDecoyPairPntrs;
	        e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::DiscriminantScore,
            &discriminantPairPntrs,
            useMonotonePairQValues,
            useTargetKeyStratifiedQValues,
            localRtBinSecondsForQValues
            ); ree;

        int discriminantTargetCount = 0;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
            discriminantCandidates,
            targetFdrThreshold,
            &discriminantTargetCount
            ); ree;

        if (isTimsRun && m_pythiaParameters.timsHighEvidenceFilterEnabled) {
            constexpr double fdrQValThreshold = 0.5;

            QVector<CandidateScores*> neuralNetFilteredCandidates = candidateScoresTargetsAndDecoysNeuralNet;
            QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> neuralNetFilteredPairPntrs =
                targetDecoyCandidateScorePairsPntrs;
            e = applyTimsHighEvidenceGlobalQValueFilter(
                &neuralNetFilteredCandidates,
                &neuralNetFilteredPairPntrs,
                QValueSettertron::QValueScoreType::NNClassifierScore,
                m_pythiaParameters
                ); ree;
            e = filterByQValueThreshold(
                fdrQValThreshold,
                &neuralNetFilteredCandidates,
                true
                ); ree;
            int neuralNetFilteredTargetCount = 0;
            e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                neuralNetFilteredCandidates,
                targetFdrThreshold,
                &neuralNetFilteredTargetCount
                ); ree;

            QVector<CandidateScores*> discriminantFilteredCandidates = candidateScoresForDiscriminantFallback;
            QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> discriminantFilteredPairPntrs =
                discriminantFallbackTargetDecoyPairPntrs.isEmpty()
                    ? targetDecoyCandidateScorePairsPntrs
                    : discriminantFallbackTargetDecoyPairPntrs;
            e = applyTimsHighEvidenceGlobalQValueFilter(
                &discriminantFilteredCandidates,
                &discriminantFilteredPairPntrs,
                QValueSettertron::QValueScoreType::DiscriminantScore,
                m_pythiaParameters
                ); ree;
            e = filterByQValueThreshold(
                fdrQValThreshold,
                &discriminantFilteredCandidates,
                false
                ); ree;
            int discriminantFilteredTargetCount = 0;
            e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                discriminantFilteredCandidates,
                targetFdrThreshold,
                &discriminantFilteredTargetCount
                ); ree;

            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMS scorer selection after high-evidence filter"
                     << "neural-net targets" << neuralNetFilteredTargetCount
                     << "discriminant targets" << discriminantFilteredTargetCount
                     << "pre-filter neural-net targets" << neuralNetTargetCount
                     << "pre-filter discriminant targets" << discriminantTargetCount;

            if (discriminantFilteredTargetCount > neuralNetFilteredTargetCount) {
                qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                         << "TIMS neural-net fallback:"
                         << "using discriminant q-values after high-evidence filter";

                *candidateScoreClassifier = discriminantFilteredCandidates;
                if (usedDiscriminantFallback != nullptr) {
                    *usedDiscriminantFallback = true;
                }

                ERR_RETURN
            }

            neuralNetFilteredCandidates = candidateScoresTargetsAndDecoysNeuralNet;
            neuralNetFilteredPairPntrs = targetDecoyCandidateScorePairsPntrs;
            e = applyTimsHighEvidenceGlobalQValueFilter(
                &neuralNetFilteredCandidates,
                &neuralNetFilteredPairPntrs,
                QValueSettertron::QValueScoreType::NNClassifierScore,
                m_pythiaParameters
                ); ree;
            e = filterByQValueThreshold(
                fdrQValThreshold,
                &neuralNetFilteredCandidates,
                true
                ); ree;

            *candidateScoreClassifier = neuralNetFilteredCandidates;
            if (usedDiscriminantFallback != nullptr) {
                *usedDiscriminantFallback = false;
            }

            ERR_RETURN
        }

        if (discriminantTargetCount > neuralNetTargetCount) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMS neural-net fallback:"
                     << "using discriminant q-values"
                     << "linear targets" << discriminantTargetCount
                     << "neural-net targets" << neuralNetTargetCount;

            *candidateScoreClassifier = discriminantCandidates;
            if (usedDiscriminantFallback != nullptr) {
                *usedDiscriminantFallback = true;
            }

	            if (m_pythiaParameters.writeFullCandidateDebug) {
	                e = writeCandidateScoresDebug(
	                    *candidateScoreClassifier,
	                    QStringLiteral(".discriminant_candidates.tsv"),
	                    m_outputFolderPath,
                    msReaderPointerAcc
	                    ); ree;
	            }

                if (isTimsRun && m_pythiaParameters.timsHighEvidenceFilterEnabled) {
                e = applyTimsHighEvidenceGlobalQValueFilter(
                    candidateScoreClassifier,
                    &discriminantPairPntrs,
                    QValueSettertron::QValueScoreType::DiscriminantScore,
                    m_pythiaParameters
                    ); ree;
                }
                else if (isTimsRun) {
                    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                             << "TIMS high-evidence global q-value filter skipped";
                }

	            constexpr double fdrQValThreshold = 0.5;
                if (useMonotonePairQValues) {
                    e = filterByQValueThreshold(
                        fdrQValThreshold,
                        candidateScoreClassifier,
                        false
                        ); ree;
                }
                else {
                    e = serialFilterByValue(
                        fdrQValThreshold,
                        candidateScoreClassifier,
                        false
                        ); ree;
                }

            ERR_RETURN
	        }

	        e = QValueSettertron::setQValueForCandidates(
	            QValueSettertron::QValueScoreType::NNClassifierScore,
	            &targetDecoyCandidateScorePairsPntrs,
                useMonotonePairQValues,
                useTargetKeyStratifiedQValues,
                localRtBinSecondsForQValues
	            ); ree;
	    }

	    *candidateScoreClassifier = candidateScoresTargetsAndDecoysNeuralNet;

	    if (m_pythiaParameters.writeFullCandidateDebug) {
	        e = writeCandidateScoresDebug(
            candidateScoresTargetsAndDecoysNeuralNet,
            QStringLiteral(".nn_candidates.tsv"),
            m_outputFolderPath,
            msReaderPointerAcc
	            ); ree;
	    }

        if (isTimsRun && m_pythiaParameters.timsHighEvidenceFilterEnabled) {
            e = applyTimsHighEvidenceGlobalQValueFilter(
                candidateScoreClassifier,
                &targetDecoyCandidateScorePairsPntrs,
                QValueSettertron::QValueScoreType::NNClassifierScore,
                m_pythiaParameters
                ); ree;
        }
        else if (isTimsRun) {
            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "TIMS high-evidence global q-value filter skipped";
        }

	    constexpr double fdrQValThreshold = 0.5;
        if (useMonotonePairQValues) {
            e = filterByQValueThreshold(
                fdrQValThreshold,
                candidateScoreClassifier,
                true
                ); ree;
        }
        else {
            e = serialFilterByValue(
                fdrQValThreshold,
                candidateScoreClassifier,
                true
                ); ree;
        }

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
        if (is_nn_decoy(candidateScores->at(i))) {
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

        if (is_nn_decoy(cs)) {
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

        if (is_nn_decoy(cs)) {
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
            if (is_nn_decoy(candidateScores)) {
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
