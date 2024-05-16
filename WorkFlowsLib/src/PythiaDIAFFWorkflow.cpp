//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAFFWorkflow.h"

#include "CandidateScores.h"
#include "ClassifierWeightsManager.h"
#include "DiscriminantScoretron.h"
#include "EigenUtils.h"
#include "FastaFileToPeptidesListWorkFlow.h"
#include "FastaReader.h"
#include "FDRCLassifierNeuralNet.h"
#include "FragLibReader.h"
#include "InterferingCandidatesEliminatomatic.h"
#include "MsFrame.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
#include "PeptideStringWithMods.h"
#include "ProteinDigestomatic.h"
#include "QValueSettertron.h"
#include "TargetDecoyCandidatePair.h"
#include "TurboXIC.h"

#include <QElapsedTimer>


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

    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(&m_targetDecoyPairPntrs); ree

    ERR_RETURN
}

namespace {

    PeptideSequenceWithModsChargeAndTargetKey buildPeptideSequenceWithModsChargeAndTargetKey(const CandidateScores* candidateScores) {
        return candidateScores->targetDecoyCandidatePair->peptideStringWithMods()
             + "|"
             + QString::number(candidateScores->targetDecoyCandidatePair->charge())
             + "|"
             + candidateScores->targetKey;
    }

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

        QHash<PeptideStringWithMods, QVector<CandidateScoresEntry>> pepStrWModsVsCandScoresEntries;

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

                if (cse.mzTargetKey == cs->targetKey) {

                    switch (cs->targetDecoyCandidatePair->charge()) {
                        case 1:
                            if(!MathUtils::tSame(cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG], cs->featuresArray[CandidateScores::Features::CosineSimSum100])
                                && cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG] > 0.01) {
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_1] = cse.cosineSimSum100;
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_1] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                                break;
                            }

                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_OG] = cs->featuresArray[CandidateScores::Features::CosineSimSum100];
                            break;

                        case 2:
                            if(!MathUtils::tSame(cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG], cs->featuresArray[CandidateScores::Features::CosineSimSum100])
                               && cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG] > 0.01) {
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_1] = cse.cosineSimSum100;
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_1] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                                break;
                            }

                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_OG] = cs->featuresArray[CandidateScores::Features::CosineSimSum100];
                            break;
                        case 3:
                            if(!MathUtils::tSame(cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG], cs->featuresArray[CandidateScores::Features::CosineSimSum100])
                               && cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG] > 0.01) {
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_1] = cse.cosineSimSum100;
                                cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_1] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                                break;
                            }

                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_OG] = cs->featuresArray[CandidateScores::Features::CosineSimSum100];
                            break;
                        case 4:
                            if(!MathUtils::tSame(cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG], cs->featuresArray[CandidateScores::Features::CosineSimSum100])
                               && cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_OG] > 0.01) {
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
                        if (cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_2] > 0.01) {
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_3] = cse.cosineSimSum100;
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_3] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                            break;
                        }

                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge1_2] = cse.cosineSimSum100;
                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge1_2] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                        break;
                    case 2:
                        if (cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_2] > 0.01) {
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_3] = cse.cosineSimSum100;
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_3] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                            break;
                        }

                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge2_2] = cse.cosineSimSum100;
                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge2_2] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                        break;
                    case 3:
                        if (cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2] > 0.01) {
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_3] = cse.cosineSimSum100;
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_3] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                            break;
                        }

                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge3_2] = cse.cosineSimSum100;
                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge3_2] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                        break;
                    case 4:
                        if (cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_2] > 0.01) {
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_3] = cse.cosineSimSum100;
                            cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_3] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                            break;
                        }

                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdCosineSimSumCharge4_2] = cse.cosineSimSum100;
                        cs->featuresArray[CandidateScores::Features::AltTargetKeyIdTimeDeltaCharge4_2] = std::abs(cs->scanTime - cse.scanTime) / cs->scanTime;
                        break;
                }

            }
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;

    m_pythiaParameters.print();

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath); ree;

    e = m_targetDecoyCandidatePairScoretron.init(
            m_pythiaParameters,
            &msReaderPointerAcc
            ); ree;

     QVector<CandidateScores*> candidateScoresTrainings;
     e = buildCalibration(
             &msReaderPointerAcc,
             &candidateScoresTrainings
             ); ree;

//     if (m_msCalibratomatic.isInit()) {
//
//         e = optimizeParameters(
//                 candidateScoresTrainings,
//                 &msReaderPointerAcc
//         ); ree;
//     }
//
//     int targetCountBelowFDRThreshold;
//     e = mainAnalysis(
//             &msReaderPointerAcc,
//             &targetCountBelowFDRThreshold
//             ); ree;
//
//     QVector<CandidateScores*> candidateScoresTargetsAndDecoys;
//     e = buildCandidateScoresPtrs(&candidateScoresTargetsAndDecoys); ree;
//
// //#define TROUBLESHOOT_MISSING
// #ifdef TROUBLESHOOT_MISSING
//
//     //try shuffling library so it's not in alphabetical order.  Check first to see if it is in alpha order.
//
//     for (CandidateScores *cs : candidateScoresTargetsAndDecoys) {
//         if (cs->targetDecoyCandidatePair->peptideStringWithMods() == "YSQAVPAVTEGPIPEVLK" && cs->targetDecoyCandidatePair->charge() == 2) {
//             qDebug() << cs->targetDecoyCandidatePair->peptideStringWithMods() << cs->isDecoy;
//             qDebug() << cs->featuresArray[CandidateScores::Features::CosineSimSum100] << cs->discriminantScore << cs->scanTime << cs->qValue;
//             qDebug() << cs->featuresArray.mid(CandidateScores::Features::CosineSimToAnchor1, 12);
//             qDebug() << "^^^^^^^^^^";
//             einfo;
//         }
//     }
// #endif
//
//     QVector<CandidateScores*> candidateScoresTargetsAndDecoys50PercentFDRFiltered;
//     e = filterScoredCandidatesTo50PercentFDR(
//             &candidateScoresTargetsAndDecoys,
//             &candidateScoresTargetsAndDecoys50PercentFDRFiltered
//     ); ree;
//     qDebug() << "Analyzing" << candidateScoresTargetsAndDecoys50PercentFDRFiltered.size() << "for filtering";
//
//     e = InterferingCandidatesEliminatomatic::removeInterferingCandidates(
//             m_pythiaParameters.ionsSharedToReject,
//             m_pythiaParameters.mzMinMS2,
//             m_pythiaParameters.mzMaxMS2,
//             &candidateScoresTargetsAndDecoys50PercentFDRFiltered
//             ); ree;
//     qDebug() << candidateScoresTargetsAndDecoys50PercentFDRFiltered.size() << "after filtering";
//
//     e = populateAltIdTargetKeys(&candidateScoresTargetsAndDecoys50PercentFDRFiltered); ree;
//
// //#define WRITE_CANDIDATE_SCORES
// #ifdef WRITE_CANDIDATE_SCORES
//     e = updateProteinGroupAnnotation(
//             m_fastaUri,
//             &candidateScoresTargetsAndDecoys50PercentFDRFiltered
//     ); ree;
//
//     QVector<CandidateScoresReaderRow> candidateScoresToWrite;
//     std::transform(
//             candidateScoresTargetsAndDecoys50PercentFDRFiltered.begin(),
//             candidateScoresTargetsAndDecoys50PercentFDRFiltered.end(),
//             std::back_inserter(candidateScoresToWrite),
//             [](const CandidateScores *cs){return CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);}
//     );
//
//     e = ParquetReader::write(
//             candidateScoresToWrite,
//             msReaderPointerAcc.ptr->filePath() + ".candidateScores"
//             ); ree;
// #endif
//
//     QVector<CandidateScores*> candidateScoreClassifierPntrs;
//     if (!m_pythiaParameters.bypassNeuralNet) {
//
//         const int seedFirstTry = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST;
//         e = applyNeuralNetClassifier(
//                 candidateScoresTargetsAndDecoys50PercentFDRFiltered,
//                 seedFirstTry,
//                 &candidateScoreClassifierPntrs
//         ); ree;
//
//         if (candidateScoreClassifierPntrs.size() <= targetCountBelowFDRThreshold) {
//
//             std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
//             std::shuffle(
//                     candidateScoresTargetsAndDecoys50PercentFDRFiltered.begin(),
//                     candidateScoresTargetsAndDecoys50PercentFDRFiltered.end(),
//                     rng
//             );
//
//             const int seedSecondTry = S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST + 111;
//             candidateScoreClassifierPntrs.clear();
//             e = applyNeuralNetClassifier(
//                     candidateScoresTargetsAndDecoys50PercentFDRFiltered,
//                     seedSecondTry,
//                     &candidateScoreClassifierPntrs
//             ); ree;
//
//             if (candidateScoreClassifierPntrs.size() <= targetCountBelowFDRThreshold){
//                 QVector<CandidateScores*> candidateScoresPntrs;
//                 buildCandidateScoresPtrs(&candidateScoresPntrs);
//
//                 std::sort(candidateScoresPntrs.rbegin(), candidateScoresPntrs.rend(), [](CandidateScores *l, CandidateScores *r){
//                     return l->discriminantScore < r->discriminantScore;
//                 });
//
//                 candidateScoresPntrs.resize(targetCountBelowFDRThreshold);
//                 candidateScoreClassifierPntrs = candidateScoresPntrs;
//             }
//
//         }
//     }
//     else {
//         candidateScoreClassifierPntrs = candidateScoresTargetsAndDecoys50PercentFDRFiltered;
//     }
//
//     if (!m_pythiaParameters.reportDecoys) {
//
//         const auto terminatorLogicFDRFilter
//             = [&](CandidateScores *cs){return cs->qValue > m_pythiaParameters.percentFDR / 100.0;};
//
//         const auto terminator = std::remove_if(
//                 candidateScoreClassifierPntrs.begin(),
//                 candidateScoreClassifierPntrs.end(),
//                 terminatorLogicFDRFilter
//                 );
//
//         candidateScoreClassifierPntrs.erase(terminator, candidateScoreClassifierPntrs.end());
//     }
//
//     qDebug() << "Updating" << candidateScoreClassifierPntrs.size() << "PSMs";
//     e = updateProteinGroupAnnotation(
//             m_fastaUri,
//             &candidateScoreClassifierPntrs
//             ); ree;
//
// #define CALC_ENTRAP
// #ifdef CALC_ENTRAP
//     int counter = 0;
//     int decoys = 0;
//     int entrap = 0;
//     for (CandidateScores *cs : candidateScoreClassifierPntrs) {
//         counter++;
//
//         if (cs->proteinGroup.contains("_ARATH") && !cs->proteinGroup.contains("_HUMAN") && !cs->isDecoy) {
//             entrap++;
//         }
//
//         if (cs->isDecoy) {
//             decoys++;
//         }
//         if (decoys/static_cast<double>(counter) >= 0.01) {
//             break;
//         }
//
//     }
//     qDebug() << "Counter:" << counter << "Decoys:" <<  decoys << "Entrap:" << entrap << "Entrap%" << entrap / (double)counter;
// #endif
//
//     QVector<CandidateScoresReaderRow> candidateScoreReaderRows;
//     std::transform(
//             candidateScoreClassifierPntrs.begin(),
//             candidateScoreClassifierPntrs.end(),
//             std::back_inserter(candidateScoreReaderRows),
//             [](const CandidateScores *cs){return CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);}
//             );
//
//     const QString resultsFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
//     e = ParquetReader::write(candidateScoreReaderRows, resultsFilePath); ree;

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

    void filterByCosineSimSum100GreaterThan80(
            float ms2CorrMin,
            QVector<CandidateScores*> *candidateScores
            ) {

        const auto terminatorLogic = [ms2CorrMin](CandidateScores *cs){
            return cs->featuresArray[CandidateScores::Features::CosineSimSum100GreaterThan80] < ms2CorrMin;
        };

        const auto terminator = std::remove_if(
                candidateScores->begin(),
                candidateScores->end(),
                terminatorLogic
                );

        candidateScores->erase(terminator, candidateScores->end());

    }

    Err buildMsCalibrationReaderRows(
            const MSLevelEnum &msLevel,
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
        const auto msCalibrationReaderRowsInsertLogic = [msLevel](CandidateScores *cs){

            MsCalibarationReaderRow row;
            row.peptideStringWithMods = cs->targetDecoyCandidatePair->peptideStringWithMods();
            row.iRTPredicted = static_cast<float>(cs->targetDecoyCandidatePair->iRt());
            row.scanTime = cs->scanTime;
            row.scanNumber = cs->scanNumber;

            if (msLevel == MSLevelEnum::MS2) {
                row.mzSearchedVec = cs->featuresArray.mid(CandidateScores::Features::MzSearched1, top6);
                row.mzFoundMeanVec = cs->featuresArray.mid(CandidateScores::Features::MzFoundMean1, top6);
                row.mzFoundStDevVec = cs->featuresArray.mid(CandidateScores::Features::MzFoundStDev1, top6);
                row.intensityFoundMaxVec = cs->featuresArray.mid(CandidateScores::Features::IntensityFoundMax1, top6);
            }
            else {
                row.mzSearchedVec = {cs->targetDecoyCandidatePair->mz()};
                row.mzFoundMeanVec = {cs->featuresArray[CandidateScores::Features::Ms1MzMeanFound100]};
                row.mzFoundStDevVec = {cs->featuresArray[CandidateScores::Features::Ms1MzStDevFound100]};
                row.intensityFoundMaxVec = {cs->featuresArray[CandidateScores::Features::Ms1IntensityFound100]};
            }

            return row;
        };

        std::transform(
                candidateScoresFiltered.begin(),
                candidateScoresFiltered.end(),
                std::back_inserter(*msCalibrationReaderRows),
                msCalibrationReaderRowsInsertLogic
        );

        ERR_RETURN
    }

    void sortCandidatePointers(QVector<CandidateScores*> *candidateScoresPntrs) {

        std::sort(candidateScoresPntrs->rbegin(), candidateScoresPntrs->rend(), [](CandidateScores *l, CandidateScores *r){
            return l->discriminantScore < r->discriminantScore;
        });

    }

    void filterMs1CandidateRowsByCorr(QVector<CandidateScores*> *candidateScoresMS1Cal) {

        const double cosineSimSumMS1Min = 0.8;

        const auto terminatorLogic = [cosineSimSumMS1Min](const CandidateScores *cs){
            return cs->featuresArray[CandidateScores::Features::CosineSim100MS1] <= cosineSimSumMS1Min;
        };

        const auto terminator = std::remove_if(
                candidateScoresMS1Cal->begin(),
                candidateScoresMS1Cal->end(),
                terminatorLogic
                );

        candidateScoresMS1Cal->erase(terminator, candidateScoresMS1Cal->end());
    }

}//namespace
Err PythiaDIAFFWorkflow::buildCalibration(
        MsReaderPointerAcc *msReaderPointerAcc,
        QVector<CandidateScores*> *candidateScoresForTrainings
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    constexpr bool useExtendedScores = false;
    constexpr bool useNeuralNetworkScores = false;
    constexpr int minTrainingCountTranche = 50;
    constexpr int topNMS2IonsCalibration = 6;

    QPair<MzMin, MzMax> precursorMzMinMax;
    e = msReaderPointerAcc->ptr->getHiLoMzPrecursors(&precursorMzMinMax); ree;

    const auto sizePerTranche = static_cast<double>(m_pythiaParameters.trancheSizeMax);
    const int numberOfTranches = std::max(static_cast<int>(m_targetDecoyPairPntrs.size() / sizePerTranche), 1);
    qDebug() << "Number of tranches for calibration:" << numberOfTranches
             << "target count:" << m_targetDecoyPairPntrs.size()
             << "sizePerTranche:" << static_cast<int>(sizePerTranche);

    QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePointersTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            m_targetDecoyPairPntrs,
            numberOfTranches,
            &targetDecoyCandidatePointersTranched
    ); ree;

    int batchCounter = 0;

    QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsTopScores;
    for (const QVector<TargetDecoyCandidatePair*> &tdcp : targetDecoyCandidatePointersTranched) {

        QElapsedTimer etBatch;
        etBatch.start();

        QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsBatch = targetDecoyCandidatePairsTopScores;
        targetDecoyCandidatePairsBatch.append(tdcp);

        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
        e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
                targetDecoyCandidatePairsBatch,
                msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos(),
                &mzTargetKeyVsTargetDecoyCandidatePointers
                ); ree;

        m_candidateScores.clear();
        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2IonsCalibration,
                m_msCalibratomatic,
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                &m_candidateScores
                ); ree

        QVector<CandidateScores*> candidateScoresVecBatchPntrs;
        e = buildCandidateScoresPtrs(&candidateScoresVecBatchPntrs); ree;
        e = buildPeptideKeyVsTargetDecoyCandidateScoresPntrs(candidateScoresVecBatchPntrs); ree;

        QVector<float> weights;
        e = DiscriminantScoretron::trainLDAClassifier(
                m_peptideKeyVsTargetDecoyCandidateScoresPntrs.values(),
                useExtendedScores,
                useNeuralNetworkScores,
                &weights
                ); ree;

        e = DiscriminantScoretron::applyWeights(
            weights,
            useExtendedScores,
            useNeuralNetworkScores,
            &candidateScoresVecBatchPntrs
            ); ree;

        e = QValueSettertron::setQValueForCandidates(
             QValueSettertron::QValueScoreType::DiscriminantScore,
             &candidateScoresVecBatchPntrs
         ); ree;

        QMap<QString, int> fdrVsCounts;
        e = FDRCLassifierNeuralNet::outputFDRResults(
            candidateScoresVecBatchPntrs,
            true,
            &fdrVsCounts
        ); ree;

        const QString fdrKey = "10";

        e = honeIRTCalibration(
            &candidateScoresVecBatchPntrs,
            fdrVsCounts.value(fdrKey)
            ); ree;

        const int fdrCutoffCount = fdrVsCounts.value(fdrKey);
        if (fdrCutoffCount >= minTrainingCountTranche) {
            // candidateScoresVecScoredPntrs.resize(std::min(idealTrainingCountAtGivenFDR * 2, candidateScoresVecScoredPntrs.size()));
            // *candidateScoresForTrainings = candidateScoresVecBatchPntrs;
            qDebug() << "mzStDevMS1:" << m_msCalibratomatic.mzStDevMS1();
            qDebug() << "mzStDevMS2:" << m_msCalibratomatic.mzStDevMS2();
            qDebug() << "scanTimeWindowStDev x 3:" << m_msCalibratomatic.scanTimeStDev(S_GLOBAL_SETTINGS.STDEV_MULTIPLIER);

            break;
        }



        qDebug() << "Processed batch" << ++batchCounter << "of" << targetDecoyCandidatePointersTranched.size() << etBatch.elapsed() << "mSec";
        break; einfo;
    }



    ERR_RETURN
}

namespace {

    void filterTargetDecoyPairPointersByPrecursorMzRange(
        float precursorMzMin,
        float precursorMzMax,
        QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairs
        ) {

        const auto terminatorLogic = [precursorMzMin, precursorMzMax](const TargetDecoyCandidatePair *tdcp) {
            const float mzPrecursorTargetDecoyPair = tdcp->mz();
            return !(precursorMzMin <= mzPrecursorTargetDecoyPair && mzPrecursorTargetDecoyPair <= precursorMzMax);
        };

        const auto terminator = std::remove_if(
            targetDecoyCandidatePairs->begin(),
            targetDecoyCandidatePairs->end(),
            terminatorLogic
            );

        targetDecoyCandidatePairs->erase(terminator, targetDecoyCandidatePairs->end());
    }

} //namespace
Err PythiaDIAFFWorkflow::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
        const QVector<MsScanInfo> &msScanInfos,
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

        const float precursorMzMin = msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower;
        const float precursorMzMax = msScanInfo.precursorTargetMz + msScanInfo.isoWindowUpper;

        QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsFiltered = targetDecoyCandidatePairs;
        filterTargetDecoyPairPointersByPrecursorMzRange(precursorMzMin, precursorMzMax, &targetDecoyCandidatePairsFiltered);

        mzTargetKeyVsTargetDecoyCandidatePointers->insert(msScanInfo.targetKey(), targetDecoyCandidatePairsFiltered);

        if (m_pythiaParameters.verbosity > 1) {
            qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyCandidatePairsFiltered.size() << "targetDecoyPairs found";
        }
    }

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::recalibrateMs1Points(
        const QVector<CandidateScores*> &candidateScoresVecBatchPntrsResized,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(diaTargetFrames->isEmpty()); ree;
    e = ErrorUtils::isFalse(scanNumberVsScanTimeMS1->isEmpty()); ree;

    QVector<CandidateScores*> candidateScoresMS1Cal = candidateScoresVecBatchPntrsResized;
    qDebug() << candidateScoresMS1Cal.size() << "Precursors for MS1 calibration!";

    filterMs1CandidateRowsByCorr(&candidateScoresMS1Cal);
    
    const int recalibrationPointCountMin = 30;
    if (candidateScoresMS1Cal.size() < recalibrationPointCountMin) {
        ERR_RETURN
    }

    QVector<MsCalibarationReaderRow> msCalibrationReaderRowsMS1;
    e = buildMsCalibrationReaderRows(
            MSLevelEnum::MS1,
            candidateScoresMS1Cal,
            &msCalibrationReaderRowsMS1
    ); ree;

    if (msCalibrationReaderRowsMS1.size() < recalibrationPointCountMin) {
        ERR_RETURN
    }

    e = m_msCalibratomatic.setMassCalibrationCoeffs(
        msCalibrationReaderRowsMS1,
        MSLevelEnum::MS1
        ); ree;

    e = recalibrateMzVals(
            MSLevelEnum::MS1,
            diaTargetFrames,
            scanNumberVsScanTimeMS1
    ); ree;

    ERR_RETURN
}

namespace {

    Err recalibrationLogic(
        MsCalibratomatic msCalibratomatic,
        const QMap<ScanNumber, ScanPoints*> &scanPoints
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msCalibratomatic.isInitCalMS2() && msCalibratomatic.isInitCalMS1()); ree;
        e = ErrorUtils::isNotEmpty(scanPoints); ree;

        e = msCalibratomatic.recalibrateScanPoints(
            MSLevelEnum::MS2,
            scanPoints
            ); ree;

        e = msCalibratomatic.recalibrateScanPoints(
            MSLevelEnum::MS1,
            scanPoints
            ); ree;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::recalibrateMzVals(
        const MSLevelEnum &msLevel,
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msCalibratomatic.isInitCalMS2()); ree;
    e = ErrorUtils::isFalse(diaTargetFrames->isEmpty()); ree;
    e = ErrorUtils::isFalse(scanNumberVsScanTimeMS1->isEmpty()); ree;

    qDebug() << "Recalibrating mz vals";

    QElapsedTimer et;
    et.start();

    if (msLevel == MSLevelEnum::MS2) {

#define RECAL_PARALLEL
#ifdef RECAL_PARALLEL
        const auto reCalBinder = std::bind(
                recalibrationLogic,
                m_msCalibratomatic,
                std::placeholders::_1
        );

        const QList<QMap<ScanNumber, ScanPoints *>> &tandemPointsVec = diaTargetFrames->values();

        QFuture<Err> futures = QtConcurrent::mapped(
                tandemPointsVec,
                reCalBinder
        );
        futures.waitForFinished();

        for (const Err &err: futures) {
            e = err; ree;
        }
#else
        const QList<MzTargetKey> &diaTargetFrameKeys = diaTargetFrames->keys();
        for (const MzTargetKey &k: diaTargetFrameKeys) {
            e = recalibrationLogic(m_msCalibratomatic, diaTargetFrames->value(k)); ree;
        }
#endif
    }

    if (msLevel == MSLevelEnum::MS1) {
        qDebug() << "Recalibrating MS1 mz vals frame";
        e = m_msCalibratomatic.recalibrateScanPoints(
            MSLevelEnum::MS1,
            scanNumberVsScanTimeMS1
            ); ree;
    }

    qDebug() << "Points recalibrated in" << et.elapsed() << "mSec";

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
//                {1.5, 2.0},
//                {2.0, 2.0},
//                {2.5,  2.0},
//                {3.5,  2.0},
                {4.5, 2.0},
                {5.5, 2.0},
                {6.5, 2.0},
                {7.5, 2.0},
                {8.5, 2.0},
                {10.0, 2.0}
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
    e = ErrorUtils::isNotEmpty(candidateScoresTrainings); ree;

    // const int topNMS2IonsOptimization = 12;
    // qDebug() << "Using top:" << topNMS2IonsOptimization << "fragments for optimization";
    //
    // qDebug() << "Optimization selection fraction" << candidateScoresTrainings.size();
    //
    // QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    // for (CandidateScores *cs : candidateScoresTrainings) {
    //     mzTargetKeyVsTargetDecoyCandidatePointers[cs->targetKey].push_back(cs->targetDecoyCandidatePair);
    // }
    //
    // QVector<PythiaParameters> pythiaParametersExperiments;
    // e = buildDOE(
    //         m_pythiaParameters,
    //         m_msCalibratomatic.mzStDevMS2(),
    //         m_msCalibratomatic.scanTimeStDev(),
    //         &pythiaParametersExperiments
    //         ); ree;
    //
    // const bool useExtendedScores = true;
    // const bool useNeuralNetworkScores = false;
    //
    // QVector<DOEResult> results;
    // for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {
    //
    //     qDebug() << "ppmTol" << pythiaParams.ms2ExtractionWidthPPM;
    //
    //     e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(pythiaParams); ree;
    //     qDebug() << "STarting opt";
    //
    //     m_candidateScores.clear();
    //
    //     e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
    //             topNMS2IonsOptimization,
    //             m_msCalibratomatic,
    //             &mzTargetKeyVsTargetDecoyCandidatePointers,
    //             &m_candidateScores
    //     ); ree;
    //
    //     QVector<CandidateScores*> candidateScoresPntrs;
    //     e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;
    //     e = DiscriminantScoretron::setDiscriminantScoreForCandidates(
    //             useExtendedScores,
    //             useNeuralNetworkScores,
    //             &candidateScoresPntrs
    //     ); ree;
    //
    //     e = QValueSettertron::setQValueForCandidates(
    //             QValueSettertron::QValueScoreType::DiscriminantScore,
    //             &candidateScoresPntrs
    //             ); ree;
    //
    //     QMap<QString, int> fdrVsCount;
    //     e = FDRCLassifierNeuralNet::outputFDRResults(
    //             m_candidateScores,
    //             true,
    //             &fdrVsCount
    //     ); ree;
    //     qDebug() << "Ending opt";
    //
    //     const double fdrThreshold = 0.1;
    //     int targetCountAboveFDRQValueThreshold;
    //     e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
    //             m_candidateScores,
    //             fdrThreshold,
    //             &targetCountAboveFDRQValueThreshold
    //             ); ree;
    //
    //     DOEResult res;
    //     res.ppm = pythiaParams.ms2ExtractionWidthPPM;
    //     res.fdrCount = targetCountAboveFDRQValueThreshold;
    //     results.push_back(res);
    // }
    //
    // e = getTopFrequencyParameters(
    //         &results,
    //         &m_pythiaParameters.ms2ExtractionWidthPPM
    //         ); ree;
    //
    // m_pythiaParameters.ms1ExtractionWidthPPM = m_pythiaParameters.ms2ExtractionWidthPPM;
    //
    // e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(m_pythiaParameters); ree;
    //
    // qDebug() << "Optimal ppm setting:" << m_pythiaParameters.ms2ExtractionWidthPPM;
    // qDebug() << "Optimal scanTimeWindow setting:" << m_msCalibratomatic.scanTimeStDev(m_pythiaParameters.scanTimeWindowStDevs);


    ERR_RETURN
}

Err PythiaDIAFFWorkflow::mainAnalysis(
        MsReaderPointerAcc *msReaderPointerAcc,
        int *targetCountBelowFDRThresholdOnePercent
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

    // m_candidateScores.clear();
    //
    // const bool useExtendedScores = true;
    // const bool useNeuralNetworkScores = false;
    //
    // const int topNMs2IonsMainAnalysis = std::max(
    //         m_minTopNMs2Ions,
    //         static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions))
    //         );
    //
    // qDebug() << "Using top:" << topNMs2IonsMainAnalysis << "fragments for main analysis";
    //
    // const double selectionFractionBypass = -1.0;
    //
    // QElapsedTimer et;
    // et.start();
    //
    // QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    // e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
    //         msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos(),
    //         selectionFractionBypass,
    //         &mzTargetKeyVsTargetDecoyCandidatePointers
    // ); ree;
    //
    // e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
    //         topNMs2IonsMainAnalysis,
    //         msReaderPointerAcc->ptr->scanTimeMinMax(),
    //         m_msCalibratomatic,
    //         &mzTargetKeyVsTargetDecoyCandidatePointers,
    //         &m_candidateScores
    // ); ree;
    // qDebug() << "Targets scored" << et.restart() << "mSec";
    //
    // QVector<CandidateScores*> candidateScoresPntrs;
    // e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;
    // e = DiscriminantScoretron::setDiscriminantScoreForCandidates(
    //         useExtendedScores,
    //         useNeuralNetworkScores,
    //         &candidateScoresPntrs
    // ); ree;
    // qDebug() << "Targets discriminated" << et.restart() << "mSec";
    //
    // e = QValueSettertron::setQValueForCandidates(
    //         QValueSettertron::QValueScoreType::DiscriminantScore,
    //         &candidateScoresPntrs
    //         ); ree;
    // qDebug() << "Targets qval'd" << et.restart() << "mSec";
    //
    // QMap<QString, int> fdrVsCount;
    // e = FDRCLassifierNeuralNet::outputFDRResults(
    //         m_candidateScores,
    //         true,
    //         &fdrVsCount
    // ); ree;
    // qDebug() << "Targets resulted" << et.restart() << "mSec";
    //
    // const double fdrThresholdOnePercent = 0.01;
    // e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
    //         m_candidateScores,
    //         fdrThresholdOnePercent,
    //         targetCountBelowFDRThresholdOnePercent
    // ); ree;
    // qDebug() << "Targets counted" << et.restart() << "mSec";
    //
    // std::sort(candidateScoresPntrs.rbegin(), candidateScoresPntrs.rend(), [](CandidateScores *l, CandidateScores *r){
    //     return l->discriminantScore < r->discriminantScore;
    // });
    // qDebug() << "Targets sorted" << et.restart() << "mSec";

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::buildCandidateScoresPtrs(QVector<CandidateScores*> *candidateScoresPntrs) {

    ERR_INIT
    e = buildCandidateScoresPtrs(
            m_candidateScores,
            candidateScoresPntrs
            ); ree;
    ERR_RETURN
}

Err PythiaDIAFFWorkflow::buildCandidateScoresPtrs(
        QVector<CandidateScores> &candidateScores,
        QVector<CandidateScores*> *candidateScoresPntrs
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;

    std::transform(
            candidateScores.begin(),
            candidateScores.end(),
            std::back_inserter(*candidateScoresPntrs),
            [](CandidateScores &cs){return &cs;}
    );

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

        const int batchSize = std::min(200, std::max(1, static_cast<int>(karnnNNTargetsNorm.size() / 100.0)));

        qDebug() << "Batch Size:" << batchSize;

        QVector<QVector<float>> xData;
        QVector<float> yData;
        for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
            xData.push_back(kt.scoreVec);
            yData.push_back(kt.isDecoy ? 1.0 : 0.0);
        }

        const int baggingSize = 8;
        const float learningRate = 0.003;
        const int epochs = 3;
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
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(m_candidateScores); ree;

    candidateScoreClassifier->clear();

    const int totalCount = candidateScoresTargetsAndDecoys50PercentFDRFiltered.size();
    const int decoyCount = static_cast<int>(std::count_if(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered.begin(),
            candidateScoresTargetsAndDecoys50PercentFDRFiltered.end(),
            [](const CandidateScores* cs){return cs->isDecoy;}
            ));

    qDebug() << "target vs decoy count" << totalCount - decoyCount << decoyCount
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

    QMap<PeptideString, QVector<FastaEntry>> peptideStringWithModsVsFastaEntriesLeucinesReplaced;
    for (auto it = peptideStringWithModsVsFastaEntries.begin(); it != peptideStringWithModsVsFastaEntries.end(); it++) {
        PeptideString peptideSeqReplacedLeucines = it.key().removeUniModChars();
        peptideSeqReplacedLeucines = PeptideStringWithMods(peptideSeqReplacedLeucines.replace('L', 'J').replace('I', 'J'));
        peptideStringWithModsVsFastaEntriesLeucinesReplaced[peptideSeqReplacedLeucines].append(it.value());
    }

    for (int i = 0; i < candidateScores->size(); i++) {

        CandidateScores *cs = (*candidateScores)[i];

        PeptideString peptideSeqReplacedLeucines = cs->targetDecoyCandidatePair->peptideStringWithMods().removeUniModChars();
        peptideSeqReplacedLeucines = PeptideStringWithMods(peptideSeqReplacedLeucines.replace('L', 'J').replace('I', 'J'));

        const QVector<FastaEntry> &fastaEntries = peptideStringWithModsVsFastaEntriesLeucinesReplaced.value(peptideSeqReplacedLeucines);

        QStringList fastaDescriptions;
        std::transform(
                fastaEntries.begin(),
                fastaEntries.end(),
                std::back_inserter(fastaDescriptions),
                [](const FastaEntry &fe){return fe.fastaDescription;}
                );

        cs->proteinGroup = fastaDescriptions.join(';');
        cs->proteinGroup.replace(",", " ");

    }

    ERR_RETURN
}

Err PythiaDIAFFWorkflow::buildPeptideKeyVsTargetDecoyCandidateScoresPntrs(const QVector<CandidateScores*> &candidateScores) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;

    m_peptideKeyVsTargetDecoyCandidateScoresPntrs.clear();

    for (CandidateScores *cs: candidateScores) {

        const PeptideSequenceWithModsChargeAndTargetKey peptideSequenceWithModsChargeAndTargetKey
                                                                    = buildPeptideSequenceWithModsChargeAndTargetKey(cs);

        if (cs->isDecoy) {
            m_peptideKeyVsTargetDecoyCandidateScoresPntrs[peptideSequenceWithModsChargeAndTargetKey].second = cs;
            continue;
        }
        m_peptideKeyVsTargetDecoyCandidateScoresPntrs[peptideSequenceWithModsChargeAndTargetKey].first = cs;
    }

    // for (auto &dd : m_peptideKeyVsTargetDecoyCandidateScoresPntrs) {
    //     qDebug() << dd.first << dd.second << "SDLFDKSJL";
    // }


    const bool allTargetsMatchedWithDecoy = std::all_of(
            m_peptideKeyVsTargetDecoyCandidateScoresPntrs.begin(),
            m_peptideKeyVsTargetDecoyCandidateScoresPntrs.end(),
            [](const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &pr){
                return pr.first != nullptr && pr.second != nullptr;
            }
    );
    e = ErrorUtils::isTrue(allTargetsMatchedWithDecoy); ree;

    ERR_RETURN
}


Err PythiaDIAFFWorkflow::honeIRTCalibration(
    QVector<CandidateScores*> *candidateScoresVecScoredPntrs,
    int topNCandidates
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScoresVecScoredPntrs->isEmpty()); ree;

    std::sort(
        candidateScoresVecScoredPntrs->rbegin(),
        candidateScoresVecScoredPntrs->rend(),
        [](CandidateScores *l, CandidateScores *r){return l->discriminantScore < r->discriminantScore;}
        );

    QVector<CandidateScores*> candidateScoresVecBatchPntrsResized = *candidateScoresVecScoredPntrs;
    candidateScoresVecBatchPntrsResized.resize(std::max(m_minTrainingCountTranche, topNCandidates));

    qDebug() << "Using" << candidateScoresVecBatchPntrsResized.size() << "for iRT Estimation";

    e = InterferingCandidatesEliminatomatic::removeInterferingCandidates(
            m_pythiaParameters.ionsSharedToReject,
            m_pythiaParameters.mzMinMS2,
            m_pythiaParameters.mzMaxMS2,
            &candidateScoresVecBatchPntrsResized
    ); ree;
    qDebug() << candidateScoresVecBatchPntrsResized.size() << "after filtering";

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = buildMsCalibrationReaderRows(
            MSLevelEnum::MS2,
            candidateScoresVecBatchPntrsResized,
            &msCalibrationReaderRows
    ); ree;

    e = m_msCalibratomatic.buildRTMapper(msCalibrationReaderRows); ree;
    qDebug() << "scanTimeWindowStDev x" << S_GLOBAL_SETTINGS.STDEV_MULTIPLIER
             <<":" << m_msCalibratomatic.scanTimeStDev(S_GLOBAL_SETTINGS.STDEV_MULTIPLIER);

    ERR_RETURN
}
