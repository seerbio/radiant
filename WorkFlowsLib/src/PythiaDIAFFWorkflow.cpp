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

     if (m_msCalibratomatic.isInitRT()) {
         e = optimizeParameters(candidateScoresTrainings); ree;
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

     e = InterferingCandidatesEliminatomatic::removeInterferingCandidates(
             m_pythiaParameters.ionsSharedToReject,
             m_pythiaParameters.mzMinMS2,
             m_pythiaParameters.mzMaxMS2,
             &candidateScoresTargetsAndDecoys50PercentFDRFiltered
             ); ree;
     qDebug() << candidateScoresTargetsAndDecoys50PercentFDRFiltered.size() << "after filtering";

     e = populateAltIdTargetKeys(&candidateScoresTargetsAndDecoys50PercentFDRFiltered); ree;

// #define WRITE_CANDIDATE_SCORES
#ifdef WRITE_CANDIDATE_SCORES
    // e = updateProteinGroupAnnotation(
    //         m_fastaUri,
    //         &candidateScoresTargetsAndDecoys
    // ); ree;

    QVector<CandidateScoresReaderRow> candidateScoresToWrite;
    std::transform(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered.begin(),
            candidateScoresTargetsAndDecoys50PercentFDRFiltered.end(),
            std::back_inserter(candidateScoresToWrite),
            [](const CandidateScores *cs){return CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);}
    );

    e = ParquetReader::write(
            candidateScoresToWrite,
            msReaderPointerAcc.ptr->filePath() + ".candidateScoresNew"
            ); ree;
#endif

    QVector<CandidateScores*> candidateScoreClassifierPntrs;
    if (!m_pythiaParameters.bypassNeuralNet) {

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
    }
    else {
        candidateScoreClassifierPntrs = candidateScoresTargetsAndDecoys50PercentFDRFiltered;
        std::sort(
            candidateScoreClassifierPntrs.rbegin(),
            candidateScoreClassifierPntrs.rend(),
            [](CandidateScores *l, CandidateScores *r){return l->discriminantScore < r->discriminantScore;}
            );
    }

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

    qDebug() << "Updating" << candidateScoreClassifierPntrs.size() << "PSMs";
    e = updateProteinGroupAnnotation(
            m_fastaUri,
            &candidateScoreClassifierPntrs
            ); ree;

#define CALC_ENTRAP
#ifdef CALC_ENTRAP
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
    qDebug() << "Counter:" << counter << "Decoys:" <<  decoys << "Entrap:" << entrap << "Entrap%" << entrap / (double)counter;
#endif

    QVector<CandidateScoresReaderRow> candidateScoreReaderRows;
    std::transform(
            candidateScoreClassifierPntrs.begin(),
            candidateScoreClassifierPntrs.end(),
            std::back_inserter(candidateScoreReaderRows),
            [](const CandidateScores *cs){return CandidateScoresReaderRow::buildCandidateScoresReaderRow(cs);}
            );

    const QString resultsFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
    e = ParquetReader::write(candidateScoreReaderRows, resultsFilePath); ree;

    ERR_RETURN
}

namespace {

    Err buildMzTargetKeyVsTurboXicPntrs(
        const QVector<MsScanInfo> &uniqueMsScanInfos,
        const QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> &diaTargetFrames,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        QMap<MzTargetKey, TurboXIC*> *mzTargetKeyVsTurboXicPntr
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(uniqueMsScanInfos); ree;
        e = ErrorUtils::isNotEmpty(diaTargetFrames); ree;
        e = ErrorUtils::isNotEmpty(scanNumberVsScanTime); ree;

        for (const MsScanInfo &msScanInfo : uniqueMsScanInfos) {

            const MzTargetKey mzTargetKey = msScanInfo.targetKey();
            e = ErrorUtils::contains(mzTargetKey, diaTargetFrames); ree;

            const QMap<ScanNumber, ScanPoints*> &scanNumberVsScanPoints = diaTargetFrames.value(mzTargetKey);

            MsFrame msFrame;
            e = msFrame.init(scanNumberVsScanPoints, scanNumberVsScanTime); ree;

            auto turboXic = new TurboXIC();
            e = turboXic->init(msFrame.frameIndexVsScanPoints()); ree;
            mzTargetKeyVsTurboXicPntr->insert(mzTargetKey, turboXic);
        }

        ERR_RETURN
    }

    void filterDuplicateCandidateScoresByDiscriminantScore(QVector<CandidateScores*> *candidateScores) {

        QMap<QString, CandidateScores*> keyVsCandidatesFoundBest;

        // for (auto *cs : *candidateScores) {
        //     qDebug() << cs->targetDecoyCandidatePair->peptideStringWithMods() << cs->targetDecoyCandidatePair->charge() << cs->isDecoy << cs->targetKey << "SDLFKJSL";
        // }

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
            const MSLevelEnum &msLevel,
            const QVector<CandidateScores*> &_candidateScores,
            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(_candidateScores); ree;

        qDebug() << _candidateScores.size() << "Found for recalibartion";

        QVector<CandidateScores*> candidateScoresFiltered = _candidateScores;
        filterDuplicateCandidateScoresByDiscriminantScore(&candidateScoresFiltered);
        qDebug() << candidateScoresFiltered.size() << "Found for recalibartion after duplicates filtered";

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

    void sortCandidatePointersDiscScoreDesc(QVector<CandidateScores*> *candidateScoresPntrs) {
        std::sort(candidateScoresPntrs->rbegin(), candidateScoresPntrs->rend(), [](const CandidateScores *l, const CandidateScores *r){
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
        const MsReaderPointerAcc *msReaderPointerAcc,
        QVector<CandidateScores*> *candidateScoresForTrainings
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    const auto sizePerTranche = static_cast<double>(m_pythiaParameters.trancheSizeMax);
    const int numberOfTranches = std::max(static_cast<int>(m_targetDecoyPairPntrs.size() / sizePerTranche), 1);
    qDebug() << "Number of tranches for calibration:" << numberOfTranches
             << "target count:" << m_targetDecoyPairPntrs.size()
             << "sizePerTranche:" << static_cast<int>(sizePerTranche);

    QVector<MsScanInfo> uniqueMsScanInfos = msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    constexpr int maxTrainingCount = 24;
    const int mzTargetKeysHalf = static_cast<int>(std::round(uniqueMsScanInfos.size() / 2.0));

    const int resizeCount = std::min(maxTrainingCount, mzTargetKeysHalf);
    uniqueMsScanInfos = uniqueMsScanInfos.mid(uniqueMsScanInfos.size() - resizeCount);

    QMap<MzTargetKey, TurboXIC*> mzTargetKeyVsTurboXicPntrs;
    e = buildMzTargetKeyVsTurboXicPntrs(
        uniqueMsScanInfos,
        *m_targetDecoyCandidatePairScoretron.diaTargetFrames(),
        msReaderPointerAcc->ptr->getScanNumberVsScanTime(),
        &mzTargetKeyVsTurboXicPntrs
        ); ree;

    QVector<QVector<TargetDecoyCandidatePair*>> targetDecoyCandidatePointersTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            m_targetDecoyPairPntrs,
            numberOfTranches,
            &targetDecoyCandidatePointersTranched
            ); ree;

    QVector<float> scanTimeStDevs;
    QVector<float> ms2PPMStDevs;

    int batchCounter = 0;
    for (const QVector<TargetDecoyCandidatePair*> &tdcp : targetDecoyCandidatePointersTranched) {

        constexpr bool useExtendedScores = false;
        constexpr bool useNeuralNetworkScores = false;
        constexpr int topNMS2IonsCalibration = 6;

        QElapsedTimer etBatch;
        etBatch.start();

        QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsBatch = m_targetDecoyCandidatePairsTopScores;
        targetDecoyCandidatePairsBatch.append(tdcp);

        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
        e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
                targetDecoyCandidatePairsBatch,
                uniqueMsScanInfos,
                &mzTargetKeyVsTargetDecoyCandidatePointers
                ); ree;

        constexpr float minPeakCountCalibration = 4.9;
        m_candidateScores.clear();
        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2IonsCalibration,
                m_msCalibratomatic,
                minPeakCountCalibration,
                mzTargetKeyVsTurboXicPntrs,
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                &m_candidateScores
                ); ree

        QVector<CandidateScores*> candidateScoresVecBatchPntrs;
        QMap<int, int> fdrVsCounts;
        e = processBatch(
            useExtendedScores,
            useNeuralNetworkScores,
            &candidateScoresVecBatchPntrs,
            &fdrVsCounts
            ); ree;

        constexpr int fdrKey = 5;
        constexpr int fdrKeyMassCal = 2;

        e = honeIRTAndMassCalibration(
            mzTargetKeyVsTurboXicPntrs.keys().toVector(),
            &candidateScoresVecBatchPntrs,
            fdrVsCounts.value(fdrKey),
            fdrVsCounts.value(fdrKeyMassCal)
            ); ree;

        scanTimeStDevs.push_back(m_msCalibratomatic.scanTimeStDev());
        ms2PPMStDevs.push_back(m_msCalibratomatic.mzStDevMS2());

        qDebug() << "********* Processed batch" << ++batchCounter << "of" << targetDecoyCandidatePointersTranched.size() << etBatch.elapsed() << "mSec ********";

        constexpr int targetTrainingCountCalibration = 1000;
        if (fdrVsCounts.value(fdrKey) > targetTrainingCountCalibration || &tdcp == targetDecoyCandidatePointersTranched.cend() - 1) {

            constexpr int fdrKey50PercentFDR = 50;
            candidateScoresVecBatchPntrs.resize(std::min(candidateScoresVecBatchPntrs.size(), fdrVsCounts.value(fdrKey50PercentFDR)));

            QVector<CandidateScores*> candidateScoresVecBatchPntrsRecal = candidateScoresVecBatchPntrs;
            candidateScoresVecBatchPntrs.resize(std::min(candidateScoresVecBatchPntrs.size(), fdrVsCounts.value(fdrKeyMassCal)));

            QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
            e = buildMsCalibrationReaderRows(
                    MSLevelEnum::MS2,
                    candidateScoresVecBatchPntrsRecal,
                    &msCalibrationReaderRows
            ); ree;

            e = m_msCalibratomatic.setMassCalibrationCoeffs(
                msCalibrationReaderRows,
                MSLevelEnum::MS2
                ); ree;

            scanTimeStDevs.push_back(m_msCalibratomatic.scanTimeStDev());
            ms2PPMStDevs.push_back(m_msCalibratomatic.mzStDevMS2());

            constexpr int minVecSize = 3;
            std::sort(scanTimeStDevs.begin(), scanTimeStDevs.end());
            if (scanTimeStDevs.size() >= minVecSize) {
                scanTimeStDevs.pop_front();
                scanTimeStDevs.pop_back();

            }

            std::sort(ms2PPMStDevs.begin(), ms2PPMStDevs.end());
            if (ms2PPMStDevs.size() >= minVecSize) {
                ms2PPMStDevs.pop_front();
                ms2PPMStDevs.pop_back();
            }

            qDebug() << scanTimeStDevs;
            qDebug() << "ScanTimeWindow Mean|Median|Min"
                     << MathUtils::mean(scanTimeStDevs)
                     << MathUtils::median(scanTimeStDevs)
                     << *std::min({scanTimeStDevs.begin(), scanTimeStDevs.end()});

            qDebug() << "Ms2 ppm Mean|Median" << MathUtils::mean(ms2PPMStDevs) << MathUtils::median(ms2PPMStDevs);

            m_msCalibratomatic.setScanTimeStDev(scanTimeStDevs.front());
            m_msCalibratomatic.setMzStDevMS2(MathUtils::mean(ms2PPMStDevs));

            e = recalibrateMzVals(
                    MSLevelEnum::MS2,
                    m_targetDecoyCandidatePairScoretron.diaTargetFrames(),
                    m_targetDecoyCandidatePairScoretron.ms1ScanNumberVsScanPoints()
            ); ree;


            *candidateScoresForTrainings = candidateScoresVecBatchPntrs;

            break;
        }

    }

    for (TurboXIC* turboXic : mzTargetKeyVsTurboXicPntrs) {delete turboXic;}

    ERR_RETURN
}

namespace {

    Err buildProcessBatchPointers(
        const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &candidateScorePairs,
        QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> *featuresArrayTargetVsDecoy,
        QVector<CandidateScores*> *candidateScoresesPntrs,
        QVector<FeaturesArray*> *featuresArrayPntrs,
        QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> *featuresArrayTargetVsDecoyPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isFalse(featuresArrayTargetVsDecoy->isEmpty()); ree;

        for (int i = 0; i < featuresArrayTargetVsDecoy->size(); i++) {

            const QPair<CandidateScoresTarget*, CandidateScoresDecoy*> &candidatePair = candidateScorePairs.at(i);
            QPair<FeaturesArrayTargets, FeaturesArrayDecoys> &featuresArrayPair = (*featuresArrayTargetVsDecoy)[i];
            featuresArrayTargetVsDecoyPntrs->push_back({&featuresArrayPair.first, &featuresArrayPair.second});

            candidateScoresesPntrs->append({candidatePair.first, candidatePair.second});
            featuresArrayPntrs->append({&featuresArrayPair.first, &featuresArrayPair.second});
        }

        ERR_RETURN
    }

    Err setCandidateScoresDiscriminantScore(
        const QVector<float> &discriminantScores,
        QVector<CandidateScores*> *candidateScoresesPntrs
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(discriminantScores); ree;
        e = ErrorUtils::isEqual(discriminantScores.size(), candidateScoresesPntrs->size()); ree;

        for (int i = 0; i < discriminantScores.size(); i++) {
            CandidateScores* cs = (*candidateScoresesPntrs)[i];
            cs->discriminantScore = discriminantScores.at(i);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::processBatch(
    bool useExtendedScores,
    bool useNeuralNetworkScores,
    QVector<CandidateScores*> *candidateScoresVecBatchPntrs,
    QMap<int, int> *fdrVsCounts
    ) {

    ERR_INIT

    e = buildCandidateScoresPtrs(candidateScoresVecBatchPntrs); ree;
    e = buildPeptideKeyVsTargetDecoyCandidateScoresPntrs(*candidateScoresVecBatchPntrs); ree;

    QVector<QPair<FeaturesArrayTargets, FeaturesArrayDecoys>> featuresArrayTargetVsDecoy;
    e = DiscriminantScoretron::convertScoreCandidatesToFeaturesArrays(
        m_peptideKeyVsTargetDecoyCandidateScoresPntrs.values().toVector(),
        useExtendedScores,
        useNeuralNetworkScores,
        &featuresArrayTargetVsDecoy
        ); ree;

    e = ErrorUtils::isEqual(
        m_peptideKeyVsTargetDecoyCandidateScoresPntrs.size(),
        featuresArrayTargetVsDecoy.size()
        ); ree;

    const QVector<QPair<CandidateScoresTarget*, CandidateScoresDecoy*>> &candidateScorePairs
                                                        = m_peptideKeyVsTargetDecoyCandidateScoresPntrs.values().toVector();

    QVector<CandidateScores*> candidateScoresesPntrs;
    QVector<FeaturesArray*> featuresArrayPntrs;
    QVector<QPair<FeaturesArrayTargets*, FeaturesArrayDecoys*>> featuresArrayTargetVsDecoyPntrs;
    e = buildProcessBatchPointers(
        candidateScorePairs,
        &featuresArrayTargetVsDecoy,
        &candidateScoresesPntrs,
        &featuresArrayPntrs,
        &featuresArrayTargetVsDecoyPntrs
        ); ree;

    QVector<float> weights;
    e = DiscriminantScoretron::trainLDAClassifier(
            featuresArrayTargetVsDecoyPntrs,
            m_pythiaParameters.threadCount,
            &weights
            ); ree;

    QVector<float> discriminantScores;
    e = DiscriminantScoretron::applyWeights(
        weights,
        m_pythiaParameters.threadCount,
        featuresArrayPntrs,
        &discriminantScores
        ); ree;

    e = setCandidateScoresDiscriminantScore(
        discriminantScores,
        &candidateScoresesPntrs
        ); ree;

    e = QValueSettertron::setQValueForCandidates(
        QValueSettertron::QValueScoreType::DiscriminantScore,
        candidateScoresVecBatchPntrs
        ); ree;

    e = FDRCLassifierNeuralNet::outputFDRResults(
        *candidateScoresVecBatchPntrs,
        true,
        fdrVsCounts
        ); ree;

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

Err PythiaDIAFFWorkflow::recalibrateMs1Points(const QVector<CandidateScores*> &candidateScoresVecBatchPntrsResized) {

    ERR_INIT

    e = ErrorUtils::isFalse(m_targetDecoyCandidatePairScoretron.diaTargetFrames()->isEmpty()); ree;
    e = ErrorUtils::isFalse(m_targetDecoyCandidatePairScoretron.ms1ScanNumberVsScanPoints()->isEmpty()); ree;

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
            m_targetDecoyCandidatePairScoretron.diaTargetFrames(),
            m_targetDecoyCandidatePairScoretron.ms1ScanNumberVsScanPoints()
    ); ree;

    ERR_RETURN
}

namespace {

    Err recalibrationLogic(
        MsCalibratomatic msCalibratomatic,
        const QMap<ScanNumber, ScanPoints*> &scanPoints
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msCalibratomatic.isInitCalMS2()); ree;
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

    e = ErrorUtils::isTrue(m_msCalibratomatic.isInitCalMS2() || m_msCalibratomatic.isInitCalMS1()); ree;
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

        const QList<QMap<ScanNumber, ScanPoints*>> &tandemPointsVec = diaTargetFrames->values();

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

        constexpr int numberOfExperiments = 9;
        double runningPPM = mzPPMStDev;

        PythiaParameters paramsInitial = pythiaParameters;
        paramsInitial.ms2ExtractionWidthPPM = runningPPM;
        pythiaParametersExperiments->push_back(paramsInitial);

        for (int exp = 0; exp < numberOfExperiments; exp++) {

            constexpr double alpha = 1.2;

            PythiaParameters params = pythiaParameters;
            runningPPM *= alpha;
            params.ms2ExtractionWidthPPM = runningPPM;
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
Err PythiaDIAFFWorkflow::optimizeParameters(const QVector<CandidateScores*> &candidateScoresTrainings) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;
    e = ErrorUtils::isNotEmpty(candidateScoresTrainings); ree;

    constexpr int topNMS2IonsOptimization = 12;
    qDebug() << "Using top:" << topNMS2IonsOptimization << "fragments for optimization";

    qDebug() << "Optimization selection fraction" << candidateScoresTrainings.size();

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    for (CandidateScores *cs : candidateScoresTrainings) {
        mzTargetKeyVsTargetDecoyCandidatePointers[cs->targetKey].push_back(cs->targetDecoyCandidatePair);
    }

    QVector<PythiaParameters> pythiaParametersExperiments;
    e = buildDOE(
            m_pythiaParameters,
            m_msCalibratomatic.mzStDevMS2(),
            m_msCalibratomatic.scanTimeStDev(),
            &pythiaParametersExperiments
            ); ree;

    QVector<DOEResult> results;
    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {

        constexpr bool useExtendedScores = true;
        constexpr bool useNeuralNetworkScores = false;

        qDebug() << "ppmTol" << pythiaParams.ms2ExtractionWidthPPM;

        e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(pythiaParams); ree;
        qDebug() << "STarting opt";

        QMap<MzTargetKey, TurboXIC*> nullToBuildTurboXICInParallelLoop;

        constexpr float minPeakCountOptimization = 2.9;
        constexpr int topNMS2Ions = 12;
        m_candidateScores.clear();
        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2Ions,
                m_msCalibratomatic,
                minPeakCountOptimization,
                nullToBuildTurboXICInParallelLoop,
                &mzTargetKeyVsTargetDecoyCandidatePointers,
                &m_candidateScores
                ); ree

        QVector<CandidateScores*> candidateScoresVecBatchPntrs;
        QMap<int, int> fdrVsCounts;
        e = processBatch(
            useExtendedScores,
            useNeuralNetworkScores,
            &candidateScoresVecBatchPntrs,
            &fdrVsCounts
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

    m_pythiaParameters.ms1ExtractionWidthPPM = m_pythiaParameters.ms2ExtractionWidthPPM;
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

    constexpr bool useExtendedScores = true;
    constexpr bool useNeuralNetworkScores = false;

    const int topNMs2IonsMainAnalysis = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions))
            );

    qDebug() << "Using top:" << topNMs2IonsMainAnalysis << "fragments for main analysis";

    QElapsedTimer et;
    et.start();

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            m_targetDecoyPairPntrs,
            msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos(),
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

    QMap<MzTargetKey, TurboXIC*> nullToBuildTurboXICInParallelLoop;

    constexpr float minPeakCountCalibration = 2.9;
    m_candidateScores.clear();
    e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
            topNMs2IonsMainAnalysis,
            m_msCalibratomatic,
            minPeakCountCalibration,
            nullToBuildTurboXICInParallelLoop,
            &mzTargetKeyVsTargetDecoyCandidatePointers,
            &m_candidateScores
            ); ree
    qDebug() << "Targets scored" << et.restart() << "mSec";

    QVector<CandidateScores*> candidateScoresVecBatchPntrs;
    QMap<int, int> fdrVsCounts;
    e = processBatch(
        useExtendedScores,
        useNeuralNetworkScores,
        &candidateScoresVecBatchPntrs,
        &fdrVsCounts
        ); ree;

    qDebug() << "Targets resulted" << et.restart() << "mSec";

    constexpr double fdrThresholdOnePercent = 0.01;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
            m_candidateScores,
            fdrThresholdOnePercent,
            targetCountBelowFDRThresholdOnePercent
    ); ree;
    qDebug() << "Targets counted" << et.restart() << "mSec";

    sortCandidatePointersDiscScoreDesc(&candidateScoresVecBatchPntrs);
    qDebug() << "Targets sorted" << et.restart() << "mSec";

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
            int threadCount,
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

        constexpr int baggingSize = 8;
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

        // explore the overfitting w/ the frequencies, RBF kernel for LDA, training on all candidates instead of top integration
        // msreader using turboxic as the point storage

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
            m_pythiaParameters.threadCount,
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

Err PythiaDIAFFWorkflow::honeIRTAndMassCalibration(
    const QVector<MzTargetKey> &mzTargetKeysToRecal,
    QVector<CandidateScores*> *candidateScoresVecScoredPntrs,
    int topNCandidates,
    int topCandidatesMass
    ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScoresVecScoredPntrs->isEmpty()); ree;

    sortCandidatePointersDiscScoreDesc(candidateScoresVecScoredPntrs);

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

    for (const CandidateScores *cs : candidateScoresVecBatchPntrsResized) {

        if (m_entered.value(cs->targetDecoyCandidatePair)) {
            continue;
        }

        m_targetDecoyCandidatePairsTopScores.push_back(cs->targetDecoyCandidatePair);
        m_entered.insert(cs->targetDecoyCandidatePair, true);
    }

    e = m_msCalibratomatic.buildRTMapper(msCalibrationReaderRows); ree;
    qDebug() << "----- scanTimeWindowStDev x" << S_GLOBAL_SETTINGS.STDEV_MULTIPLIER
             <<":" << m_msCalibratomatic.scanTimeStDev(S_GLOBAL_SETTINGS.STDEV_MULTIPLIER);

    constexpr int ms2MassRecalCountMin = 200;
    if (topCandidatesMass > ms2MassRecalCountMin) {

        msCalibrationReaderRows.resize(topCandidatesMass);

        e = recalibrateMs1Points(candidateScoresVecBatchPntrsResized); ree;

        e = m_msCalibratomatic.setMassCalibrationCoeffs(
            msCalibrationReaderRows,
            MSLevelEnum::MS2
            ); ree;

        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames = m_targetDecoyCandidatePairScoretron.diaTargetFrames();
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> selectDIATargetFrames;
        for (auto it = diaTargetFrames->begin(); it != diaTargetFrames->end(); ++it) {

            const MzTargetKey &mzTargetKey = it.key();
            const QMap<ScanNumber, ScanPoints*> &spp = it.value();

            if (!mzTargetKeysToRecal.contains(mzTargetKey)) {
                continue;
            }
            selectDIATargetFrames.insert(mzTargetKey, spp);
        }

        e = recalibrateMzVals(
                MSLevelEnum::MS2,
                &selectDIATargetFrames,
                m_targetDecoyCandidatePairScoretron.ms1ScanNumberVsScanPoints()
        ); ree;
    }

    ERR_RETURN
}
