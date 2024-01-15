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
#include "ProteinDigestomatic.h"
#include "QValueSettertron.h"
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

    QVector<CandidateScores*> candidateScoresTrainings;
    e = buildCalibration(
            &msReaderPointerAcc,
            &diaTargetFrames,
            &scanNumberVsScanPointsMS1,
            &candidateScoresTrainings
            ); ree;

    if (m_msCalibratomatic.isInit()) {

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

//#define TROUBLESHOOT_MISSING
#ifdef TROUBLESHOOT_MISSING
    for (CandidateScores *cs : candidateScoresTargetsAndDecoys) {
        if (cs->targetDecoyCandidatePair->peptideStringWithMods() == "EAQGNSSAGVEAAEQRPVEDGER" && cs->targetDecoyCandidatePair->charge() == 3) {
            qDebug() << cs->targetDecoyCandidatePair->peptideStringWithMods() << cs->isDecoy;
            qDebug() << cs->featuresArray[CandidateScores::Features::CosineSimSum100] << cs->discriminantScore << cs->classifierScore << cs->qValue;
            qDebug() << cs->featuresArray.mid(CandidateScores::Features::CosineSimToAnchor1, 12);
            qDebug() << "^^^^^^^^^^";
            einfo;
        }
    }
#endif

    QVector<CandidateScores*> candidateScoresTargetsAndDecoys50PercentFDRFiltered;
    e = filterScoredCandidatesTo50PercentFDR(
            &candidateScoresTargetsAndDecoys,
            &candidateScoresTargetsAndDecoys50PercentFDRFiltered
    ); ree;

#ifdef TROUBLESHOOT_MISSING
    for (CandidateScores *cs : candidateScoresTargetsAndDecoys50PercentFDRFiltered) {
        if (cs->targetDecoyCandidatePair->peptideStringWithMods() == "EAQGNSSAGVEAAEQRPVEDGER" && cs->targetDecoyCandidatePair->charge() == 3) {
            qDebug() << cs->targetDecoyCandidatePair->peptideStringWithMods() << cs->isDecoy;
            qDebug() << cs->featuresArray[CandidateScores::Features::CosineSimSum100] << cs->discriminantScore << cs->classifierScore << cs->qValue;
            qDebug() << cs->featuresArray.mid(CandidateScores::Features::CosineSimToAnchor1, 12);
            qDebug() << "^^^^^^^^^^";
            einfo;
        }
    }
#endif

    qDebug() << "Analyzing" << candidateScoresTargetsAndDecoys50PercentFDRFiltered.size() << "for filtering";

    e = InterferingCandidatesEliminatomatic::removeInterferingCandidates(
            m_pythiaParameters.ionsSharedToReject,
            m_pythiaParameters.mzMinMS2,
            m_pythiaParameters.mzMaxMS2,
            &candidateScoresTargetsAndDecoys50PercentFDRFiltered
            ); ree;
    qDebug() << candidateScoresTargetsAndDecoys50PercentFDRFiltered.size() << "after filtering";

    e = populateAltIdTargetKeys(&candidateScoresTargetsAndDecoys50PercentFDRFiltered); ree;

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
            row.scanNumber = cs->scanNumber;
            row.mzSearchedVec = cs->featuresArray.mid(CandidateScores::Features::MzSearched1, top6);
            row.mzFoundMeanVec = cs->featuresArray.mid(CandidateScores::Features::MzFoundMean1, top6);
            row.mzFoundStDevVec = cs->featuresArray.mid(CandidateScores::Features::MzFoundStDev1, top6);
            row.intensityFoundMaxVec = cs->featuresArray.mid(CandidateScores::Features::IntensityFoundMax1, top6);

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

        QVector<CandidateScores*> candidateScoresPntrs;
        e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;

        e = DiscriminantScoretron::setDiscriminantScoreForCandidates(
                useExtendedScores,
                useNeuralNetworkScores,
                &candidateScoresPntrs
        ); ree;

        e = QValueSettertron::setQValueForCandidates(
                QValueSettertron::QValueScoreType::DiscriminantScore,
                &candidateScoresPntrs
                ); ree;

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

        e = InterferingCandidatesEliminatomatic::removeInterferingCandidates(
                m_pythiaParameters.ionsSharedToReject,
                m_pythiaParameters.mzMinMS2,
                m_pythiaParameters.mzMaxMS2,
                &candidateScoresPntrs
                ); ree;
        qDebug() << candidateScoresPntrs.size() << "after filtering";

        QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
        e = buildMsCalibrationReaderRows(candidateScoresPntrs, &msCalibrationReaderRows); ree;
        topScoresOfTranches.clear();
        for (CandidateScores* p : candidateScoresPntrs) {
            topScoresOfTranches[p->targetKey].push_back(p->targetDecoyCandidatePair);
        }

        const int numberOfStDevs = 3;
        e = m_msCalibratomatic.initRtOnly(msCalibrationReaderRows); ree;
        qDebug() << "scanTimeWindowStDev x" << numberOfStDevs <<":" << m_msCalibratomatic.scanTimeStDev(numberOfStDevs);

        const int minMzTrainingSize = 200;
        const QString twoPercenFDRKey = "2";
        if (fdrVsCount.value(twoPercenFDRKey) >= minMzTrainingSize) {

            msCalibrationReaderRows.resize(fdrVsCount.value(twoPercenFDRKey));

//#define TIGHT_IRT_CAL
#ifdef TIGHT_IRT_CAL

            //TODO Test this en masse for stats

                e = m_msCalibratomatic.initRtOnly(msCalibrationReaderRows); ree;
                qDebug() << "scanTimeWindowStDev x" << numberOfStDevs <<":" << m_msCalibratomatic.scanTimeStDev(numberOfStDevs);
#endif
            e = m_msCalibratomatic.initMzOnly(msCalibrationReaderRows); ree;
            e = recalibrateMzVals(
                    diaTargetFrames,
                    scanNumberVsScanTimeMS1
            ); ree;
        }

        if (minTrainingCount >= idealTrainingCountAtGivenFDR) {

            QVector<CandidateScores*> candidateScoresPntrsOpt;
            e = buildCandidateScoresPtrs(&candidateScoresPntrsOpt); ree;
            std::sort(candidateScoresPntrsOpt.rbegin(), candidateScoresPntrsOpt.rend(), [](CandidateScores *l, CandidateScores *r){
                return l->discriminantScore < r->discriminantScore;
            });

            candidateScoresPntrsOpt.resize(idealTrainingCountAtGivenFDR * 2);

            *candidateScoresForTrainings = candidateScoresPntrsOpt;
            qDebug() << "scanTimeWindowStDev x 3:" << m_msCalibratomatic.scanTimeStDev(numberOfStDevs);

//#define WRITE_CALIBRATION_ROWS_TS
#ifdef WRITE_CALIBRATION_ROWS_TS
            qDebug() << "ACHTUNG, ACHTUNG, ACHTUNG!!!! WRITING CAL FILE IN:"; einfo;
            const QString filename = QStringLiteral("/home/anichols/Desktop/Data/ConfigFiles/cal2.prq");
            e = ParquetReader::write(msCalibrationReaderRows, filename); ree;
#endif
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

    Err recalibrationLogic(
        MsCalibratomatic msCalibratomatic,
        const QMap<ScanNumber, ScanPoints*> &scanPoints
        ) {

        ERR_INIT

        e = ErrorUtils::isTrue(msCalibratomatic.isInit()); ree;
        e = ErrorUtils::isNotEmpty(scanPoints); ree;

        e = msCalibratomatic.recalibrateScanPoints(scanPoints); ree;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::recalibrateMzVals(
        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrames,
        QMap<ScanNumber, ScanPoints> *scanNumberVsScanTimeMS1
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_msCalibratomatic.isInit()); ree;

    e = ErrorUtils::isFalse(diaTargetFrames->isEmpty()); ree;
    e = ErrorUtils::isFalse(scanNumberVsScanTimeMS1->isEmpty()); ree;

    qDebug() << "Recalibrating mz vals";

    QElapsedTimer et;
    et.start();

    const QList<MzTargetKey> &diaTargetFrameKeys = diaTargetFrames->keys();

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

    for (const Err &err : futures) {
        e = err; ree;
    }
#else
    for (const MzTargetKey &k: diaTargetFrameKeys) {
        e = recalibrationLogic(m_msCalibratomatic, diaTargetFrames->value(k)); ree;
    }
#endif

    qDebug() << "Recalibrating MS1 mz vals frame";
    e = m_msCalibratomatic.recalibrateScanPoints(scanNumberVsScanTimeMS1); ree;

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
                {2.0, 2.0},
                {2.5,  2.0},
                {3.5,  2.0},
                {4.5, 2.0},
                {5.5, 2.0}
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

        QVector<CandidateScores*> candidateScoresPntrs;
        e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;
        e = DiscriminantScoretron::setDiscriminantScoreForCandidates(
                useExtendedScores,
                useNeuralNetworkScores,
                &candidateScoresPntrs
        ); ree;

        e = QValueSettertron::setQValueForCandidates(
                QValueSettertron::QValueScoreType::DiscriminantScore,
                &candidateScoresPntrs
                ); ree;

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

    e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
            topNMs2IonsMainAnalysis,
            msReaderPointerAcc->ptr->scanTimeMinMax(),
            m_msCalibratomatic,
            &mzTargetKeyVsTargetDecoyCandidatePointers,
            &m_candidateScores
    ); ree;
    qDebug() << "Targets scored" << et.restart() << "mSec";

    QVector<CandidateScores*> candidateScoresPntrs;
    e = buildCandidateScoresPtrs(&candidateScoresPntrs); ree;
    e = DiscriminantScoretron::setDiscriminantScoreForCandidates(
            useExtendedScores,
            useNeuralNetworkScores,
            &candidateScoresPntrs
    ); ree;
    qDebug() << "Targets discriminated" << et.restart() << "mSec";

    e = QValueSettertron::setQValueForCandidates(
            QValueSettertron::QValueScoreType::DiscriminantScore,
            &candidateScoresPntrs
            ); ree;
    qDebug() << "Targets qval'd" << et.restart() << "mSec";

    QMap<QString, int> fdrVsCount;
    e = FDRCLassifierNeuralNet::outputFDRResults(
            m_candidateScores,
            true,
            &fdrVsCount
    ); ree;
    qDebug() << "Targets resulted" << et.restart() << "mSec";

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
            double fdrCutOff,
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
            m_pythiaParameters.reportDecoys,
            m_pythiaParameters.percentFDR,
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
        cs->proteinGroup.replace(",", " ");

    }

    ERR_RETURN
}
