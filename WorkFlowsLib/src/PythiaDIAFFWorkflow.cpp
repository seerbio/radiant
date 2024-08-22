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
#include "InterferingCandidatesEliminatomatic.h"
#include "MsFrame.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
#include "PeptideStringWithMods.h"
#include "QuanFileBuilder.h"
#include "QuanTransitionRefinertron.h"
#include "QValueSettertron.h"
#include "SequenceSubstringSearchomatic.h"
#include "SpectrumCentricMzTargetFrameSearch.h"
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

    PeptideSequenceWithModsChargeAndTargetKey buildPeptideSequenceWithModsChargeAndTargetKey(const CandidateScores* candidateScores) {
        return candidateScores->targetDecoyCandidatePair->peptideStringWithMods()
             + "|"
             + QString::number(candidateScores->targetDecoyCandidatePair->charge())
             + "|"
             + candidateScores->targetKey;
    }

    Err filterScoredCandidatesTo60PercentFDR(
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
            if (constexpr double fdrThreshold = 0.6; csp->qValue >= fdrThreshold && !csp->isDecoy) {
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
    MsReaderPointerAcc msReaderPointerAcc;

    e = msReaderPointerAcc.openFile(msDataFilePath); ree;
    msReaderPointerAcc.ptr->printSize();

    e = m_targetDecoyCandidatePairScoretron.init(
            m_pythiaParameters,
            &msReaderPointerAcc
            ); ree;

    e = buildCalibration(&msReaderPointerAcc); ree;

    if (m_msCalibratomatic.isInitRT()) {
        e = optimizeParameters(&msReaderPointerAcc); ree;
        // m_pythiaParameters.ms2ExtractionWidthPPM = 6.75;
        // m_targetDecoyCandidatePairScoretron.setPythiaParameters(m_pythiaParameters);
    }

    int targetCountBelowFDRThreshold;
    e = mainAnalysis(
        &msReaderPointerAcc,
        &targetCountBelowFDRThreshold
        ); ree;

    QVector<CandidateScores*> candidateScoresTargetsAndDecoys;
    e = buildCandidateScoresPtrs(&candidateScoresTargetsAndDecoys); ree;

    QVector<CandidateScores*> candidateScoresTargetsAndDecoys50PercentFDRFiltered;
    e = filterScoredCandidatesTo60PercentFDR(
        m_pythiaParameters.minMs2FragCount,
        &candidateScoresTargetsAndDecoys,
        &candidateScoresTargetsAndDecoys50PercentFDRFiltered
        ); ree;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Analyzing" << candidateScoresTargetsAndDecoys50PercentFDRFiltered.size() << "for filtering";

    // qDebug() << "Starting spectrum centric search";
    // QElapsedTimer et;
    // et.start();
    // e = spectrumCentricSearch(
    //     candidateScoresTargetsAndDecoys50PercentFDRFiltered,
    //     m_msCalibratomatic,
    //     &msReaderPointerAcc
    //     ); ree;
    // qDebug() << "Spectrum centric searched finsihed in" << et.elapsed() << "mSec";

    e = populateAltIdTargetKeys(&candidateScoresTargetsAndDecoys50PercentFDRFiltered); ree;

    QVector<CandidateScores*> candidateScoreClassifierPntrs;
    e = applyNeuralNetClassifier(
            candidateScoresTargetsAndDecoys50PercentFDRFiltered,
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
        m_targetDecoyCandidatePairScoretron.mzTargetKeyVsMsFramePntr(),
        quanFilePath,
        static_cast<float>(m_pythiaParameters.ms2ExtractionWidthPPM)
        ); ree;

    ERR_RETURN
}

namespace {

    struct TurboXICLoadInput {
        QMap<ScanNumber, ScanPoints*> scanNumberVsScanPoints;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        MzTargetKey mzTargetKey;
    };

    std::tuple<Err, MzTargetKey, TurboXIC*> buildTurboXICLogic(const TurboXICLoadInput &turboXicLoadInput) {

        ERR_INIT

        MsFrame msFrame;
        e = msFrame.init(
            turboXicLoadInput.scanNumberVsScanPoints,
            turboXicLoadInput.scanNumberVsScanTime
            ); rtee;

        auto turboXic = new TurboXIC();
        e = turboXic->init(msFrame.frameIndexVsScanPoints()); rtee;

        return {e, turboXicLoadInput.mzTargetKey, turboXic};
    }

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

        QVector<TurboXICLoadInput> turboXicLoadInputs;
        for (const MsScanInfo &msScanInfo : uniqueMsScanInfos) {

            const MzTargetKey mzTargetKey = msScanInfo.targetKey();
            e = ErrorUtils::contains(mzTargetKey, diaTargetFrames); ree;

            TurboXICLoadInput turboXicLoadInput;
            turboXicLoadInput.scanNumberVsScanPoints = diaTargetFrames.value(mzTargetKey);
            turboXicLoadInput.scanNumberVsScanTime = scanNumberVsScanTime;
            turboXicLoadInput.mzTargetKey = mzTargetKey;

            turboXicLoadInputs.push_back(turboXicLoadInput);
        }

        QFuture<std::tuple<Err, MzTargetKey, TurboXIC*>> futures = QtConcurrent::mapped(
            turboXicLoadInputs,
            buildTurboXICLogic
            );
        futures.waitForFinished();

        for (const std::tuple<Err, MzTargetKey, TurboXIC*> &result : futures) {
            e = std::get<0>(result); ree;
            mzTargetKeyVsTurboXicPntr->insert(std::get<1>(result), std::get<2>(result));
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
            int verbosity,
            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(_candidateScores); ree;

        QVector<CandidateScores*> candidateScoresFiltered = _candidateScores;
        filterDuplicateCandidateScoresByDiscriminantScore(&candidateScoresFiltered);

        if(verbosity > 0) {
            qDebug() << _candidateScores.size() << "Found for recalibartion";
            qDebug() << candidateScoresFiltered.size() << "Found for recalibartion after duplicates filtered";
        }


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

        constexpr double cosineSimSumMS1Min = 0.9;
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

    Err buildUniqueMsScanInfosForProcessing(
        const QVector<MsScanInfo> &uniqueMsScanInfos,
        int numberOfUniqueScanInfosForCalibration,
        int offset,
        QVector<MsScanInfo> *uniqueMsScanInfosCalibration
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(uniqueMsScanInfos); ree;
        e = ErrorUtils::isAboveThreshold(numberOfUniqueScanInfosForCalibration, 8, ErrorUtilsParam::IncludeThreshold); ree;

        const double halfSizeScanInfos = uniqueMsScanInfos.size() / 2.0;
        const int incrementSize = std::max(1, static_cast<int>(halfSizeScanInfos / numberOfUniqueScanInfosForCalibration));

        for (int i = 0; i < uniqueMsScanInfos.size(); i += incrementSize) {

            const int distanceFromCenterRight
                    = std::min(static_cast<int>(halfSizeScanInfos) + i + offset, uniqueMsScanInfos.size() - 1);

            const int distanceFromCenterLeft = std::max(static_cast<int>(halfSizeScanInfos) - i - offset, 0);

            if (i == 0) {
                uniqueMsScanInfosCalibration->push_back(uniqueMsScanInfos.at(distanceFromCenterRight));
                continue;
            }

            uniqueMsScanInfosCalibration->push_back(uniqueMsScanInfos.at(distanceFromCenterRight));
            uniqueMsScanInfosCalibration->push_back(uniqueMsScanInfos.at(distanceFromCenterLeft));

            if (uniqueMsScanInfosCalibration->size() >= numberOfUniqueScanInfosForCalibration) {
                break;
            }

        }

        uniqueMsScanInfosCalibration->resize(std::min(
            numberOfUniqueScanInfosForCalibration,
            uniqueMsScanInfosCalibration->size()
            ));

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::buildCalibration(const MsReaderPointerAcc *msReaderPointerAcc) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    constexpr int topNMS2IonsCalibration = 6;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Using top:" << topNMS2IonsCalibration << "fragments for calibration";

    const auto sizePerTranche = static_cast<double>(m_pythiaParameters.trancheSizeMax);
    const int numberOfTranches = std::max(static_cast<int>(m_targetDecoyPairPntrs.size() / sizePerTranche), 1);
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "Number of tranches for calibration:" << numberOfTranches
             << "target count:" << m_targetDecoyPairPntrs.size()
             << "sizePerTranche:" << static_cast<int>(sizePerTranche);

    const QVector<MsScanInfo> uniqueMsScanInfos = msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    constexpr int maxUniqueScanInfosTrainingCount = 16;
    constexpr int offset = 0;
    QVector<MsScanInfo> uniqueMsScanInfosCalibration;
    e = buildUniqueMsScanInfosForProcessing(
        uniqueMsScanInfos,
        maxUniqueScanInfosTrainingCount,
        offset,
        &uniqueMsScanInfosCalibration
        ); ree;

    QMap<MzTargetKey, TurboXIC*> mzTargetKeyVsTurboXicPntrs;
    e = buildMzTargetKeyVsTurboXicPntrs(
        uniqueMsScanInfosCalibration,
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

        QElapsedTimer etBatch;
        etBatch.start();

        QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsBatch = m_targetDecoyCandidatePairsTopScores;
        targetDecoyCandidatePairsBatch.append(tdcp);

        QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
        e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
                targetDecoyCandidatePairsBatch,
                uniqueMsScanInfosCalibration,
                &mzTargetKeyVsTargetDecoyCandidatePointers
                ); ree;

        constexpr float minPeakCountCalibration = 4.9;
        m_candidateScores.clear();
        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2IonsCalibration,
                m_msCalibratomatic,
                minPeakCountCalibration,
                maxUniqueScanInfosTrainingCount,
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
        constexpr int fdrKeyMassCalMS2 = 2;
        constexpr int fdrKeyMassCalMS1 = 5;

        const QVector<MzTargetKey> targetKeysToRecalibrateWhileBuildingCalibration = mzTargetKeyVsTurboXicPntrs.keys().toVector();
        e = honeIRTAndMassCalibration(
            targetKeysToRecalibrateWhileBuildingCalibration,
            &candidateScoresVecBatchPntrs,
            fdrVsCounts.value(fdrKey),
            fdrVsCounts.value(fdrKeyMassCalMS2)
            ); ree;

        scanTimeStDevs.push_back(m_msCalibratomatic.scanTimeStDev());
        ms2PPMStDevs.push_back(m_msCalibratomatic.mzStDevMS2());

        QString fdrString;
        e = FDRCLassifierNeuralNet::outPutFDRCounts(fdrVsCounts, &fdrString); ree;

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                << "Processed batch" << ++batchCounter
                << "of"
                << targetDecoyCandidatePointersTranched.size()
                << etBatch.elapsed()
                << "mSec"
                << qPrintable(fdrString);

        if (fdrVsCounts.value(fdrKey) > m_pythiaParameters.calibrationTrainingVolume
            || &tdcp == targetDecoyCandidatePointersTranched.cend() - 1
            ) {

            constexpr int fdrKey50PercentFDR = 50;
            candidateScoresVecBatchPntrs.resize(std::min(candidateScoresVecBatchPntrs.size(), fdrVsCounts.value(fdrKey50PercentFDR)));

            QVector<CandidateScores*> candidateScoresVecBatchPntrsRecal = candidateScoresVecBatchPntrs;
            candidateScoresVecBatchPntrsRecal.resize(std::min(candidateScoresVecBatchPntrsRecal.size(), fdrVsCounts.value(fdrKeyMassCalMS2)));

            QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
            e = buildMsCalibrationReaderRows(
                    MSLevelEnum::MS2,
                    candidateScoresVecBatchPntrsRecal,
                    m_pythiaParameters.verbosity,
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

            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                     << "ScanTimeWindow Mean|Median|Min"
                     << MathUtils::mean(scanTimeStDevs)
                     << MathUtils::median(scanTimeStDevs)
                     << *std::min({scanTimeStDevs.begin(), scanTimeStDevs.end()});

            qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Ms2 ppm Mean|Median" << MathUtils::mean(ms2PPMStDevs) << MathUtils::median(ms2PPMStDevs);

            m_msCalibratomatic.setScanTimeStDev(scanTimeStDevs.front());
            m_msCalibratomatic.setMzStDevMS2(MathUtils::mean(ms2PPMStDevs));

            e = recalibrateMzVals(
                    MSLevelEnum::MS2,
                    m_targetDecoyCandidatePairScoretron.diaTargetFrames(),
                    m_targetDecoyCandidatePairScoretron.ms1ScanNumberVsScanPoints()
                    ); ree;

            candidateScoresVecBatchPntrsRecal = candidateScoresVecBatchPntrs;
            candidateScoresVecBatchPntrsRecal.resize(std::min(candidateScoresVecBatchPntrsRecal.size(), fdrVsCounts.value(fdrKeyMassCalMS1)));
            filterMs1CandidateRowsByCorr(&candidateScoresVecBatchPntrsRecal);
            constexpr int recalibrationPointCountMin = 400;
            qDebug() << candidateScoresVecBatchPntrsRecal.size() << "found for MS1 Recalibration";
            if (candidateScoresVecBatchPntrsRecal.size() < recalibrationPointCountMin) {
                qWarning() << "Skipping MS1 recalibration.  Not enough points found";
                ERR_RETURN
            }

            QVector<MsCalibarationReaderRow> msCalibrationReaderRowsMS1;
            e = buildMsCalibrationReaderRows(
                    MSLevelEnum::MS1,
                    candidateScoresVecBatchPntrsRecal,
                    m_pythiaParameters.verbosity,
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

            e = m_targetDecoyCandidatePairScoretron.reloadTurboXICMS1(); ree;

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

    const auto terminatorLogic = [&](CandidateScores *cs) {
        return cs->featuresArray[CandidateScores::Features::CosineSimSum100] < m_pythiaParameters.minMs2FragCount;
    };
    const auto terminator = std::remove_if(candidateScoresVecBatchPntrs->begin(), candidateScoresVecBatchPntrs->end(), terminatorLogic);
    candidateScoresVecBatchPntrs->erase(terminator, candidateScoresVecBatchPntrs->end());

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
            m_pythiaParameters.verbosity,
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
        m_pythiaParameters.verbosity,
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

        const float precursorMzMin = msScanInfo.precursorTargetMz - (msScanInfo.isoWindowLower + m_pythiaParameters.precursorExtractionWindowThomsons);
        const float precursorMzMax = msScanInfo.precursorTargetMz + (msScanInfo.isoWindowLower + m_pythiaParameters.precursorExtractionWindowThomsons);

        QVector<TargetDecoyCandidatePair*> targetDecoyCandidatePairsFiltered = targetDecoyCandidatePairs;
        filterTargetDecoyPairPointersByPrecursorMzRange(precursorMzMin, precursorMzMax, &targetDecoyCandidatePairsFiltered);

        mzTargetKeyVsTargetDecoyCandidatePointers->insert(msScanInfo.targetKey(), targetDecoyCandidatePairsFiltered);

        if (m_pythiaParameters.verbosity > 0) {
            qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyCandidatePairsFiltered.size() << "targetDecoyPairs found";
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

        e = ErrorUtils::isTrue(msCalibratomatic.isInitCalMS2()); ree;
        e = ErrorUtils::isNotEmpty(scanPoints); ree;

        e = msCalibratomatic.recalibrateScanPoints(
            MSLevelEnum::MS2,
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

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Recalibrating mz vals";

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

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Points recalibrated in" << et.elapsed() << "mSec";

    ERR_RETURN
}

namespace {

    struct DOEResult {
        double ppm = -1.0;
        int fdrCount = -1;
    };

    Err buildDOE(
            const PythiaParameters &pythiaParameters,
            double mzPPMStDev,
            double scanTimeStDev,
            int verbosity,
            QVector<PythiaParameters> *pythiaParametersExperiments
            ) {

        ERR_INIT

        e = ErrorUtils::isTrue(mzPPMStDev > 0.0); ree;
        e = ErrorUtils::isTrue(scanTimeStDev > 0.0); ree;
        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;

        const QVector<double> ppmList = {
            5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 23, 26, 30, 35, 40, 50
        };

        for (double ppm : ppmList ) {
            PythiaParameters params = pythiaParameters;
            params.ms2ExtractionWidthPPM = ppm;
            pythiaParametersExperiments->push_back(params);
        }

        if (verbosity > 0) {
            qDebug() << "Testing PPM Values";
            for (const PythiaParameters &pp : *pythiaParametersExperiments) {
                qDebug() << "ppmTol" << pp.ms2ExtractionWidthPPM;
            }
        }


        ERR_RETURN
    }

    Err getTopFrequencyParameters(
            const QVector<DOEResult> &results,
            int verbosity,
            double *ppmSetting
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(results); ree;

        Eigen::MatrixX<double> xyMat(results.size() + 1, 2);
        xyMat.setZero();
        for (int row = 0; row < results.size(); row++) {
            const DOEResult &doeResult = results.at(row);
            xyMat.coeffRef(row + 1, 0) = doeResult.ppm;
            xyMat.coeffRef(row + 1, 1) = static_cast<double>(doeResult.fdrCount);
            if (verbosity > 0) {
                qDebug() << doeResult.ppm << doeResult.fdrCount << xyMat.coeff(row, 1);

            }
        }

        constexpr int polynomialOrder = 2;

        QVector<double> coeffs;
        EigenUtils::fitPolynomialQRDecomposition(xyMat, polynomialOrder, &coeffs);

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "PPM Coeffs" << coeffs;

        QVector<double> xPoints = {results.first().ppm};
        while (xPoints.back() < results.back().ppm) {
            constexpr double incrementVal = 0.25;
            xPoints.push_back(xPoints.back() + incrementVal);
        }

        QVector<double> yPoints;
        for (double x : xPoints) {
            double y = 0.0;
            for (int i = 0; i < coeffs.size(); i++) {
                y += coeffs.at(i) * std::pow(x, i);
            }
            if (verbosity > 0) {
                qDebug() << x << y;

            }
            yPoints.push_back(y);
        }

        *ppmSetting = xPoints.at(std::max_element(yPoints.begin(), yPoints.end()) - yPoints.begin());

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::optimizeParameters(MsReaderPointerAcc *msReaderPointerAcc) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairScoretron.isInit()); ree;

    constexpr int topNMS2IonsOptimization = 12;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Using top:" << topNMS2IonsOptimization << "fragments for optimization";

    constexpr int optimizationMultiplicationFactor = 6;
    const auto sizePerTranche = static_cast<double>(m_pythiaParameters.trancheSizeMax * optimizationMultiplicationFactor);
    const int numberOfTranches = std::max(static_cast<int>(m_targetDecoyPairPntrs.size() / sizePerTranche), 1);
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
             << "target count:" << m_targetDecoyPairPntrs.size()
             << "sizePerTranche:" << static_cast<int>(sizePerTranche);

    const QVector<MsScanInfo> uniqueMsScanInfos = msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    constexpr int maxUniqueScanInfosTrainingCount = 16;
    constexpr int offset = 1;
    QVector<MsScanInfo> uniqueMsScanInfosCalibration;
    e = buildUniqueMsScanInfosForProcessing(
        uniqueMsScanInfos,
        maxUniqueScanInfosTrainingCount,
        offset,
        &uniqueMsScanInfosCalibration
        ); ree;

    QMap<MzTargetKey, TurboXIC*> mzTargetKeyVsTurboXicPntrs;
    e = buildMzTargetKeyVsTurboXicPntrs(
        uniqueMsScanInfosCalibration,
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

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            targetDecoyCandidatePointersTranched.first(),
            uniqueMsScanInfosCalibration,
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

    QVector<PythiaParameters> pythiaParametersExperiments;
    e = buildDOE(
            m_pythiaParameters,
            m_msCalibratomatic.mzStDevMS2(),
            m_msCalibratomatic.scanTimeStDev(),
            m_pythiaParameters.verbosity,
            &pythiaParametersExperiments
            ); ree;

    int bestResultCount = 0;
    int countSinceLastHigh = 0;
    QVector<DOEResult> results;
    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {

        constexpr bool useExtendedScores = true;
        constexpr bool useNeuralNetworkScores = false;

        e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(pythiaParams); ree;

        constexpr float minPeakCountOptimization = 2.9;
        constexpr int topNMS2Ions = 12;
        m_candidateScores.clear();
        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2Ions,
                m_msCalibratomatic,
                minPeakCountOptimization,
                maxUniqueScanInfosTrainingCount,
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

        QString fdrString;
        e = FDRCLassifierNeuralNet::outPutFDRCounts(fdrVsCounts, &fdrString); ree;

        qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
                 << "ppmTol"
                 << pythiaParams.ms2ExtractionWidthPPM
                 << "Finished"
                 << qPrintable(fdrString);

        constexpr double fdrThresholdAccuracyOptimization = 0.1;
        int targetCountAboveFDRQValueThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                m_candidateScores,
                fdrThresholdAccuracyOptimization,
                &targetCountAboveFDRQValueThreshold
                ); ree;

        DOEResult res;
        res.ppm = pythiaParams.ms2ExtractionWidthPPM;
        res.fdrCount = targetCountAboveFDRQValueThreshold;
        results.push_back(res);

        if (res.fdrCount >= bestResultCount) {
            bestResultCount = res.fdrCount;
            countSinceLastHigh = 0;
            continue;
        }

        constexpr int maxCountsSinceLastHigh = 3;
        if (++countSinceLastHigh >= maxCountsSinceLastHigh) {
            break;
        }
    }

    e = getTopFrequencyParameters(
            results,
            m_pythiaParameters.verbosity,
            &m_pythiaParameters.ms2ExtractionWidthPPM
            ); ree;

    m_pythiaParameters.ms1ExtractionWidthPPM = m_pythiaParameters.ms2ExtractionWidthPPM;
    e = m_targetDecoyCandidatePairScoretron.setPythiaParameters(m_pythiaParameters); ree;

    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Optimal ppm setting:" << m_pythiaParameters.ms2ExtractionWidthPPM;
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Optimal scanTimeWindow setting:" << m_msCalibratomatic.scanTimeStDev(m_pythiaParameters.scanTimeWindowStDevs) << "minutes";

    for (TurboXIC* turboXic : mzTargetKeyVsTurboXicPntrs) {delete turboXic;}

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
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed())
            << "Using top:"
            << topNMs2IonsMainAnalysis
            << "fragments for main analysis";

    QElapsedTimer et;
    et.start();

    const QVector<MsScanInfo> uniqueMsScanInfos = msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos();

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            m_targetDecoyPairPntrs,
            uniqueMsScanInfos,
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

    QMap<MzTargetKey, TurboXIC*> nullToBuildTurboXICInParallelLoop;

    const int threadCount = std::min(uniqueMsScanInfos.size(), m_pythiaParameters.threadCount);

    constexpr float minPeakCount = 2.9;
    m_candidateScores.clear();
    e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
            topNMs2IonsMainAnalysis,
            m_msCalibratomatic,
            minPeakCount,
            threadCount,
            nullToBuildTurboXICInParallelLoop,
            &mzTargetKeyVsTargetDecoyCandidatePointers,
            &m_candidateScores
            ); ree
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Targets scored" << et.restart() << "mSec";

    QVector<CandidateScores*> candidateScoresVecBatchPntrs;
    QMap<int, int> fdrVsCounts;
    e = processBatch(
        useExtendedScores,
        useNeuralNetworkScores,
        &candidateScoresVecBatchPntrs,
        &fdrVsCounts
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

    sortCandidatePointersDiscScoreDesc(&candidateScoresVecBatchPntrs);
    qDebug() << qPrintable(S_GLOBAL_TIMER.elapsed()) << "Targets sorted" << et.restart() << "mSec";

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
            const QVector<CandidateScores*> &candidateScoresTargetsAndDecoysFDRFiltered,
            QVector<KarnnNNTarget> *karnnNNTargetsNorm
    ){

        ERR_INIT

        e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoysFDRFiltered); ree;

        QVector<KarnnNNTarget> karnnNNTargets;
        karnnNNTargets.reserve(candidateScoresTargetsAndDecoysFDRFiltered.size());
        for (int i = 0; i < candidateScoresTargetsAndDecoysFDRFiltered.size(); i++) {
            const CandidateScores *cs = candidateScoresTargetsAndDecoysFDRFiltered.at(i);
            KarnnNNTarget karnnNnTarget;
            karnnNnTarget.seq = cs->targetDecoyCandidatePair->peptideStringWithMods();
            karnnNnTarget.isDecoy = cs->isDecoy;
            karnnNnTarget.index = i;
            karnnNnTarget.scoreVec = cs->featuresArray;

            for (int j = CandidateScores::Features::ColumnApexIndexRatiosToAnchor1;
                j <= CandidateScores::Features::ColumnApexIndexRatiosToAnchor12; j++) {
                karnnNnTarget.scoreVec[j] = 0.0f;
            }
            for (int j = CandidateScores::Features::MzSearched1;
                    j <= CandidateScores::Features::TheoIntensity12; j++) {
                                karnnNnTarget.scoreVec[j] = 0.0f;
                    }
            for (int j = CandidateScores::Features::MzAccuracy1;
                j <= CandidateScores::Features::MzAccuracy12; j++) {
                        karnnNnTarget.scoreVec[j] = 0.0f;
                }
            for (int j = CandidateScores::Features::ShadowsIntensityRatio1;
                j <= CandidateScores::Features::ShadowsIntensityRatio12; j++) {
                karnnNnTarget.scoreVec[j] = 0.0f;
                }

            // for (int j = CandidateScores::Features::CosineSimToAnchor7;
            //     j <= CandidateScores::Features::CosineSimToAnchor12; j++) {
            //     karnnNnTarget.scoreVec[j] = 0.0f;
            //     }

            // for (int j = CandidateScores::Features::CosineSimToAnchor7;
            //     j <= CandidateScores::Features::CosineSimToAnchor12; j++) {
            //     karnnNnTarget.scoreVec[j] = 0.0f;
            //     }




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

    e = InterferingCandidatesEliminatomatic::removeInterferingCandidates(
            m_pythiaParameters.ionsSharedToReject,
            m_pythiaParameters.mzMinMS2,
            m_pythiaParameters.mzMaxMS2,
            &candidateScoresVecBatchPntrsResized
            ); ree;

    if (m_pythiaParameters.verbosity > 0) {
        qDebug() << "Using" << candidateScoresVecBatchPntrsResized.size() << "for iRT Estimation";
        qDebug() << candidateScoresVecBatchPntrsResized.size() << "after filtering";
    }

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = buildMsCalibrationReaderRows(
            MSLevelEnum::MS2,
            candidateScoresVecBatchPntrsResized,
            m_pythiaParameters.verbosity,
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

    if (m_pythiaParameters.verbosity > 0) {
        qDebug() << "----- scanTimeWindowStDev x" << m_pythiaParameters.scanTimeWindowStDevs
                 <<":" << m_msCalibratomatic.scanTimeStDev(m_pythiaParameters.scanTimeWindowStDevs);
    }

    constexpr int ms2MassRecalCountMin = 200;
    if (topCandidatesMass > ms2MassRecalCountMin) {

        msCalibrationReaderRows.resize(topCandidatesMass);

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

    }

    ERR_RETURN
}

namespace {

    struct SpectrumCentricParallelInput {
        MzTargetKey mzTargetKey;
        QVector<CandidateScores*> candidateScoresPntrs;
        QMap<ScanNumber, ScanPoints*> diaTargetFrame;
        PythiaParameters pythiaParameters;
        MsCalibratomatic msCalibratomatic;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    };

    Err buildSpectrumCentricParallelInput(
        const QMap<MzTargetKey, QVector<CandidateScores*>> &mzTargetKeyVsCandidateScoresPointers,
        const QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> &diaTargetFrames,
        const MsCalibratomatic &msCalibratomatic,
        const PythiaParameters &pythiaParameters,
        const QMap<ScanNumber, ScanTime> &scanNumberVsScanTime,
        QVector<SpectrumCentricParallelInput> *spectrumCentricParallelInputs
        ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(mzTargetKeyVsCandidateScoresPointers); ree;
        e = ErrorUtils::isNotEmpty(diaTargetFrames); ree;

        for(auto it = mzTargetKeyVsCandidateScoresPointers.begin(); it != mzTargetKeyVsCandidateScoresPointers.end(); ++it) {
            SpectrumCentricParallelInput spectrumCentricParallelInput;
            spectrumCentricParallelInput.mzTargetKey = it.key();
            spectrumCentricParallelInput.candidateScoresPntrs = it.value();
            spectrumCentricParallelInput.msCalibratomatic = msCalibratomatic;
            spectrumCentricParallelInput.pythiaParameters = pythiaParameters;
            spectrumCentricParallelInput.scanNumberVsScanTime = scanNumberVsScanTime;

            e = ErrorUtils::contains(spectrumCentricParallelInput.mzTargetKey, diaTargetFrames); ree;
            spectrumCentricParallelInput.diaTargetFrame = diaTargetFrames.value(spectrumCentricParallelInput.mzTargetKey);
            spectrumCentricParallelInputs->push_back(spectrumCentricParallelInput);
        }

        ERR_RETURN
    }

    QPair<Err, QVector<QPair<CandidateScores*, DeconvolvotronResult>>> spectrumCentricParallelLogic(const SpectrumCentricParallelInput &scpi) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scpi.candidateScoresPntrs); rree;
        e = ErrorUtils::isNotEmpty(scpi.diaTargetFrame); rree;
        e = ErrorUtils::isNotEmpty(scpi.mzTargetKey); rree;

        QElapsedTimer et;
        et.start();

        SpectrumCentricMzTargetFrameSearch spectrumCentricMzTargetFrameSearch;
        e = spectrumCentricMzTargetFrameSearch.init(
            scpi.pythiaParameters,
            scpi.msCalibratomatic,
            scpi.diaTargetFrame,
            scpi.candidateScoresPntrs,
            scpi.scanNumberVsScanTime
            ); rree;

        QVector<QPair<CandidateScores*, DeconvolvotronResult>> candidateScoresPntrsVsScore;
        e = spectrumCentricMzTargetFrameSearch.assignIdsToScans(&candidateScoresPntrsVsScore); rree;

        qDebug() << "Processed MzTargetKey" << scpi.mzTargetKey << et.elapsed() << "mSec";

        return {e, candidateScoresPntrsVsScore};
    }

}//namespace
Err PythiaDIAFFWorkflow::spectrumCentricSearch(
    const QVector<CandidateScores*> &candidateScoresPntrs,
    const MsCalibratomatic &msCalibratomatic,
    const MsReaderPointerAcc *msReaderPointerAcc
    ) const {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_targetDecoyPairPntrs); ree;

    QMap<MzTargetKey, QVector<CandidateScores*>> mzTargetKeyVsCandidateScoresPntrs;
    for (CandidateScores *cs : candidateScoresPntrs) {
        mzTargetKeyVsCandidateScoresPntrs[cs->targetKey].push_back(cs);
    }

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    e = msReaderPointerAcc->ptr->collateMS2MzTargetFrames(&diaTargetFrames); ree;

    QVector<SpectrumCentricParallelInput> spectrumCentricParallelInputs;
    e = buildSpectrumCentricParallelInput(
        mzTargetKeyVsCandidateScoresPntrs,
        diaTargetFrames,
        msCalibratomatic,
        m_pythiaParameters,
        msReaderPointerAcc->ptr->getScanNumberVsScanTime(),
        &spectrumCentricParallelInputs
        ); ree;

#define PARALLEL_DCON
#ifdef PARALLEL_DCON
    QFuture<QPair<Err, QVector<QPair<CandidateScores*, DeconvolvotronResult>>>> futures = QtConcurrent::mapped(
        spectrumCentricParallelInputs,
        spectrumCentricParallelLogic
        ); ree;
    futures.waitForFinished();

    for (const QPair<Err, QVector<QPair<CandidateScores*, DeconvolvotronResult>>> &res : futures) {
        e = res.first; ree;
        const QVector<QPair<CandidateScores*, DeconvolvotronResult>> &prs = res.second;
    }
#else
    for (const SpectrumCentricParallelInput &inp : spectrumCentricParallelInputs) {
        QPair<Err, QVector<QPair<CandidateScores*, DeconvolvotronResult>>> res = spectrumCentricParallelLogic(inp); ree;
        e = res.first; ree;
        qDebug() << res.second.size() << "SDLFKDJSL";
    }
#endif

    ERR_RETURN
}
