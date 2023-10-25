//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "MsReaderParquet.h"

//#include "CandidateClassifier.h"
#include "ClassifierWeightsManager.h"
//#include "EigenUtils.h"
//#include "ErrorUtils.h"
//#include "FastaFileToPeptidesListWorkFlow.h"
//#include "FastaReader.h"
#include "FDRCLassifierNeuralNet.h"
//#include "FeatureFinderHillBuilder.h"
#include "MathUtils.h"
//#include "MS2ChargeDeconvolvotron.h"
//#include "MS2DataExtractomatic.h"
//#include "MsCalibrationReader.h"
//#include "MsFrameScoretron.h"
//#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
//#include "PeakIntegratomatic.h"
//#include "TandemSpectraDeconvolvotron.h"
#include "TandemFragmentPredictotron.h"
#include "TargetDecoyCandidatePairScoretron.h"
#include "TurboXIC.h"

#include <QtConcurrent/QtConcurrent>


PythiaDIAWorkflow::PythiaDIAWorkflow()
: m_minTopNMs2Ions(6)
, m_byIonsOnly(true)
{}

Err PythiaDIAWorkflow::init(
        const PythiaParameters &pythiaParameters,
        const QString &fragLibUri,
        const QString &fastaUri
        ) {

    ERR_INIT

    pythiaParameters.print();

    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::fileExists(fragLibUri); ree;

    if (!fastaUri.isEmpty()) {
        e = ErrorUtils::fileExists(fastaUri); ree;
        m_fastaUri = fastaUri;
    }

    m_pythiaParameters = pythiaParameters;
    m_fragLibUri = fragLibUri;

    e = m_targetDecoyCandidatePairManager.init(m_pythiaParameters, m_fragLibUri); ree;

    ERR_RETURN
}

namespace{

    Err filterScoreCandidatesByFDR(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            double qValueThreshold,
            bool filterDecoys,
            QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;
        e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

//        scoredCandidatesTargetsFDRThresholded->clear();
//
//        for (const ScoredCandidate &sc : scoredCandidatesAll) {
//
//            if (filterDecoys && sc.isDecoy) {
//                continue;
//            }
//
//            else if (sc.qValue > qValueThreshold) {
//                continue;
//            }
//
//            scoredCandidatesTargetsFDRThresholded->push_back(sc);
//        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::processFile(const QString &_msDataFilePath) {

    ERR_INIT

    m_msScanInfos.clear();
    m_scanTimeMinMax = {-1.0, -1.0};

    QString msDataFilePath = _msDataFilePath;

//#define USE_FILE_CACHING
#ifdef USE_FILE_CACHING
    {
        const QString msDataFilePathCached = msDataFilePath + S_GLOBAL_SETTINGS.DOT_CACHED_FILE_EXTENSION;
        const bool cacheFileExists = QFileInfo::exists(msDataFilePathCached);
        qDebug() << "Using cached msdatafile" << cacheFileExists;

        if (cacheFileExists) {
            msDataFilePath = msDataFilePathCached;
        }
        else {
            MsReaderPointerAcc msReaderPointerAcc;
            e = msReaderPointerAcc.openFile(msDataFilePath); ree;
            e = deisotopeScans(&msReaderPointerAcc); ree;
            e = MsReaderParquet::writeMsReaderToParquet(
                    msDataFilePathCached,
                    msReaderPointerAcc.ptr
            ); ree;
            msDataFilePath = msDataFilePathCached;
        }
    }
#endif

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath); ree;
    m_msScanInfos = msReaderPointerAcc.ptr->getUniqueTandemMsScanInfos();
    m_scanTimeMinMax = msReaderPointerAcc.ptr->scanTimeMinMax();

    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> diaTargetFrame;
    e = msReaderPointerAcc.ptr->collateTandemPrecursorTargetsDIA(
            &diaTargetFrame
            ); ree;

    const int msLevel = 1;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
    e = msReaderPointerAcc.ptr->getScanPoints(msLevel, &scanNumberVsScanTimeMS1); ree;
    
    TargetDecoyCandidatePairScoretron targetDecoyCandidatePairScoretron;
    e = targetDecoyCandidatePairScoretron.init(
            m_pythiaParameters,
            scanNumberVsScanTimeMS1,
            &msReaderPointerAcc,
            &diaTargetFrame,
            &m_targetDecoyCandidatePairManager
            );

#ifndef USE_FILE_CACHING
//    if (m_pythiaParameters.deisotopeScans) {
//        e = deisotopeScans(&msReaderPointerAcc); ree;
//    }
#endif

    e = buildCalibration(&targetDecoyCandidatePairScoretron); ree;

//#define BYPASS_OPTI
#ifndef BYPASS_OPTI
    e = optimizeParameters(&targetDecoyCandidatePairScoretron); ree;
#else
    //Pythia Main Library
//    m_pythiaParameters.ms2ExtractionWidthPPM = 12.2715;
//    m_pythiaParameters.scanTimeWindowMinutes = 1.79397;
//    m_pythiaParameters.cosineSimToAnchorThreshold = 0.9;

    //Entrapment libarary
    m_pythiaParameters.ms2ExtractionWidthPPM = 11.945;
    m_pythiaParameters.scanTimeWindowMinutes = 3.60323;
    m_pythiaParameters.cosineSimToAnchorThreshold = 0.9;
#endif

////#define BYPASS_MAIN_ANALYSIS
//#ifndef BYPASS_MAIN_ANALYSIS
//
//    QVector<ScoredCandidate> scoredCandidatesAll;
//    e = mainAnalysis(
//            &msReaderPointerAcc,
//            &scoredCandidatesAll
//            ); ree;
//
//    const double fdrThreshold = 0.2;
//
//    QVector<ScoredCandidate> scoredCandidatesTargetsFDRThresholded;
//    e = MS2DataExtractomatic::filterScoreCandidatesByFDR(
//            scoredCandidatesAll,
//            fdrThreshold,
//            true,
//            &scoredCandidatesTargetsFDRThresholded
//    ); ree;
//
//    QVector<ScoredCandidate> scoredCandidatesAllUpdated;
//    e = removeInterferingCandidates(
//            &msReaderPointerAcc,
//            scoredCandidatesTargetsFDRThresholded,
//            scoredCandidatesAll,
//            &scoredCandidatesAllUpdated
//            ); ree;
//
//    qDebug() << "scoredCandidatesAll size" << scoredCandidatesAll.size();
//    qDebug() << "ScoredCandidatesAllUpdated size" << scoredCandidatesAllUpdated.size();
//    qDebug() << "ScoredCandidatesThresholded size" << scoredCandidatesTargetsFDRThresholded.size();
//
////#define WRITE_FILTERED_PSM
//#ifdef WRITE_FILTERED_PSM
//    qDebug() << scoredCandidatesAllUpdated.size() << "Candidates written";
//    e = updateProteinGroupAnnotation(
//            "/home/anichols/Downloads/human_plasma_arath_entrapment.fasta", //TODO make this proper input
//            &scoredCandidatesAllUpdated
//    ); ree;
//    e = ParquetReader::write(scoredCandidatesAll, "scoredCandidatesAll.parquet"); ree;
//    e = ParquetReader::write(scoredCandidatesAllUpdated, "scoredCandidatesAllUpdated.parquet"); ree;
//#endif
//
//#else
//    QVector<ScoredCandidate> scoredCandidatesAllUpdated;
//    e = ParquetReader::read(
//            "/home/anichols/Repositories/Builds/PythiaDIACpp/bin/scoredCandidatesAllUpdated.parquet",
//            &scoredCandidatesAllUpdated
//            ); ree;
//
//    QVector<ScoredCandidate> scoredCandidatesAll;
//    e = ParquetReader::read(
//            "/home/anichols/Repositories/Builds/PythiaDIACpp/bin/scoredCandidatesAll.parquet",
//            &scoredCandidatesAll
//    ); ree;
//#endif
//
//#define USE_NEURAL_NET_CLASSIFIER
//#ifdef USE_NEURAL_NET_CLASSIFIER
//    QVector<ScoredCandidate> scoredCandidatesClassifierUpdated;
//    e = applyNeuralNetClassifier(
//            scoredCandidatesAll,
//            scoredCandidatesAllUpdated,
//            msReaderPointerAcc.ptr->scanTimeMinMax(),
//            &scoredCandidatesClassifierUpdated
//            ); ree;
//
//    QVector<ScoredCandidate> scoredCandidatesAllThresholded;
//    filterScoreCandidatesByFDR(
//            scoredCandidatesClassifierUpdated,
//            0.01,
//            true,
//            &scoredCandidatesAllThresholded
//    );
//    scoredCandidatesClassifierUpdated = scoredCandidatesAllThresholded;
//
//    e = updateProteinGroupAnnotation(
//            m_fastaUri,
//            &scoredCandidatesClassifierUpdated
//            ); ree;
//
//    const QString resultsFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
//    e = ParquetReader::write(scoredCandidatesClassifierUpdated, resultsFilePath); ree;
//#else
//    QVector<ScoredCandidate> scoredCandidatesAllThresholded;
//    filterScoreCandidatesByFDR(
//            scoredCandidatesAllUpdated,
//            0.01,
//            true,
//            &scoredCandidatesAllThresholded
//    );
//    scoredCandidatesAllUpdated = scoredCandidatesAllThresholded;
//
//    e = updateProteinGroupAnnotation(
//            "/home/anichols/Downloads/human_plasma_arath_entrapment.fasta", //TODO make this proper input
//            &scoredCandidatesAllUpdated
//            ); ree;
//
//    const QString resultsFilePath = msReaderParquet.filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
//    e = ParquetReader::write(scoredCandidatesAllUpdated, resultsFilePath); ree;
//#endif

    ERR_RETURN
}

namespace {

    QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> deisotopeLogic(const QVector<QPair<ScanNumber, ScanPoints>> &scanPointPairs, double ppmTol) {

        ERR_INIT

//        const QString &chargeModelFilePath
//                = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");
//
//        const QString &monoModelFilePath
//                = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");
//
//        MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
//        e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);

        QVector<QPair<ScanNumber, ScanPoints>> deisotopedScanPoints;
//        for (const QPair<ScanNumber, ScanPoints> &pr : scanPointPairs) {
//
//            if (pr.first % 1000 == 0) {
//                qDebug() << "Deisotoping" << pr.first;
//            }
//
//            ScanPoints scanPointsIterDeisotoped;
//            e = ms2ChargeDeconvolvotron.deisotopeScanPoints(pr.second, &scanPointsIterDeisotoped); rree;
//
//            deisotopedScanPoints.push_back({pr.first, scanPointsIterDeisotoped});
//        }

        return {e, deisotopedScanPoints};
    }

}//namespace
Err PythiaDIAWorkflow::deisotopeScans(MsReaderPointerAcc *msReaderPointerAcc) {

    ERR_INIT

    qDebug() << "Deisotoping start";

    QElapsedTimer et;
    et.start();

    const auto deisotopeLogicBinder = std::bind(
            deisotopeLogic,
            std::placeholders::_1,
            m_pythiaParameters.ms2ExtractionWidthPPM
    );

//    QVector<QVector<QPair<ScanNumber, ScanPoints>>> scanPointsTranched;
//    ParallelUtils::trancheMapForParallelization(
//            msReaderPointerAcc->ptr->getScanPoints(),
//            ParallelUtils::numberOfAvailableSystemProcessors(),
//            &scanPointsTranched
//    );
//
//    QFuture<QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>>> futures
//            = QtConcurrent::mapped(scanPointsTranched, deisotopeLogicBinder);
//    futures.waitForFinished();
//
//    QMap<ScanNumber, ScanPoints> deisotopedScanPoints;
//    for (const QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> &res : futures) {
//
//        e = res.first; ree;
//        for (const QPair<ScanNumber, ScanPoints> &r : res.second) {
//            const ScanNumber scanNumber = r.first;
//            const ScanPoints &scanPoints = r.second;
//            deisotopedScanPoints.insert(scanNumber, scanPoints);
//        }
//    }
//
//    e = msReaderPointerAcc->ptr->setScanPoints(deisotopedScanPoints); ree;

    qDebug() << "Scans deisotoped in" << et.elapsed() << "mSec";

    ERR_RETURN
}

namespace {

    struct TargetDecoyPairTargetKey {
        TargetDecoyCandidatePair* targetDecoyCandidatePair;
        QPair<ScoresTargets, ScoresDecoys> scoresTargetVsScoresDecoys;
        UniqueMsInfoScanKey uniqueMsInfoScanKey;
    };

    Err buildClassifierInput(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            int theoMzIonsSize,
            const QPair<double, double> &scanTimeMinMax,
            QVector<QPair<ScoresTargets, ScoresDecoys>> *targetsScoresVsDecoyScores,
            QVector<TargetDecoyPairTargetKey> *targetDecoyPairTargetKeys
            ) {

        ERR_INIT

        const int theoMzIonsSizeMin = 6;

        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairPntrs); ree;
        e = ErrorUtils::isTrue(scanTimeMinMax.second > scanTimeMinMax.first); ree;
        e = ErrorUtils::isTrue(theoMzIonsSize >= theoMzIonsSizeMin); ree;

        targetsScoresVsDecoyScores->clear();

        for (TargetDecoyCandidatePair *tdc : targetDecoyCandidatePairPntrs) {

            const QList<UniqueMsInfoScanKey> &uniqueInfoScanKeys = tdc->uniqueInfoScanKeyVsScoresTarget()->keys();
            for (const UniqueMsInfoScanKey &key : uniqueInfoScanKeys) {

                const bool decoyContainsTargetKey = tdc->uniqueInfoScanKeyVsScoresDecoy()->contains(key);
                if (!decoyContainsTargetKey) {
                    qDebug() << "Decoy scores not found for" << key << tdc->peptideStringWithMods();
                    qDebug() << "Keys" << uniqueInfoScanKeys << "Mz" << tdc->mz();
                }
                e = ErrorUtils::isTrue(decoyContainsTargetKey); ree;

                const ScoresTargets scoresTargets = FDRCLassifierNeuralNet::buildScoreVector(
                        tdc->uniqueInfoScanKeyVsScoresTarget()->value(key),
                        useExtendedScores,
                        useNeuralNetworkScores,
                        theoMzIonsSize,
                        scanTimeMinMax
                        ); ree;

                const ScoresDecoys scoresDecoys = FDRCLassifierNeuralNet::buildScoreVector(
                        tdc->uniqueInfoScanKeyVsScoresDecoy()->value(key),
                        useExtendedScores,
                        useNeuralNetworkScores,
                        theoMzIonsSize,
                        scanTimeMinMax
                        ); ree;

                TargetDecoyPairTargetKey targetDecoyPairTargetKey;
                targetDecoyPairTargetKey.targetDecoyCandidatePair = tdc;
                targetDecoyPairTargetKey.scoresTargetVsScoresDecoys = {scoresTargets, scoresDecoys};
                targetDecoyPairTargetKey.uniqueMsInfoScanKey = key;

                targetDecoyPairTargetKeys->push_back(targetDecoyPairTargetKey);
                targetsScoresVsDecoyScores->push_back({scoresTargets, scoresDecoys});
            }

        }

        ERR_RETURN
    }

    Err getBestFDRFraction(
            const QMap<QString, int> &fdrResults,
            int minTrainingThreshold,
            double *bestFDRFraction
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(fdrResults); ree;

        const QStringList fdrFractions = {
                "1",
                "2",
                "5",
                "10",
                "20"
        };

        for (const QString &fdrStr : fdrFractions) {
            if (fdrResults.value(fdrStr) > minTrainingThreshold) {

                double fdrPercent;
                e = ErrorUtils::toDouble(fdrStr, &fdrPercent); ree;

                *bestFDRFraction = fdrPercent / 100.0;
                ERR_RETURN
            }
        }

        double fdrPercent;
        e = ErrorUtils::toDouble(fdrResults.lastKey(), &fdrPercent); ree;

        *bestFDRFraction = fdrPercent;

        ERR_RETURN
    }

    Err buildMsCalibrationReaderRows(
            const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointersFDRThresholded,
            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredTargetDecoyPointersFDRThresholded); ree;

        qDebug() << scoredTargetDecoyPointersFDRThresholded.size() << "Found for recalibartion";

        const auto msCalibrationReaderRowsInsertLogic = [](TargetDecoyCandidatePair* tdp){

            CandidateScores *candidateScoresTargetBestDiscriminantScore = tdp->candidateScoresBestDiscriminantScorePtrTarget();

            MsCalibarationReaderRow row;
            row.peptideStringWithMods = tdp->peptideStringWithMods();
            row.iRTPredicted = static_cast<float>(tdp->iRt());
            row.scanTime = candidateScoresTargetBestDiscriminantScore->scanTime;
            row.mzSearchedVec = candidateScoresTargetBestDiscriminantScore->mzSearchedVec;
            row.mzFoundMeanVec = candidateScoresTargetBestDiscriminantScore->mzFoundMeanVec;
            row.mzFoundStDevVec = candidateScoresTargetBestDiscriminantScore->mzFoundStDevVec;
            return row;
        };

        std::transform(
                scoredTargetDecoyPointersFDRThresholded.begin(),
                scoredTargetDecoyPointersFDRThresholded.end(),
                std::back_inserter(*msCalibrationReaderRows),
                msCalibrationReaderRowsInsertLogic
        );

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::buildCalibration(TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron){

    ERR_INIT

    e = ErrorUtils::isTrue(targetDecoyCandidatePairScoretron->isInit()); ree;

    const double calibrationTrainingFraction = 0.2;
    const bool useExtendedScores = false;
    const bool useNeuralNetworkScores = false;

    const int minTrainingCount = 100;

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;
    QMap<QString, int> fdrVsCount;
    e = setTargetDecoyCandidateScores(
            targetDecoyCandidatePairScoretron,
            m_minTopNMs2Ions,
            calibrationTrainingFraction,
            useExtendedScores,
            useNeuralNetworkScores,
            &scoredTargetDecoyPointers,
            &fdrVsCount
            ); ree;

    double fallBackFDR;
    e = getBestFDRFraction(fdrVsCount, minTrainingCount, &fallBackFDR); ree;
    qDebug() << "Fallback FDR" << fallBackFDR  << "Count" << fdrVsCount.value(QString::number(static_cast<int>(fallBackFDR * 100)));

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointersFDRThresholded;
    e = FDRCLassifierNeuralNet::filterScoreCandidatesByFDR(
            scoredTargetDecoyPointers,
            fallBackFDR,
            &scoredTargetDecoyPointersFDRThresholded
        ); ree;

    QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
    e = buildMsCalibrationReaderRows(
            scoredTargetDecoyPointersFDRThresholded,
            &msCalibrationReaderRows
            ); ree;

//#define WRITE_CALIBRATION
#ifdef WRITE_CALIBRATION
    const QString resultsFilePath = msReaderParquet->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_CAL_FILE_EXTENSION;
    e = ParquetReader::write(msCalibrationReaderRows, resultsFilePath); ree;
    e = m_msCalibratomatic.init(resultsFilePath); ree;
#else
    e = m_msCalibratomatic.init(msCalibrationReaderRows); ree;
#endif

    ERR_RETURN
}

Err PythiaDIAWorkflow::setTargetDecoyCandidateScores(
        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
        int topNMS2Ions,
        double calibrationTrainingFraction,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers,
        QMap<QString, int> *fdrVsCount
        ) {

    ERR_INIT

    e = targetDecoyCandidatePairScoretron->scoreTargetDecoyPairs(
            topNMS2Ions,
            calibrationTrainingFraction,
            m_msCalibratomatic,
            scoredTargetDecoyPointers
    ); ree;

    e = setDiscriminantScoreForCandidates(
            *scoredTargetDecoyPointers,
            useExtendedScores,
            useNeuralNetworkScores,
            topNMS2Ions
    ); ree;

    e = setQValueForCandidates(*scoredTargetDecoyPointers); ree;

    e = FDRCLassifierNeuralNet::outputFDRResults(
            *scoredTargetDecoyPointers,
            true,
            fdrVsCount
    ); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::setDiscriminantScoreForCandidates(
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int theoMzIonsSize
        ) {

    ERR_INIT

    QVector<QPair<ScoresTargets, ScoresDecoys>> targetsScoresVsDecoyScores;
    QVector<TargetDecoyPairTargetKey> targetDecoyPairTargetKeys;
    e = buildClassifierInput(
            targetDecoyCandidatePairPntrs,
            useExtendedScores,
            useNeuralNetworkScores,
            theoMzIonsSize,
            m_scanTimeMinMax,
            &targetsScoresVsDecoyScores,
            &targetDecoyPairTargetKeys
    ); ree;

    QVector<QVector<double>> A;
    QVector<double> b;
    e = ClassifierWeightsManager::buildDataClassifier1(
            targetsScoresVsDecoyScores,
            &A,
            &b
    ); ree;

    QVector<double> weights;
    e = ClassifierWeightsManager::fitWeights(A, b, &weights); ree;

    qDebug() << "Weights:" << weights;
    qDebug() << "b:" << b;

    //TODO change fitweights out of for loop and do a matrix calc
    for(const TargetDecoyPairTargetKey &tdp : targetDecoyPairTargetKeys) {

        QVector<double> discScoreTargets;
        e = ClassifierWeightsManager::applyWeights({tdp.scoresTargetVsScoresDecoys.first}, weights, &discScoreTargets); ree;

        QVector<double> discScoreDecoys;
        e = ClassifierWeightsManager::applyWeights({tdp.scoresTargetVsScoresDecoys.second}, weights, &discScoreDecoys); ree;

        TargetDecoyCandidatePair *tdcp = tdp.targetDecoyCandidatePair;

        (*tdcp->uniqueInfoScanKeyVsScoresTarget())[tdp.uniqueMsInfoScanKey].discriminateScore = discScoreTargets.front();
        (*tdcp->uniqueInfoScanKeyVsScoresDecoy())[tdp.uniqueMsInfoScanKey].discriminateScore = discScoreDecoys.front();

    }

    ERR_RETURN
}

namespace {

    Err buildsetQValueForCandidatesInputs(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
            QMap<PeptideSequenceChargeKey, double> *identifierVsTargets,
            QMap<PeptideSequenceChargeKey, double> *identifierVsDecoys
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairPntrs); ree;
        identifierVsTargets->clear();
        identifierVsDecoys->clear();

        for (TargetDecoyCandidatePair *tdcp : targetDecoyCandidatePairPntrs) {

            const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
                    tdcp->peptideStringWithMods(),
                    tdcp->charge()
            );

            const double discriminantScoreTarget
                    = tdcp->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore;

            const double discriminantScoreDecoy
                    = tdcp->candidateScoresBestDiscriminantScorePtrDecoy()->discriminateScore;

            identifierVsTargets->insert(peptideSequenceChargeKey, discriminantScoreTarget);
            identifierVsDecoys->insert(peptideSequenceChargeKey, discriminantScoreDecoy);

        }

        ERR_RETURN
    }

    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
            const QMap<PeptideSequenceChargeKey, double> &identifierVsQValue,
            const QMap<PeptideSequenceChargeKey, double> &identifierVsDecoyRatio
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairPntrs); ree;

        for (TargetDecoyCandidatePair *tdcp : targetDecoyCandidatePairPntrs) {

            const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
                    tdcp->peptideStringWithMods(),
                    tdcp->charge()
            );

            e = ErrorUtils::isTrue(identifierVsQValue.contains(peptideSequenceChargeKey)); ree;
            e = ErrorUtils::isTrue(identifierVsDecoyRatio.contains(peptideSequenceChargeKey)); ree;

            tdcp->candidateScoresBestDiscriminantScorePtrTarget()->qValue = identifierVsQValue.value(peptideSequenceChargeKey);
            tdcp->candidateScoresBestDiscriminantScorePtrTarget()->decoyRatio = identifierVsDecoyRatio.value(peptideSequenceChargeKey);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::setQValueForCandidates(const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs){

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairPntrs); ree;

    QMap<PeptideSequenceChargeKey, double> identifierVsTargets;
    QMap<PeptideSequenceChargeKey, double> identifierVsDecoys;
    e = buildsetQValueForCandidatesInputs(
            targetDecoyCandidatePairPntrs,
            &identifierVsTargets,
            &identifierVsDecoys
            ); ree;

    QMap<PeptideSequenceChargeKey, double> identifierVsQValue;
    QMap<PeptideSequenceChargeKey, double> identifierVsDecoyRatio;
    e = MathUtils::calculateQValue(
            identifierVsTargets,
            identifierVsDecoys,
            &identifierVsQValue,
            &identifierVsDecoyRatio
            ); ree;

    e = setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            targetDecoyCandidatePairPntrs,
            identifierVsQValue,
            identifierVsDecoyRatio
            ); ree;

    ERR_RETURN
}

namespace {

    struct DOEResult {
        double ppm = -1.0;
        double scanTimeStDev = -1.0;
        double cosineSimAnchor = -1.0;
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

        const QVector<QVector<double>> experiments = {
                {1.5, 1.0,  2.0},
                {3.5, 1.0,  2.0},
                {1.5,  3.0,  2.0},
                {3.5,  3.0,  2.0},
                {1.5,  2.0, 1.0},
                {3.5,  2.0, 1.0},
                {1.5,  2.0,  3.0},
                {3.5,  2.0,  3.0},
                {2.5, 1.0, 1.0},
                {2.5,  3.0, 1.0},
                {2.5, 1.0,  3.0},
                {2.5,  3.0,  3.0},
                {2.5,  2.0,  2.0}
        };

        for (const QVector<double> &exp : experiments) {

            PythiaParameters params = pythiaParameters;
            params.ms2ExtractionWidthPPM = mzPPMStDev * exp.at(0);
            params.scanTimeWindowMinutes = scanTimeStDev * exp.at(1);

            switch (static_cast<int>(exp.at(2))) {
                case 1:
                    params.cosineSimToAnchorThreshold = 0.9;
                    break;
                case 2:
                    params.cosineSimToAnchorThreshold = 0.935;
                    break;
                case 3:
                    params.cosineSimToAnchorThreshold = 0.97;
                    break;
                default:
                    params.cosineSimToAnchorThreshold = pythiaParameters.cosineSimToAnchorThreshold;
            }

            pythiaParametersExperiments->push_back(params);
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

        e = ErrorUtils::toDouble(countsVector.front().first, value); ree;

        ERR_RETURN
    }

    Err getTopFrequencyParameters(
            QVector<DOEResult> *results,
            double *ppmSetting,
            double *scanTimeWidthSetting,
            double *cosineSimSetting
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(*results); ree;

        std::sort(
                results->rbegin(),
                results->rend(),
                [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
        );

        for (const DOEResult &r : *results) {
            qDebug() << r.ppm << r.scanTimeStDev << r.cosineSimAnchor << r.fdrCount;
        }

        *ppmSetting = results->front().ppm;
        *scanTimeWidthSetting = results->front().scanTimeStDev;
        *cosineSimSetting = results->front().cosineSimAnchor;

//        const int topNResults = 5;
//        results->resize(topNResults);
//
//        QMap<QString, int> ppmCounts;
//        QMap<QString, int> scanTimeWidthCounts;
//        QMap<QString, int> cosineSimCounts;
//        for (const DOEResult &r : *results) {
//
//            const QString ppmString = QString::number(r.ppm);
//            const QString scanTimeWidthString = QString::number(r.scanTimeStDev);
//            const QString cosineSimString = QString::number(r.cosineSimAnchor);
//
//            ppmCounts[ppmString]++;
//            scanTimeWidthCounts[scanTimeWidthString]++;
//            cosineSimCounts[cosineSimString]++;
//        }

//        const QVector<QPair<QString, int>> ppmCountsVector = ParallelUtils::convertMapToVectorPairs(ppmCounts);
//        const QVector<QPair<QString, int>> scanTimeWidthCountsVector = ParallelUtils::convertMapToVectorPairs(scanTimeWidthCounts);
//        const QVector<QPair<QString, int>> cosineSimCountsVector = ParallelUtils::convertMapToVectorPairs(cosineSimCounts);
//
//        e = findMostFrequentValue(ppmCountsVector, ppmSetting); ree;
//        e = findMostFrequentValue(scanTimeWidthCountsVector, scanTimeWidthSetting); ree;
//        e = findMostFrequentValue(cosineSimCountsVector, cosineSimSetting); ree;

        ERR_RETURN
    }

}//namespace
Err PythiaDIAWorkflow::optimizeParameters(TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron) {

    ERR_INIT

    const int topNMs2IonsOptimization = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
    );

    qDebug() << "Using top:" << topNMs2IonsOptimization << "fragments for optimization";

    const double selectionFractionValue = 0.05;
    const double fdrThreshold = 0.01;

    QVector<PythiaParameters> pythiaParametersExperiments;
    e = buildDOE(
            m_pythiaParameters,
            m_msCalibratomatic.mzStDev(),
            m_msCalibratomatic.scanTimeStDev(),
            &pythiaParametersExperiments
            ); ree;

    const bool useExtendedScores = true;
    const bool useNeuralNetworkScores = false;
    const int minTrainingCount = 100;

    QVector<DOEResult> results;
    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {

        e = targetDecoyCandidatePairScoretron->setPythiaParameters(pythiaParams); ree;

        QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;
        QMap<QString, int> fdrVsCount;
        e = setTargetDecoyCandidateScores(
                targetDecoyCandidatePairScoretron,
                m_minTopNMs2Ions,
                selectionFractionValue,
                useExtendedScores,
                useNeuralNetworkScores,
                &scoredTargetDecoyPointers,
                &fdrVsCount
                ); ree;

        double fallBackFDR;
        e = getBestFDRFraction(fdrVsCount, minTrainingCount, &fallBackFDR); ree;
        qDebug() << "Fallback FDR" << fallBackFDR  << "Count" << fdrVsCount.value(QString::number(static_cast<int>(fallBackFDR * 100)));


       int targetCountAboveFDRQValueThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                scoredTargetDecoyPointers,
                fdrThreshold,
                &targetCountAboveFDRQValueThreshold
                ); ree;

        DOEResult res;
        res.ppm = pythiaParams.ms2ExtractionWidthPPM;
        res.scanTimeStDev = pythiaParams.scanTimeWindowMinutes;
        res.cosineSimAnchor = pythiaParams.cosineSimToAnchorThreshold;
        res.fdrCount = targetCountAboveFDRQValueThreshold;
        results.push_back(res);
    }

    e = getTopFrequencyParameters(
            &results,
            &m_pythiaParameters.ms2ExtractionWidthPPM,
            &m_pythiaParameters.scanTimeWindowMinutes,
            &m_pythiaParameters.cosineSimToAnchorThreshold
            ); ree;

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

Err PythiaDIAWorkflow::mainAnalysis(
        MsReaderPointerAcc *msReaderPointerAcc,
        QVector<ScoredCandidate> *scoredCandidatesAll
        ) {

    ERR_INIT

//    m_pythiaParameters.print();
//
//    const double selectionFractionBypassValue = -1.0;
//
//    const int topNMs2IonsMainAnalysis = std::max(
//            m_minTopNMs2Ions,
//            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions))
//    );
//
//    qDebug() << "Using top:" << topNMs2IonsMainAnalysis << "fragments for main analysis";
//
//    QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> uniqueInfoScanKeyVsCandidatePeptides;
//    e = buildCandidates(
//            topNMs2IonsMainAnalysis,
//            selectionFractionBypassValue,
//            &uniqueInfoScanKeyVsCandidatePeptides
//    ); ree;
//
//    const bool useExtendedScores = true;
//    const bool useNeuralNetworkScores = false;
//
//    MS2DataExtractomatic ms2DataExtractomatic;
//    e = ms2DataExtractomatic.init(
//            m_pythiaParameters,
//            topNMs2IonsMainAnalysis,
//            useExtendedScores,
//            useNeuralNetworkScores,
//            msReaderPointerAcc,
//            m_msCalibratomatic
//            ); ree;
//
//    e = ms2DataExtractomatic.extractMS2ForCandidates(
//            uniqueInfoScanKeyVsCandidatePeptides,
//            scoredCandidatesAll
//    ); ree;

    ERR_RETURN
}

namespace {

    Err buildPeptideStringWithModsVsScoreCandidatesDecoys(
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            QMap<PeptideStringWithMods, ScoredCandidate> *peptideStringWithModsVsScoreCandidatesDecoys
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;

//        for (const ScoredCandidate &sc : scoredCandidatesAll) {
//
////            if (!sc.isDecoy) {
////                continue;
////            }
////
////            const QStringList splitPepStrWithMods = sc.peptideStringWithMods.split(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP);
////            e = ErrorUtils::isEqual(splitPepStrWithMods.size(), 2); ree;
////
////            const PeptideStringWithMods &peptideStringWithMods = splitPepStrWithMods.front();
////            peptideStringWithModsVsScoreCandidatesDecoys->insert(peptideStringWithMods, sc);
//
//        }

        e = ErrorUtils::isNotEmpty(*peptideStringWithModsVsScoreCandidatesDecoys); ree

        ERR_RETURN
    }

    Err buildMs2IonsFromScoredCandidate(
            const ScoredCandidate &scoredCandidate,
            QVector<MS2Ion> *ms2Ions
    ) {

        ERR_INIT

//        if (scoredCandidate.theoIntensityVec.isEmpty()) {
//            qDebug() << scoredCandidate.peptideStringWithMods;
//            qDebug() << scoredCandidate.mzSearchedVec;
//            qDebug() << scoredCandidate.intensityFoundMaxVec;
//        }
//
//        e = ErrorUtils::isNotEmpty(scoredCandidate.theoIntensityVec); ree;
//        e = ErrorUtils::isEqual(scoredCandidate.theoIntensityVec.size(), scoredCandidate.mzSearchedVec.size()); ree;
//
//        for (int i = 0; i < scoredCandidate.theoIntensityVec.size(); i++) {
//            const double intensity = scoredCandidate.theoIntensityVec.at(i);
//            const double mz = scoredCandidate.mzSearchedVec.at(i);
//
//            if (MathUtils::tZero(mz) || MathUtils::tZero(intensity)) {
//                continue;
//            }
//
//            MS2Ion ms2Ion;
//            ms2Ion.mz = mz;
//            ms2Ion.intensity = intensity;
//
//            ms2Ions->push_back(ms2Ion);
//        }

        ERR_RETURN
    }

    Err buildTandemDeconvolutionInput(
            const QVector<ScoredCandidate> &scoredCandidatesTargetsFDRThresholded,
            const QVector<ScoredCandidate> &scoredCandidatesAll,
            QMap<ScanNumber, QMap<PeptideStringWithMods, QVector<MS2Ion>>> *scanNumberVsTandemPredictions
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesTargetsFDRThresholded); ree;
        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;

//        QMap<PeptideStringWithMods, ScoredCandidate> peptideStringWithModsVsScoreCandidatesDecoys;
//        e = buildPeptideStringWithModsVsScoreCandidatesDecoys(
//                scoredCandidatesAll,
//                &peptideStringWithModsVsScoreCandidatesDecoys
//        ); ree;
//
//        for (const ScoredCandidate &sc : scoredCandidatesTargetsFDRThresholded) {
//            e = ErrorUtils::isTrue(peptideStringWithModsVsScoreCandidatesDecoys.contains(sc.peptideStringWithMods)); ree;
//
//            const ScoredCandidate &scoredCandidateDecoy
//                    = peptideStringWithModsVsScoreCandidatesDecoys.value(sc.peptideStringWithMods);
//
//            QVector<MS2Ion> ms2IonsTarget;
//            e = buildMs2IonsFromScoredCandidate(sc, &ms2IonsTarget); ree;
//            (*scanNumberVsTandemPredictions)[sc.scanNumber].insert(sc.peptideStringWithMods, ms2IonsTarget);
//
//            if(scoredCandidateDecoy.mzSearchedVec.isEmpty()) {
//                continue;
//            }
//
//            QVector<MS2Ion> ms2IonsDecoy;
//            e = buildMs2IonsFromScoredCandidate(scoredCandidateDecoy, &ms2IonsDecoy); ree;
//            (*scanNumberVsTandemPredictions)[sc.scanNumber].insert(scoredCandidateDecoy.peptideStringWithMods, ms2IonsDecoy);
//        }

        ERR_RETURN
    }

//    struct DeconVol {
//        ScanNumber scanNumber = -1;
//        ScanPoints scanPoints;
//        QMap<PeptideStringWithMods, QVector<MS2Ion>> tandemPredictions;
//    };
//
//    struct DeconResult {
//        Err e = eNoError;
//        ScanNumber scanNumber = -1;
//        QMap<PeptideStringWithMods, TandemDeconvolverResult> tandemDeconvolverResult;
//    };
//
//    DeconResult tandemDeconvolutionLogic(
//            const DeconVol &deconVol,
//            const PythiaParameters &params
//    ) {
//
//        ERR_INIT
//
//        const int maxIteration= 20;
//        const double stopTol = 0.000000001;
//        const double pValThreshold = 0.05;
//
//        DeconResult deconResult;
//
//        TandemSpectraDeconvolvotron tandemSpectraDeconvolvotron;
//        e = tandemSpectraDeconvolvotron.init(
//                S_GLOBAL_SETTINGS.HASHING_PRECISION,
//                params.mzMaxDataStructure,
//                params.ms2ExtractionWidthPPM,
//                maxIteration,
//                stopTol,
//                pValThreshold
//        );
//        if (e != eNoError){
//            return deconResult;
//        }
//
//        QMap<PeptideStringWithMods, TandemDeconvolverResult> pepSeqVsWeight;
//        e = tandemSpectraDeconvolvotron.deconvolveTandemSpectra(
//                deconVol.scanPoints,
//                deconVol.tandemPredictions,
//                &pepSeqVsWeight
//        );
//        if (e != eNoError){
//            return deconResult;
//        }
//
//        deconResult.e = e;
//        deconResult.scanNumber = deconVol.scanNumber;
//        deconResult.tandemDeconvolverResult = pepSeqVsWeight;
//
//        return deconResult;
//    }
//
//    Err deconvolveTandemSpectra(
//            const QMap<ScanNumber, QMap<PeptideStringWithMods, QVector<MS2Ion>>> &scanNumberVsTandemPredictions,
//            const PythiaParameters &params,
//            MsReaderPointerAcc *msReaderPointerAcc,
//            QMap<ScanNumber, QMap<PeptideStringWithMods, TandemDeconvolverResult>> *scanNumberVsTandemDeconvolverResult
//    ) {
//
//        ERR_INIT
//        e = ErrorUtils::isNotEmpty(scanNumberVsTandemPredictions); ree;
//
//        QElapsedTimer et;
//        et.start();
//
//        const auto insertLogic = [msReaderPointerAcc, scanNumberVsTandemPredictions](const ScanNumber sn){
//
//            DeconVol deconVol;
//            deconVol.scanNumber = sn;
//            deconVol.tandemPredictions = scanNumberVsTandemPredictions.value(sn);
//            msReaderPointerAcc->ptr->getScanPoints(sn, &deconVol.scanPoints);
//
//            return deconVol;
//        };
//
//        QVector<DeconVol> scanPointsVsTandemPredictions;
//        const QList<ScanNumber> &scanNumbers = scanNumberVsTandemPredictions.keys();
//        std::transform(
//                scanNumbers.begin(),
//                scanNumbers.end(),
//                std::back_inserter(scanPointsVsTandemPredictions),
//                insertLogic
//        );
//
//#define PARALLEL_DECONVOLVE
//#ifdef PARALLEL_DECONVOLVE
//        const auto deconvolutionLogicBinder = std::bind(
//            tandemDeconvolutionLogic,
//            std::placeholders::_1,
//            params
//        );
//
//        QFuture<DeconResult> futures = QtConcurrent::mapped(
//                scanPointsVsTandemPredictions,
//                deconvolutionLogicBinder
//                );
//        futures.waitForFinished();
//
//        for (const DeconResult &result : futures) {
//            e = result.e; ree;
//            (*scanNumberVsTandemDeconvolverResult)[result.scanNumber] = result.tandemDeconvolverResult;
//        }
//#else
//        for (const DeconVol &deconVol : scanPointsVsTandemPredictions) {
//
//            const DeconResult result = tandemDeconvolutionLogic(
//                    deconVol,
//                    params
//            );
//
//            e = result.e; ree;
//            (*scanNumberVsTandemDeconvolverResult)[result.scanNumber] = result.tandemDeconvolverResult;
//        }
//#endif
//
//        qDebug() << "Deconvolution finished in" << et.elapsed() << "mSec";
//
//        ERR_RETURN
//    }
//
//    void filterScoredCandidatesByWeightAndPVal(
//            QVector<ScoredCandidate> *scoredCandidatesAllUpdated,
//            double pValThreshold) {
//
//        const double weightThreshold = 0.0;
//        const auto terminatorLogic = [weightThreshold, pValThreshold](const ScoredCandidate &s){
//            return s.matrixWeight < weightThreshold || s.matrixPValue > pValThreshold;
//        };
//
//        const auto terminator = std::remove_if(
//                scoredCandidatesAllUpdated->begin(),
//                scoredCandidatesAllUpdated->end(),
//                terminatorLogic
//        );
//        scoredCandidatesAllUpdated->erase(terminator, scoredCandidatesAllUpdated->end());
//    }

}//namespace
Err PythiaDIAWorkflow::removeInterferingCandidates(
        MsReaderPointerAcc *msReaderPointerAcc,
        const QVector<ScoredCandidate> &scoredCandidatesTargetsFDRThresholded,
        const QVector<ScoredCandidate> &scoredCandidatesAll,
        QVector<ScoredCandidate> *scoredCandidatesAllUpdated
        ) {

    ERR_INIT

    qDebug() << "Starting interference removal of shared tandem fragments";

    QMap<ScanNumber, QMap<PeptideStringWithMods, QVector<MS2Ion>>> scanNumberVsTandemPredictions;
    e = buildTandemDeconvolutionInput(
            scoredCandidatesTargetsFDRThresholded,
            scoredCandidatesAll,
            &scanNumberVsTandemPredictions
    ); ree;

//    QMap<ScanNumber, QMap<PeptideStringWithMods, TandemDeconvolverResult>> scanNumberVsTandemDeconvolverResult;
//    e = deconvolveTandemSpectra(
//            scanNumberVsTandemPredictions,
//            m_pythiaParameters,
//            msReaderPointerAcc,
//            &scanNumberVsTandemDeconvolverResult
//    ); ree;
//
//    for (const ScoredCandidate &sc : scoredCandidatesAll) {
//
//        const QMap<PeptideStringWithMods, TandemDeconvolverResult> &tandemResult
//                = scanNumberVsTandemDeconvolverResult.value(sc.scanNumber);
//
//        const TandemDeconvolverResult &tandemDeconvolverResult = tandemResult.value(sc.peptideStringWithMods);
//
//        ScoredCandidate scNew = sc;
//        scNew.matrixWeight = tandemDeconvolverResult.discScore;
//        scNew.matrixPValue = tandemDeconvolverResult.pVal;
//        scNew.matrixError = tandemDeconvolverResult.frameError;
//        scNew.scanNumberCandidateCount = tandemDeconvolverResult.scanNumberCandidateCount;
//
//        scoredCandidatesAllUpdated->push_back(scNew);
//    }
//
//    filterScoredCandidatesByWeightAndPVal(
//            scoredCandidatesAllUpdated,
//            m_pythiaParameters.pValThreshold
//            );
//
//    QMap<QString, int> fdrCountResultsUnused;
//    e = MS2DataExtractomatic::outputFDRResults(*scoredCandidatesAllUpdated, false, &fdrCountResultsUnused); ree;

    ERR_RETURN
}

namespace {

//    void filterTargetsOut(QVector<ScoredCandidate> *scoredCandidatesAll) {
//
//        const auto terminatorLogic = [](const ScoredCandidate &sc){
//            return sc.isDecoy == 0;
//        };
//
//        const auto terminator
//            = std::remove_if(scoredCandidatesAll->begin(), scoredCandidatesAll->end(), terminatorLogic);
//
//        scoredCandidatesAll->erase(terminator, scoredCandidatesAll->end());
//    }
//
//    Err getTopDecoys(
//            const QVector<ScoredCandidate> &scoredCandidatesAll,
//            int decoyCount,
//            QVector<ScoredCandidate> *decoys
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;
//
//        QVector<ScoredCandidate> scoredCandidatesDecoys = scoredCandidatesAll;
//        filterTargetsOut(&scoredCandidatesDecoys);
//
//        const auto sortLogicCosineSimSumDesc = [](const ScoredCandidate &l, const ScoredCandidate &r){
//            return l.cosineSimSum100 < r.cosineSimSum100;
//        };
//
//        std::sort(scoredCandidatesDecoys.rbegin(), scoredCandidatesDecoys.rend(), sortLogicCosineSimSumDesc);
//        scoredCandidatesDecoys.resize(decoyCount);
//
//        *decoys = scoredCandidatesDecoys;
//
//        ERR_RETURN
//    }

}//namespace
Err PythiaDIAWorkflow::applyNeuralNetClassifier(
        const QVector<ScoredCandidate> &scoredCandidatesAll,
        const QVector<ScoredCandidate> &scoredCandidatesCulled,
        const QPair<double, double> &scanTimeMinMax,
        QVector<ScoredCandidate> *scoredCandidatesClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;
    e = ErrorUtils::isNotEmpty(scoredCandidatesCulled); ree;

//    QVector<ScoredCandidate> topDecoys;
//    e = getTopDecoys(scoredCandidatesAll, scoredCandidatesCulled.size(), &topDecoys); ree
//
//    QVector<ScoredCandidate> trainingData = scoredCandidatesCulled;
//    trainingData.append(topDecoys);
//
//    const int epochs = 10;
//    const int baggingSize = 6;
//    const double batchFraction = 0.01;
//    const double learningRate = 0.001;
//    FDRCLassifierNeuralNet fdrClassifierNeuralNet;
//    e = fdrClassifierNeuralNet.init(
//            m_pythiaParameters.topNMs2Ions,
//            epochs,
//            baggingSize,
//            batchFraction,
//            learningRate,
//            scanTimeMinMax
//            ); ree;
//
//    e = fdrClassifierNeuralNet.exec(
//            trainingData,
//            scoredCandidatesClassifier
//            ); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::updateProteinGroupAnnotation(
        const QString &fastaFilePath,
        QVector<ScoredCandidate> *scoredCandidates
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(*scoredCandidates); ree;

//    FastaReader fastaReader;
//    e = fastaReader.parseFastaFile(fastaFilePath); ree;
//
//    QVector<PeptideSequence> peptideSequences;
//    QMap<PeptideStringWithMods, QVector<FastaEntry>> peptideStringWithModsVsFastaEntries;
//    e = FastaFileToPeptidesListWorkFlow::digestFastaEntries(
//            m_pythiaParameters,
//            fastaReader.fastaEntries(),
//            &peptideSequences,
//            &peptideStringWithModsVsFastaEntries
//            ); ree;
//
//    QMap<PeptideStringWithMods, QVector<FastaEntry>> peptideStringWithModsVsFastaEntriesLeucinesReplaced;
//    for (auto it = peptideStringWithModsVsFastaEntries.begin(); it != peptideStringWithModsVsFastaEntries.end(); it++) {
//        QString peptideSeqReplacedLeucines = it.key();
//        peptideSeqReplacedLeucines = peptideSeqReplacedLeucines.replace('L', 'J').replace('I', 'J');
//        peptideStringWithModsVsFastaEntriesLeucinesReplaced.insert(peptideSeqReplacedLeucines, it.value());
//    }
//
//    for (int i = 0; i < scoredCandidates->size(); i++) {
//
//        ScoredCandidate &sc = (*scoredCandidates)[i];
//
//        QString peptideSeqReplacedLeucines = sc.peptideStringWithMods;
//        peptideSeqReplacedLeucines = peptideSeqReplacedLeucines.replace('L', 'J').replace('I', 'J');
//
//        const QVector<FastaEntry> &fastaEntries = peptideStringWithModsVsFastaEntriesLeucinesReplaced.value(peptideSeqReplacedLeucines);
//
//        QStringList fastaDescriptions;
//        std::transform(
//                fastaEntries.begin(),
//                fastaEntries.end(),
//                std::back_inserter(fastaDescriptions),
//                [](const FastaEntry &fe){return fe.fastaDescription;}
//                );
//
//        sc.proteinGroup = fastaDescriptions.join(';');
//    }

    ERR_RETURN
}
