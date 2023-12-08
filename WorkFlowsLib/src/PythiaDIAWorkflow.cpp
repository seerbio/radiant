//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "MsReaderParquet.h"

//#include "CandidateClassifier.h"
//#include "ClassifierWeightsManager.h"
//#include "EigenUtils.h"
//#include "ErrorUtils.h"
//#include "FastaFileToPeptidesListWorkFlow.h"
//#include "FastaReader.h"
//#include "FDRCLassifierNeuralNet.h"
////#include "FeatureFinderHillBuilder.h"
//#include "MathUtils.h"
//#include "MS2ChargeDeconvolvotron.h"
////#include "MS2DataExtractomatic.h"
//#include "MsCalibrationReader.h"
////#include "MsFrameScoretron.h"
//#include "MsReaderPointerAcc.h"
//#include "ParallelUtils.h"
////#include "PeakIntegratomatic.h"
//#include "TandemSpectraDeconvolvotron.h"
//#include "TandemFragmentPredictotron.h"
//#include "TargetDecoyCandidatePairScoretron.h"
//#include "TurboXIC.h"
//
//#include <QtConcurrent/QtConcurrent>
//
//#include <iostream>
//
//
//PythiaDIAWorkflow::PythiaDIAWorkflow()
//: m_minTopNMs2Ions(6)
//{}
//
//Err PythiaDIAWorkflow::init(
//        const PythiaParameters &pythiaParameters,
//        const QString &fragLibUri,
//        const QString &fastaUri
//        ) {
//
//    ERR_INIT
//
//    pythiaParameters.print();
//
//    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
//    e = ErrorUtils::fileExists(fragLibUri); ree;
//
//    if (!fastaUri.isEmpty()) {
//        e = ErrorUtils::fileExists(fastaUri); ree;
//        m_fastaUri = fastaUri;
//    }
//
//    m_pythiaParameters = pythiaParameters;
//    m_fragLibUri = fragLibUri;
//
//    e = m_targetDecoyCandidatePairManager.init(m_pythiaParameters, m_fragLibUri); ree;
//
//    ERR_RETURN
//}
//
//Err PythiaDIAWorkflow::processFile(const QString &_msDataFilePath) {
//
//    ERR_INIT
//
//    m_msScanInfos.clear();
//    m_scanTimeMinMax = {-1.0, -1.0};
//
//    QString msDataFilePath = _msDataFilePath;
//
//    MsReaderPointerAcc msReaderPointerAcc;
//    e = msReaderPointerAcc.openFile(msDataFilePath); ree;
//
//    m_msScanInfos = msReaderPointerAcc.ptr->getUniqueTandemMsScanInfos();
//    m_scanTimeMinMax = msReaderPointerAcc.ptr->scanTimeMinMax();
//
//    if (m_pythiaParameters.deisotopeScans) {
//        e = deisotopeScans(&msReaderPointerAcc); ree;
//    }
//
//    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
//    e = msReaderPointerAcc.ptr->collateTandemPrecursorTargetsDIA(
//            &diaTargetFrames
//            ); ree;
//
//    const int msLevel = 1;
//    QMap<ScanNumber, ScanPoints*> scanNumberVsScanTimeMS1;
//    e = msReaderPointerAcc.ptr->getScanPoints(msLevel, &scanNumberVsScanTimeMS1); ree;
//
//    TargetDecoyCandidatePairScoretron targetDecoyCandidatePairScoretron;
//    e = targetDecoyCandidatePairScoretron.init(
//            m_pythiaParameters,
//            scanNumberVsScanTimeMS1,
//            &msReaderPointerAcc,
//            &diaTargetFrames,
//            &m_targetDecoyCandidatePairManager
//            ); ree;
//
//    const double calibrationTrainingFraction = 0.2;
//    e = buildCalibration(
//            calibrationTrainingFraction,
//            false,
//            &targetDecoyCandidatePairScoretron
//            ); ree;
//
//    e = recalibrateMzVals(
//            &diaTargetFrames,
//            &scanNumberVsScanTimeMS1,
//            &targetDecoyCandidatePairScoretron,
//            &msReaderPointerAcc
//            ); ree;
//
//    e = optimizeParameters(&targetDecoyCandidatePairScoretron); ree;
//
////    e = buildCalibration(
////            calibrationTrainingFraction,
////            true,
////            &targetDecoyCandidatePairScoretron
////            ); ree;
////
////    e = recalibrateMzVals(
////            &diaTargetFrames,
////            &scanNumberVsScanTimeMS1,
////            &targetDecoyCandidatePairScoretron,
////            &msReaderPointerAcc
////            ); ree;
//
//    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;
//    int psmCountOnePercentFDR;
//    e = mainAnalysis(
//            &targetDecoyCandidatePairScoretron,
//            &scoredTargetDecoyPointers,
//            &psmCountOnePercentFDR
//            ); ree;
//
//    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointersFDRThresholded;
//    const double fdrThreshold = 0.5;
//    e = FDRCLassifierNeuralNet::filterScoreCandidatesByFDR(
//            scoredTargetDecoyPointers,
//            fdrThreshold,
//            &scoredTargetDecoyPointersFDRThresholded
//    ); ree;
//
//    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointersFDRFiltered;
//
////#define REMOVE_INTERFERENCE
//#ifdef REMOVE_INTERFERENCE
//    e = removeInterferingCandidates(
//            &msReaderPointerAcc,
//            scoredTargetDecoyPointersFDRThresholded,
//            &scoredTargetDecoyPointersFDRFiltered
//            ); ree;
//#else
//    scoredTargetDecoyPointersFDRFiltered = scoredTargetDecoyPointersFDRThresholded;
//#endif
//
//    QVector<CandidateScores> scoredCandidatesClassifierUpdated;
//    e = applyNeuralNetClassifier(
//            scoredTargetDecoyPointers,
//            scoredTargetDecoyPointersFDRFiltered,
//            msReaderPointerAcc.ptr->scanTimeMinMax(),
//            psmCountOnePercentFDR,
//            m_pythiaParameters.reportDecoys,
//            &scoredCandidatesClassifierUpdated
//            ); ree;
//
//    qDebug() << "Updating" << scoredCandidatesClassifierUpdated.size() << "PSMs";
//
//    e = updateProteinGroupAnnotation(
//            m_fastaUri,
//            &scoredCandidatesClassifierUpdated
//            ); ree;
//
//    const QString resultsFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
//    e = ParquetReader::write(scoredCandidatesClassifierUpdated, resultsFilePath); ree;
//
//    ERR_RETURN
//}
//
//namespace {
//
//    QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> deisotopeLogic(const QVector<QPair<ScanNumber, ScanPoints>> &scanPointPairs, double ppmTol) {
//
//        ERR_INIT
//
//        const QString &chargeModelFilePath
//                = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");
//
//        const QString &monoModelFilePath
//                = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");
//
//        MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
//        e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);
//
//        QVector<QPair<ScanNumber, ScanPoints>> deisotopedScanPoints;
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
//
//        return {e, deisotopedScanPoints};
//    }
//
//}//namespace
//Err PythiaDIAWorkflow::deisotopeScans(MsReaderPointerAcc *msReaderPointerAcc) {
//
//    ERR_INIT
//
//    qDebug() << "Deisotoping start";
//
//    QElapsedTimer et;
//    et.start();
//
//    const auto deisotopeLogicBinder = std::bind(
//            deisotopeLogic,
//            std::placeholders::_1,
//            m_pythiaParameters.ms2ExtractionWidthPPM
//    );
//
//    QVector<QVector<QPair<ScanNumber, ScanPoints>>> scanPointsTranched;
//    ParallelUtils::trancheMapForParallelization(
//            msReaderPointerAcc->ptr->getScanPoints(),
//            ParallelUtils::numberOfAvailableSystemProcessors(),
//            &scanPointsTranched
//    );
//
//    QMap<ScanNumber, ScanPoints> deisotopedScanPoints;
//
////#define PARALLEL_DEISO
//#ifdef PARALLEL_DEISO
//    QFuture<QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>>> futures
//            = QtConcurrent::mapped(scanPointsTranched, deisotopeLogicBinder);
//    futures.waitForFinished();
//
//    for (const QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> &res : futures) {
//        e = res.first; ree;
//        for (const QPair<ScanNumber, ScanPoints> &r : res.second) {
//            const ScanNumber scanNumber = r.first;
//            const ScanPoints &scanPoints = r.second;
//            deisotopedScanPoints.insert(scanNumber, scanPoints);
//        }
//    }
//#else
//    for (const QVector<QPair<ScanNumber, ScanPoints>> &spPairs :  scanPointsTranched) {
//
//        QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> deisotopedScan = deisotopeLogic(
//                spPairs,
//                m_pythiaParameters.ms2ExtractionWidthPPM
//                );
//
//        e = deisotopedScan.first; ree;
//        const QVector<QPair<ScanNumber, ScanPoints>> &res = deisotopedScan.second;
//        for (const QPair<ScanNumber, ScanPoints> &r : res) {
//            deisotopedScanPoints.insert(r.first, r.second);
//        }
//    }
//#endif
//
//    e = msReaderPointerAcc->ptr->setScanPoints(deisotopedScanPoints); ree;
//
//    qDebug() << "Scans deisotoped in" << et.elapsed() << "mSec";
//
//    ERR_RETURN
//}
//
//namespace {
//
//    struct TargetDecoyPairTargetKey {
//        TargetDecoyCandidatePair* targetDecoyCandidatePair;
//        QPair<ScoresTargets, ScoresDecoys> scoresTargetVsScoresDecoys;
//        MzTargetKey uniqueMsInfoScanKey;
//    };
//
//    Err buildClassifierInput(
//            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
//            bool useExtendedScores,
//            bool useNeuralNetworkScores,
//            int theoMzIonsSize,
//            const QPair<double, double> &scanTimeMinMax,
//            QVector<QPair<ScoresTargets, ScoresDecoys>> *targetsScoresVsDecoyScores,
//            QVector<TargetDecoyPairTargetKey> *targetDecoyPairTargetKeys
//            ) {
//
//        ERR_INIT
//
//        const int theoMzIonsSizeMin = 6;
//        const double cosineSimSum100MinForTraining = 0.9;
//
//        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairPntrs); ree;
//        e = ErrorUtils::isTrue(scanTimeMinMax.second > scanTimeMinMax.first); ree;
//        e = ErrorUtils::isTrue(theoMzIonsSize >= theoMzIonsSizeMin); ree;
//
//        targetsScoresVsDecoyScores->clear();
//
//        for (TargetDecoyCandidatePair *tdc : targetDecoyCandidatePairPntrs) {
//
//            const QList<MzTargetKey> &uniqueInfoScanKeys = tdc->uniqueInfoScanKeyVsScoresTarget()->keys();
//            for (const MzTargetKey &key : uniqueInfoScanKeys) {
//
//                const bool decoyContainsTargetKey = tdc->uniqueInfoScanKeyVsScoresDecoy()->contains(key);
//                if (!decoyContainsTargetKey) {
//                    qDebug() << "Decoy scores not found for" << key << tdc->peptideStringWithMods();
//                    qDebug() << "Keys" << uniqueInfoScanKeys << "Mz" << tdc->mz();
//                    continue;
//                }
//                e = ErrorUtils::isTrue(decoyContainsTargetKey); ree;
//
//                const CandidateScores &candidateScores = tdc->uniqueInfoScanKeyVsScoresTarget()->value(key);
//                if (candidateScores.cosineSimSum100 < cosineSimSum100MinForTraining) {
//                    continue;
//                }
//
//                ScoresTargets scoresTargets;
//                e = FDRCLassifierNeuralNet::buildScoreVector(
//                        candidateScores,
//                        useExtendedScores,
//                        useNeuralNetworkScores,
//                        theoMzIonsSize,
//                        scanTimeMinMax,
//                        &scoresTargets
//                        ); ree;
//
//                ScoresDecoys scoresDecoys;
//                e = FDRCLassifierNeuralNet::buildScoreVector(
//                        tdc->uniqueInfoScanKeyVsScoresDecoy()->value(key),
//                        useExtendedScores,
//                        useNeuralNetworkScores,
//                        theoMzIonsSize,
//                        scanTimeMinMax,
//                        &scoresDecoys
//                        ); ree;
//
//                TargetDecoyPairTargetKey targetDecoyPairTargetKey;
//                targetDecoyPairTargetKey.targetDecoyCandidatePair = tdc;
//                targetDecoyPairTargetKey.scoresTargetVsScoresDecoys = {scoresTargets, scoresDecoys};
//                targetDecoyPairTargetKey.uniqueMsInfoScanKey = key;
//
//                targetDecoyPairTargetKeys->push_back(targetDecoyPairTargetKey);
//                targetsScoresVsDecoyScores->push_back({scoresTargets, scoresDecoys});
//            }
//
//        }
//
//        ERR_RETURN
//    }
//
//    Err getBestFDRFraction(
//            const QMap<QString, int> &fdrResults,
//            int minTrainingThreshold,
//            double *bestFDRFraction
//    ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(fdrResults); ree;
//
//        const QStringList fdrFractions = {
//                "1",
//                "2",
//                "5",
//                "10"
//        };
//
//        for (const QString &fdrStr : fdrFractions) {
//            if (fdrResults.value(fdrStr) > minTrainingThreshold) {
//
//                double fdrPercent;
//                e = ErrorUtils::toDouble(fdrStr, &fdrPercent); ree;
//
//                *bestFDRFraction = fdrPercent / 100.0;
//                ERR_RETURN
//            }
//        }
//
//        double fdrPercent;
//        e = ErrorUtils::toDouble(fdrFractions.back(), &fdrPercent); ree;
//
//        *bestFDRFraction = fdrPercent;
//
//        ERR_RETURN
//    }
//
//    Err buildMsCalibrationReaderRows(
//            const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointersFDRThresholded,
//            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
//    ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(scoredTargetDecoyPointersFDRThresholded); ree;
//
//        qDebug() << scoredTargetDecoyPointersFDRThresholded.size() << "Found for recalibartion";
//
//        const auto msCalibrationReaderRowsInsertLogic = [](TargetDecoyCandidatePair* tdp){
//
//            CandidateScores *candidateScoresTargetBestDiscriminantScore = tdp->candidateScoresBestDiscriminantScorePtrTarget();
//
//            MsCalibarationReaderRow row;
//            row.peptideStringWithMods = tdp->peptideStringWithMods();
//            row.iRTPredicted = static_cast<float>(tdp->iRt());
//            row.scanTime = candidateScoresTargetBestDiscriminantScore->scanTime;
//            row.mzSearchedVec = candidateScoresTargetBestDiscriminantScore->mzSearchedVec;
//            row.mzFoundMeanVec = candidateScoresTargetBestDiscriminantScore->mzFoundMeanVec;
//            row.mzFoundStDevVec = candidateScoresTargetBestDiscriminantScore->mzFoundStDevVec;
//            return row;
//        };
//
//        std::transform(
//                scoredTargetDecoyPointersFDRThresholded.begin(),
//                scoredTargetDecoyPointersFDRThresholded.end(),
//                std::back_inserter(*msCalibrationReaderRows),
//                msCalibrationReaderRowsInsertLogic
//        );
//
//        ERR_RETURN
//    }
//
//}//namespace
//Err PythiaDIAWorkflow::buildCalibration(
//        double calibrationTrainingFraction,
//        bool useExtendedScores,
//        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron
//        ){
//
//    ERR_INIT
//
//    e = ErrorUtils::isTrue(targetDecoyCandidatePairScoretron->isInit()); ree;
//
//    const bool useNeuralNetworkScores = false;
//    const int minTrainingCountTranche = 50;
//
//    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;
//
//    double calFractionIter = calibrationTrainingFraction;
//
//    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
//            m_pythiaParameters.mzMinDataStructure,
//            m_pythiaParameters.mzMaxDataStructure,
//            calFractionIter,
//            &scoredTargetDecoyPointers
//            ); ree;
//
//    const double sizePerTranche = 5000.0;
//    const int trancheSize = std::max(static_cast<int>(scoredTargetDecoyPointers.size() / sizePerTranche), 1);
//
//    QVector<QVector<TargetDecoyCandidatePair*>> scoredTargetDecoyPointersTranched;
//    e = ParallelUtils::trancheVectorForParallelization(
//            scoredTargetDecoyPointers,
//            trancheSize,
//            &scoredTargetDecoyPointersTranched
//            ); ree;
//
//    const double fdrThreshold = 0.1;
//    for (int i = 0; i < scoredTargetDecoyPointersTranched.size(); i++) {
//
//        QVector<TargetDecoyCandidatePair*> tdcp;
//        for (int j = 0; j <= i; j++) {
//            tdcp.append(scoredTargetDecoyPointersTranched.at(j));
//        }
//
//        QMap<QString, int> fdrVsCount;
//        e = setTargetDecoyCandidateScores(
//                targetDecoyCandidatePairScoretron,
//                m_minTopNMs2Ions,
//                useExtendedScores,
//                useNeuralNetworkScores,
//                &tdcp,
//                &fdrVsCount
//                ); ree;
//
//        int targetCountBelowFDRThreshold;
//        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
//                tdcp,
//                fdrThreshold,
//                &targetCountBelowFDRThreshold
//        ); ree;
//
//        std::sort(tdcp.rbegin(), tdcp.rend(), [](TargetDecoyCandidatePair *l, TargetDecoyCandidatePair *r){
//            return l->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore < r->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore;
//        });
//
//        const int minTrainingCount = std::max(minTrainingCountTranche, targetCountBelowFDRThreshold);
//        qDebug() << "Training RT count 10% FDR:" << minTrainingCount;
//
//        tdcp.resize(minTrainingCount);
//        QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
//        e = buildMsCalibrationReaderRows(
//                tdcp,
//                &msCalibrationReaderRows
//                ); ree;
//
////#define WRITE_CALIBRATION
//#ifdef WRITE_CALIBRATION
//        const QString resultsFilePath = "testing" + S_GLOBAL_SETTINGS.DOT_PYTHIA_CAL_FILE_EXTENSION;
//        e = ParquetReader::write(msCalibrationReaderRows, resultsFilePath); ree;
//        e = m_msCalibratomatic.init(resultsFilePath); ree;
//#else
//        e = m_msCalibratomatic.init(msCalibrationReaderRows); ree;
//        m_pythiaParameters.scanTimeWindowMinutes = m_msCalibratomatic.scanTimeStDev() * 3.0;
//        qDebug() << "Setting scanTimeWindowMinutes to:" << m_pythiaParameters.scanTimeWindowMinutes;
//        e = targetDecoyCandidatePairScoretron->setPythiaParameters(m_pythiaParameters); ree;
//
//#endif
//        if (minTrainingCount > 1000) {
//            break;
//        }
//    }
//
////#define WRITE_TARGET_DECOYS
//#ifdef WRITE_TARGET_DECOYS
//    if(useExtendedScores) {
//        QVector<CandidateScores> candidateScores;
//        for (TargetDecoyCandidatePair* ptr : scoredTargetDecoyPointers) {
//            candidateScores.push_back(*ptr->candidateScoresBestDiscriminantScorePtrTarget());
//        }
//        e = ParquetReader::write(candidateScores, "scoresCalibrationAsstroll.prq"); ree;
//    }
//#endif
//
//    ERR_RETURN
//}
//
//Err PythiaDIAWorkflow::recalibrateMzVals(
//        QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> *diaTargetFrame,
//        QMap<ScanNumber, ScanPoints*> *scanNumberVsScanTimeMS1,
//        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
//        MsReaderPointerAcc *msReaderPointerAcc
//        ) {
//
//    ERR_INIT
//
//    e = ErrorUtils::isTrue(m_msCalibratomatic.isInit()); ree;
//    e = ErrorUtils::isFalse(diaTargetFrame->isEmpty()); ree;
//    e = ErrorUtils::isFalse(scanNumberVsScanTimeMS1->isEmpty()); ree;
//
//    qDebug() << "Recalibrating mz vals";
//
//    const QList<MzTargetKey> &diaTargetFrameKeys = diaTargetFrame->keys();
//    for (const MzTargetKey &k: diaTargetFrameKeys) {
////        qDebug() << "Recalibrating mz vals frame" << k;
//        e = m_msCalibratomatic.recalibrateScanPoints(diaTargetFrame->value(k)); ree;
//    }
//
//    qDebug() << "Recalibrating MS1 mz vals frame";
//    e = m_msCalibratomatic.recalibrateScanPoints(*scanNumberVsScanTimeMS1); ree;
//
//    e = targetDecoyCandidatePairScoretron->init(
//            m_pythiaParameters,
//            *scanNumberVsScanTimeMS1,
//            msReaderPointerAcc,
//            diaTargetFrame,
//            &m_targetDecoyCandidatePairManager
//    ); ree;
//
//    ERR_RETURN
//}
//
//Err PythiaDIAWorkflow::setTargetDecoyCandidateScores(
//        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
//        int topNMS2Ions,
//        bool useExtendedScores,
//        bool useNeuralNetworkScores,
//        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers,
//        QMap<QString, int> *fdrVsCount
//        ) {
//
//    ERR_INIT
//
//    e = targetDecoyCandidatePairScoretron->scoreTargetDecoyPairs(
//            topNMS2Ions,
//            m_msCalibratomatic,
//            scoredTargetDecoyPointers
//    ); ree;
//
//    e = setDiscriminantScoreForCandidates(
//            *scoredTargetDecoyPointers,
//            useExtendedScores,
//            useNeuralNetworkScores,
//            topNMS2Ions
//    ); ree;
//
//    e = setQValueForCandidates(*scoredTargetDecoyPointers); ree;
//
//    e = FDRCLassifierNeuralNet::outputFDRResults(
//            *scoredTargetDecoyPointers,
//            true,
//            fdrVsCount
//    ); ree;
//
//    ERR_RETURN
//}
//
//Err PythiaDIAWorkflow::setDiscriminantScoreForCandidates(
//        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
//        bool useExtendedScores,
//        bool useNeuralNetworkScores,
//        int theoMzIonsSize
//        ) {
//
//    ERR_INIT
//
//    QVector<QPair<ScoresTargets, ScoresDecoys>> targetsScoresVsDecoyScores;
//    QVector<TargetDecoyPairTargetKey> targetDecoyPairTargetKeys;
//    e = buildClassifierInput(
//            targetDecoyCandidatePairPntrs,
//            useExtendedScores,
//            useNeuralNetworkScores,
//            theoMzIonsSize,
//            m_scanTimeMinMax,
//            &targetsScoresVsDecoyScores,
//            &targetDecoyPairTargetKeys
//    ); ree;
//
//    QVector<QVector<double>> A;
//    QVector<double> b;
//    e = ClassifierWeightsManager::buildDataClassifier1(
//            targetsScoresVsDecoyScores,
//            &A,
//            &b
//    ); ree;
//
//    QVector<double> weights;
//    e = ClassifierWeightsManager::fitWeights(A, b, &weights); ree;
//
//    qDebug() << "Weights:" << weights;
//    qDebug() << "b:" << b;
//
//
//    QVector<ScoresTargets> scoresTargets;
//    QVector<ScoresDecoys> scoresDecoys;
//    for(const TargetDecoyPairTargetKey &tdp : targetDecoyPairTargetKeys) {
//        scoresTargets.push_back(tdp.scoresTargetVsScoresDecoys.first);
//        scoresDecoys.push_back(tdp.scoresTargetVsScoresDecoys.second);
//    }
//
//    QVector<double> discScoreTargets;
//    e = ClassifierWeightsManager::applyWeights(scoresTargets, weights, &discScoreTargets); ree;
//
//    QVector<double> discScoreDecoys;
//    e = ClassifierWeightsManager::applyWeights(scoresDecoys, weights, &discScoreDecoys); ree;
//
//    e = ErrorUtils::isEqual(scoresTargets.size(), scoresDecoys.size()); ree;
//    e = ErrorUtils::isEqual(scoresTargets.size(), targetDecoyPairTargetKeys.size()); ree;
//
//    for (int i = 0; i < targetDecoyPairTargetKeys.size(); i++) {
//
//        const TargetDecoyPairTargetKey &tdptk = targetDecoyPairTargetKeys.at(i);
//        TargetDecoyCandidatePair *tdcp = tdptk.targetDecoyCandidatePair;
//
//        (*tdcp->uniqueInfoScanKeyVsScoresTarget())[tdptk.uniqueMsInfoScanKey].discriminateScore = discScoreTargets.at(i);
//        (*tdcp->uniqueInfoScanKeyVsScoresDecoy())[tdptk.uniqueMsInfoScanKey].discriminateScore = discScoreDecoys.at(i);
//    }
//
//    ERR_RETURN
//}
//
//namespace {
//
//    Err buildsetQValueForCandidatesInputs(
//            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
//            QMap<PeptideSequenceChargeKey, double> *identifierVsTargets,
//            QMap<PeptideSequenceChargeKey, double> *identifierVsDecoys
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairPntrs); ree;
//        identifierVsTargets->clear();
//        identifierVsDecoys->clear();
//
//        for (TargetDecoyCandidatePair *tdcp : targetDecoyCandidatePairPntrs) {
//
//            const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
//                    tdcp->peptideStringWithMods(),
//                    tdcp->charge()
//            );
//
//            const double discriminantScoreTarget
//                    = tdcp->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore;
//
//            const double discriminantScoreDecoy
//                    = tdcp->candidateScoresBestDiscriminantScorePtrDecoy()->discriminateScore;
//
//            identifierVsTargets->insert(peptideSequenceChargeKey, discriminantScoreTarget);
//            identifierVsDecoys->insert(peptideSequenceChargeKey, discriminantScoreDecoy);
//
//        }
//
//        ERR_RETURN
//    }
//
//    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
//            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs,
//            const QMap<PeptideSequenceChargeKey, double> &identifierVsQValue,
//            const QMap<PeptideSequenceChargeKey, double> &identifierVsDecoyRatio
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairPntrs); ree;
//
//        for (TargetDecoyCandidatePair *tdcp : targetDecoyCandidatePairPntrs) {
//
//            const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
//                    tdcp->peptideStringWithMods(),
//                    tdcp->charge()
//            );
//
//            e = ErrorUtils::isTrue(identifierVsQValue.contains(peptideSequenceChargeKey)); ree;
//            e = ErrorUtils::isTrue(identifierVsDecoyRatio.contains(peptideSequenceChargeKey)); ree;
//
//            tdcp->candidateScoresBestDiscriminantScorePtrTarget()->qValue = identifierVsQValue.value(peptideSequenceChargeKey);
//            tdcp->candidateScoresBestDiscriminantScorePtrTarget()->decoyRatio = identifierVsDecoyRatio.value(peptideSequenceChargeKey);
//        }
//
//        ERR_RETURN
//    }
//
//}//namespace
//Err PythiaDIAWorkflow::setQValueForCandidates(const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairPntrs){
//
//    ERR_INIT
//
//    e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairPntrs); ree;
//
//    QMap<PeptideSequenceChargeKey, double> identifierVsTargets;
//    QMap<PeptideSequenceChargeKey, double> identifierVsDecoys;
//    e = buildsetQValueForCandidatesInputs(
//            targetDecoyCandidatePairPntrs,
//            &identifierVsTargets,
//            &identifierVsDecoys
//            ); ree;
//
//    QMap<PeptideSequenceChargeKey, double> identifierVsQValue;
//    QMap<PeptideSequenceChargeKey, double> identifierVsDecoyRatio;
//    e = MathUtils::calculateQValue(
//            identifierVsTargets,
//            identifierVsDecoys,
//            &identifierVsQValue,
//            &identifierVsDecoyRatio
//            ); ree;
//
//    e = setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
//            targetDecoyCandidatePairPntrs,
//            identifierVsQValue,
//            identifierVsDecoyRatio
//            ); ree;
//
//    ERR_RETURN
//}
//
//namespace {
//
//    Err buildsetQValueForCandidateScoresInputs(
//            QVector<CandidateScores> *candidateScores,
//            QMap<PeptideSequenceChargeKey, double> *identifierVsTargets,
//            QMap<PeptideSequenceChargeKey, double> *identifierVsDecoys
//    ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
//        identifierVsTargets->clear();
//        identifierVsDecoys->clear();
//
//        for (const CandidateScores &cs : *candidateScores) {
//
//            const QString decoyToString = cs.isDecoy ? "_1" : "_0";
//            const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
//                    cs.peptideStringWithMods,
//                    cs.charge
//            ) + decoyToString;
//
//            const double classifierScore = cs.classifierScore;
//
//            if (cs.isDecoy) {
//                identifierVsDecoys->insert(peptideSequenceChargeKey, classifierScore);
//                continue;
//            }
//
//            identifierVsTargets->insert(peptideSequenceChargeKey, classifierScore);
//
//        }
//
//        ERR_RETURN
//    }
//
//    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
//            const QMap<PeptideSequenceChargeKey, double> &identifierVsQValue,
//            const QMap<PeptideSequenceChargeKey, double> &identifierVsDecoyRatio,
//            QVector<CandidateScores> *candidateScores
//    ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
//
//        for (CandidateScores &cs : *candidateScores) {
//
//            if (cs.isDecoy) {
//                continue;
//            }
//
//            const QString decoyToString = cs.isDecoy ? "_1" : "_0";
//            const PeptideSequenceChargeKey peptideSequenceChargeKey = TandemFragmentPredictotron::buildPeptideSequenceChargeKey(
//                    cs.peptideStringWithMods,
//                    cs.charge
//            ) + decoyToString;
//
//            e = ErrorUtils::isTrue(identifierVsQValue.contains(peptideSequenceChargeKey)); ree;
//            e = ErrorUtils::isTrue(identifierVsDecoyRatio.contains(peptideSequenceChargeKey)); ree;
//
//            cs.qValue = identifierVsQValue.value(peptideSequenceChargeKey);
//            cs.decoyRatio = identifierVsDecoyRatio.value(peptideSequenceChargeKey);
//        }
//
//        ERR_RETURN
//    }
//
//}//namespace
//Err PythiaDIAWorkflow::setQValueForCandidates(QVector<CandidateScores> *candidateScores) {
//
//    ERR_INIT
//
//    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
//
//    QMap<PeptideSequenceChargeKey, double> identifierVsTargets;
//    QMap<PeptideSequenceChargeKey, double> identifierVsDecoys;
//    e = buildsetQValueForCandidateScoresInputs(
//            candidateScores,
//            &identifierVsTargets,
//            &identifierVsDecoys
//    ); ree;
//
//    QMap<PeptideSequenceChargeKey, double> identifierVsQValue;
//    QMap<PeptideSequenceChargeKey, double> identifierVsDecoyRatio;
//    e = MathUtils::calculateQValue(
//            identifierVsTargets,
//            identifierVsDecoys,
//            &identifierVsQValue,
//            &identifierVsDecoyRatio
//    ); ree;
//
//    e = setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
//            identifierVsQValue,
//            identifierVsDecoyRatio,
//            candidateScores
//    ); ree;
//
//    ERR_RETURN
//}
//
//namespace {
//
//    struct DOEResult {
//        double ppm = -1.0;
//        double scanTimeStDev = -1.0;
////        double cosineSimAnchor = -1.0;
//        int fdrCount = -1;
//    };
//
//    Err buildDOE(
//            const PythiaParameters &pythiaParameters,
//            double mzPPMStDev,
//            double scanTimeStDev,
//            QVector<PythiaParameters> *pythiaParametersExperiments
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isTrue(mzPPMStDev > 0.0); ree;
//        e = ErrorUtils::isTrue(scanTimeStDev > 0.0); ree;
//        e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
//
////        const QVector<QVector<double>> experiments = {
////                {1.5, 1.0,  2.0},
////                {3.5, 1.0,  2.0},
////                {1.5,  3.0,  2.0},
////                {3.5,  3.0,  2.0},
////                {1.5,  2.0, 1.0},
////                {3.5,  2.0, 1.0},
////                {1.5,  2.0,  3.0},
////                {3.5,  2.0,  3.0},
////                {2.5, 1.0, 1.0},
////                {2.5,  3.0, 1.0},
////                {2.5, 1.0,  3.0},
////                {2.5,  3.0,  3.0},
////                {2.5,  2.0,  2.0}
////        };
//
////        const QVector<QVector<double>> experiments = {
////                {2.5, 1.0},
////                {3.5, 1.0},
////                {2.5,  3.0},
////                {3.5,  3.0},
////                {2.5,  2.0},
////                {3.5,  2.0}
////        };
//
//        const QVector<QVector<double>> experiments = {
//                {2.0, 2.0},
//                {2.5,  2.0},
//                {3.5,  2.0},
//                {4.5, 2.0},
//                {5.5, 2.0}
//        };
//
//        for (const QVector<double> &exp : experiments) {
//
//            PythiaParameters params = pythiaParameters;
//            params.ms2ExtractionWidthPPM = mzPPMStDev * exp.at(0);
//            params.scanTimeWindowMinutes = scanTimeStDev * exp.at(1);
//
////            switch (static_cast<int>(exp.at(2))) {
////                case 1:
////                    params.cosineSimToAnchorThreshold = 0.9;
////                    break;
////                case 2:
////                    params.cosineSimToAnchorThreshold = 0.935;
////                    break;
////                case 3:
////                    params.cosineSimToAnchorThreshold = 0.97;
////                    break;
////                default:
////                    params.cosineSimToAnchorThreshold = pythiaParameters.cosineSimToAnchorThreshold;
////            }
//
//            pythiaParametersExperiments->push_back(params);
//        }
//
//        for (const PythiaParameters &pp : *pythiaParametersExperiments) {
//            qDebug() << "ppmTol" << pp.ms2ExtractionWidthPPM
//                     << "scanTimeWinMin" << pp.scanTimeWindowMinutes
//                     << "cosineSimSumAnchThreshold" << pp.cosineSimToAnchorThreshold;
//        }
//
//        ERR_RETURN
//    }
//
//
//    Err findMostFrequentValue(
//            const QVector<QPair<QString, int>> &countsVector,
//            double *value
//            ) {
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(countsVector); ree;
//
//        if (countsVector.size() == 1) {
//            e = ErrorUtils::toDouble(countsVector.front().first, value); ree;
//            ERR_RETURN
//        }
//
//        QVector<QPair<QString, int>> cv = countsVector;
//
//        std::sort(cv.rbegin(), cv.rend(), [](const QPair<QString, int> &l, const QPair<QString, int> &r){
//            return l.second < r.second;
//        });
//
//        if (cv.at(0).second == cv.at(1).second) {
//            double v1;
//            e = ErrorUtils::toDouble(cv.at(0).first, &v1); ree;
//
//            double v2;
//            e = ErrorUtils::toDouble(cv.at(1).first, &v2); ree;
//
//            * value = (v1 + v2) / 2;
//            ERR_RETURN
//        }
//
//        e = ErrorUtils::toDouble(cv.front().first, value); ree;
//
//        ERR_RETURN
//    }
//
//    Err getTopFrequencyParameters(
//            QVector<DOEResult> *results,
//            double *ppmSetting
////            double *scanTimeWidthSetting
////            double *cosineSimSetting
//            ) {
//
//        ERR_INIT
//        e = ErrorUtils::isNotEmpty(*results); ree;
//
//        std::sort(
//                results->rbegin(),
//                results->rend(),
//                [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
//                );
//
//        for (const DOEResult &r : *results) {
////            qDebug() << r.ppm << r.scanTimeStDev << r.cosineSimAnchor << r.fdrCount;
//            qDebug() << r.ppm << r.scanTimeStDev << r.fdrCount;
//        }
//
//        *ppmSetting = results->front().ppm;
////        *scanTimeWidthSetting = results->front().scanTimeStDev;
//
////        const int topNResults = 3;
////        results->resize(topNResults);
////
////        QMap<QString, int> ppmCounts;
////        QMap<QString, int> scanTimeWidthCounts;
//////        QMap<QString, int> cosineSimCounts;
////        for (const DOEResult &r : *results) {
////
////            const QString ppmString = QString::number(r.ppm);
////            const QString scanTimeWidthString = QString::number(r.scanTimeStDev);
//////            const QString cosineSimString = QString::number(r.cosineSimAnchor);
////
////            ppmCounts[ppmString]++;
////            scanTimeWidthCounts[scanTimeWidthString]++;
//////            cosineSimCounts[cosineSimString]++;
////        }
////
////        const QVector<QPair<QString, int>> ppmCountsVector = ParallelUtils::convertMapToVectorPairs(ppmCounts);
////        const QVector<QPair<QString, int>> scanTimeWidthCountsVector = ParallelUtils::convertMapToVectorPairs(scanTimeWidthCounts);
//////        const QVector<QPair<QString, int>> cosineSimCountsVector = ParallelUtils::convertMapToVectorPairs(cosineSimCounts);
////
////        e = findMostFrequentValue(ppmCountsVector, ppmSetting); ree;
////        e = findMostFrequentValue(scanTimeWidthCountsVector, scanTimeWidthSetting); ree;
//////        e = findMostFrequentValue(cosineSimCountsVector, cosineSimSetting); ree;
//
//        ERR_RETURN
//    }
//
//}//namespace
//Err PythiaDIAWorkflow::optimizeParameters(TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron) {
//
//    ERR_INIT
//
//    const int topNMs2IonsOptimization = std::max(
//            m_minTopNMs2Ions,
//            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions / 2.0))
//    );
//
//    qDebug() << "Using top:" << topNMs2IonsOptimization << "fragments for optimization";
//
//    const double selectionFractionValue = 0.01;
//
//    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;
//
//    double selectionFractionValueIter = selectionFractionValue;
//    while (selectionFractionValueIter < 1.0) {
//
//        scoredTargetDecoyPointers.clear();
//
//        e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
//                m_pythiaParameters.mzMinDataStructure,
//                m_pythiaParameters.mzMaxDataStructure,
//                selectionFractionValueIter,
//                &scoredTargetDecoyPointers
//        ); ree;
//
//        if (scoredTargetDecoyPointers.size() < m_pythiaParameters.trancheSizeMax) {
//            selectionFractionValueIter *= 2;
//            continue;
//        }
//
//        break;
//    }
//
//    const double fdrThreshold = 0.01;
//
//    QVector<PythiaParameters> pythiaParametersExperiments;
//    e = buildDOE(
//            m_pythiaParameters,
//            m_msCalibratomatic.mzStDev(),
//            m_msCalibratomatic.scanTimeStDev(),
//            &pythiaParametersExperiments
//            ); ree;
//
//    const bool useExtendedScores = true;
//    const bool useNeuralNetworkScores = false;
//    const int minTrainingCount = 100;
//
//    QVector<DOEResult> results;
//    for (const PythiaParameters &pythiaParams : pythiaParametersExperiments) {
//
//        qDebug() << "ppmTol" << pythiaParams.ms2ExtractionWidthPPM
//                 << "scanTimeWinMin" << pythiaParams.scanTimeWindowMinutes
//                 << "cosineSimSumAnchThreshold" << pythiaParams.cosineSimToAnchorThreshold;
//
//        e = targetDecoyCandidatePairScoretron->setPythiaParameters(pythiaParams); ree;
//        qDebug() << "STarting opt";
//        QMap<QString, int> fdrVsCount;
//        e = setTargetDecoyCandidateScores(
//                targetDecoyCandidatePairScoretron,
//                topNMs2IonsOptimization,
//                useExtendedScores,
//                useNeuralNetworkScores,
//                &scoredTargetDecoyPointers,
//                &fdrVsCount
//                ); ree;
//        qDebug() << "Ending opt";
//
//        double fallBackFDR;
//        e = getBestFDRFraction(fdrVsCount, minTrainingCount, &fallBackFDR); ree;
//        qDebug() << "Fallback FDR" << fallBackFDR  << "Count" << fdrVsCount.value(QString::number(static_cast<int>(fallBackFDR * 100)));
//
//
//       int targetCountAboveFDRQValueThreshold;
//        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
//                scoredTargetDecoyPointers,
//                fdrThreshold,
//                &targetCountAboveFDRQValueThreshold
//                ); ree;
//
//        DOEResult res;
//        res.ppm = pythiaParams.ms2ExtractionWidthPPM;
////        res.scanTimeStDev = pythiaParams.scanTimeWindowMinutes;
////        res.cosineSimAnchor = pythiaParams.cosineSimToAnchorThreshold;
//        res.fdrCount = targetCountAboveFDRQValueThreshold;
//        results.push_back(res);
//    }
//
//    e = getTopFrequencyParameters(
//            &results,
//            &m_pythiaParameters.ms2ExtractionWidthPPM
////            &m_pythiaParameters.scanTimeWindowMinutes
////            &m_pythiaParameters.cosineSimToAnchorThreshold
//            ); ree;
//
//    e = targetDecoyCandidatePairScoretron->setPythiaParameters(m_pythiaParameters); ree;
//
//    qDebug() << "Optimal ppm setting:" << m_pythiaParameters.ms2ExtractionWidthPPM;
//    qDebug() << "Optimal scanTimeWindow setting:" << m_pythiaParameters.scanTimeWindowMinutes;
//    qDebug() << "Optimal cosineSimSum setting:" << m_pythiaParameters.cosineSimToAnchorThreshold;
//
////    //TODO replace this with response surface derived DOE parameters when you figure out how to do it.
//////    const DOEResult bestParametersFDR = *std::max_element(
//////            results.rbegin(),
//////            results.rend(),
//////            [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
//////    );
//////
//////    m_pythiaParameters.ms2ExtractionWidthPPM = bestParametersFDR.ppm;
//////    m_pythiaParameters.scanTimeWindowMinutes = bestParametersFDR.scanTimeStDev;
//////    m_pythiaParameters.cosineSimToAnchorThreshold = bestParametersFDR.cosineSimAnchor;
////    ///////////////////////////////////////////////////////////////////////////////////////////////////
//
//    ERR_RETURN
//}
//
//Err PythiaDIAWorkflow::mainAnalysis(
//        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
//        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers,
//        int *psmCountOnePercentFDR
//        ) {
//
//    ERR_INIT
//
//    e = ErrorUtils::isTrue(targetDecoyCandidatePairScoretron->isInit()); ree;
//
//    m_pythiaParameters.print();
//
//
//    const bool useExtendedScores = true;
//    const bool useNeuralNetworkScores = true;
//
//    const int topNMs2IonsMainAnalysis = std::max(
//            m_minTopNMs2Ions,
//            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions))
//    );
//
//    qDebug() << "Using top:" << topNMs2IonsMainAnalysis << "fragments for main analysis";
//
//    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
//            m_pythiaParameters.mzMinDataStructure,
//            m_pythiaParameters.mzMaxDataStructure,
//            scoredTargetDecoyPointers
//    ); ree;
//
//    QMap<QString, int> fdrVsCount;
//    e = setTargetDecoyCandidateScores(
//            targetDecoyCandidatePairScoretron,
//            topNMs2IonsMainAnalysis,
//            useExtendedScores,
//            useNeuralNetworkScores,
//            scoredTargetDecoyPointers,
//            &fdrVsCount
//    ); ree;
//
//    *psmCountOnePercentFDR = fdrVsCount.value("1");
//
//    ERR_RETURN
//}
//
//namespace {
//
//    double getMinDiscScoreOfTargets(QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairs) {
//
//        TargetDecoyCandidatePair* discScoreMinElement = *std::min_element(
//                targetDecoyCandidatePairs->begin(),
//                targetDecoyCandidatePairs->end(),
//                [](TargetDecoyCandidatePair *l, TargetDecoyCandidatePair *r){return
//                        l->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore <
//                        r->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore;
//                }
//        );
//
//        return discScoreMinElement->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore;
//    }
//
//    Err getDecoyPointersAboveDiscScore(
//            double discScoreMin,
//            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
//            QVector<TargetDecoyCandidatePair*> *decoyCandidatePairs
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairs); ree;
//
//        *decoyCandidatePairs = targetDecoyCandidatePairs;
//
//        const auto terminatorLogic = [discScoreMin](TargetDecoyCandidatePair *tdcp){
//            return tdcp->candidateScoresBestDiscriminantScorePtrDecoy()->discriminateScore < discScoreMin;
//        };
//
//        const auto terminator = std::remove_if(
//                decoyCandidatePairs->begin(),
//                decoyCandidatePairs->end(),
//                terminatorLogic
//                );
//
//        decoyCandidatePairs->erase(terminator, decoyCandidatePairs->end());
//
//        ERR_RETURN
//    }
//
//    Err buildMs2IonsFromScoredCandidate(
//            bool getDecoyMS2Ions,
//            TargetDecoyCandidatePair* targetDecoyCandidatePair,
//            QVector<MS2Ion> *ms2Ions
//    ) {
//
//        ERR_INIT
//
//        QVector<double> theoIntensityVec;
//        QVector<double> mzSearchedVec;
//
//        if (getDecoyMS2Ions) {
//            theoIntensityVec = targetDecoyCandidatePair->candidateScoresBestDiscriminantScorePtrDecoy()->theoIntensityVec;
//            mzSearchedVec = targetDecoyCandidatePair->candidateScoresBestDiscriminantScorePtrDecoy()->mzSearchedVec;
//        }
//        else {
//            theoIntensityVec = targetDecoyCandidatePair->candidateScoresBestDiscriminantScorePtrTarget()->theoIntensityVec;
//            mzSearchedVec = targetDecoyCandidatePair->candidateScoresBestDiscriminantScorePtrTarget()->mzSearchedVec;
//        }
//
//        e = ErrorUtils::isNotEmpty(theoIntensityVec); ree;
//        e = ErrorUtils::isEqual(theoIntensityVec.size(), mzSearchedVec.size()); ree;
//
//        for (int i = 0; i < theoIntensityVec.size(); i++) {
//            const double intensity = theoIntensityVec.at(i);
//            const double mz = mzSearchedVec.at(i);
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
//
//        ERR_RETURN
//    }
//
//    Err buildTandemDeconvolutionInput(
//            const QVector<TargetDecoyCandidatePair*> &scoredTargetsPointersFDRThresholded,
//            const QVector<TargetDecoyCandidatePair*> &scoredDecoysPointersFDRThresholded,
//            QMap<ScanNumber, QMap<TargetDecoyCandidatePair*, QVector<MS2Ion>>> *scanNumberVsTandemPredictions
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(scoredTargetsPointersFDRThresholded); ree;
//        e = ErrorUtils::isNotEmpty(scoredDecoysPointersFDRThresholded); ree;
//
//        for (TargetDecoyCandidatePair *tdcp : scoredTargetsPointersFDRThresholded) {
//
//            const CandidateScores &bestCandidateScoresTarget = *tdcp->candidateScoresBestDiscriminantScorePtrTarget();
//
//            QVector<MS2Ion> ms2IonsTarget;
//            e = buildMs2IonsFromScoredCandidate(false, tdcp, &ms2IonsTarget); ree;
//            (*scanNumberVsTandemPredictions)[bestCandidateScoresTarget.scanNumber].insert(tdcp, ms2IonsTarget);
//        }
//
//        for (TargetDecoyCandidatePair *tdcp : scoredDecoysPointersFDRThresholded) {
//
//            const CandidateScores &bestCandidateScoresTarget = *tdcp->candidateScoresBestDiscriminantScorePtrDecoy();
//
//            QVector<MS2Ion> ms2IonsDecoy;
//            e = buildMs2IonsFromScoredCandidate(true, tdcp, &ms2IonsDecoy); ree;
//            (*scanNumberVsTandemPredictions)[bestCandidateScoresTarget.scanNumber].insert(tdcp, ms2IonsDecoy);
//        }
//
//        ERR_RETURN
//    }
//
//    struct DeconVol {
//        ScanNumber scanNumber = -1;
//        ScanPoints scanPoints;
//        QMap<TargetDecoyCandidatePair*, QVector<MS2Ion>> tandemPredictions;
//    };
//
//    struct DeconResult {
//        Err e = eNoError;
//        ScanNumber scanNumber = -1;
//        QMap<TargetDecoyCandidatePair*, TandemDeconvolverResult> tandemDeconvolverResult;
//    };
//
//    Err tandemDeconvolutionLogic(
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
//        ); ree;
//
//        QMap<TargetDecoyCandidatePair*, TandemDeconvolverResult> result;
//        e = tandemSpectraDeconvolvotron.deconvolveTandemSpectra(
//                deconVol.scanPoints,
//                deconVol.tandemPredictions,
//                &result
//        ); ree;
//
//        for (auto it = result.begin(); it != result.end(); it++) {
//
//            TargetDecoyCandidatePair* tdcp = it.key();
//            const TandemDeconvolverResult &tdr = it.value();
//
//            CandidateScores *candidateScoresTarget = tdcp->candidateScoresBestDiscriminantScorePtrTarget();
//            candidateScoresTarget->matrixPValue = tdr.pVal;
//            candidateScoresTarget->matrixError = tdr.frameError;
//            candidateScoresTarget->matrixWeight = tdr.weight;
//
//            CandidateScores *candidateScoresDecoy = tdcp->candidateScoresBestDiscriminantScorePtrDecoy();
//            candidateScoresDecoy->matrixPValue = tdr.pVal;
//            candidateScoresDecoy->matrixError = tdr.frameError;
//            candidateScoresDecoy->matrixWeight = tdr.weight;
//        }
//
//        ERR_RETURN
//    }
//
//    Err deconvolveTandemSpectra(
//            const QMap<ScanNumber, QMap<TargetDecoyCandidatePair*, QVector<MS2Ion>>> &scanNumberVsTandemPredictions,
//            const PythiaParameters &params,
//            MsReaderPointerAcc *msReaderPointerAcc
//            ) {
//
//        ERR_INIT
//        e = ErrorUtils::isNotEmpty(scanNumberVsTandemPredictions); ree;
//
//        QElapsedTimer et;
////        et.start();
////
////        const auto insertLogic = [msReaderPointerAcc, scanNumberVsTandemPredictions](const ScanNumber sn){
////
////            DeconVol deconVol;
////            deconVol.scanNumber = sn;
////            deconVol.tandemPredictions = scanNumberVsTandemPredictions.value(sn);
////            QPair<Err, ScanPoints*> scanPointsResult = msReaderPointerAcc->ptr->getScanPoints(sn);
////
////            deconVol.scanPoints = scanPointsResult.second;
////
////            return deconVol;
////        };
////
////        QVector<DeconVol> scanPointsVsTandemPredictions;
////        const QList<ScanNumber> &scanNumbers = scanNumberVsTandemPredictions.keys();
////        std::transform(
////                scanNumbers.begin(),
////                scanNumbers.end(),
////                std::back_inserter(scanPointsVsTandemPredictions),
////                insertLogic
////        );
////
////#define PARALLEL_DECONVOLVE
////#ifdef PARALLEL_DECONVOLVE
////        const auto deconvolutionLogicBinder = std::bind(
////            tandemDeconvolutionLogic,
////            std::placeholders::_1,
////            params
////        );
////
////        QFuture<Err> futures = QtConcurrent::mapped(
////                scanPointsVsTandemPredictions,
////                deconvolutionLogicBinder
////                );
////        futures.waitForFinished();
////
////        for (const Err &result : futures) {
////            e = result; ree;
////        }
////#else
////        for (const DeconVol &deconVol : scanPointsVsTandemPredictions) {
////
////            const DeconResult result = tandemDeconvolutionLogic(
////                    deconVol,
////                    params
////            );
////
////            e = result.e; ree;
////            (*scanNumberVsTandemDeconvolverResult)[result.scanNumber] = result.tandemDeconvolverResult;
////        }
////#endif
//
//        qDebug() << "Deconvolution finished in" << et.elapsed() << "mSec";
//
//        ERR_RETURN
//    }
//
//    void filterScoredCandidatesByWeightAndPVal(
//            double pValThreshold,
//            bool filterByDecoyScores,
//            QVector<TargetDecoyCandidatePair*> *TargetDecoyCandidatePairs
//        ) {
//
//        const double weightThreshold = 0.0;
//        const auto terminatorLogic = [weightThreshold, pValThreshold, filterByDecoyScores](TargetDecoyCandidatePair *tdcp){
//
//            CandidateScores *candidateScores;
//            if (filterByDecoyScores) {
//                candidateScores = tdcp->candidateScoresBestDiscriminantScorePtrTarget();
//            }
//            else {
//                candidateScores = tdcp->candidateScoresBestDiscriminantScorePtrDecoy();
//            }
//
//            return candidateScores->matrixWeight < weightThreshold || candidateScores->matrixPValue > pValThreshold;
//        };
//
//        const auto terminator = std::remove_if(
//                TargetDecoyCandidatePairs->begin(),
//                TargetDecoyCandidatePairs->end(),
//                terminatorLogic
//        );
//        TargetDecoyCandidatePairs->erase(terminator, TargetDecoyCandidatePairs->end());
//    }
//
//}//namespace
//Err PythiaDIAWorkflow::removeInterferingCandidates(
//        MsReaderPointerAcc *msReaderPointerAcc,
//        const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointers,
//        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointersUpdatedTargets
//        ) {
//
//    ERR_INIT
//
//    qDebug() << "Starting interference removal of shared tandem fragments";
//
//    const double fdrThreshold = 0.1;
//    e = FDRCLassifierNeuralNet::filterScoreCandidatesByFDR(
//            scoredTargetDecoyPointers,
//            fdrThreshold,
//            scoredTargetDecoyPointersUpdatedTargets
//            ); ree;
//
//    const double targetsDiscScoreMin = getMinDiscScoreOfTargets(scoredTargetDecoyPointersUpdatedTargets);
//
//    QVector<TargetDecoyCandidatePair*> scoredDecoysPointers50PercentFDR;
//    e = getDecoyPointersAboveDiscScore(
//            targetsDiscScoreMin,
//            scoredTargetDecoyPointers,
//            &scoredDecoysPointers50PercentFDR
//            ); ree;
//
//    QMap<ScanNumber, QMap<TargetDecoyCandidatePair*, QVector<MS2Ion>>> scanNumberVsTandemPredictions;
//    e = buildTandemDeconvolutionInput(
//            *scoredTargetDecoyPointersUpdatedTargets,
//            scoredDecoysPointers50PercentFDR,
//            &scanNumberVsTandemPredictions
//            ); ree;
//
//    e = deconvolveTandemSpectra(
//            scanNumberVsTandemPredictions,
//            m_pythiaParameters,
//            msReaderPointerAcc
//            ); ree;
//
//    filterScoredCandidatesByWeightAndPVal(
//            m_pythiaParameters.pValThreshold,
//            false,
//            scoredTargetDecoyPointersUpdatedTargets
//            );
//
//    QMap<QString, int> fdrVsCount;
//    FDRCLassifierNeuralNet::outputFDRResults(*scoredTargetDecoyPointersUpdatedTargets, true, &fdrVsCount);
//
//    ERR_RETURN
//}
//
//namespace {
//
//    void filterDecoys(
//            double discriminantScoreMin,
//            QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointersFiltered
//            ) {
//
//        const auto terminatorLogic = [discriminantScoreMin](TargetDecoyCandidatePair* p){
//            return p->candidateScoresBestDiscriminantScorePtrDecoy()->discriminateScore < discriminantScoreMin;
//        };
//
//        const auto terminator = std::remove_if(
//                scoredTargetDecoyPointersFiltered->begin(),
//                scoredTargetDecoyPointersFiltered->end(),
//                terminatorLogic
//                );
//
//        scoredTargetDecoyPointersFiltered->erase(terminator, scoredTargetDecoyPointersFiltered->end());
//    }
//
//    Err minMaxScaleScores(
//            const QVector<KarnnNNTarget> &karnnNNTargets,
//            QVector<KarnnNNTarget> *karnnNNTargetsNorm
//            ) {
//
//        ERR_INIT
//
//        e = ErrorUtils::isNotEmpty(karnnNNTargets); ree;
//
//        QVector<QVector<double>> vecs;
//        std::transform(
//                karnnNNTargets.begin(),
//                karnnNNTargets.end(),
//                std::back_inserter(vecs),
//                [](const KarnnNNTarget &kt){return kt.scoreVec;}
//                );
//
//        Eigen::MatrixX<double> mat = EigenUtils::convertQVectorsToEigenMatrix(vecs);
//        EigenUtils::minMaxScaleMatrix(&mat);
//
//        const QVector<QVector<double>> vecsNorm = EigenUtils::convertEigenMatrixToQVectors(mat);
//
//        e = ErrorUtils::isEqual(vecsNorm.size(), karnnNNTargets.size()); ree;
//
//        for (int i = 0; i < vecsNorm.size(); i++) {
//            KarnnNNTarget ktNew = karnnNNTargets.at(i);
//            ktNew.scoreVec = vecsNorm.at(i);
//            karnnNNTargetsNorm->push_back(ktNew);
//        }
//
//        ERR_RETURN
//    }
//
//}//namespace
//Err PythiaDIAWorkflow::applyNeuralNetClassifier(
//        const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointers,
//        const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointersFDRFiltered,
//        const QPair<double, double> &scanTimeMinMax,
//        int psmCountOnePercentFDR,
//        bool reportDecoys,
//        QVector<CandidateScores> *candidateScoreClassifier
//        ) {
//
//    ERR_INIT
//
//    e = ErrorUtils::isNotEmpty(scoredTargetDecoyPointers); ree;
//    e = ErrorUtils::isNotEmpty(scoredTargetDecoyPointersFDRFiltered); ree;
//
//    QVector<CandidateScores> candidateScoresTargetsAndDecoys;
//    for (TargetDecoyCandidatePair *tdcp : scoredTargetDecoyPointers) {
//
//        QMap<MzTargetKey, CandidateScores> *uniqueInfoScanKeyVsScoresTarget = tdcp->uniqueInfoScanKeyVsScoresTarget();
//        for (auto it = uniqueInfoScanKeyVsScoresTarget->begin(); it != uniqueInfoScanKeyVsScoresTarget->end(); it++) {
//            candidateScoresTargetsAndDecoys.push_back(it.value());
//        }
//
//        QMap<MzTargetKey, CandidateScores> *uniqueInfoScanKeyVsScoresDecoy = tdcp->uniqueInfoScanKeyVsScoresDecoy();
//        for (auto it = uniqueInfoScanKeyVsScoresDecoy->begin(); it != uniqueInfoScanKeyVsScoresDecoy->end(); it++) {
//            candidateScoresTargetsAndDecoys.push_back(it.value());
//        }
//    }
//
//    std::sort(
//            candidateScoresTargetsAndDecoys.rbegin(),
//            candidateScoresTargetsAndDecoys.rend(),
//            [](const CandidateScores &l, const CandidateScores &r){return l.discriminateScore < r.discriminateScore;}
//            );
//
//    TargetDecoyCandidatePair* discriminantScoreMinElement = *std::min_element(
//            scoredTargetDecoyPointersFDRFiltered.begin(),
//            scoredTargetDecoyPointersFDRFiltered.end(),
//            [](TargetDecoyCandidatePair *l, TargetDecoyCandidatePair *r){
//                return l->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore
//                        < r->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore;
//            });
//
//
//
//    int decoyCounter = 0;
//    int counter = 0;
//    for (const CandidateScores &cs : candidateScoresTargetsAndDecoys) {
//
//        counter++;
//        if (cs.isDecoy) {
//            decoyCounter++;
//        }
//
//        if (cs.discriminateScore < discriminantScoreMinElement->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore) {
//            break;
//        }
//    }
//
//    candidateScoresTargetsAndDecoys.resize(counter);
//
////#define BYPASS_NEURAL_NET
//#ifdef BYPASS_NEURAL_NET
//    std::sort(
//            candidateScoresTargetsAndDecoysShuffled.rbegin(),
//            candidateScoresTargetsAndDecoysShuffled.rend(),
//            [](const CandidateScores &l, const CandidateScores &r){return l.discriminateScore < r.discriminateScore;}
//            );
//
//    int counter = 0;
//    int decoyCounter = 0;
//    for (const CandidateScores &cs : candidateScoresTargetsAndDecoysShuffled) {
//
//        if (cs.isDecoy) {
//            decoyCounter++;
//        }
//
//        const double qVal = static_cast<double>(decoyCounter) / ++counter;
//        if (qVal > 0.01) {
//            break;
//        }
//
//        candidateScoreClassifier->push_back(cs);
//    }
//#else
//    std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
//    std::shuffle(candidateScoresTargetsAndDecoys.begin(), candidateScoresTargetsAndDecoys.end(),rng);
//
//    qDebug() << "target vs decoy count" << counter - decoyCounter << decoyCounter << "total" << counter;
//
//    QVector<KarnnNNTarget> karnnNNTargets;
//    for (int i = 0; i < candidateScoresTargetsAndDecoys.size(); i++) {
//        const CandidateScores &cs = candidateScoresTargetsAndDecoys.at(i);
//        KarnnNNTarget karnnNnTarget;
//        karnnNnTarget.seq = cs.peptideStringWithMods;
//        karnnNnTarget.isDecoy = cs.isDecoy;
//        karnnNnTarget.index = i;
//
//        e = FDRCLassifierNeuralNet::buildScoreVector(
//                cs,
//                true,
//                true,
//                m_pythiaParameters.topNMs2Ions,
//                scanTimeMinMax,
//                &karnnNnTarget.scoreVec
//                ); ree;
//        karnnNNTargets.push_back(karnnNnTarget);
//    }
//
//    QVector<KarnnNNTarget> karnnNNTargetsNorm;
//    e = minMaxScaleScores(karnnNNTargets, &karnnNNTargetsNorm); ree;
//
//    const int baggingSize = 12;
//    const int hiddenLayers = 5;
//    const float learningRate = 0.003;
//    const int epochs = 1;
//
//    const int batchSize = std::min(50, std::max(1, static_cast<int>(karnnNNTargetsNorm.size() / 100.0)));
//    const int maxIters = 1;
//    qDebug() << "Batch Size:" << batchSize << "MaxIters:" << maxIters;
//
//    QVector<QVector<float>> xData;
//    QVector<float> yData;
//    for (const KarnnNNTarget &kt : karnnNNTargetsNorm) {
//        xData.push_back(QVector<float>(kt.scoreVec.begin(), kt.scoreVec.end()));
//        yData.push_back(kt.isDecoy ? 1.0 : 0.0);
//    }
//
//    int cycles = 0;
//    counter = 0;
//
//    candidateScoreClassifier->clear();
//
//    FDRCLassifierNeuralNet fdrcLassifierNeuralNet;
//    e = fdrcLassifierNeuralNet.init(
//            epochs,
//            baggingSize,
//            batchSize,
//            learningRate
//            ); ree;
//
//    QVector<float> predictions;
//    e = fdrcLassifierNeuralNet.exec(xData, yData, &predictions); ree;
//
//    for (int i = 0; i < yData.size(); i++) {
//        karnnNNTargetsNorm[i].nnScore = predictions.at(i);
//    }
//
//    std::sort(
//            karnnNNTargetsNorm.begin(),
//            karnnNNTargetsNorm.end(),
//            [](const KarnnNNTarget &l, const KarnnNNTarget &r){return l.nnScore < r.nnScore;}
//    );
//
//    counter = 0;
//    int falsePositives = 0;
//    for (const KarnnNNTarget &rp : karnnNNTargetsNorm) {
//
//        CandidateScores candidateScoresNew = candidateScoresTargetsAndDecoys.at(rp.index);
//        candidateScoresNew.classifierScore = rp.nnScore;
//        candidateScoreClassifier->push_back(candidateScoresNew);
//
//        ++counter;
//
//        if (rp.nnScore > 0.5 || (falsePositives / static_cast<double>(counter)) > 0.008) {
//            if (!reportDecoys) {
//                break;
//            }
//            else {
//                counter--;
//                continue;
//            }
//        }
//
//        if (rp.isDecoy){
//            falsePositives++;
//        }
//    }
//
//    qDebug() << "False Pos" << falsePositives << "Total" << counter << "FDR 0.5 nnScore cuttoff" << falsePositives / (counter + 0.0);
//
//    e = setQValueForCandidates(candidateScoreClassifier); ree
//#endif
//
//    ERR_RETURN
//}
//
//Err PythiaDIAWorkflow::updateProteinGroupAnnotation(
//        const QString &fastaFilePath,
//        QVector<CandidateScores> *candidateScores
//) {
//
//    ERR_INIT
//
//    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
//
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
//    for (int i = 0; i < candidateScores->size(); i++) {
//
//        CandidateScores &cs = (*candidateScores)[i];
//
//        QString peptideSeqReplacedLeucines = cs.peptideStringWithMods;
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
//        cs.proteinGroup = fastaDescriptions.join(';');
//
//    }
//
//    ERR_RETURN
//}
