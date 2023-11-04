//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAWorkflow.h"

#include "MsReaderParquet.h"

//#include "CandidateClassifier.h"
#include "ClassifierWeightsManager.h"
//#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "FastaFileToPeptidesListWorkFlow.h"
#include "FastaReader.h"
#include "FDRCLassifierNeuralNet.h"
//#include "FeatureFinderHillBuilder.h"
#include "MathUtils.h"
#include "MS2ChargeDeconvolvotron.h"
//#include "MS2DataExtractomatic.h"
#include "MsCalibrationReader.h"
//#include "MsFrameScoretron.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"
//#include "PeakIntegratomatic.h"
#include "TandemSpectraDeconvolvotron.h"
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

Err PythiaDIAWorkflow::processFile(const QString &_msDataFilePath) {

    ERR_INIT

    m_msScanInfos.clear();
    m_scanTimeMinMax = {-1.0, -1.0};

    QString msDataFilePath = _msDataFilePath;

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath); ree;

    m_msScanInfos = msReaderPointerAcc.ptr->getUniqueTandemMsScanInfos();
    m_scanTimeMinMax = msReaderPointerAcc.ptr->scanTimeMinMax();

    if (m_pythiaParameters.deisotopeScans) {
        e = deisotopeScans(&msReaderPointerAcc); ree;
    }

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

    e = buildCalibration(
            0.2,
            false,
            &targetDecoyCandidatePairScoretron
            ); ree;

    e = optimizeParameters(&targetDecoyCandidatePairScoretron); ree;

    e = buildCalibration(
            0.2,
            true,
            &targetDecoyCandidatePairScoretron
            ); ree;

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;
    e = mainAnalysis(
            &targetDecoyCandidatePairScoretron,
            &scoredTargetDecoyPointers
            ); ree;

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointersFDRThresholded;
    const double fdrThreshold = 0.01;
    e = FDRCLassifierNeuralNet::filterScoreCandidatesByFDR(
            scoredTargetDecoyPointers,
            fdrThreshold,
            &scoredTargetDecoyPointersFDRThresholded
    ); ree;

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointersUpdatedTargets;
    e = removeInterferingCandidates(
            &msReaderPointerAcc,
            scoredTargetDecoyPointers,
            &scoredTargetDecoyPointersUpdatedTargets
            ); ree;

    e = updateProteinGroupAnnotation(
            m_fastaUri,
            &scoredTargetDecoyPointersUpdatedTargets
            ); ree;

    std::sort(scoredTargetDecoyPointersUpdatedTargets.rbegin(), scoredTargetDecoyPointersUpdatedTargets.rend(),
              [](TargetDecoyCandidatePair*l, TargetDecoyCandidatePair*r){
        return l->candidateScoresBestDiscriminantScorePtrTarget()->cosineSimSum100
        < r->candidateScoresBestDiscriminantScorePtrTarget()->cosineSimSum100;
    });

    QVector<CandidateScores> candidateScoresTargetsAndDecoys;
    for (TargetDecoyCandidatePair *tdcp : scoredTargetDecoyPointersUpdatedTargets) {
        candidateScoresTargetsAndDecoys.push_back(*tdcp->candidateScoresBestDiscriminantScorePtrTarget());
        candidateScoresTargetsAndDecoys.push_back(*tdcp->candidateScoresBestDiscriminantScorePtrDecoy());
    }

#define USE_NEURAL_NET_CLASSIFIER
#ifdef USE_NEURAL_NET_CLASSIFIER
    QVector<CandidateScores> scoredCandidatesClassifierUpdated;
    e = applyNeuralNetClassifier(
            candidateScoresTargetsAndDecoys,
            msReaderPointerAcc.ptr->scanTimeMinMax(),
            &scoredCandidatesClassifierUpdated
            ); ree;
#endif

    const QString resultsFilePath = msReaderPointerAcc.ptr->filePath() + S_GLOBAL_SETTINGS.DOT_PYTHIA_DIA_FILE_EXTENSION;
    e = ParquetReader::write(candidateScoresTargetsAndDecoys, resultsFilePath); ree;

    ERR_RETURN
}

namespace {

    QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> deisotopeLogic(const QVector<QPair<ScanNumber, ScanPoints>> &scanPointPairs, double ppmTol) {

        ERR_INIT

        const QString &chargeModelFilePath
                = QDir(qApp->applicationDirPath()).filePath("MS2_Charge_Model.json");

        const QString &monoModelFilePath
                = QDir(qApp->applicationDirPath()).filePath("MS2_Mono_Model.json");

        MS2ChargeDeconvolvotron ms2ChargeDeconvolvotron;
        e = ms2ChargeDeconvolvotron.init(chargeModelFilePath, monoModelFilePath, ppmTol);

        QVector<QPair<ScanNumber, ScanPoints>> deisotopedScanPoints;
        for (const QPair<ScanNumber, ScanPoints> &pr : scanPointPairs) {

            if (pr.first % 1000 == 0) {
                qDebug() << "Deisotoping" << pr.first;
            }

            ScanPoints scanPointsIterDeisotoped;
            e = ms2ChargeDeconvolvotron.deisotopeScanPoints(pr.second, &scanPointsIterDeisotoped); rree;

            deisotopedScanPoints.push_back({pr.first, scanPointsIterDeisotoped});
        }

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

    QVector<QVector<QPair<ScanNumber, ScanPoints>>> scanPointsTranched;
    ParallelUtils::trancheMapForParallelization(
            msReaderPointerAcc->ptr->getScanPoints(),
            ParallelUtils::numberOfAvailableSystemProcessors(),
            &scanPointsTranched
    );

    QMap<ScanNumber, ScanPoints> deisotopedScanPoints;

//#define PARALLEL_DEISO
#ifdef PARALLEL_DEISO
    QFuture<QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>>> futures
            = QtConcurrent::mapped(scanPointsTranched, deisotopeLogicBinder);
    futures.waitForFinished();

    for (const QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> &res : futures) {
        e = res.first; ree;
        for (const QPair<ScanNumber, ScanPoints> &r : res.second) {
            const ScanNumber scanNumber = r.first;
            const ScanPoints &scanPoints = r.second;
            deisotopedScanPoints.insert(scanNumber, scanPoints);
        }
    }
#else
    for (const QVector<QPair<ScanNumber, ScanPoints>> &spPairs :  scanPointsTranched) {

        QPair<Err, QVector<QPair<ScanNumber, ScanPoints>>> deisotopedScan = deisotopeLogic(
                spPairs,
                m_pythiaParameters.ms2ExtractionWidthPPM
                );

        e = deisotopedScan.first; ree;
        const QVector<QPair<ScanNumber, ScanPoints>> &res = deisotopedScan.second;
        for (const QPair<ScanNumber, ScanPoints> &r : res) {
            deisotopedScanPoints.insert(r.first, r.second);
        }
    }
#endif

    e = msReaderPointerAcc->ptr->setScanPoints(deisotopedScanPoints); ree;

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
        const double cosineSimSum100MinForTraining = 0.9;

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
                    continue;
                }
                e = ErrorUtils::isTrue(decoyContainsTargetKey); ree;

                const CandidateScores &candidateScores = tdc->uniqueInfoScanKeyVsScoresTarget()->value(key);
                if (candidateScores.cosineSimSum100 < cosineSimSum100MinForTraining) {
                    continue;
                }

                const ScoresTargets scoresTargets = FDRCLassifierNeuralNet::buildScoreVector(
                        candidateScores,
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
Err PythiaDIAWorkflow::buildCalibration(
        double calibrationTrainingFraction,
        bool useExtendedScores,
        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron
        ){

    ERR_INIT

    e = ErrorUtils::isTrue(targetDecoyCandidatePairScoretron->isInit()); ree;

    const bool useNeuralNetworkScores = false;
    const int minTrainingCount = 100;

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;

    double calFractionIter = calibrationTrainingFraction;
    while (calFractionIter < 1.0) {

        scoredTargetDecoyPointers.clear();

        e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
                m_pythiaParameters.mzMinDataStructure,
                m_pythiaParameters.mzMaxDataStructure,
                calFractionIter,
                &scoredTargetDecoyPointers
        ); ree;

        if (scoredTargetDecoyPointers.size() < m_pythiaParameters.trancheSizeMax) {
            calFractionIter *= 2;
            continue;
        }

        break;
    }

    QMap<QString, int> fdrVsCount;
    e = setTargetDecoyCandidateScores(
            targetDecoyCandidatePairScoretron,
            m_minTopNMs2Ions,
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
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers,
        QMap<QString, int> *fdrVsCount
        ) {

    ERR_INIT

    e = targetDecoyCandidatePairScoretron->scoreTargetDecoyPairs(
            topNMS2Ions,
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


    QVector<ScoresTargets> scoresTargets;
    QVector<ScoresDecoys> scoresDecoys;
    for(const TargetDecoyPairTargetKey &tdp : targetDecoyPairTargetKeys) {
        scoresTargets.push_back(tdp.scoresTargetVsScoresDecoys.first);
        scoresDecoys.push_back(tdp.scoresTargetVsScoresDecoys.second);
    }

    QVector<double> discScoreTargets;
    e = ClassifierWeightsManager::applyWeights(scoresTargets, weights, &discScoreTargets); ree;

    QVector<double> discScoreDecoys;
    e = ClassifierWeightsManager::applyWeights(scoresDecoys, weights, &discScoreDecoys); ree;

    e = ErrorUtils::isEqual(scoresTargets.size(), scoresDecoys.size()); ree;
    e = ErrorUtils::isEqual(scoresTargets.size(), targetDecoyPairTargetKeys.size()); ree;

    for (int i = 0; i < targetDecoyPairTargetKeys.size(); i++) {

        const TargetDecoyPairTargetKey &tdptk = targetDecoyPairTargetKeys.at(i);
        TargetDecoyCandidatePair *tdcp = tdptk.targetDecoyCandidatePair;

        (*tdcp->uniqueInfoScanKeyVsScoresTarget())[tdptk.uniqueMsInfoScanKey].discriminateScore = discScoreTargets.at(i);
        (*tdcp->uniqueInfoScanKeyVsScoresDecoy())[tdptk.uniqueMsInfoScanKey].discriminateScore = discScoreDecoys.at(i);
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

        const QVector<QVector<double>> experiments = {
                {2.5, 1.0},
                {3.5, 1.0},
                {2.5,  3.0},
                {3.5,  3.0},
                {2.5,  2.0},
                {3.5,  2.0}
        };

        for (const QVector<double> &exp : experiments) {

            PythiaParameters params = pythiaParameters;
            params.ms2ExtractionWidthPPM = mzPPMStDev * exp.at(0);
            params.scanTimeWindowMinutes = scanTimeStDev * exp.at(1);

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
            qDebug() << "ppmTol" << pp.ms2ExtractionWidthPPM
                     << "scanTimeWinMin" << pp.scanTimeWindowMinutes
                     << "cosineSimSumAnchThreshold" << pp.cosineSimToAnchorThreshold;
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

        e = ErrorUtils::toDouble(cv.front().first, value); ree;

        ERR_RETURN
    }

    Err getTopFrequencyParameters(
            QVector<DOEResult> *results,
            double *ppmSetting,
            double *scanTimeWidthSetting
//            double *cosineSimSetting
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(*results); ree;

        std::sort(
                results->rbegin(),
                results->rend(),
                [](const DOEResult &l, const DOEResult &r){return l.fdrCount < r.fdrCount;}
                );

        for (const DOEResult &r : *results) {
//            qDebug() << r.ppm << r.scanTimeStDev << r.cosineSimAnchor << r.fdrCount;
            qDebug() << r.ppm << r.scanTimeStDev << r.fdrCount;
        }

        const int topNResults = 3;
        results->resize(topNResults);

        QMap<QString, int> ppmCounts;
        QMap<QString, int> scanTimeWidthCounts;
//        QMap<QString, int> cosineSimCounts;
        for (const DOEResult &r : *results) {

            const QString ppmString = QString::number(r.ppm);
            const QString scanTimeWidthString = QString::number(r.scanTimeStDev);
//            const QString cosineSimString = QString::number(r.cosineSimAnchor);

            ppmCounts[ppmString]++;
            scanTimeWidthCounts[scanTimeWidthString]++;
//            cosineSimCounts[cosineSimString]++;
        }

        const QVector<QPair<QString, int>> ppmCountsVector = ParallelUtils::convertMapToVectorPairs(ppmCounts);
        const QVector<QPair<QString, int>> scanTimeWidthCountsVector = ParallelUtils::convertMapToVectorPairs(scanTimeWidthCounts);
//        const QVector<QPair<QString, int>> cosineSimCountsVector = ParallelUtils::convertMapToVectorPairs(cosineSimCounts);

        e = findMostFrequentValue(ppmCountsVector, ppmSetting); ree;
        e = findMostFrequentValue(scanTimeWidthCountsVector, scanTimeWidthSetting); ree;
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

    QVector<TargetDecoyCandidatePair*> scoredTargetDecoyPointers;

    double selectionFractionValueIter = selectionFractionValue;
    while (selectionFractionValueIter < 1.0) {

        scoredTargetDecoyPointers.clear();

        e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
                m_pythiaParameters.mzMinDataStructure,
                m_pythiaParameters.mzMaxDataStructure,
                selectionFractionValueIter,
                &scoredTargetDecoyPointers
        ); ree;

        if (scoredTargetDecoyPointers.size() < m_pythiaParameters.trancheSizeMax) {
            selectionFractionValueIter *= 2;
            continue;
        }

        break;
    }

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

        qDebug() << "ppmTol" << pythiaParams.ms2ExtractionWidthPPM
                 << "scanTimeWinMin" << pythiaParams.scanTimeWindowMinutes
                 << "cosineSimSumAnchThreshold" << pythiaParams.cosineSimToAnchorThreshold;

        e = targetDecoyCandidatePairScoretron->setPythiaParameters(pythiaParams); ree;

        QMap<QString, int> fdrVsCount;
        e = setTargetDecoyCandidateScores(
                targetDecoyCandidatePairScoretron,
                topNMs2IonsOptimization,
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
//        res.cosineSimAnchor = pythiaParams.cosineSimToAnchorThreshold;
        res.fdrCount = targetCountAboveFDRQValueThreshold;
        results.push_back(res);
    }

    e = getTopFrequencyParameters(
            &results,
            &m_pythiaParameters.ms2ExtractionWidthPPM,
            &m_pythiaParameters.scanTimeWindowMinutes
//            &m_pythiaParameters.cosineSimToAnchorThreshold
            ); ree;

    e = targetDecoyCandidatePairScoretron->setPythiaParameters(m_pythiaParameters); ree;

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
        TargetDecoyCandidatePairScoretron *targetDecoyCandidatePairScoretron,
        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(targetDecoyCandidatePairScoretron->isInit()); ree;

    m_pythiaParameters.print();


    const bool useExtendedScores = true;
    const bool useNeuralNetworkScores = false;

    const int topNMs2IonsMainAnalysis = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_pythiaParameters.topNMs2Ions))
    );

    qDebug() << "Using top:" << topNMs2IonsMainAnalysis << "fragments for main analysis";

    e = m_targetDecoyCandidatePairManager.getTargetDecoyCandidatePairPointers(
            m_pythiaParameters.mzMinDataStructure,
            m_pythiaParameters.mzMaxDataStructure,
            scoredTargetDecoyPointers
    ); ree;

    QMap<QString, int> fdrVsCount;
    e = setTargetDecoyCandidateScores(
            targetDecoyCandidatePairScoretron,
            topNMs2IonsMainAnalysis,
            useExtendedScores,
            useNeuralNetworkScores,
            scoredTargetDecoyPointers,
            &fdrVsCount
    ); ree;

    ERR_RETURN
}

namespace {

    double getMinDiscScoreOfTargets(QVector<TargetDecoyCandidatePair*> *targetDecoyCandidatePairs) {

        TargetDecoyCandidatePair* discScoreMinElement = *std::min_element(
                targetDecoyCandidatePairs->begin(),
                targetDecoyCandidatePairs->end(),
                [](TargetDecoyCandidatePair *l, TargetDecoyCandidatePair *r){return
                        l->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore <
                        r->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore;
                }
        );

        return discScoreMinElement->candidateScoresBestDiscriminantScorePtrTarget()->discriminateScore;
    }

    Err getDecoyPointersAboveDiscScore(
            double discScoreMin,
            const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
            QVector<TargetDecoyCandidatePair*> *decoyCandidatePairs
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(targetDecoyCandidatePairs); ree;

        *decoyCandidatePairs = targetDecoyCandidatePairs;

        const auto terminatorLogic = [discScoreMin](TargetDecoyCandidatePair *tdcp){
            return tdcp->candidateScoresBestDiscriminantScorePtrDecoy()->discriminateScore < discScoreMin;
        };

        const auto terminator = std::remove_if(
                decoyCandidatePairs->begin(),
                decoyCandidatePairs->end(),
                terminatorLogic
                );

        decoyCandidatePairs->erase(terminator, decoyCandidatePairs->end());

        ERR_RETURN
    }

    Err buildMs2IonsFromScoredCandidate(
            bool getDecoyMS2Ions,
            TargetDecoyCandidatePair* targetDecoyCandidatePair,
            QVector<MS2Ion> *ms2Ions
    ) {

        ERR_INIT

        QVector<double> theoIntensityVec;
        QVector<double> mzSearchedVec;

        if (getDecoyMS2Ions) {
            theoIntensityVec = targetDecoyCandidatePair->candidateScoresBestDiscriminantScorePtrDecoy()->theoIntensityVec;
            mzSearchedVec = targetDecoyCandidatePair->candidateScoresBestDiscriminantScorePtrDecoy()->mzSearchedVec;
        }
        else {
            theoIntensityVec = targetDecoyCandidatePair->candidateScoresBestDiscriminantScorePtrTarget()->theoIntensityVec;
            mzSearchedVec = targetDecoyCandidatePair->candidateScoresBestDiscriminantScorePtrTarget()->mzSearchedVec;
        }

        e = ErrorUtils::isNotEmpty(theoIntensityVec); ree;
        e = ErrorUtils::isEqual(theoIntensityVec.size(), mzSearchedVec.size()); ree;

        for (int i = 0; i < theoIntensityVec.size(); i++) {
            const double intensity = theoIntensityVec.at(i);
            const double mz = mzSearchedVec.at(i);

            if (MathUtils::tZero(mz) || MathUtils::tZero(intensity)) {
                continue;
            }

            MS2Ion ms2Ion;
            ms2Ion.mz = mz;
            ms2Ion.intensity = intensity;

            ms2Ions->push_back(ms2Ion);
        }

        ERR_RETURN
    }

    Err buildTandemDeconvolutionInput(
            const QVector<TargetDecoyCandidatePair*> &scoredTargetsPointersFDRThresholded,
            const QVector<TargetDecoyCandidatePair*> &scoredDecoysPointersFDRThresholded,
            QMap<ScanNumber, QMap<TargetDecoyCandidatePair*, QVector<MS2Ion>>> *scanNumberVsTandemPredictions
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredTargetsPointersFDRThresholded); ree;
        e = ErrorUtils::isNotEmpty(scoredDecoysPointersFDRThresholded); ree;

        for (TargetDecoyCandidatePair *tdcp : scoredTargetsPointersFDRThresholded) {

            const CandidateScores &bestCandidateScoresTarget = *tdcp->candidateScoresBestDiscriminantScorePtrTarget();

            QVector<MS2Ion> ms2IonsTarget;
            e = buildMs2IonsFromScoredCandidate(false, tdcp, &ms2IonsTarget); ree;
            (*scanNumberVsTandemPredictions)[bestCandidateScoresTarget.scanNumber].insert(tdcp, ms2IonsTarget);
        }

        for (TargetDecoyCandidatePair *tdcp : scoredDecoysPointersFDRThresholded) {

            const CandidateScores &bestCandidateScoresTarget = *tdcp->candidateScoresBestDiscriminantScorePtrDecoy();

            QVector<MS2Ion> ms2IonsDecoy;
            e = buildMs2IonsFromScoredCandidate(true, tdcp, &ms2IonsDecoy); ree;
            (*scanNumberVsTandemPredictions)[bestCandidateScoresTarget.scanNumber].insert(tdcp, ms2IonsDecoy);
        }

        ERR_RETURN
    }

    struct DeconVol {
        ScanNumber scanNumber = -1;
        ScanPoints scanPoints;
        QMap<TargetDecoyCandidatePair*, QVector<MS2Ion>> tandemPredictions;
    };

    struct DeconResult {
        Err e = eNoError;
        ScanNumber scanNumber = -1;
        QMap<TargetDecoyCandidatePair*, TandemDeconvolverResult> tandemDeconvolverResult;
    };

    Err tandemDeconvolutionLogic(
            const DeconVol &deconVol,
            const PythiaParameters &params
    ) {

        ERR_INIT

        const int maxIteration= 20;
        const double stopTol = 0.000000001;
        const double pValThreshold = 0.05;

        DeconResult deconResult;

        TandemSpectraDeconvolvotron tandemSpectraDeconvolvotron;
        e = tandemSpectraDeconvolvotron.init(
                S_GLOBAL_SETTINGS.HASHING_PRECISION,
                params.mzMaxDataStructure,
                params.ms2ExtractionWidthPPM,
                maxIteration,
                stopTol,
                pValThreshold
        ); ree;

        QMap<TargetDecoyCandidatePair*, TandemDeconvolverResult> result;
        e = tandemSpectraDeconvolvotron.deconvolveTandemSpectra(
                deconVol.scanPoints,
                deconVol.tandemPredictions,
                &result
        ); ree;

        for (auto it = result.begin(); it != result.end(); it++) {

            TargetDecoyCandidatePair* tdcp = it.key();
            const TandemDeconvolverResult &tdr = it.value();

            CandidateScores *candidateScoresTarget = tdcp->candidateScoresBestDiscriminantScorePtrTarget();
            candidateScoresTarget->matrixPValue = tdr.pVal;
            candidateScoresTarget->matrixError = tdr.frameError;
            candidateScoresTarget->matrixWeight = tdr.weight;

            CandidateScores *candidateScoresDecoy = tdcp->candidateScoresBestDiscriminantScorePtrDecoy();
            candidateScoresDecoy->matrixPValue = tdr.pVal;
            candidateScoresDecoy->matrixError = tdr.frameError;
            candidateScoresDecoy->matrixWeight = tdr.weight;
        }

        ERR_RETURN
    }

    Err deconvolveTandemSpectra(
            const QMap<ScanNumber, QMap<TargetDecoyCandidatePair*, QVector<MS2Ion>>> &scanNumberVsTandemPredictions,
            const PythiaParameters &params,
            MsReaderPointerAcc *msReaderPointerAcc
            ) {

        ERR_INIT
        e = ErrorUtils::isNotEmpty(scanNumberVsTandemPredictions); ree;

        QElapsedTimer et;
        et.start();

        const auto insertLogic = [msReaderPointerAcc, scanNumberVsTandemPredictions](const ScanNumber sn){

            DeconVol deconVol;
            deconVol.scanNumber = sn;
            deconVol.tandemPredictions = scanNumberVsTandemPredictions.value(sn);
            msReaderPointerAcc->ptr->getScanPoints(sn, &deconVol.scanPoints);

            return deconVol;
        };

        QVector<DeconVol> scanPointsVsTandemPredictions;
        const QList<ScanNumber> &scanNumbers = scanNumberVsTandemPredictions.keys();
        std::transform(
                scanNumbers.begin(),
                scanNumbers.end(),
                std::back_inserter(scanPointsVsTandemPredictions),
                insertLogic
        );

#define PARALLEL_DECONVOLVE
#ifdef PARALLEL_DECONVOLVE
        const auto deconvolutionLogicBinder = std::bind(
            tandemDeconvolutionLogic,
            std::placeholders::_1,
            params
        );

        QFuture<Err> futures = QtConcurrent::mapped(
                scanPointsVsTandemPredictions,
                deconvolutionLogicBinder
                );
        futures.waitForFinished();

        for (const Err &result : futures) {
            e = result; ree;
        }
#else
        for (const DeconVol &deconVol : scanPointsVsTandemPredictions) {

            const DeconResult result = tandemDeconvolutionLogic(
                    deconVol,
                    params
            );

            e = result.e; ree;
            (*scanNumberVsTandemDeconvolverResult)[result.scanNumber] = result.tandemDeconvolverResult;
        }
#endif

        qDebug() << "Deconvolution finished in" << et.elapsed() << "mSec";

        ERR_RETURN
    }

    void filterScoredCandidatesByWeightAndPVal(
            double pValThreshold,
            bool filterByDecoyScores,
            QVector<TargetDecoyCandidatePair*> *TargetDecoyCandidatePairs
        ) {

        const double weightThreshold = 0.0;
        const auto terminatorLogic = [weightThreshold, pValThreshold, filterByDecoyScores](TargetDecoyCandidatePair *tdcp){

            CandidateScores *candidateScores;
            if (filterByDecoyScores) {
                candidateScores = tdcp->candidateScoresBestDiscriminantScorePtrTarget();
            }
            else {
                candidateScores = tdcp->candidateScoresBestDiscriminantScorePtrDecoy();
            }

            return candidateScores->matrixWeight < weightThreshold || candidateScores->matrixPValue > pValThreshold;
        };

        const auto terminator = std::remove_if(
                TargetDecoyCandidatePairs->begin(),
                TargetDecoyCandidatePairs->end(),
                terminatorLogic
        );
        TargetDecoyCandidatePairs->erase(terminator, TargetDecoyCandidatePairs->end());
    }

}//namespace
Err PythiaDIAWorkflow::removeInterferingCandidates(
        MsReaderPointerAcc *msReaderPointerAcc,
        const QVector<TargetDecoyCandidatePair*> &scoredTargetDecoyPointers,
        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointersUpdatedTargets
        ) {

    ERR_INIT

    qDebug() << "Starting interference removal of shared tandem fragments";

    const double fdrThreshold = 0.2;
    e = FDRCLassifierNeuralNet::filterScoreCandidatesByFDR(
            scoredTargetDecoyPointers,
            fdrThreshold,
            scoredTargetDecoyPointersUpdatedTargets
            ); ree;

    const double targetsDiscScoreMin = getMinDiscScoreOfTargets(scoredTargetDecoyPointersUpdatedTargets);

    QVector<TargetDecoyCandidatePair*> scoredDecoysPointers50PercentFDR;
    e = getDecoyPointersAboveDiscScore(
            targetsDiscScoreMin,
            scoredTargetDecoyPointers,
            &scoredDecoysPointers50PercentFDR
            ); ree;

    QMap<ScanNumber, QMap<TargetDecoyCandidatePair*, QVector<MS2Ion>>> scanNumberVsTandemPredictions;
    e = buildTandemDeconvolutionInput(
            *scoredTargetDecoyPointersUpdatedTargets,
            scoredDecoysPointers50PercentFDR,
            &scanNumberVsTandemPredictions
            ); ree;

    e = deconvolveTandemSpectra(
            scanNumberVsTandemPredictions,
            m_pythiaParameters,
            msReaderPointerAcc
            ); ree;

    filterScoredCandidatesByWeightAndPVal(
            m_pythiaParameters.pValThreshold,
            false,
            scoredTargetDecoyPointersUpdatedTargets
            );

    QMap<QString, int> fdrVsCount;
    FDRCLassifierNeuralNet::outputFDRResults(*scoredTargetDecoyPointersUpdatedTargets, true, &fdrVsCount);

    ERR_RETURN
}

Err PythiaDIAWorkflow::applyNeuralNetClassifier(
        const QVector<CandidateScores> &candidateScoresTargetsAndDecoys,
        const QPair<double, double> &scanTimeMinMax,
        QVector<CandidateScores> *candidateScoreClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoys); ree;

    const long decoyCount = std::count_if(
            candidateScoresTargetsAndDecoys.begin(),
            candidateScoresTargetsAndDecoys.end(), [](const CandidateScores &c){return c.isDecoy;}
            );

    const long targetCount = std::count_if(
            candidateScoresTargetsAndDecoys.begin(),
            candidateScoresTargetsAndDecoys.end(), [](const CandidateScores &c){return !c.isDecoy;}
            );

    QVector<CandidateScores> candidateScoresTargetsAndDecoysShuffled = candidateScoresTargetsAndDecoys;
    std::mt19937 rng(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
    std::shuffle(candidateScoresTargetsAndDecoysShuffled.begin(), candidateScoresTargetsAndDecoysShuffled.end(),rng);

    qDebug() << "target vs decoy count" << targetCount << decoyCount;

    const int epochs = 5;
    const int baggingSize = 6;
    const double batchFraction = 0.01;
    const double learningRate = 0.001;
    FDRCLassifierNeuralNet fdrClassifierNeuralNet;
    e = fdrClassifierNeuralNet.init(
            m_pythiaParameters.topNMs2Ions,
            epochs,
            baggingSize,
            batchFraction,
            learningRate,
            scanTimeMinMax
            ); ree;

    e = fdrClassifierNeuralNet.exec(
            candidateScoresTargetsAndDecoysShuffled,
            candidateScoreClassifier
            ); ree;

    ERR_RETURN
}

Err PythiaDIAWorkflow::updateProteinGroupAnnotation(
        const QString &fastaFilePath,
        QVector<TargetDecoyCandidatePair*> *scoredTargetDecoyPointers
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(*scoredTargetDecoyPointers); ree;

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

    for (int i = 0; i < scoredTargetDecoyPointers->size(); i++) {

        TargetDecoyCandidatePair *tdcp = scoredTargetDecoyPointers->at(i);

        QString peptideSeqReplacedLeucines = tdcp->peptideStringWithMods();
        peptideSeqReplacedLeucines = peptideSeqReplacedLeucines.replace('L', 'J').replace('I', 'J');

        const QVector<FastaEntry> &fastaEntries = peptideStringWithModsVsFastaEntriesLeucinesReplaced.value(peptideSeqReplacedLeucines);

        QStringList fastaDescriptions;
        std::transform(
                fastaEntries.begin(),
                fastaEntries.end(),
                std::back_inserter(fastaDescriptions),
                [](const FastaEntry &fe){return fe.fastaDescription;}
                );

        CandidateScores *bestScoreTarget = tdcp->candidateScoresBestDiscriminantScorePtrTarget();
        CandidateScores *bestScoreDecoy = tdcp->candidateScoresBestDiscriminantScorePtrDecoy();

        bestScoreTarget->proteinGroup = fastaDescriptions.join(';');
        bestScoreTarget->isDecoy = false;
        bestScoreTarget->iRTPredicted = tdcp->iRt();
        bestScoreTarget->targetKey = tdcp->bestDiscriminateScoreKeyTarget();

        bestScoreDecoy->proteinGroup = "DECOY";
        bestScoreDecoy->isDecoy = true;
        bestScoreDecoy->iRTPredicted = tdcp->iRt();
        bestScoreDecoy->peptideStringWithMods = tdcp->peptideStringWithMods() + "_decoy";
        bestScoreDecoy->targetKey = tdcp->bestDiscriminateScoreKeyDecoy();
    }

    ERR_RETURN
}
