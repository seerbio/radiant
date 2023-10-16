//
// Created by anichols on 9/30/23.
//

#include "MS2DataExtractomatic.h"

#include "CandidateProcessertron.h"
#include "ClassifierWeightsManager.h"
#include "FDRCLassifierNeuralNet.h"
#include "FragLibReader.h"
#include "MsFrameScoretron.h"

#include <QtConcurrent/QtConcurrent>


MS2DataExtractomatic::MS2DataExtractomatic()
: m_topNMs2Ions(-1)
, m_useExtendedScores(false)
, m_useNeuralNetworkScores(false)
{}

Err MS2DataExtractomatic::init(
        const PythiaParameters &pythiaParameters,
        int topNMs2Ions,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        MsReaderPointerAcc *msReaderPointerAcc
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(msReaderPointerAcc->ptr->isInit()); ree;
    e = ErrorUtils::isTrue(pythiaParameters.isValid()); ree;
    e = ErrorUtils::isTrue(topNMs2Ions > 0); ree;

    m_pythiaParameters = pythiaParameters;
    m_msReaderPointerAcc = msReaderPointerAcc;
    m_topNMs2Ions = topNMs2Ions;
    m_useExtendedScores = useExtendedScores;
    m_useNeuralNetworkScores = useNeuralNetworkScores;

    ERR_RETURN
}

Err MS2DataExtractomatic::init(
        const PythiaParameters &pythiaParameters,
        int topNMs2Ions,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        MsReaderPointerAcc *msReaderPointerAcc,
        const MsCalibratomatic &msCalibratomatic
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(msCalibratomatic.isInit()); ree;

    e = init(
            pythiaParameters,
            topNMs2Ions,
            useExtendedScores,
            useNeuralNetworkScores,
            msReaderPointerAcc
            ); ree;

    m_msCalibratomatic = msCalibratomatic;

    ERR_RETURN
}

Err MS2DataExtractomatic::extractMS2ForCandidates(
        const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptide,
        QVector<ScoredCandidate> *scoredCandidatesAll
        ) {

    ERR_INIT

    QVector<ScoredCandidate> scoredCandidates;
    e = extractTargetDecoyData(
            uniqueInfoScanKeyVsCandidatePeptide,
            &scoredCandidates
    ); ree;

    e = buildScoredCandidatesFDR(
            scoredCandidates,
            m_msReaderPointerAcc->ptr->scanTimeMinMax(),
            scoredCandidatesAll
    ); ree;

    ERR_RETURN
}

namespace {

    class ParallelProcessingInput {

    public:
        ParallelProcessingInput() = default;
        ~ParallelProcessingInput() = default;

        UniqueMsInfoScanKey uniqueMsInfoScanKey;
        PythiaParameters pythiaParameters;
        MsCalibratomatic msCalibratomatic;
        QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
        QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
        QMap<ScanNumber , ScanPoints> *scanPoints = nullptr;
        QMap<PeptideStringWithMods, CandidatePeptide> peptideStringWithModsVsCandidatePeptide;
        int topNMS2Ions = -1;
    };

    QPair<Err, QVector<ScoredCandidate>> parallelProciessingLogic(const ParallelProcessingInput &ppi) {

        ERR_INIT

        QElapsedTimer et;
        et.start();

        MsFrameScoretron msFrameScoretron;

        e = msFrameScoretron.init(
                ppi.uniqueMsInfoScanKey,
                ppi.pythiaParameters,
                ppi.topNMS2Ions,
                *ppi.scanPoints,
                ppi.scanNumberVsScanTimeMS1,
                ppi.peptideStringWithModsVsCandidatePeptide,
                ppi.scanNumberVsScanTime,
                ppi.msCalibratomatic
        ); rree;

        QVector<ScoredCandidate> scoredCandidates;
        e = msFrameScoretron.scoreFrameCandidates(&scoredCandidates); rree;

        qDebug() << "Processed" << ppi.uniqueMsInfoScanKey << "Count" << scoredCandidates.size() << et.elapsed() << "mSec";

        return {e, scoredCandidates};
    }

}//namespace
Err MS2DataExtractomatic::extractTargetDecoyData(
        const QMap<UniqueMsInfoScanKey, QMap<PeptideStringWithMods, CandidatePeptide>> &uniqueInfoScanKeyVsCandidatePeptideCalibration,
        QVector<ScoredCandidate> *combinedResults
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(uniqueInfoScanKeyVsCandidatePeptideCalibration); ree;

    combinedResults->clear();

    QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> uniqueInfoScanKeyVsScanPoints;
    QMap<ScanNumber, ScanTime> scanNumberVsScanTime;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
    e = buildUniqueMsInfoScanKeyVsScanPoints(
            m_msReaderPointerAcc,
            &uniqueInfoScanKeyVsScanPoints,
            &scanNumberVsScanTime,
            &scanNumberVsScanTimeMS1
    ); ree;

    QVector<ParallelProcessingInput> parallelProcessingInputs;
    for (const UniqueMsInfoScanKey &uniqueMsInfoScanKey : uniqueInfoScanKeyVsCandidatePeptideCalibration.keys()) {
        ParallelProcessingInput ppi;
        ppi.uniqueMsInfoScanKey = uniqueMsInfoScanKey;
        ppi.scanNumberVsScanTime = scanNumberVsScanTime;
        ppi.scanNumberVsScanTimeMS1 = scanNumberVsScanTimeMS1;
        ppi.msCalibratomatic = m_msCalibratomatic;
        ppi.pythiaParameters = m_pythiaParameters;
        ppi.scanPoints = &uniqueInfoScanKeyVsScanPoints[uniqueMsInfoScanKey];
        ppi.peptideStringWithModsVsCandidatePeptide = uniqueInfoScanKeyVsCandidatePeptideCalibration[uniqueMsInfoScanKey];
        ppi.topNMS2Ions = m_topNMs2Ions;

        parallelProcessingInputs.push_back(ppi);
    }

#define PROCESS_PARALLEL
#ifdef PROCESS_PARALLEL
    QFuture<QPair<Err, QVector<ScoredCandidate>>> futures = QtConcurrent::mapped(
            parallelProcessingInputs,
            parallelProciessingLogic
    );
    futures.waitForFinished();

    for (const QPair<Err, QVector<ScoredCandidate>> &res : futures) {

        e = res.first; ree;
        combinedResults->append(res.second);
    }
#else
    for (const ParallelProcessingInput &inp : parallelProcessingInputs) {

        if (inp.uniqueMsInfoScanKey != "695066") {
            continue;
        }

        const QPair<Err, QVector<ScoredCandidate>> res = parallelProciessingLogic(inp);
        e = res.first; ree;
        combinedResults->append(res.second);
    }
#endif

    ERR_RETURN
}

namespace {

    Err buildClassifierInput(
            const QVector<ScoredCandidate> &scoredCandidates,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            int theoMzIonsSize,
            const QPair<double, double> &scanTimeMinMax,
            QVector<QVector<double>> *targetScoresVector,
            QVector<QVector<double>> *decoyScoresVector
    ) {

        ERR_INIT

        QMap<PeptideStringWithMods, ScoredCandidate> targets;
        QMap<PeptideStringWithMods, ScoredCandidate> decoys;
        for (const ScoredCandidate &sc : scoredCandidates) {

            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    sc.peptideStringWithMods,
                    sc.targetKey,
                    sc.charge
            );

            if(sc.isDecoy) {
                decoys.insert(key, sc);
            }
            else {
                targets.insert(key, sc);
            }
        }

        e = ErrorUtils::isEqual(targets.size(), decoys.size()); ree;

        for (auto it = targets.begin(); it != targets.end(); it++) {

            const ScoredCandidate &scTarget = it.value();
            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    scTarget.peptideStringWithMods,
                    scTarget.targetKey,
                    scTarget.charge
            );

            e = ErrorUtils::isTrue(decoys.contains(key)); ree;

            const ScoredCandidate scDecoy = decoys.value(key);

            const QVector<double> targetScores = FDRCLassifierNeuralNet::buildScoreVector(
                    scTarget,
                    useExtendedScores,
                    useNeuralNetworkScores,
                    theoMzIonsSize,
                    scanTimeMinMax
            );

            const QVector<double> decoyScores = FDRCLassifierNeuralNet::buildScoreVector(
                    scDecoy,
                    useExtendedScores,
                    useNeuralNetworkScores,
                    theoMzIonsSize,
                    scanTimeMinMax
            );

            targetScoresVector->push_back(targetScores);
            decoyScoresVector->push_back(decoyScores);
        }

        e = ErrorUtils::isEqual(targetScoresVector->size(), decoyScoresVector->size()); ree;

        ERR_RETURN
    }

}//namespace
Err MS2DataExtractomatic::buildScoredCandidatesFDR(
        const QVector<ScoredCandidate> &scoredCandidates,
        const QPair<double, double> &scanTimeMinMax,
        QVector<ScoredCandidate> *scoredCandidatesAll
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scoredCandidates); ree;

    scoredCandidatesAll->clear();

    QVector<QVector<double>> targetScoresVector;
    QVector<QVector<double>> decoyScoresVector;
    e = buildClassifierInput(
            scoredCandidates,
            m_useExtendedScores,
            m_useNeuralNetworkScores,
            m_topNMs2Ions,
            scanTimeMinMax,
            &targetScoresVector,
            &decoyScoresVector
    ); ree;

    QVector<QVector<double>> A;
    QVector<double> b;
    e = ClassifierWeightsManager::buildDataClassifier1(
            targetScoresVector,
            decoyScoresVector,
            &A,
            &b
    ); ree;

    QVector<double> weights;
    weights.reserve(b.size());

    e = fitWeightsLogic(
            scoredCandidates,
            scanTimeMinMax,
            A,
            b,
            &weights,
            scoredCandidatesAll
    ); ree;

    int psmCountTenPercentFDR;
    e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
            *scoredCandidatesAll,
            0.1,
            &psmCountTenPercentFDR
    ); ree;

    qDebug() << "Adjusted weights:" << weights;
    qDebug() << "*******";
    qDebug() << "Averages:" << b;

    QMap<QString, int> unused;
    e = MS2DataExtractomatic::outputFDRResults(*scoredCandidatesAll, &unused); ree;

    ERR_RETURN
}

Err MS2DataExtractomatic::buildUniqueMsInfoScanKeyVsScanPoints(
        MsReaderPointerAcc *msReaderPointerAcc,
        QMap<UniqueMsInfoScanKey, QMap<ScanNumber, ScanPoints>> *diaTargetFrames,
        QMap<ScanNumber, ScanTime> *scanNumberVsScanTime,
        QMap<ScanNumber, ScanPoints > *scanNumberVsScanTimeMS1
) {

    ERR_INIT

    e = msReaderPointerAcc->ptr->collateTandemPrecursorTargetsDIA(diaTargetFrames); ree

    const QMap<ScanNumber, MsScanInfo> msScanInfos = msReaderPointerAcc->ptr->getMsScanInfos();
    for (auto it = msScanInfos.begin(); it != msScanInfos.end(); it++) {
        scanNumberVsScanTime->insert(it.key(), it.value().scanTime);
    }

    const int msLevel = 1;
    e = msReaderPointerAcc->ptr->getScanPoints(msLevel, scanNumberVsScanTimeMS1); ree;

    ERR_RETURN
}

namespace {

    Err separateScoredTargetDecoys(
            const QVector<ScoredCandidate> &scoredCandidatesCalibration,
            const QVector<double> &weights,
            bool useExtendedScores,
            bool useNeuralNetworkScores,
            const int theoMzIonsSize,
            const QPair<double, double> &scanTimeMinMax,
            QMap<PeptideStringWithMods, ScoredCandidate> *peptideStringWithModsVsScoredCandidateTargets,
            QMap<PeptideStringWithMods, ScoredCandidate> *peptideStringWithModsVsScoredCandidateDecoys,
            QMap<PeptideStringWithMods, double> *peptideStringWithModsVsDiscScoreTargets,
            QMap<PeptideStringWithMods, double> *peptideStringWithModsVsDiscScoreDecoys
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidatesCalibration); ree;
        e = ErrorUtils::isNotEmpty(weights); ree;

        for (const ScoredCandidate &sc : scoredCandidatesCalibration) {

            QVector<double> scores = FDRCLassifierNeuralNet::buildScoreVector(
                    sc,
                    useExtendedScores,
                    useNeuralNetworkScores,
                    theoMzIonsSize,
                    scanTimeMinMax
            );

            QVector<double> results;
            e = ClassifierWeightsManager::applyWeights({scores}, weights, &results); ree;

            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    sc.peptideStringWithMods,
                    sc.targetKey,
                    sc.charge
            );

            ScoredCandidate scoredCandidate = sc;
            scoredCandidate.discriminateScore = results.front();

            if(sc.isDecoy) {
                peptideStringWithModsVsDiscScoreDecoys->insert(key, scoredCandidate.discriminateScore);
                peptideStringWithModsVsScoredCandidateDecoys->insert(key, scoredCandidate);
                continue;
            }

            peptideStringWithModsVsDiscScoreTargets->insert(key, scoredCandidate.discriminateScore);
            peptideStringWithModsVsScoredCandidateTargets->insert(key, scoredCandidate);
        }

        ERR_RETURN
    }

    Err buildScoredCandidatesAll(
            const QMap<PeptideStringWithMods, ScoredCandidate> &peptideStringWithModsVsScoredCandidateTargets,
            const QMap<PeptideStringWithMods, ScoredCandidate> &peptideStringWithModsVsScoredCandidateDecoys,
            QVector<ScoredCandidate> *scoredCandidatesAll
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(peptideStringWithModsVsScoredCandidateTargets); ree;
        e = ErrorUtils::isNotEmpty(peptideStringWithModsVsScoredCandidateDecoys); ree;

        for (const ScoredCandidate &sc : peptideStringWithModsVsScoredCandidateTargets) {
            scoredCandidatesAll->push_back(sc);
        }

        for (const ScoredCandidate &sc : peptideStringWithModsVsScoredCandidateDecoys) {
            ScoredCandidate scNew = sc;
            scNew.peptideStringWithMods += (static_cast<QString>(S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP) + "DECOY");
            scoredCandidatesAll->push_back(scNew);
        }

        ERR_RETURN
    }

    Err assignQValuesAndDecoyRatios(
            const QMap<QString, double> &identifierVsQValue,
            const QMap<QString, double> &identifierVsDecoyRatio,
            QVector<ScoredCandidate> *scoredCandidatesQValsUpdated
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(identifierVsDecoyRatio); ree;
        e = ErrorUtils::isNotEmpty(identifierVsQValue); ree;
        e = ErrorUtils::isNotEmpty(*scoredCandidatesQValsUpdated); ree;

        for (int i = 0; i < scoredCandidatesQValsUpdated->size(); i++) {

            ScoredCandidate &sc = (*scoredCandidatesQValsUpdated)[i];

            if (sc.isDecoy) {
                continue;
            }

            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    sc.peptideStringWithMods,
                    sc.targetKey,
                    sc.charge
            );

            e = ErrorUtils::isTrue(identifierVsQValue.contains(key)); ree;
            e = ErrorUtils::isTrue(identifierVsDecoyRatio.contains(key)); ree;

            sc.qValue = identifierVsQValue.value(key);
            sc.decoyRatio = identifierVsDecoyRatio.value(key);
        }

        ERR_RETURN
    }

}//namespace
Err MS2DataExtractomatic::fitWeightsLogic(
        const QVector<ScoredCandidate> &scoredCandidatesCalibration,
        const QPair<double, double> &scanTimeMinMax,
        const QVector<QVector<double>> &A,
        const QVector<double> &b,
        QVector<double> *weights,
        QVector<ScoredCandidate> *scoredCandidatesAll
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scoredCandidatesCalibration); ree;

    e = ClassifierWeightsManager::fitWeights(A, b, weights); ree;

    QMap<PeptideStringWithMods, ScoredCandidate> peptideStringWithModsVsScoredCandidateTargets;
    QMap<PeptideStringWithMods, ScoredCandidate> peptideStringWithModsVsScoredCandidateDecoys;
    QMap<PeptideStringWithMods, double> peptideStringWithModsVsDiscScoreTargets;
    QMap<PeptideStringWithMods, double> peptideStringWithModsVsDiscScoreDecoys;
    e = separateScoredTargetDecoys(
            scoredCandidatesCalibration,
            *weights,
            m_useExtendedScores,
            m_useNeuralNetworkScores,
            m_topNMs2Ions,
            scanTimeMinMax,
            &peptideStringWithModsVsScoredCandidateTargets,
            &peptideStringWithModsVsScoredCandidateDecoys,
            &peptideStringWithModsVsDiscScoreTargets,
            &peptideStringWithModsVsDiscScoreDecoys
    ); ree;

    QMap<QString, double> identifierVsQValue;
    QMap<QString, double> identifierVsDecoyRatio;
    e = MathUtils::calculateQValue(
            peptideStringWithModsVsDiscScoreTargets,
            peptideStringWithModsVsDiscScoreDecoys,
            &identifierVsQValue,
            &identifierVsDecoyRatio
    ); ree;

    e = buildScoredCandidatesAll(
            peptideStringWithModsVsScoredCandidateTargets,
            peptideStringWithModsVsScoredCandidateDecoys,
            scoredCandidatesAll
    ); ree;

    e = assignQValuesAndDecoyRatios(
            identifierVsQValue,
            identifierVsDecoyRatio,
            scoredCandidatesAll
    ); ree;

    ERR_RETURN
}

Err MS2DataExtractomatic::outputFDRResults(
        const QVector<ScoredCandidate> &scoredCandidatesAll,
        QMap<QString, int> *fdrVsCount
        ) {

    ERR_INIT

    const QVector<double> fdrFractions = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005};
    for (double fdrThresh : fdrFractions) {
        int foundAtThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                scoredCandidatesAll,
                fdrThresh,
                &foundAtThreshold
        ); ree;
        const double fdrPercent = fdrThresh * 100;
        fdrVsCount->insert(QString::number(fdrPercent), foundAtThreshold);
        qDebug() << foundAtThreshold << "PSMs found at" << fdrPercent  << "% FDR";
    }

    ERR_RETURN
}

Err MS2DataExtractomatic::filterScoreCandidatesByFDR(
        const QVector<ScoredCandidate> &scoredCandidatesAll,
        double qValueThreshold,
        bool filterDecoys,
        QVector<ScoredCandidate> *scoredCandidatesTargetsFDRThresholded
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;
    e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

    scoredCandidatesTargetsFDRThresholded->clear();

    for (const ScoredCandidate &sc : scoredCandidatesAll) {

        if (filterDecoys && sc.isDecoy) {
            continue;
        }

        else if (sc.qValue > qValueThreshold) {
            continue;
        }

        scoredCandidatesTargetsFDRThresholded->push_back(sc);
    }

    ERR_RETURN
}