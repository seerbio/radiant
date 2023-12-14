//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAFFWorkflow.h"

#include "ClassifierWeightsManager.h"
#include "FDRCLassifierNeuralNet.h"
#include "FragLibReader.h"
#include "MsReaderPointerAcc.h"
#include "ParallelUtils.h"

#include <QElapsedTimer>

PythiaDIAFFWorkflow::PythiaDIAFFWorkflow(){}

Err PythiaDIAFFWorkflow::init(
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

Err PythiaDIAFFWorkflow::processFile(const QString &msDataFilePath) {

    ERR_INIT

    e = ErrorUtils::fileExists(msDataFilePath); ree;

    MsReaderPointerAcc msReaderPointerAcc;
    e = msReaderPointerAcc.openFile(msDataFilePath); ree;

    QMap<MzTargetKey, QMap<ScanNumber, ScanPoints*>> diaTargetFrames;
    e = msReaderPointerAcc.ptr->collateMS2MzTargetFrames(
            &diaTargetFrames
            ); ree;

    const int msLevel = 1;
    QMap<ScanNumber, ScanPoints> scanNumberVsScanTimeMS1;
    e = msReaderPointerAcc.ptr->getScanPoints(msLevel, &scanNumberVsScanTimeMS1); ree;

    e = m_targetDecoyCandidatePairScoretron.init(
            m_pythiaParameters,
            scanNumberVsScanTimeMS1,
            &msReaderPointerAcc,
            &diaTargetFrames
            ); ree;

    e = buildCalibration(&msReaderPointerAcc); ree;

    ERR_RETURN
}

namespace {

    Err buildMsCalibrationReaderRows(
            const QVector<CandidateScores> &candidateScores,
            QVector<MsCalibarationReaderRow> *msCalibrationReaderRows
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(candidateScores); ree;

        qDebug() << candidateScores.size() << "Found for recalibartion";

        const auto msCalibrationReaderRowsInsertLogic = [](const CandidateScores &cs){

            MsCalibarationReaderRow row;
            row.peptideStringWithMods = cs.peptideStringWithMods;
            row.iRTPredicted = static_cast<float>(cs.iRTPredicted);
            row.scanTime = cs.scanTime;
            row.mzSearchedVec = cs.mzSearchedVec;
            row.mzFoundMeanVec = cs.mzFoundMeanVec;
            row.mzFoundStDevVec = cs.mzFoundStDevVec;
            return row;
        };

        std::transform(
                candidateScores.begin(),
                candidateScores.end(),
                std::back_inserter(*msCalibrationReaderRows),
                msCalibrationReaderRowsInsertLogic
        );

        ERR_RETURN
    }


}//namespace
Err PythiaDIAFFWorkflow::buildCalibration(MsReaderPointerAcc *msReaderPointerAcc) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    const double calibrationTrainingFraction = -1.0;
    const bool useExtendedScores = false;
    const bool useNeuralNetworkScores = false;
    const int minTrainingCountTranche = 50;

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos(),
            calibrationTrainingFraction,
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

    const int numberOfTrances = 5;
    QVector<QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>>> mzTargetKeyVsTargetDecoyCandidatePointersTranched;
    e = ParallelUtils::trancheMapValueVectorsByKeyForParallelization(
            mzTargetKeyVsTargetDecoyCandidatePointers,
            numberOfTrances,
            &mzTargetKeyVsTargetDecoyCandidatePointersTranched
            ); ree;

    const int topNMS2IonsCalibration = 6;

    QVector<double> timeWindowStDevs;
    for (QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> &tranche : mzTargetKeyVsTargetDecoyCandidatePointersTranched) {

        QVector<CandidateScores> candidateScores;
        e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
                topNMS2IonsCalibration,
                m_msCalibratomatic,
                &tranche,
                &candidateScores
        ); ree;

        const QPair<double, double> scanTimeMinMax = msReaderPointerAcc->ptr->scanTimeMinMax();
        e = setDiscriminantScoreForCandidates(
                scanTimeMinMax,
                useExtendedScores,
                useNeuralNetworkScores,
                topNMS2IonsCalibration,
                &candidateScores
        ); ree;

        e = setQValueForCandidates(QValueScoreType::DiscriminantScore, &candidateScores); ree;

        QMap<QString, int> fdrVsCount;
        e = FDRCLassifierNeuralNet::outputFDRResults(
                candidateScores,
                true,
                &fdrVsCount
        ); ree;

        const double fdrThreshold = 0.1;
        int targetCountBelowFDRThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                candidateScores,
                fdrThreshold,
                &targetCountBelowFDRThreshold
        ); ree;

        std::sort(candidateScores.rbegin(), candidateScores.rend(), [](const CandidateScores &l, const CandidateScores &r){
            return l.discriminateScore < r.discriminateScore;
        });

        const int minTrainingCount = std::max(minTrainingCountTranche, targetCountBelowFDRThreshold);
        qDebug() << "Training RT count 10% FDR:" << minTrainingCount;

        candidateScores.resize(minTrainingCount);
        QVector<MsCalibarationReaderRow> msCalibrationReaderRows;
        e = buildMsCalibrationReaderRows(candidateScores, &msCalibrationReaderRows); ree;

        e = m_msCalibratomatic.init(msCalibrationReaderRows); ree;
        const int numberOfStDevs = 3;
        timeWindowStDevs.push_back(m_msCalibratomatic.scanTimeStDev(1));
        qDebug() << "scanTimeWindowStDev x 3:" << m_msCalibratomatic.scanTimeStDev(numberOfStDevs);

    }

    qDebug() << timeWindowStDevs << MathUtils::mean(timeWindowStDevs);
    m_msCalibratomatic.setScanTimeStDev(MathUtils::mean(timeWindowStDevs));

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
                msScanInfo.precursorTargetMz + (msScanInfo.isoWindowLower + m_pythiaParameters.precursorExtractionWindowThomsons),
                selectionFraction,
                &targetDecoyPointers
                ); ree;

        mzTargetKeyVsTargetDecoyCandidatePointers->insert(msScanInfo.targetKey(), targetDecoyPointers);
        qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyPointers.size() << "targetDecoyPairs found";
    }

    ERR_RETURN
}

namespace {

    struct TargetDecoyCandidateScores {
        CandidateScores *candidateScoresTarget;
        CandidateScores *candidateScoresDecoy;
        ScoresTargets scoresTarget;
        ScoresDecoys  scoresDecoy;
    };

    QString buildCandidateKey(const CandidateScores &candidateScores) {
        const QString decoyToString = candidateScores.isDecoy ? "_1" : "_0";
        return candidateScores.peptideStringWithMods + QString::number(candidateScores.charge) + candidateScores.targetKey + decoyToString;
    }

}//namespace
Err PythiaDIAFFWorkflow::setDiscriminantScoreForCandidates(
        const QPair<double, double> &scanTimeMinMax,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int theoMzIonsSize,
        QVector<CandidateScores> *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    QHash<QString, TargetDecoyCandidateScores> keyVsTargetDecoyCandidateScores;
    for(CandidateScores &cs : *candidateScores) {

        const QString key = cs.peptideStringWithMods + QString::number(cs.charge) + cs.targetKey;

        QVector<double> scoreVector;
        e = CandidateScores::buildScoreVector(
                cs,
                useExtendedScores,
                useNeuralNetworkScores,
                theoMzIonsSize,
                scanTimeMinMax,
                &scoreVector
        ); ree;

        if (cs.isDecoy) {
            keyVsTargetDecoyCandidateScores[key].scoresDecoy = scoreVector;
            keyVsTargetDecoyCandidateScores[key].candidateScoresDecoy = &cs;
            continue;
        }

        keyVsTargetDecoyCandidateScores[key].candidateScoresTarget = &cs;
        keyVsTargetDecoyCandidateScores[key].scoresTarget = scoreVector;
    }

    QVector<ScoresTargets> scoresTargets;
    QVector<ScoresDecoys> scoresDecoys;
    QVector<CandidateScores*> candidateScoresTargetsPtrs;
    QVector<CandidateScores*> candidateScoresDecoysPtrs;

    for (auto it = keyVsTargetDecoyCandidateScores.begin(); it != keyVsTargetDecoyCandidateScores.end(); it++) {

        const QString &key = it.key();
        const TargetDecoyCandidateScores &tdcs = it.value();

        e = ErrorUtils::isEqual(tdcs.scoresTarget.size(), tdcs.scoresDecoy.size());
        if (e != eNoError) {
            qDebug() << "target decoys not paired for key" << key;
            rrr(eValueError);
        }

        scoresTargets.push_back(tdcs.scoresTarget);
        scoresDecoys.push_back(tdcs.scoresDecoy);
        candidateScoresTargetsPtrs.push_back(tdcs.candidateScoresTarget);
        candidateScoresDecoysPtrs.push_back(tdcs.candidateScoresDecoy);
    }

    QVector<QVector<double>> A;
    QVector<double> b;

    e = ClassifierWeightsManager::buildDataClassifier1(
            scoresTargets,
            scoresDecoys,
            &A,
            &b
    ); ree;

    QVector<double> weights;
    e = ClassifierWeightsManager::fitWeights(A, b, &weights); ree;

    qDebug() << "Weights:" << weights;
    qDebug() << "b:" << b;

    QVector<double> discScoreTargets;
    e = ClassifierWeightsManager::applyWeights(scoresTargets, weights, &discScoreTargets); ree;

    QVector<double> discScoreDecoys;
    e = ClassifierWeightsManager::applyWeights(scoresDecoys, weights, &discScoreDecoys); ree;

    e = ErrorUtils::isEqual(scoresTargets.size(), scoresDecoys.size()); ree;
    e = ErrorUtils::isEqual(scoresTargets.size(), candidateScoresTargetsPtrs.size()); ree;
    e = ErrorUtils::isEqual(scoresDecoys.size(), candidateScoresDecoysPtrs.size()); ree;
    e = ErrorUtils::isEqual(discScoreTargets.size(), scoresTargets.size());
    e = ErrorUtils::isEqual(discScoreDecoys.size(), scoresDecoys.size());

    for (int i = 0; i < scoresDecoys.size(); i++) {

        CandidateScores* csTarget = candidateScoresTargetsPtrs.at(i);
        csTarget->discriminateScore = discScoreTargets.at(i);

        CandidateScores* csDecoy = candidateScoresDecoysPtrs.at(i);
        csDecoy->discriminateScore = discScoreDecoys.at(i);

    }

    ERR_RETURN

}

namespace {

    Err buildsetQValueForCandidateScoresInputs(
            const QValueScoreType &qValueScoreType,
            QVector<CandidateScores> *candidateScores,
            QHash<PeptideSequenceChargeKey, double> *identifierVsTargets,
            QHash<PeptideSequenceChargeKey, double> *identifierVsDecoys

    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
        identifierVsTargets->clear();
        identifierVsDecoys->clear();

        for (const CandidateScores &cs : *candidateScores) {


            const PeptideSequenceChargeKey peptideSequenceChargeKey = buildCandidateKey(cs);

            double classifierScore = -1.0;
            if (qValueScoreType == QValueScoreType::NNClassifierScore) {
                classifierScore = cs.classifierScore;
            }
            else {
                classifierScore = cs.discriminateScore;
            }


            if (cs.isDecoy) {
                identifierVsDecoys->insert(peptideSequenceChargeKey, classifierScore);
                continue;
            }

            identifierVsTargets->insert(peptideSequenceChargeKey, classifierScore);

        }

        ERR_RETURN
    }

    Err setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            const QHash<PeptideSequenceChargeKey, double> &identifierVsQValue,
            const QHash<PeptideSequenceChargeKey, double> &identifierVsDecoyRatio,
            QVector<CandidateScores> *candidateScores
    ) {

        ERR_INIT

        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

        for (CandidateScores &cs : *candidateScores) {

            if (cs.isDecoy) {
                continue;
            }

            const PeptideSequenceChargeKey peptideSequenceChargeKey = buildCandidateKey(cs);

            e = ErrorUtils::isTrue(identifierVsQValue.contains(peptideSequenceChargeKey)); ree;
            e = ErrorUtils::isTrue(identifierVsDecoyRatio.contains(peptideSequenceChargeKey)); ree;

            cs.qValue = identifierVsQValue.value(peptideSequenceChargeKey);
            cs.decoyRatio = identifierVsDecoyRatio.value(peptideSequenceChargeKey);
        }

        ERR_RETURN
    }

}//namespace
Err PythiaDIAFFWorkflow::setQValueForCandidates(
        const QValueScoreType &qValueScoreType,
        QVector<CandidateScores> *candidateScores
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;

    QElapsedTimer et;
    et.start();

    QHash<PeptideSequenceChargeKey, double> identifierVsTargets;
    QHash<PeptideSequenceChargeKey, double> identifierVsDecoys;
    e = buildsetQValueForCandidateScoresInputs(
            qValueScoreType,
            candidateScores,
            &identifierVsTargets,
            &identifierVsDecoys
    ); ree;

    QHash<PeptideSequenceChargeKey, double> identifierVsQValue;
    QHash<PeptideSequenceChargeKey, double> identifierVsDecoyRatio;
    e = MathUtils::calculateQValue(
            identifierVsTargets,
            identifierVsDecoys,
            &identifierVsQValue,
            &identifierVsDecoyRatio
    ); ree;

    e = setQValueAndDecoyRatioToTargetDecoyCandidatePairs(
            identifierVsQValue,
            identifierVsDecoyRatio,
            candidateScores
    ); ree;

    ERR_RETURN
}
