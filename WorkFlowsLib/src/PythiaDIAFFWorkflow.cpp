//
// Created by anichols on 2/28/23.
//

#include "PythiaDIAFFWorkflow.h"

#include "ClassifierWeightsManager.h"
#include "FDRCLassifierNeuralNet.h"
#include "FragLibReader.h"
#include "MsReaderPointerAcc.h"


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

Err PythiaDIAFFWorkflow::buildCalibration(MsReaderPointerAcc *msReaderPointerAcc) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_targetDecoyCandidatePairManager.isInit()); ree;

    const double calibrationTrainingFraction = 0.2;
    const bool useExtendedScores = false;
    const bool useNeuralNetworkScores = false;
    const int minTrainingCountTranche = 50;

    QMap<MzTargetKey, QVector<TargetDecoyCandidatePair*>> mzTargetKeyVsTargetDecoyCandidatePointers;
    e = buildUniqueInfoScanKeyVsTargetDecoyCandidatePointers(
            msReaderPointerAcc->ptr->getUniqueTandemMsScanInfos(),
            calibrationTrainingFraction,
            &mzTargetKeyVsTargetDecoyCandidatePointers
            ); ree;

    const int topNMS2IonsCalibration = 6;

    QVector<CandidateScores> candidateScores;
    e = m_targetDecoyCandidatePairScoretron.scoreTargetDecoyPairs(
            topNMS2IonsCalibration,
            m_msCalibratomatic,
            &mzTargetKeyVsTargetDecoyCandidatePointers,
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
                msScanInfo.precursorTargetMz - msScanInfo.isoWindowLower,
                msScanInfo.precursorTargetMz + msScanInfo.isoWindowLower,
                selectionFraction,
                &targetDecoyPointers
                ); ree;

        mzTargetKeyVsTargetDecoyCandidatePointers->insert(msScanInfo.targetKey(), targetDecoyPointers);
        qDebug() << "MzTargetKey" << msScanInfo.targetKey() << targetDecoyPointers.size() << "targetDecoyPairs found";
    }

    ERR_RETURN
}

namespace {

//    Err buildClassifierInput(
//        bool useExtendedScores,
//        bool useNeuralNetworkScores,
//        int theoMzIonsSize,
//        const QPair<double, double> &scanTimeMinMax,
//        QVector<CandidateScores> *candidateScores,
//        QVector<QPair<ScoresTargets, ScoresDecoys>> *targetsScoresVsDecoyScores
//        ) {
//
//        ERR_INIT
//
//        const int theoMzIonsSizeMin = 6;
//        const double cosineSimSum100MinForTraining = 0.9;
//
//        e = ErrorUtils::isFalse(candidateScores->isEmpty()); ree;
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

}//namespace
Err PythiaDIAFFWorkflow::setDiscriminantScoreForCandidates(
        const QPair<double, double> &scanTimeMinMax,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int theoMzIonsSize,
        QVector<CandidateScores> *candidateScores
        ) {

    ERR_INIT

    const auto sortLogic = [](const CandidateScores &l, const CandidateScores &r){

        if (l.peptideStringWithMods == r.peptideStringWithMods) {

            if (l.charge == r.charge) {

                if (l.targetKey == r.targetKey) {
                    return l.isDecoy < r.isDecoy;
                }

                return l.targetKey < r.targetKey;
            }

            return l.charge < r.charge;
        }

        return l.peptideStringWithMods < r.peptideStringWithMods;
    };
    std::sort(candidateScores->begin(), candidateScores->end(), sortLogic);

    for (const CandidateScores &cs : *candidateScores) {
        qDebug() << cs.peptideStringWithMods << cs.charge << cs.targetKey << cs.isDecoy;
    }


//    e = buildClassifierInput(
//            useExtendedScores,
//            useNeuralNetworkScores,
//            theoMzIonsSize,
//            scanTimeMinMax,
//            candidateScores,
//            &targetsScoresVsDecoyScores
//    ); ree;

//    QVector<CandidateScores*> candidateScoresTargetsPtrs;
//    QVector<ScoresTargets> scoresTargets;
//    QVector<CandidateScores*> candidateScoresDecoysPtrs;
//    QVector<ScoresDecoys> scoresDecoys;
//    for(CandidateScores &cs : *candidateScores) {
//
//        QVector<double> scoreVector;
//        e = FDRCLassifierNeuralNet::buildScoreVector(
//                cs,
//                useExtendedScores,
//                useNeuralNetworkScores,
//                theoMzIonsSize,
//                scanTimeMinMax,
//                &scoreVector
//        ); ree;
//
//        if (cs.isDecoy) {
//            scoresDecoys.push_back(scoreVector);
//            candidateScoresDecoysPtrs.push_back(&cs);
//            continue;
//        }
//
//        scoresTargets.push_back(scoreVector);
//        candidateScoresTargetsPtrs.push_back(&cs);
//    }
//
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
//
//    QVector<double> discScoreTargets;
//    e = ClassifierWeightsManager::applyWeights(scoresTargets, weights, &discScoreTargets); ree;
//
//    QVector<double> discScoreDecoys;
//    e = ClassifierWeightsManager::applyWeights(scoresDecoys, weights, &discScoreDecoys); ree;
//
//    e = ErrorUtils::isEqual(scoresTargets.size(), scoresDecoys.size()); ree;
//    e = ErrorUtils::isEqual(scoresTargets.size(), candidateScoresTargetsPtrs.size()); ree;
//    e = ErrorUtils::isEqual(scoresDecoys.size(), candidateScoresDecoysPtrs.size()); ree;
//    e = ErrorUtils::isEqual(discScoreTargets.size(), scoresTargets.size());
//    e = ErrorUtils::isEqual(discScoreDecoys.size(), scoresDecoys.size());
//
//    for (int i = 0; i < scoresDecoys.size(); i++) {
//
//        CandidateScores* csTarget = candidateScoresTargetsPtrs.at(i);
//        csTarget->discriminateScore = discScoreTargets.at(i);
//
//        CandidateScores* csDecoy = candidateScoresDecoysPtrs.at(i);
//        csDecoy->discriminateScore = discScoreDecoys.at(i);
//
//    }

    ERR_RETURN


}
