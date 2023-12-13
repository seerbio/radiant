//
// Created by anichols on 9/26/23.
//

#include "FDRCLassifierNeuralNet.h"

#include "BiophysicalCalcs.h"
#include "CandidateScores.h"
#include "EigenUtils.h"
#include "ParallelUtils.h"
#include "TargetDecoyCandidatePair.h"

#include <random>


FDRCLassifierNeuralNet::FDRCLassifierNeuralNet()
: m_epochs(-1)
, m_batchSize(-1)
, m_learningRate(-1.0)
, m_minTopNMs2Ions(6)
, m_isInit(false)
, m_topNMs2Ions(-1)
{}

FDRCLassifierNeuralNet::~FDRCLassifierNeuralNet() {
    for (CandidateClassifier *cc : m_candidateClassifiers) {
        delete cc;
    }
}

Err FDRCLassifierNeuralNet::init(
        int epochs,
        int baggingSize,
        int batchSize,
        double learningRate
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(epochs > 0); ree;
    e = ErrorUtils::isTrue(baggingSize >= 1); ree;
    e = ErrorUtils::isTrue(batchSize > 0); ree;
    e = ErrorUtils::isTrue(learningRate > 0 & learningRate < 1 & !MathUtils::tZero(learningRate)); ree;

    m_epochs = epochs;
    m_baggingSize = baggingSize;
    m_batchSize = batchSize;
    m_learningRate = learningRate;

    m_isInit = true;

    ERR_RETURN
}

namespace {

    Err buildNeuralNetworkInput(
            const QMap<QString, CandidateScores> &keyVsScoredCandidates,
            int topNMs2IonsFull,
            const QPair<double, double> &scanTimeMinMax,
            QVector<NeuralNetData> *trainingData
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(keyVsScoredCandidates); ree;

        const bool useExtendedScores = true;
        const bool useNeuralNetworkScores = false;

        QVector<QPair<CandidateScores, QVector<double>>> targetPairs;
        for (const CandidateScores &sc: keyVsScoredCandidates) {

            QVector<double> scoreVector;
            e = CandidateScores::buildScoreVector(
                    sc,
                    useExtendedScores,
                    useNeuralNetworkScores,
                    topNMs2IonsFull,
                    scanTimeMinMax,
                    &scoreVector
            ); ree;

            targetPairs.push_back({sc, scoreVector});
        }

        std::mt19937 g(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);
        std::shuffle(targetPairs.begin(),targetPairs.end(), g);

        for (const QPair<CandidateScores, QVector<double>> &scPair : targetPairs) {
            NeuralNetData td;
            td.peptideStringWithMods = scPair.first.peptideStringWithMods;
            td.isDecoy = scPair.first.isDecoy;
            td.discScore = scPair.first.discriminateScore;
            td.scores = scPair.second;
            td.charge = scPair.first.charge;
            td.targetKey = scPair.first.targetKey;
            trainingData->push_back(td);
        }

//#define WRITE_NN_TRAIN_DATA
#ifdef WRITE_NN_TRAIN_DATA
        const QString dataFilePath = "/home/anichols/Desktop/Data/MsData/EXP22092_2022ms0742X32_A.raw.mzML.prq.nnTrainingData";
        e = ParquetReader::write(*trainingData, dataFilePath); ree;
#endif

        ERR_RETURN
    }

    Err buildKeyVsScoredCandidates(
            const QVector<CandidateScores> &scoredCandidates,
            QMap<QString, CandidateScores> *keyVsScoredCandidateCulled
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidates); ree;

        for (const CandidateScores &sc : scoredCandidates) {
            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    sc.peptideStringWithMods,
                    sc.targetKey,
                    sc.charge
                    );
            keyVsScoredCandidateCulled->insert(key, sc);
        }

        ERR_RETURN
    }

    Err buildQValCalcInput(
        const QVector<NeuralNetData> &trainingData,
        const QVector<float> &predictions,
        QHash<QString, double> *identifierVsTarget,
        QHash<QString, double> *identifierVsDecoys
        ) {

        ERR_INIT

        e = ErrorUtils::isEqual(trainingData.size(), predictions.size());


        for (int i = 0; i < trainingData.size(); i++) {

            const NeuralNetData &nnd = trainingData.at(i);

            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    nnd.peptideStringWithMods,
                    nnd.targetKey,
                    nnd.charge
            );

            if (nnd.isDecoy) {
                identifierVsDecoys->insert(key, 1 / predictions.at(i));
                continue;
            }

            identifierVsTarget->insert(key, 1 / predictions.at(i));
        }

        ERR_RETURN
    }

    Err updateScoreCandidatesClassifier(
        const QVector<CandidateScores> &scoredCandidatesAllFullFragIons,
        const QHash<QString, double> &identifierVsTarget,
        const QHash<QString, double> &identifierVsQValue,
        const QHash<QString, double> &identifierVsDecoyRatio,
        QVector<CandidateScores> *scoredCandidatesUpdateQVal
        ) {

        ERR_INIT

        for (const CandidateScores &sc : scoredCandidatesAllFullFragIons) {

            if (sc.isDecoy) {
                continue;
            }

            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    sc.peptideStringWithMods,
                    sc.targetKey,
                    sc.charge
                    );

            if (identifierVsQValue.contains(key)) {

                CandidateScores scUpdate = sc;
                scUpdate.qValue = identifierVsQValue.value(key);
                scUpdate.decoyRatio = identifierVsDecoyRatio.value(key);

                // take reciprical because it was done in buildQValCalcInput() to calc QVal because
                // score for qVal higher is better, but for classifier lower is better.
                scUpdate.classifierScore = 1 / identifierVsTarget.value(key);

                scoredCandidatesUpdateQVal->push_back(scUpdate);
            }
        }

        ERR_RETURN
    }

    Err recalculateQValsFromNeuralNetScore(
        const QVector<NeuralNetData> &trainingData,
        const QVector<float> &predictions,
        const QVector<CandidateScores> &scoredCandidatesAllFullFragIons,
        QVector<CandidateScores> *scoredCandidatesUpdateQVal
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(trainingData); ree;
        e = ErrorUtils::isNotEmpty(predictions); ree;
        e = ErrorUtils::isNotEmpty(scoredCandidatesAllFullFragIons); ree;

        scoredCandidatesUpdateQVal->clear();

        QHash<QString, double> identifierVsTarget;
        QHash<QString, double> identifierVsDecoys;
        e = buildQValCalcInput(
                trainingData,
                predictions,
                &identifierVsTarget,
                &identifierVsDecoys
                ); ree;

        qDebug() << identifierVsTarget.size() << identifierVsDecoys.size() << "SLDKJFS";

        QHash<QString, double> identifierVsQValue;
        QHash<QString, double> identifierVsDecoyRatio;
        e = MathUtils::calculateQValue(
                identifierVsTarget,
                identifierVsDecoys,
                &identifierVsQValue,
                &identifierVsDecoyRatio
        ); ree;

        e = updateScoreCandidatesClassifier(
                scoredCandidatesAllFullFragIons,
                identifierVsTarget,
                identifierVsQValue,
                identifierVsDecoyRatio,
                scoredCandidatesUpdateQVal
                ); ree;

        int targetCount;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                *scoredCandidatesUpdateQVal,
                0.01,
                &targetCount
                );
        qDebug() << targetCount << "Targets at 1% FDR";

        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                *scoredCandidatesUpdateQVal,
                0.1,
                &targetCount
                );
        qDebug() << targetCount << "Targets at 10% FDR";

        ERR_RETURN
    }

}//namespace
Err FDRCLassifierNeuralNet::exec(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        QVector<float> *meanPredictions
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isNotEmpty(xData); ree;

    e = trainClassifier(
            xData,
            yData
            ); ree;

    e = predictBaggedClassifiers(xData, meanPredictions); ree;

//    e = recalculateQValsFromNeuralNetScore(
//        trainingData,
//        meanPredictions,
//        keyVsScoredCandidateTargetsAndDecoys.values().toVector(),
//        candidateScoreClassifier
//        ); ree;

    ERR_RETURN
}


namespace {

    Err buildNeuralNetworkNormalizedVectors(
            const QVector<NeuralNetData> &trainingData,
            QVector<QVector<float>> *allDataVecs,
            QVector<float> *yData
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(trainingData); ree;

        QVector<QVector<double>> allData;
        for (const NeuralNetData &nnd : trainingData) {
            allData.push_back(nnd.scores);
        }

        Eigen::MatrixX<double> mat = EigenUtils::convertQVectorsToEigenMatrix(allData);
        EigenUtils::minMaxScaleMatrix(&mat);

        Eigen::MatrixX<float> matFloat = mat.cast<float>();
        *allDataVecs = EigenUtils::convertEigenMatrixToQVectors(matFloat);

        for (const NeuralNetData &nnd : trainingData) {
            yData->push_back(static_cast<float>(nnd.isDecoy));
        }

        e = ErrorUtils::isNotEmpty(*yData); ree;
        e = ErrorUtils::isEqual(yData->size(), allDataVecs->size()); ree;

        ERR_RETURN
    }

}//namespace
Err FDRCLassifierNeuralNet::trainClassifier(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(xData); ree;
    e = ErrorUtils::isEqual(xData.size(), yData.size());

    e = trainBaggedNeuralNets(
            xData,
            yData
            ); ree;

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::trainBaggedNeuralNets(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(xData);
    e = ErrorUtils::isEqual(xData.size(), yData.size()); ree;

    for (int bag = 0; bag < m_baggingSize; bag++) {

//        const int firstRowSize = trainingDataVecs.front().size();
//
//        qDebug() << "Tranche" << tranche;
//        qDebug() << "X Size" << trainingDataVecs.size();
//        qDebug() << "X input size" << firstRowSize;
//        qDebug() << "Y Size" << yData.size();
//        qDebug() << "X rows are uniform" << std::all_of(
//                trainingDataVecs.begin(),
//                trainingDataVecs.end(),
//                [firstRowSize](const QVector<float> &f){return f.size() == firstRowSize;}
//                );

        auto *candidateClassifier = new CandidateClassifier();
        bool trainingCompletedNoErrors = candidateClassifier->trainCandidateClassifier(
                xData,
                yData,
                m_epochs,
                m_batchSize,
                m_learningRate,
                bag
        );
        e = ErrorUtils::isTrue(trainingCompletedNoErrors); ree;

        m_candidateClassifiers.push_back(candidateClassifier);
    }

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::predictBaggedClassifiers(
        const QVector<QVector<float>> &allDataVecs,
        QVector<float> *meanPredictions
        ) {

    ERR_INIT
    e = ErrorUtils::isNotEmpty(m_candidateClassifiers); ree;
    e = ErrorUtils::isNotEmpty(allDataVecs); ree;

    Eigen::VectorX<float> predictionsSum(allDataVecs.size());ree
    predictionsSum.setZero();

    for (CandidateClassifier *cc : m_candidateClassifiers) {

        QVector<float> predictions;
        const bool predictionOK = cc->predict(
                allDataVecs,
                &predictions
                );
        e = ErrorUtils::isTrue(predictionOK); ree;

        const Eigen::VectorX<float> tranchePred = EigenUtils::convertQVectorToEigenVector(predictions);
        predictionsSum += tranchePred;
    }

    predictionsSum /= static_cast<float>(m_baggingSize);

    *meanPredictions = EigenUtils::convertEigenVectorToQVector(predictionsSum);

    ERR_RETURN
}

QString FDRCLassifierNeuralNet::buildTargetDecoyKey(
        const PeptideStringWithMods &peptideStringWithMods,
        const MzTargetKey &uniqueMsInfoScanKey,
        Charge charge
        ) {

    const QString &sep = S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP;
    return peptideStringWithMods + sep + QString::number(charge) + sep + uniqueMsInfoScanKey;
}

Err FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
        const QVector<CandidateScores> &candidateScores,
        double qValueThreshold,
        int *targetCountBelowFDRThreshold
) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;
    e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

    const auto countLogic = [qValueThreshold](const CandidateScores &cs){
        return cs.qValue < qValueThreshold;
    };

    *targetCountBelowFDRThreshold
            = static_cast<int>(std::count_if(candidateScores.begin(), candidateScores.end(), countLogic));

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::outputFDRResults(
        const QVector<CandidateScores> &candidateScores,
        bool verbose,
        QMap<QString, int> *fdrVsCount
        ) {

    ERR_INIT

    const QVector<double> fdrFractions = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005};
    for (double fdrThresh : fdrFractions) {
        int foundAtThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                candidateScores,
                fdrThresh,
                &foundAtThreshold
        ); ree;
        const double fdrPercent = fdrThresh * 100;
        fdrVsCount->insert(QString::number(fdrPercent), foundAtThreshold);

        if (!verbose) {
            continue;
        }

        qDebug() << foundAtThreshold << "PSMs found at" << fdrPercent  << "% FDR";
    }

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::filterScoreCandidatesByFDR(
        const QVector<TargetDecoyCandidatePair *> &_targetDecoyCandidatePairs,
        double qValueThreshold,
        QVector<TargetDecoyCandidatePair *> *targetDecoyCandidatePairsFDRThresholded
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(_targetDecoyCandidatePairs); ree;
    e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

    *targetDecoyCandidatePairsFDRThresholded = _targetDecoyCandidatePairs;

//    const auto terminatorLogic = [qValueThreshold](TargetDecoyCandidatePair *tdp){
//        return tdp->candidateScoresBestDiscriminantScorePtrTarget()->qValue > qValueThreshold;
//    };

//    const auto terminator = std::remove_if(
//            targetDecoyCandidatePairsFDRThresholded->begin(),
//            targetDecoyCandidatePairsFDRThresholded->end(),
//            terminatorLogic
//            );
//
//    targetDecoyCandidatePairsFDRThresholded->erase(terminator, targetDecoyCandidatePairsFDRThresholded->end());

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::filterScoreCandidatesByFDR(
        const QVector<CandidateScores> &candidateScores,
        double qValueThreshold,
        QVector<CandidateScores> *candidateScoresFDRThresholded
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(candidateScores); ree;
    e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

    *candidateScoresFDRThresholded = candidateScores;

    const auto terminatorLogic = [qValueThreshold](const CandidateScores &cs){
        return cs.qValue > qValueThreshold;
    };

    const auto terminator = std::remove_if(
            candidateScoresFDRThresholded->begin(),
            candidateScoresFDRThresholded->end(),
            terminatorLogic
    );

    candidateScoresFDRThresholded->erase(terminator, candidateScoresFDRThresholded->end());


    ERR_RETURN
}