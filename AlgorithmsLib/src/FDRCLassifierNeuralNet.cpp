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
, m_batchFraction(-1.0)
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
        int topNMs2Ions,
        int epochs,
        int baggingSize,
        double batchFraction,
        double learningRate,
        const QPair<double, double> &scanTimeMinMax
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(topNMs2Ions >= m_minTopNMs2Ions); ree;
    e = ErrorUtils::isTrue(epochs > 0); ree;
    e = ErrorUtils::isTrue(baggingSize >= 1); ree;
    e = ErrorUtils::isTrue(batchFraction > 0 & batchFraction < 1 & !MathUtils::tZero(batchFraction)); ree;
    e = ErrorUtils::isTrue(learningRate > 0 & learningRate < 1 & !MathUtils::tZero(learningRate)); ree;

    m_topNMs2Ions = topNMs2Ions;
    m_epochs = epochs;
    m_baggingSize = baggingSize;
    m_batchFraction = batchFraction;
    m_learningRate = learningRate;
    m_scanTimeMinMax = scanTimeMinMax;

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
            e = FDRCLassifierNeuralNet::buildScoreVector(
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

#define WRITE_NN_TRAIN_DATA
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
        QMap<QString, double> *identifierVsTarget,
        QMap<QString, double> *identifierVsDecoys
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
        const QMap<QString, double> &identifierVsTarget,
        const QMap<QString, double> &identifierVsQValue,
        const QMap<QString, double> &identifierVsDecoyRatio,
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

        QMap<QString, double> identifierVsTarget;
        QMap<QString, double> identifierVsDecoys;
        e = buildQValCalcInput(
                trainingData,
                predictions,
                &identifierVsTarget,
                &identifierVsDecoys
                ); ree;

        qDebug() << identifierVsTarget.size() << identifierVsDecoys.size() << "SLDKJFS";

        QMap<QString, double> identifierVsQValue;
        QMap<QString, double> identifierVsDecoyRatio;
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
        const QVector<CandidateScores> &candidateScoresTargetsAndDecoys,
        QVector<CandidateScores> *candidateScoreClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isNotEmpty(candidateScoresTargetsAndDecoys); ree;
    candidateScoreClassifier->clear();

    QMap<QString, CandidateScores> keyVsScoredCandidateTargetsAndDecoys;
    e = buildKeyVsScoredCandidates(
            candidateScoresTargetsAndDecoys,
            &keyVsScoredCandidateTargetsAndDecoys
            ); ree;

    QVector<QVector<float>> allDataVecs;
    QVector<NeuralNetData> trainingData;
    e = buildNeuralNetworkInput(
            keyVsScoredCandidateTargetsAndDecoys,
            m_topNMs2Ions,
            m_scanTimeMinMax,
            &trainingData
    ); ree;

    e = trainClassifier(
            keyVsScoredCandidateTargetsAndDecoys,
            &allDataVecs,
            &trainingData
            ); ree;

    QVector<float> meanPredictions;
    e = predictBaggedClassifiers(allDataVecs, &meanPredictions); ree;

    e = recalculateQValsFromNeuralNetScore(
        trainingData,
        meanPredictions,
        keyVsScoredCandidateTargetsAndDecoys.values().toVector(),
        candidateScoreClassifier
        ); ree;

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
        const QMap<QString, CandidateScores> &keyVsScoredCandidateCulled,
        QVector<QVector<float>> *allDataVecs,
        QVector<NeuralNetData> *trainingData
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(keyVsScoredCandidateCulled); ree;

    QVector<float> yData;
    e = buildNeuralNetworkNormalizedVectors(
            *trainingData,
            allDataVecs,
            &yData
            ); ree;

    QVector<QVector<QVector<float>>> allTrainingVecsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            *allDataVecs,
            m_baggingSize,
            &allTrainingVecsTranched
            ); ree;

    QVector<QVector<float>> yDataTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            yData,
            m_baggingSize,
            &yDataTranched
            ); ree;

    e = trainBaggedNeuralNets(
            allTrainingVecsTranched,
            yDataTranched
            ); ree;

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::trainBaggedNeuralNets(
        const QVector<QVector<QVector<float>>> &trainingDataVecsTranched,
        const QVector<QVector<float>> &yTrainingDataTranched
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(trainingDataVecsTranched);
    e = ErrorUtils::isEqual(trainingDataVecsTranched.size(), yTrainingDataTranched.size()); ree;

    for (int tranche = 0; tranche < trainingDataVecsTranched.size(); tranche++) {

        const QVector<QVector<float>> &trainingDataVecs = trainingDataVecsTranched.at(tranche);
        const QVector<float> &yData = yTrainingDataTranched.at(tranche);
        const int firstRowSize = trainingDataVecs.front().size();

        qDebug() << "Tranche" << tranche;
        qDebug() << "X Size" << trainingDataVecs.size();
        qDebug() << "X input size" << firstRowSize;
        qDebug() << "Y Size" << yData.size();
        qDebug() << "X rows are uniform" << std::all_of(
                trainingDataVecs.begin(),
                trainingDataVecs.end(),
                [firstRowSize](const QVector<float> &f){return f.size() == firstRowSize;}
                );

        auto *candidateClassifier = new CandidateClassifier();
        bool trainingCompletedNoErrors = candidateClassifier->trainCandidateClassifier(
                trainingDataVecs,
                yData,
                m_epochs,
                m_batchFraction,
                m_learningRate
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
        const UniqueMsInfoScanKey &uniqueMsInfoScanKey,
        Charge charge
        ) {

    const QString &sep = S_GLOBAL_SETTINGS.MODIFICATION_INTERNAL_SEP;
    return peptideStringWithMods + sep + QString::number(charge) + sep + uniqueMsInfoScanKey;
}

Err FDRCLassifierNeuralNet::buildScoreVector(
        const CandidateScores &candidateScores,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int theoMzIonsSize,
        const QPair<double, double> &scanTimeMinMax,
        QVector<double> *scoreVec
        ) {

    ERR_INIT

    QVector<double> scores = {
            static_cast<double>(std::max(0, candidateScores.allignedMaxIndexesCount)),
            std::max(candidateScores.cosineSimSum100, 0.0), //1
            std::max(candidateScores.cosineSim100MS1, 0.0), //2
            std::pow(std::max(0.0, candidateScores.cosineSimSpectrum), 3), //3
            std::pow(std::max(0.0, candidateScores.klDivSpectrum), 1/3.0) //4
    };

    if (useExtendedScores || useNeuralNetworkScores) {

        const double scanTimeDelta = candidateScores.scanTime - candidateScores.scanTimePredicted;
        scores.push_back(std::abs(scanTimeDelta)); //5

        const double scanTimeRange = std::max(
                std::numeric_limits<double>::min(),
                scanTimeMinMax.second - scanTimeMinMax.first
                );

        const double pdScanTime = std::sqrt(std::min(std::abs(scanTimeDelta), scanTimeRange) / scanTimeRange);
        scores.push_back(pdScanTime); //6

        scores.push_back(std::max(candidateScores.cosineSim100MS1Iso1, 0.0)); //7
        scores.push_back(std::max(candidateScores.cosineSim100MS1Iso2, 0.0)); //8
        scores.push_back(std::max(candidateScores.cosineSimSum45, 0.0)); //9
        scores.push_back(std::max(candidateScores.cosineSimSum20, 0.0)); //10
        scores.push_back(std::max(candidateScores.cosineSim45MS1, 0.0)); //11
        scores.push_back(std::max(candidateScores.cosineSim20MS1, 0.0)); //12

        const double pepLength = (-10.0 + candidateScores.peptideStringWithMods.size()) / 10.0;
        scores.push_back(pepLength); //13
        scores.push_back(candidateScores.theoFragmentCount); //14

        const QVector<double> cosineSimToAnchors
                = extractScoresFromVecFeatures(candidateScores.cosineSimToAnchorVec, theoMzIonsSize);
        scores.append(cosineSimToAnchors); //28-33

        const QVector<double> topHalfCosineSimScores = cosineSimToAnchors.mid(0, theoMzIonsSize / 2);
        const double topHalfCosineSimScoresSum = std::accumulate(topHalfCosineSimScores.begin(), topHalfCosineSimScores.end(), 0.0);
        scores.push_back(topHalfCosineSimScoresSum); //34

        const QVector<double> bottomHalfCosineSimScores = cosineSimToAnchors.mid(theoMzIonsSize / 2, theoMzIonsSize / 2);
        const double bottomHalfCosineSimScoresSum = std::accumulate(bottomHalfCosineSimScores.begin(), bottomHalfCosineSimScores.end(), 0.0);
        scores.push_back(bottomHalfCosineSimScoresSum); //35

        const double topBottomRatio
            = std::log(std::max(1.0, topHalfCosineSimScoresSum) / (topHalfCosineSimScoresSum + bottomHalfCosineSimScoresSum + 1.0));
        scores.push_back(topBottomRatio); //36

        const QVector<double> theoApexIntensity
                = extractScoresFromVecFeatures(candidateScores.theoIntensityVec, theoMzIonsSize);
        scores.append(theoApexIntensity);

        const QVector<double> intensityFoundMaxVec
                = extractScoresFromVecFeatures(candidateScores.intensityFoundMaxVec, theoMzIonsSize);
        const double totalIntensity = std::accumulate(intensityFoundMaxVec.begin(), intensityFoundMaxVec.end(), 0.0001);
        const double totalIntensityLog = std::log(std::max(totalIntensity, std::numeric_limits<double>::min()));
        scores.push_back(totalIntensityLog);

        if (useNeuralNetworkScores) {

            for (double intzFound : intensityFoundMaxVec) {
                scores.push_back(intzFound / totalIntensity);
            }

            const int topSix = 6;
            QVector<double> ppmVec;
            for (int i = 0; i < topSix; i++) {

                const double mzSearched = candidateScores.mzSearchedVec.at(i);
                if (i >= candidateScores.mzFoundMeanVec.size()) {
                    break;
                }

                const double mzFound = candidateScores.mzFoundMeanVec[i];
                const double ppm = cosineSimToAnchors.at(i) * std::abs(mzFound - mzSearched) / mzSearched;
                ppmVec.push_back(ppm);
            }
        }

        const double charge = -2.0 + candidateScores.charge;
        scores.push_back(charge);
    }

    if (useNeuralNetworkScores) {

        scores.push_back(std::max(0.0, candidateScores.cosineSimSpectrum));
        scores.push_back(candidateScores.discriminateScore);
        scores.push_back(std::pow(std::max(0.0, candidateScores.klDivSpectrum), 3));

        QMap<QChar, int> aminoAcidCounts = {
                {'A', 0},
                {'C', 0},
                {'D', 0},
                {'E', 0},
                {'F', 0},
                {'G', 0},
                {'H', 0},
                {'I', 0},
                {'K', 0},
                {'L', 0},
                {'M', 0},
                {'N', 0},
                {'P', 0},
                {'Q', 0},
                {'R', 0},
                {'S', 0},
                {'T', 0},
                {'V', 0},
                {'W', 0},
                {'Y', 0}
        };

        for (const QChar aminoAcid : candidateScores.peptideStringWithMods) {

            if (!aminoAcidCounts.contains(aminoAcid)) {
                qDebug() << candidateScores.peptideStringWithMods << "missing amino acid";
            }

            e = ErrorUtils::isTrue(aminoAcidCounts.contains(aminoAcid)); ree;
            aminoAcidCounts[aminoAcid]++;
        }

        for (int cnt : aminoAcidCounts.values()) {
            scores.push_back(static_cast<double>(cnt));
        }

        const double mz = BiophysicalCalcs::calculateThomsonFromMass(
                candidateScores.mass,
                candidateScores.charge
                );
        scores.push_back((mz - 600.0) * 0.002);

        scores.push_back(candidateScores.peakShapeRatio1);
        scores.push_back(candidateScores.peakShapeRatio2);
        scores.push_back(candidateScores.peakShapeRatio3);
        scores.push_back(candidateScores.cosineSimSumBottom6);
        scores.push_back(candidateScores.cosineSimSumBottom6 / candidateScores.theoFragmentCount);

    }

    *scoreVec = scores;

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
        const QVector<TargetDecoyCandidatePair *> &targetDecoyCandidatePair,
        double qValueThreshold,
        int *targetCountBelowFDRThreshold
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(targetDecoyCandidatePair); ree;
    e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

    const auto countLogic = [qValueThreshold](TargetDecoyCandidatePair *tdcp){
        return tdcp->candidateScoresBestDiscriminantScorePtrTarget()->qValue < qValueThreshold;
    };

    *targetCountBelowFDRThreshold
            = static_cast<int>(std::count_if(targetDecoyCandidatePair.begin(), targetDecoyCandidatePair.end(), countLogic));

    ERR_RETURN
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
        const QVector<TargetDecoyCandidatePair*> &targetDecoyCandidatePairs,
        bool verbose,
        QMap<QString, int> *fdrVsCount
        ) {

    ERR_INIT

    const QVector<double> fdrFractions = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005};
    for (double fdrThresh : fdrFractions) {
        int foundAtThreshold;
        e = FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
                targetDecoyCandidatePairs,
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

    const auto terminatorLogic = [qValueThreshold](TargetDecoyCandidatePair *tdp){
        return tdp->candidateScoresBestDiscriminantScorePtrTarget()->qValue > qValueThreshold;
    };

    const auto terminator = std::remove_if(
            targetDecoyCandidatePairsFDRThresholded->begin(),
            targetDecoyCandidatePairsFDRThresholded->end(),
            terminatorLogic
            );

    targetDecoyCandidatePairsFDRThresholded->erase(terminator, targetDecoyCandidatePairsFDRThresholded->end());

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