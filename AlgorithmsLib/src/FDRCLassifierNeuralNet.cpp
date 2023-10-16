//
// Created by anichols on 9/26/23.
//

#include "FDRCLassifierNeuralNet.h"

#include "ParallelUtils.h"

#include "EigenUtils.h"

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
            const QMap<QString, ScoredCandidate> &keyVsScoredCandidateCulled,
            int topNMs2IonsFull,
            const QPair<double, double> &scanTimeMinMax,
            QVector<NeuralNetData> *trainingData
    ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(keyVsScoredCandidateCulled); ree;

        const double qValueMin = 0.01;
        const double cosineSimSumMin = 0.0;
        const bool useExtendedScores = true;
        const bool useNeuralNetworkScores = false;

        QVector<QPair<ScoredCandidate, QVector<double>>> targetPairs;
        QVector<QPair<ScoredCandidate, QVector<double>>> decoyPairs;
        QVector<QPair<ScoredCandidate, QVector<double>>> testPairs;
        for (const ScoredCandidate &sc: keyVsScoredCandidateCulled) {

            const QVector<double> scoreVector = FDRCLassifierNeuralNet::buildScoreVector(
                    sc,
                    useExtendedScores,
                    useNeuralNetworkScores,
                    topNMs2IonsFull,
                    scanTimeMinMax
            );

            if (sc.isDecoy && sc.cosineSimSum100 > cosineSimSumMin) {
                decoyPairs.push_back({sc, scoreVector});
                continue;
            }

            if (sc.qValue > qValueMin && !sc.isDecoy) {
                testPairs.push_back({sc, scoreVector});
                continue;
            }

            targetPairs.push_back({sc, scoreVector});
        }

        // TODO decide whether to shuffle this rather than sort.
        using PR = QPair<ScoredCandidate, QVector<double>>;
        std::sort(
                decoyPairs.rbegin(),
                decoyPairs.rend(),
                [](const PR &l, const PR &r){return l.first.discriminateScore < r.first.discriminateScore;}
        );

        const int decoySize = std::min(decoyPairs.size(), targetPairs.size());
        decoyPairs.resize(decoySize);

        qDebug() << "target count training" << targetPairs.size() << "decoy count training" << decoyPairs.size();

        QVector<QPair<ScoredCandidate, QVector<double>>> targetDecoyPairs;
        targetDecoyPairs.append(targetPairs);
        targetDecoyPairs.append(decoyPairs);

        std::random_device rd;
        std::mt19937 g(S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST);

        std::shuffle(targetDecoyPairs.begin(), targetDecoyPairs.end(), g);

        for (const QPair<ScoredCandidate, QVector<double>> &scPair : targetDecoyPairs) {
            NeuralNetData td;
            td.peptideStringWithMods = scPair.first.peptideStringWithMods;
            td.isDecoy = scPair.first.isDecoy;
            td.scores = scPair.second;
            td.isTest = false;
            td.charge = scPair.first.charge;
            td.targetKey = scPair.first.targetKey;
            trainingData->push_back(td);
        }

        for (const QPair<ScoredCandidate, QVector<double>> &scPair : testPairs) {
            NeuralNetData td;
            td.peptideStringWithMods = scPair.first.peptideStringWithMods;
            td.isDecoy = scPair.first.isDecoy;
            td.scores = scPair.second;
            td.isTest = true;
            td.charge = scPair.first.charge;
            td.targetKey = scPair.first.targetKey;
            trainingData->push_back(td);
        }

//#define WRITE_NN_TRAIN_DATA
#ifdef WRITE_NN_TRAIN_DATA
        const QString dataFilePath = "/home/anichols/Desktop/Testing/LatestStuff/trainingData.parquet";
        e = ParquetReader::write(*trainingData, dataFilePath); ree;
#endif

        ERR_RETURN
    }

    Err buildKeyVsScoredCandidateCulled(
            const QVector<ScoredCandidate> &scoredCandidateCulled,
            QMap<QString, ScoredCandidate> *keyVsScoredCandidateCulled
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(scoredCandidateCulled); ree;

        for (const ScoredCandidate &sc : scoredCandidateCulled) {
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
        const QVector<ScoredCandidate> &scoredCandidatesAllFullFragIons,
        const QMap<QString, double> &identifierVsTarget,
        const QMap<QString, double> &identifierVsQValue,
        const QMap<QString, double> &identifierVsDecoyRatio,
        QVector<ScoredCandidate> *scoredCandidatesUpdateQVal
        ) {

        ERR_INIT

        for (const ScoredCandidate &sc : scoredCandidatesAllFullFragIons) {

            if (sc.isDecoy) {
                continue;
            }

            const QString key = FDRCLassifierNeuralNet::buildTargetDecoyKey(
                    sc.peptideStringWithMods,
                    sc.targetKey,
                    sc.charge
                    );

            if (identifierVsQValue.contains(key)) {

                ScoredCandidate scUpdate = sc;
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
        const QVector<ScoredCandidate> &scoredCandidatesAllFullFragIons,
        QVector<ScoredCandidate> *scoredCandidatesUpdateQVal
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
        const QVector<ScoredCandidate> &trainingDataTargetsAndDecoys,
        QVector<ScoredCandidate> *scoredCandidatesClassifier
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isNotEmpty(trainingDataTargetsAndDecoys); ree;

    const int topNMs2IonsFull = std::max(
            m_minTopNMs2Ions,
            static_cast<int>(std::round(m_topNMs2Ions / 2))
    );

    QMap<QString, ScoredCandidate> keyVsScoredCandidateTargetsAndDecoys;
    e = buildKeyVsScoredCandidateCulled(
            trainingDataTargetsAndDecoys,
            &keyVsScoredCandidateTargetsAndDecoys
            ); ree;

    QVector<QVector<float>> allDataVecs;
    QVector<NeuralNetData> trainingData;
    e = buildNeuralNetworkInput(
            keyVsScoredCandidateTargetsAndDecoys,
            topNMs2IonsFull,
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
        scoredCandidatesClassifier
        ); ree;

    ERR_RETURN
}


namespace {

    Err separateTrainingDataFromAllData(
            const QVector<NeuralNetData> &trainingData,
            const QVector<QVector<float>> &allDataVecs,
            QVector<QVector<float>> *trainingDataVecs
            ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(trainingData); ree;
        e = ErrorUtils::isNotEmpty(allDataVecs); ree;
        e = ErrorUtils::isEqual(trainingData.size(), allDataVecs.size()); ree;

        for (int i = 0; i < trainingData.size(); i++) {

            const bool isTest = trainingData.at(i).isTest;

            if (isTest) {
                continue;
            }

            trainingDataVecs->push_back(allDataVecs.at(i));
        }

        ERR_RETURN
    }

    Err buildNeuralNetworkNormalizedVectors(
            const QVector<NeuralNetData> &trainingData,
            QVector<QVector<float>> *allDataVecs,
            QVector<QVector<float>> *trainingDataVecs,
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

        e = separateTrainingDataFromAllData(
                trainingData,
                *allDataVecs,
                trainingDataVecs
                ); ree;

        for (const NeuralNetData &nnd : trainingData) {

            if (nnd.isTest) {
                continue;
            }

            yData->push_back(static_cast<float>(nnd.isDecoy));
        }

        e = ErrorUtils::isNotEmpty(*yData); ree;
        e = ErrorUtils::isEqual(yData->size(), trainingDataVecs->size()); ree;

        ERR_RETURN
    }

}//namespace
Err FDRCLassifierNeuralNet::trainClassifier(
        const QMap<QString, ScoredCandidate> &keyVsScoredCandidateCulled,
        QVector<QVector<float>> *allDataVecs,
        QVector<NeuralNetData> *trainingData
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(keyVsScoredCandidateCulled); ree;

    QVector<QVector<float>> trainingDataVecs;
    QVector<float> yData;
    e = buildNeuralNetworkNormalizedVectors(
            *trainingData,
            allDataVecs,
            &trainingDataVecs,
            &yData
            ); ree;

    QVector<QVector<QVector<float>>> allTrainingVecsTranched;
    e = ParallelUtils::trancheVectorForParallelization(
            trainingDataVecs,
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

QVector<double> FDRCLassifierNeuralNet::buildScoreVector(
        const ScoredCandidate &scoreCandidate,
        bool useExtendedScores,
        bool useNeuralNetworkScores,
        int theoMzIonsSize,
        const QPair<double, double> &scanTimeMinMax

) {

    QVector<double> scores = {
            std::max(scoreCandidate.cosineSimSum100, 0.0), //1
            std::max(scoreCandidate.cosineSim100MS1, 0.0), //2
            std::pow(std::max(0.0, scoreCandidate.cosineSimSpectrum), 3), //3
            std::pow(std::max(0.0, scoreCandidate.klDivSpectrum), 1/3.0) //4
    };

    if (useExtendedScores || useNeuralNetworkScores) {

        scores.push_back(std::abs(scoreCandidate.scanTime - scoreCandidate.scanTimePredicted)); //5

        const double scanTimeRange = std::max(
                std::numeric_limits<double>::min(),
                scanTimeMinMax.second - scanTimeMinMax.first
                );
        const double scanTimeDelta = scoreCandidate.scanTime - scoreCandidate.scanTimePredicted;
        const double pdScanTime = std::sqrt(std::min(std::abs(scanTimeDelta), scanTimeRange) / scanTimeRange);
        scores.push_back(pdScanTime); //6
        scores.push_back(std::max(scoreCandidate.cosineSim100MS1Iso1, 0.0)); //8
        scores.push_back(std::max(scoreCandidate.cosineSim100MS1Iso2, 0.0)); //8
        scores.push_back(std::max(scoreCandidate.cosineSimSum45, 0.0)); //8
        scores.push_back(std::max(scoreCandidate.cosineSimSum20, 0.0)); //8
        scores.push_back(std::max(scoreCandidate.cosineSim45MS1, 0.0)); //8
        scores.push_back(std::max(scoreCandidate.cosineSim20MS1, 0.0)); //8
        scores.push_back(scoreCandidate.peptideStringWithMods.size()); //7
        scores.push_back(scoreCandidate.theoFragmentCount); //7
//        scores.push_back(scoreCandidate.peakShapeRatio1);
//        scores.push_back(scoreCandidate.peakShapeRatio2);
//        scores.push_back(scoreCandidate.peakShapeRatio3);
//        scores.push_back(std::max(scoreCandidate.b2Corr, 0.0));
//        scores.push_back(std::max(scoreCandidate.b3Corr, 0.0));

        const QVector<double> cosineSimToAnchors
                = extractScoresFromVecFeatures(scoreCandidate.cosineSimToAnchorVec, theoMzIonsSize);
        scores.append(cosineSimToAnchors); //11-16

        const QVector<double> topHalfCosineSimScores = cosineSimToAnchors.mid(0, theoMzIonsSize / 2);
        const double topHalfCosineSimScoresSum = std::accumulate(topHalfCosineSimScores.begin(), topHalfCosineSimScores.end(), 0.0);
        scores.push_back(topHalfCosineSimScoresSum); //17

        const QVector<double> bottomHalfCosineSimScores = cosineSimToAnchors.mid(theoMzIonsSize / 2, theoMzIonsSize / 2);
        const double bottomHalfCosineSimScoresSum = std::accumulate(bottomHalfCosineSimScores.begin(), bottomHalfCosineSimScores.end(), 0.0);
        scores.push_back(bottomHalfCosineSimScoresSum); //18

        const double topBottomRatio
            = std::log(std::max(1.0, topHalfCosineSimScoresSum) / (topHalfCosineSimScoresSum + bottomHalfCosineSimScoresSum + 1.0));
        scores.push_back(topBottomRatio);

        const QVector<double> theoApexIntensity
                = extractScoresFromVecFeatures(scoreCandidate.theoIntensityVec, theoMzIonsSize);
        scores.append(theoApexIntensity);

        const QVector<double> intensityFoundMaxVec
                = extractScoresFromVecFeatures(scoreCandidate.intensityFoundMaxVec, theoMzIonsSize);
        const double totalIntensity = std::accumulate(intensityFoundMaxVec.begin(), intensityFoundMaxVec.end(), 0.0);
        const double totalIntensityLog = std::log(std::max(totalIntensity, std::numeric_limits<double>::min()));
        scores.append(totalIntensityLog);

        const double maxIntensity = std::max(
                *std::max(intensityFoundMaxVec.begin(), intensityFoundMaxVec.end()),
                1.0
        );
        QVector<double> intensityFoundMaxVecNorm;
        std::transform(
                intensityFoundMaxVec.begin(),
                intensityFoundMaxVec.end(),
                std::back_inserter(intensityFoundMaxVecNorm),
                [maxIntensity](double d){return d / maxIntensity;}
        );
        scores.append(intensityFoundMaxVecNorm);

        const QVector<double> mzStDev
                = extractScoresFromVecFeatures(scoreCandidate.mzFoundStDevVec, theoMzIonsSize);
        scores.append(mzStDev);

    }

    if (useNeuralNetworkScores) {

        scores.push_back(std::max(scoreCandidate.b2Corr, 0.0));
        scores.push_back(std::max(scoreCandidate.b3Corr, 0.0));
        scores.push_back(std::max(scoreCandidate.b2b3CosineSimSum, 0.0));

        scores.push_back(std::max(0.0, scoreCandidate.cosineSimSpectrum));
        scores.push_back(scoreCandidate.discriminateScore);
        scores.push_back(std::pow(std::max(0.0, scoreCandidate.klDivSpectrum), 3));

//        QVector<double> ppmVec;
//        for (int i = 0; i < scoreCandidate.mzSearchedVec.size(); i++) {
//
//            const double mzSearched = scoreCandidate.mzSearchedVec.at(i);
//            if (i >= scoreCandidate.mzFoundMeanVec.size()) {
//                break;
//            }
//
//            const double mzFound = scoreCandidate.mzFoundMeanVec[i];
//
//            const double ppm = 1e6 * (mzFound - mzSearched) / mzSearched;
//            ppmVec.push_back(std::min(ppm, 100.0));
//        }
//        const QVector<double> ppmMz = extractScoresFromVecFeatures(ppmVec, theoMzIonsSize);
//        scores.append(ppmMz);

        scores.push_back(scoreCandidate.matrixPValue);
        scores.push_back(scoreCandidate.matrixWeight);
        scores.push_back(scoreCandidate.matrixError);
        scores.push_back(scoreCandidate.scanNumberCandidateCount);

    }

    return scores;
}

Err FDRCLassifierNeuralNet::countScoreCandidatesByFDR(
        const QVector<ScoredCandidate> &scoredCandidatesAll,
        double qValueThreshold,
        int *targetCountBelowFDRThreshold
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scoredCandidatesAll); ree;
    e = ErrorUtils::isTrue(qValueThreshold > 0.0); ree;

    const auto countLogic = [qValueThreshold](const ScoredCandidate &sc){
        return !sc.isDecoy && sc.qValue < qValueThreshold;
    };

    *targetCountBelowFDRThreshold
            = static_cast<int>(std::count_if(scoredCandidatesAll.begin(), scoredCandidatesAll.end(), countLogic));

    ERR_RETURN
}
