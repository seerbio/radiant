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
, m_isInit(false)
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

Err FDRCLassifierNeuralNet::exec(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int seed,
        QVector<float> *meanPredictions
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isNotEmpty(xData); ree;

    e = trainClassifier(
            xData,
            yData,
            seed
            ); ree;

    e = predictBaggedClassifiers(xData, meanPredictions); ree;

    ERR_RETURN
}

Err FDRCLassifierNeuralNet::trainClassifier(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int seed
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(xData); ree;
    e = ErrorUtils::isEqual(xData.size(), yData.size());

    e = trainBaggedNeuralNets(
            xData,
            yData,
            seed
            ); ree;

    ERR_RETURN
}

namespace {

    struct CandidateClassifierParallelInput {
        CandidateClassifier *candidateClassifier;
        QVector<QVector<float>> xData;
        QVector<float> yData;
        int epochs = -1;
        int batchSize = -1;
        double learningRate = -1.0;
        int bag = -1;
    };

    Err trainingLogic(const CandidateClassifierParallelInput &input) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(input.xData); ree;
        e = ErrorUtils::isEqual(input.xData.size(), input.yData.size()); ree;

        bool trainingCompletedNoErrors = input.candidateClassifier->trainCandidateClassifier(
                input.xData,
                input.yData,
                input.epochs,
                input.batchSize,
                input.learningRate,
                input.bag
        );
        e = ErrorUtils::isTrue(trainingCompletedNoErrors); ree;

        ERR_RETURN
    }

}//namespace
Err FDRCLassifierNeuralNet::trainBaggedNeuralNets(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int seed
        ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(xData);
    e = ErrorUtils::isEqual(xData.size(), yData.size()); ree;

    QElapsedTimer et;
    et.start();

    QVector<CandidateClassifierParallelInput> parallelInputs;
    for (int bag = 0; bag < m_baggingSize; bag++) {

        auto *candidateClassifier = new CandidateClassifier();
        m_candidateClassifiers.push_back(candidateClassifier);

        CandidateClassifierParallelInput ccpi;
        ccpi.candidateClassifier = candidateClassifier;
        ccpi.xData = xData;
        ccpi.yData = yData;
        ccpi.epochs = m_epochs;
        ccpi.batchSize = m_batchSize;
        ccpi.learningRate = m_learningRate;
        ccpi.bag = bag + S_GLOBAL_SETTINGS.NUMBER_OF_THE_BEAST + seed;

        parallelInputs.push_back(ccpi);
    }

//#define PARALLEL_TRAIN_NN
#ifdef PARALLEL_TRAIN_NN
    QFuture<Err> futures = QtConcurrent::mapped(
            parallelInputs,
            trainingLogic
            );
    futures.waitForFinished();

    for (const Err &result : futures) {
        e = result; ree;
    }
#else
    for (const CandidateClassifierParallelInput &ccpi : parallelInputs) {
        e = trainingLogic(ccpi); ree;
    }
#endif

    qDebug() << "Neural nets trained in:" << et.elapsed() << "mSec";
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
