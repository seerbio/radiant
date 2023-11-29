//
// Created by anichols on 11/13/23.
//

#include "KarnnNeuralNet.h"


#include "ErrorUtils.h"
#include "ParallelUtils.h"

#include "cranium.h"
#include "optimizer.h"


KarnnNeuralNet::KarnnNeuralNet()
: m_regularization(0.000001)
{}

Err KarnnNeuralNet::init(
        int baggingCount,
        int hiddenLayerCount,
        float learningRate
        ) {

    ERR_INIT

    e = ErrorUtils::isTrue(baggingCount > 0); ree;
    e = ErrorUtils::isTrue(hiddenLayerCount > 0); ree;
    e = ErrorUtils::isTrue(learningRate > 0); ree;

    m_bagginCount = baggingCount;
    m_hiddenLayerCount = hiddenLayerCount;
    m_learningRate = learningRate;

    ERR_RETURN
}

namespace {

    void tallyTrainingStatisitics(const QVector<bool> &labels) {

        int targetCount = 0;
        int decoyCount = 0;

        for (bool label : labels) {

            if (label == true) {
                decoyCount++;
            }
            else {
                targetCount++;
            }
        }
        qDebug() << "Targets:" << targetCount << "Decoys:" << decoyCount;

    }

    static float target_nn[2] = { 1.0, 0.0 };
    static float decoy_nn[2] = { 0.0, 1.0 };

    Err buildLabelVec(
            const QVector<bool> &labels,
            std::vector<float*> *labelVec
            ) {

        ERR_INIT

        for (bool isDecoy : labels) {

            if (isDecoy) {
                labelVec->emplace_back(target_nn);
            }
            else {
                labelVec->emplace_back(decoy_nn);
            }

        }

        ERR_RETURN
    }

    Err train(NN *net) {
        net->optimise();
        return eNoError;
    }

}//namespace
Err KarnnNeuralNet::run(
        int epochs,
        int seed,
        QVector<KarnnNNTarget> *karnnNNTargets
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(karnnNNTargets->isEmpty()); ree;

    QVector<QVector<double>> trainingData;
    QVector<bool> labels;
    for (const KarnnNNTarget &kt : *karnnNNTargets) {
        trainingData.push_back(kt.scoreVec);
        labels.push_back(kt.isDecoy);
    }

    e = ErrorUtils::isNotEmpty(trainingData); ree;
    e = ErrorUtils::isEqual(trainingData.size(), labels.size()); ree;

    tallyTrainingStatisitics(labels);

    std::vector<float*> labelVec;
    e = buildLabelVec(labels, &labelVec); ree;

    const auto transformLogic = [](const QVector<double> &dd){
        std::vector<float> ff;
        for (double d : dd) {
            ff.push_back(static_cast<float>(d));
        }
        return ff;
    };

    std::vector<std::vector<float>> dataVec;
    std::transform(
            trainingData.begin(),
            trainingData.end(),
            std::back_inserter(dataVec),
            transformLogic
            );

    std::vector<float*> dataVecPointers;
    std::transform(
            dataVec.begin(),
            dataVec.end(),
            std::back_inserter(dataVecPointers),
            [](std::vector<float> &ff){return ff.data();}
        );

    QVector<NN> nets;
    const int cols = trainingData.front().size();
    const int rows = trainingData.size();
    e = initNeuralNets(
            rows,
            cols,
            epochs,
            seed,
            dataVecPointers,
            labelVec,
            &nets
            ); ree;

    QVector<NN*> netsPointers;
    std::transform(
            nets.begin(),
            nets.end(),
            std::back_inserter(netsPointers),
            [](NN &n){return &n;}
            );

    QFuture<Err> futures = QtConcurrent::mapped(
            netsPointers,
            train
    );
    futures.waitForFinished();

    e = calculatedNNScores(cols, &nets, karnnNNTargets); ree;

//    for (int i = 0; i < nets.size(); i++) {
////        free(nets[i].data);
////        free(nets[i].classes);
////        destroyNetwork(nets[i].network);
//    }

    ERR_RETURN
}

Err KarnnNeuralNet::initNeuralNets(
        int rows,
        int cols,
        int epochs,
        int seed,
        std::vector<float*> &dataVecPointers,
        std::vector<float*> &labelVec,
        QVector<NN> *nets
        ) {

    ERR_INIT

    e = ErrorUtils::isFalse(dataVecPointers.empty()); ree;
    e = ErrorUtils::isFalse(labelVec.empty()); ree;

    size_t* hiddenSize = (size_t*)alloca(m_hiddenLayerCount * sizeof(size_t));
    Activation* hiddenActivation = (Activation*)alloca(m_hiddenLayerCount * sizeof(Activation));
    for (int i = 0; i < m_hiddenLayerCount; i++) {
        hiddenSize[i] = (m_hiddenLayerCount - i) * 5;
        hiddenActivation[i] = tanH;
    }

    DataSet* trainingDataSet = matrix::createDataSet(rows, cols, &(dataVecPointers[0]));
    DataSet* trainingClasses = matrix::createDataSet(rows, 2, &(labelVec[0]));

    for (int i = 0; i < m_bagginCount; i++) {

        NN net;

        net.network = createNetwork(cols, m_hiddenLayerCount, hiddenSize, hiddenActivation, 2, softmax, net.random);
        net.seed(i + seed);
        net.lossFunction = CROSS_ENTROPY_LOSS;
        net.batchSize = std::min(50, std::max(1, static_cast<int>(rows / 100.0)));
        net.learningRate = m_learningRate;
        net.searchTime = 0;
        net.regularizationStrength = m_regularization;
        net.B1 = 0.9;
        net.B2 = 0.999;
        net.shuffle = 1;
        net.verbose = true;
        net.data = trainingDataSet;
        net.classes = trainingClasses;
        net.maxIters = static_cast<int>((rows * epochs) / net.batchSize);

        nets->push_back(net);
    }

    ERR_RETURN
}

Err KarnnNeuralNet::calculatedNNScores(
        int cols,
        QVector<NN> *net,
        QVector<KarnnNNTarget> *karnnNNTargets
) {

    ERR_INIT

    e = ErrorUtils::isFalse(net->empty()); ree;
    e = ErrorUtils::isTrue(cols > 0);
    e = ErrorUtils::isFalse(karnnNNTargets->isEmpty()); ree;

    std::vector<std::vector<float>> dataVecs;
    for (const KarnnNNTarget &t : *karnnNNTargets) {
        std::vector<float> f;
        for (double d : t.scoreVec) {
            f.push_back(static_cast<float>(d));
        }
        dataVecs.push_back(f);
    }

    std::vector<float*> dataVecsPointers;
    std::transform(
            dataVecs.begin(),
            dataVecs.end(),
            std::back_inserter(dataVecsPointers),
            [](std::vector<float> &ff){return ff.data();}
            );

    DataSet* dataset = matrix::createDataSet(dataVecsPointers.size(), cols, &(dataVecsPointers[0]));

    for (int i = 0; i < m_bagginCount; i++) {
        forwardPassDataSet((*net)[i].network, dataset);
        Matrix* result = getOuput((*net)[i].network);
        for (int j = 0; j < result->rows; j++) {
            (*karnnNNTargets)[j].nnScore += matrix::getMatrix(result, j, 0) / m_bagginCount;
        }
    }

    ERR_RETURN
}
