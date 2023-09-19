//
// Created by anichols on 9/17/23.
//

#include "CandidateClassifier.h"

#include <torch/torch.h>
#include <torch/nn/modules/loss.h>
#include <torch/optim/adam.h>

#include <QCoreApplication>
#include <QDir>

#include <iostream>
#include <random>


struct Net : torch::nn::Module {

    torch::nn::Linear layer1{nullptr}, layer2{nullptr}, layer3{nullptr}, layer4{nullptr};
    torch::nn::BatchNorm1d batchNorm1{nullptr}, batchNorm2{nullptr}, batchNorm3{nullptr};

    Net(int input_size, int nodes, int num_classes) {
        layer1 = register_module("layer1", torch::nn::Linear(torch::nn::LinearOptions(input_size, nodes).bias(false)));
        batchNorm1 = register_module("batchNorm1", torch::nn::BatchNorm1d(nodes));
        layer2 = register_module("layer2", torch::nn::Linear(torch::nn::LinearOptions(nodes, nodes).bias(false)));
        batchNorm2 = register_module("batchNorm2", torch::nn::BatchNorm1d(nodes));
        layer3 = register_module("layer3", torch::nn::Linear(torch::nn::LinearOptions(nodes, nodes).bias(false)));
        batchNorm3 = register_module("batchNorm3", torch::nn::BatchNorm1d(nodes));
        layer4 = register_module("layer4", torch::nn::Linear(nodes, num_classes));

        std::mt19937 rng(666);
        torch::nn::init::xavier_uniform_(layer1->weight);
        torch::nn::init::xavier_uniform_(layer2->weight);
        torch::nn::init::xavier_uniform_(layer3->weight);
        torch::nn::init::xavier_uniform_(layer4->weight);
    }

    torch::Tensor forward(torch::Tensor x) {
        x = torch::relu(batchNorm1->forward(layer1->forward(x)));
        x = torch::relu(batchNorm2->forward(layer2->forward(x)));
        x = torch::relu(batchNorm3->forward(layer3->forward(x)));
        x = torch::sigmoid(layer4->forward(x));
        return x;
    }
};

class Q_DECL_HIDDEN CandidateClassifier::Private {

public:

    Private();
    ~Private();

    bool trainCandidateClassifier(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData,
            int epochsMax,
            double batchFraction,
            double learningRate
    );

private:

    Net *m_net;
    bool m_isTrained;

};

CandidateClassifier::Private::Private()
: m_net(nullptr)
, m_isTrained(false)
{}

CandidateClassifier::Private::~Private() { delete m_net;}

namespace {

    bool checkIfInputsAreValid(
            const QVector<QVector<float>> &xData,
            const QVector<float> &yData
            ) {

        const bool xyDataIsSameSize = xData.size() == yData.size();

        const int firstRowSize = xData.front().size();
        const bool xDataAllSameSize = std::all_of(
                xData.begin(),
                xData.end(),
                [firstRowSize](const QVector<float> &dd){return dd.size() == firstRowSize;}
                );

        return xDataAllSameSize && xyDataIsSameSize;
    }

} //namespace
bool CandidateClassifier::Private::trainCandidateClassifier(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int epochsMax,
        double batchFraction,
        double learningRate
        ) {

    const bool dataInputIsValid = checkIfInputsAreValid(
            xData,
            yData
            );

    if (!dataInputIsValid) {
        return false;
    }

    const int input_size = xData.front().size();
    const int nodes = input_size;
    const int num_classes = 1;
    const int batchSize = static_cast<int>(batchFraction * xData.size());

    m_net = new Net(input_size, nodes, num_classes);

    std::vector<float> flatData;
    for (const QVector<float> &innerVec : xData) {
        flatData.insert(flatData.end(), innerVec.begin(), innerVec.end());
    }
    torch::Tensor X = torch::from_blob(flatData.data(), {static_cast<long>(xData.size()), input_size}, torch::kFloat);

    std::vector<float> yStdVec = yData.toStdVector();
    torch::Tensor y = torch::from_blob(yStdVec.data(), {static_cast<long>(yStdVec.size()), 1}, torch::kFloat);

    torch::nn::BCELoss loss_function;
    torch::optim::Adam optimizer(m_net->parameters(), torch::optim::AdamOptions(learningRate));

//    const QString &weightsPath
//            = QDir(qApp->applicationDirPath()).filePath("model_weights.pt");

    std::ostringstream oss;
    float bestEpochLoss = std::numeric_limits<float>::max();
    m_net->train();

    for (int epoch = 0; epoch < epochsMax; ++epoch) {

        float batchLossSum = 0.0;
        int iters = 0;

        for (int i = 0; i < X.size(0); i += batchSize) {

            iters++;

            torch::Tensor batchX = X.index({torch::indexing::Slice(i, i + batchSize)});
            torch::Tensor batchY = y.index({torch::indexing::Slice(i, i + batchSize)});

            torch::Tensor output = m_net->forward(batchX);
            torch::Tensor loss = loss_function(output, batchY);

            const auto batchLoss = loss.item<float>();
            batchLossSum += batchLoss;

            optimizer.zero_grad();
            loss.backward();
            optimizer.step();

            std::cout << "Epoch " << epoch + 1 << ", Batch " << (i / batchSize) + 1 << ", Loss: " << loss.item<float>() << std::endl;
        }

        const float meanBatchLoss = batchLossSum / static_cast<float>(iters);

//        if (meanBatchLoss < bestEpochLoss) {
//
//            std::cout << "Updating weights" << std::endl;
//
//            bestEpochLoss = meanBatchLoss;
//            torch::serialize::OutputArchive archive;
//            m_net->save(archive);
//            archive.save_to(oss);
//        }

        std::cout << "****" << "Epoch " << epoch + 1 << ", Best loss: " << bestEpochLoss << " Mean loss " << meanBatchLoss <<  std::endl;


    }

    torch::Tensor output = m_net->forward(X);
    std::cout << output << std::endl;

    torch::serialize::OutputArchive archive;
    m_net->save(archive);
    archive.save_to(oss);

    m_net->eval();

    std::istringstream iss(oss.str());
    torch::serialize::InputArchive inputArchive;
    inputArchive.load_from(iss);
    m_net->load(inputArchive);

    torch::Tensor output2 = m_net->forward(X);
    std::cout << output2 << std::endl;

    m_isTrained = true;

    return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
//END PRIVATE
///////////////////////////////////////////////////////////////////////////////////////////

CandidateClassifier::CandidateClassifier() : d_ptr(new Private()) {}

CandidateClassifier::~CandidateClassifier() {}

bool CandidateClassifier::trainCandidateClassifier(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int epochsMax,
        double batchFraction,
        double learningRate
        ) {
    return d_ptr->trainCandidateClassifier(
            xData,
            yData,
            epochsMax,
            batchFraction,
            learningRate
            );
}
