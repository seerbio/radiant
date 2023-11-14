//
// Created by anichols on 9/17/23.
//

#include "CandidateClassifier.h"

#include <torch/torch.h>
#include <torch/nn/modules/loss.h>
#include <torch/optim/adam.h>

//NOTE: QDebug and many QT classes that rely on it have to come after pytorch or a compile error occurs.
#include <QCoreApplication>
#include <QDebug>

#include <iostream>
#include <random>

const int NUMBER_OF_THE_BEAST = 666;

struct Net : torch::nn::Module {

    torch::nn::Linear layer1{nullptr};
    torch::nn::Linear layer2{nullptr};
//    torch::nn::Linear layer3{nullptr};
    torch::nn::Linear layer4{nullptr};

    torch::nn::BatchNorm1d batchNorm1{nullptr};
    torch::nn::BatchNorm1d batchNorm2{nullptr};
//    torch::nn::BatchNorm1d batchNorm3{nullptr};

    Net(int input_size, int nodes, int num_classes) {
        layer1 = register_module("layer1", torch::nn::Linear(torch::nn::LinearOptions(input_size, nodes).bias(false)));
        batchNorm1 = register_module("batchNorm1", torch::nn::BatchNorm1d(nodes));
        layer2 = register_module("layer2", torch::nn::Linear(torch::nn::LinearOptions(nodes, nodes).bias(false)));
        batchNorm2 = register_module("batchNorm2", torch::nn::BatchNorm1d(nodes));
//        layer3 = register_module("layer3", torch::nn::Linear(torch::nn::LinearOptions(nodes, nodes).bias(false)));
//        batchNorm3 = register_module("batchNorm3", torch::nn::BatchNorm1d(nodes));
        layer4 = register_module("layer4", torch::nn::Linear(nodes, num_classes));

        torch::nn::init::xavier_uniform_(layer1->weight);
        torch::nn::init::xavier_uniform_(layer2->weight);
//        torch::nn::init::xavier_uniform_(layer3->weight);
        torch::nn::init::xavier_uniform_(layer4->weight);
        torch::nn::init::constant_(layer4->bias, 0.01);
    }

    torch::Tensor forward(torch::Tensor x) {
        x = torch::relu(batchNorm1->forward(layer1->forward(x)));
        x = torch::relu(batchNorm2->forward(layer2->forward(x)));
//        x = torch::relu(batchNorm3->forward(layer3->forward(x)));
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

    bool predict(
            const QVector<QVector<float>> &xData,
            QVector<float> *predictions
    );

    bool isTrained();

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

    QVector<float> tensorToVector(const torch::Tensor& tensor) {

        const long numel = tensor.numel();
        const float* data = tensor.data_ptr<float>();
        std::vector<float> vec(data, data + numel);

        return QVector<float>(vec.begin(), vec.end());
    }

} //namespace
bool CandidateClassifier::Private::trainCandidateClassifier(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int epochsMax,
        double batchFraction,
        double learningRate
        ) {

    torch::manual_seed(NUMBER_OF_THE_BEAST);

    const bool dataInputIsValid = checkIfInputsAreValid(
            xData,
            yData
            );

    if (!dataInputIsValid) {
        qDebug() << "Training input is not valid";
        return false;
    }

    const int input_size = xData.front().size();
    const int nodes = 20;
    const int num_classes = 1;
    const int batchSize = static_cast<int>(batchFraction * xData.size());

    m_net = new Net(input_size, nodes, num_classes);

    std::vector<float> flatData;
    for (const QVector<float> &innerVec : xData) {
        flatData.insert(flatData.end(), innerVec.begin(), innerVec.end());
    }
    torch::Tensor X = torch::from_blob(flatData.data(), {static_cast<long>(xData.size()), input_size}, torch::kFloat);

    std::vector<float> yStdVec(yData.begin(), yData.end());
    torch::Tensor y = torch::from_blob(yStdVec.data(), {static_cast<long>(yStdVec.size()), 1}, torch::kFloat);

    torch::nn::BCELoss loss_function;
    torch::optim::Adam optimizer(m_net->parameters(), torch::optim::AdamOptions(learningRate));

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

            const int tensorRows = batchY.sizes().at(0);
            if (tensorRows < 2) {
                qDebug() << "Skipping batch due to low training volume" << tensorRows;
                continue;
            }

            torch::Tensor output = m_net->forward(batchX);

//#define REGULARIZATION
#ifdef REGULARIZATION
            // Regularization strength (lambda)
            float lambda = 0.01; // Adjust as needed
            torch::Tensor l2_regularization = torch::tensor(0.0, torch::kFloat32);
            for (const auto& parameter : m_net->parameters()) {
                l2_regularization += torch::norm(parameter, 2);
            }
            torch::Tensor loss = loss_function(output, batchY) + lambda * l2_regularization;
#else
            torch::Tensor loss = loss_function(output, batchY);
#endif
            const auto batchLoss = loss.item<float>();
            batchLossSum += batchLoss;

            optimizer.zero_grad();
            loss.backward();
            optimizer.step();
        }

        const float meanBatchLoss = batchLossSum / static_cast<float>(iters);

        if (meanBatchLoss < bestEpochLoss) {

            qDebug() << "Updating weights";

            bestEpochLoss = meanBatchLoss;
        }

        qDebug() << "****" << "Epoch" << epoch + 1 << "Best loss:" << bestEpochLoss << "Mean loss" << meanBatchLoss;


    }

    m_net->eval();

    torch::Tensor output = m_net->forward(X);
    const QVector<float> outputVec = tensorToVector(output);

    int falsePos = 0;
    int falseNeg = 0;
    for (int i = 0; i < yData.size(); i++) {

        const float act = yData.at(i);
        const float pred =  outputVec.at(i);

        if (act > 0.5 && pred < 0.5) {
            falsePos++;
        }
        else if (act < 0.5 && pred >= 0.5) {
            falseNeg++;
        }

    }

    m_isTrained = true;

    return true;
}

bool CandidateClassifier::Private::isTrained() {
    return m_isTrained;
}

bool CandidateClassifier::Private::predict(
        const QVector<QVector<float>> &xData,
        QVector<float> *predictions
        ) {

    if (!isTrained()) {
        qDebug() << "Neural network is not trained";
        return false;
    }

    const int input_size = xData.front().size();
    m_net->eval();

    std::vector<float> flatData;
    for (const QVector<float> &innerVec : xData) {
        flatData.insert(flatData.end(), innerVec.begin(), innerVec.end());
    }
    torch::Tensor X = torch::from_blob(flatData.data(), {static_cast<long>(xData.size()), input_size}, torch::kFloat);

    torch::Tensor output = m_net->forward(X);
    *predictions = tensorToVector(output);

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

bool CandidateClassifier::predict(
        const QVector<QVector<float>> &xData,
        QVector<float> *predictions
        ) {
    return d_ptr->predict(xData, predictions);
}
