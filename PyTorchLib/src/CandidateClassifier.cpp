//
// Created by anichols on 9/17/23.
//

#include "CandidateClassifier.h"

#include <torch/torch.h>
#include <torch/nn/modules/loss.h>
#include <torch/optim/adam.h>
#include <omp.h>

//NOTE: QDebug and many QT classes that rely on it have to come after pytorch or a compile error occurs.
#include <QCoreApplication>
#include <QDebug>
#include <QElapsedTimer>

#include <iostream>
#include <random>

#include "ParallelUtils.h"

struct TransformerModel : torch::nn::Module {

    TransformerModel(
            int64_t input_dim,
            int64_t output_dim,
            int64_t nhead,
            int64_t dim_feedforward,
            int64_t num_encoder_layers
            )
            : transformer_encoder(torch::nn::TransformerEncoderOptions(torch::nn::TransformerEncoderLayerOptions(input_dim, nhead).dim_feedforward(dim_feedforward), num_encoder_layers))
            , fc(input_dim, output_dim) {

        register_module("transformer_encoder", transformer_encoder);
        register_module("fc", fc);
    }

    torch::Tensor forward(torch::Tensor src) {
        src = transformer_encoder->forward(src);
        return fc->forward(src);
    }

    torch::nn::TransformerEncoder transformer_encoder{nullptr};
    torch::nn::Linear fc{nullptr};
};

struct Net : torch::nn::Module {

    torch::nn::Linear layer1{nullptr};
    torch::nn::Linear layer2{nullptr};
    torch::nn::Linear layer3{nullptr};
    // torch::nn::Linear layer4{nullptr};
    // torch::nn::Linear layer5{nullptr};
    torch::nn::Linear layer6{nullptr};


    Net(int input_size, int nodes, int num_classes) {
        layer1 = register_module("layer1", torch::nn::Linear(torch::nn::LinearOptions(input_size, nodes).bias(true)));
        layer2 = register_module("layer2", torch::nn::Linear(torch::nn::LinearOptions(nodes, nodes).bias(true)));
        layer3 = register_module("layer3", torch::nn::Linear(torch::nn::LinearOptions(nodes, nodes).bias(true)));
        // layer4 = register_module("layer4", torch::nn::Linear(torch::nn::LinearOptions(nodes, nodes).bias(true)));
        // layer5 = register_module("layer5", torch::nn::Linear(torch::nn::LinearOptions(nodes, nodes).bias(true)));
        layer6 = register_module("layer6", torch::nn::Linear(nodes, num_classes));

        torch::nn::init::kaiming_uniform_(layer1->weight, torch::nn::init::calculate_gain(torch::kReLU));
        torch::nn::init::kaiming_uniform_(layer2->weight, torch::nn::init::calculate_gain(torch::kReLU));
        torch::nn::init::kaiming_uniform_(layer3->weight, torch::nn::init::calculate_gain(torch::kReLU));
        // torch::nn::init::kaiming_uniform_(layer4->weight, torch::nn::init::calculate_gain(torch::kReLU));
        // torch::nn::init::kaiming_uniform_(layer5->weight, torch::nn::init::calculate_gain(torch::kReLU));
        torch::nn::init::kaiming_uniform_(layer6->weight, torch::nn::init::calculate_gain(torch::kSigmoid));

        torch::nn::init::constant_(layer1->bias, 0.01);
        torch::nn::init::constant_(layer2->bias, 0.01);
        torch::nn::init::constant_(layer3->bias, 0.01);
        // torch::nn::init::constant_(layer4->bias, 0.01);
        // torch::nn::init::constant_(layer5->bias, 0.01);
        torch::nn::init::constant_(layer6->bias, 0.01);
    }

    torch::Tensor forward(torch::Tensor x) {
        x = torch::relu(layer1->forward(x));
        x = torch::relu(layer2->forward(x));
        x = torch::relu(layer3->forward(x));
        // x = torch::relu(layer4->forward(x));
        // x = torch::relu(layer5->forward(x));
        x = torch::sigmoid(layer6->forward(x));
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
            int batchSize,
            double learningRate,
            float focalLossGamma,
            int seed,
            double nodeFraction,
            int verbosity
    );

    bool predict(
            const QVector<QVector<float>> &xData,
            QVector<float> *predictions
    );

    bool isTrained() const;

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

        return {vec.begin(), vec.end()};
    }

    /*!
     * @brief Computes focal loss for binary classification
     * @param predictions: Model predictions (sigmoid output) in range [0,1]
     * @param targets: Ground truth labels (0 or 1)
     * @param alpha: Class balancing parameter (default 1.0)
     * @param gamma: Focusing parameter that reduces loss for well-classified examples (default 2.0)
     * @return Focal loss tensor
     *
     * Focal Loss: FL(p_t) = -α(1-p_t)^γ log(p_t)
     * where p_t = p if y=1, else 1-p
     *
     * Future enhancement: Expose alpha and gamma as configurable parameters
     */
    class FocalLoss {
    private:
        float m_alpha;
        float m_gamma;

    public:
        FocalLoss(float alpha = 1.0f, float gamma = 2.0f) : m_alpha(alpha), m_gamma(gamma) {}

        torch::Tensor operator()(const torch::Tensor& predictions, const torch::Tensor& targets) const {
            constexpr float epsilon = 1e-8f;
            torch::Tensor pred_clamped = torch::clamp(predictions, epsilon, 1.0f - epsilon);

            torch::Tensor p_t = torch::where(targets == 1, pred_clamped, 1 - pred_clamped);

            torch::Tensor focal_weight = m_alpha * torch::pow((1 - p_t), m_gamma);

            torch::Tensor ce_loss = -torch::log(p_t);

            return torch::mean(focal_weight * ce_loss);
        }
    };

} //namespace
bool CandidateClassifier::Private::trainCandidateClassifier(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int epochsMax,
        int batchSize,
        double learningRate,
        float focalLossGamma,
        int seed,
        double nodeFraction,
        int verbosity
        ) {

    torch::manual_seed(seed);
    if (torch::cuda::is_available()) {
        qDebug() << "CUDA IS AVAILABLE";
        torch::cuda::manual_seed_all(seed);
        torch::globalContext().setDeterministicCuDNN(true);
        torch::globalContext().setBenchmarkCuDNN(false);
    }

    omp_set_num_threads(1);
    torch::set_num_threads(1);
    std::srand(seed);

    torch::data::DataLoaderOptions options;
    options.enforce_ordering(true);

    const bool dataInputIsValid = checkIfInputsAreValid(
            xData,
            yData
            );

    if (!dataInputIsValid) {
        qDebug() << "Training input is not valid";
        return false;
    }

    const int input_size = xData.front().size();
    const int nodes = std::max(static_cast<int>(xData.front().size() * nodeFraction), 1);
    constexpr int num_classes = 1;

    m_net = new Net(input_size, nodes, num_classes);

    std::vector<float> flatData;
    for (const QVector<float> &innerVec : xData) {
        flatData.insert(flatData.end(), innerVec.begin(), innerVec.end());
    }
    torch::Tensor X = torch::from_blob(flatData.data(), {static_cast<long>(xData.size()), input_size}, torch::kFloat);

    std::vector<float> yStdVec(yData.begin(), yData.end());
    torch::Tensor y = torch::from_blob(yStdVec.data(), {static_cast<long>(yStdVec.size()), 1}, torch::kFloat);

    std::function<torch::Tensor(const torch::Tensor&, const torch::Tensor&)> loss_function;
    if (focalLossGamma > 0 && std::isfinite(focalLossGamma)) {
        loss_function = FocalLoss(1.0f, focalLossGamma);
    } else {
        loss_function = torch::nn::BCELoss();
    }

    torch::optim::Adam optimizer(m_net->parameters(), torch::optim::AdamOptions(learningRate));

    std::ostringstream oss;
    float bestEpochLoss = std::numeric_limits<float>::max();
    m_net->train();

    QElapsedTimer et;
    et.start();

    for (int epoch = 0; epoch < epochsMax; ++epoch) {

        float batchLossSum = 0.0;
        int iters = 0;

        for (int i = 0; i < X.size(0); i += batchSize) {

            iters++;

            torch::Tensor batchX = X.index({torch::indexing::Slice(i, i + batchSize)});
            torch::Tensor batchY = y.index({torch::indexing::Slice(i, i + batchSize)});

            const int tensorRows = static_cast<int>(batchY.sizes().at(0));
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
            if (verbosity > 0) {
                qDebug() << "Updating weights";
            }

            bestEpochLoss = meanBatchLoss;
        }

        if (verbosity > 0) {
            qDebug() << "****" << "Epoch" << epoch + 1 << "Best loss:" << bestEpochLoss << "Mean loss" << meanBatchLoss << et.restart() << "mSec";

        }
    }

    m_net->eval();

    m_isTrained = true;

    return true;
}

bool CandidateClassifier::Private::isTrained() const {
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

CandidateClassifier::CandidateClassifier()
: d_ptr(QScopedPointer<Private>(new Private))
{}

CandidateClassifier::~CandidateClassifier() {}

bool CandidateClassifier::trainCandidateClassifier(
        const QVector<QVector<float>> &xData,
        const QVector<float> &yData,
        int epochsMax,
        int batchSize,
        double learningRate,
        int seed,
        double nodeFraction,
        float focalLossGamma,
        int verbosity
        ) const {

    QVector<QVector<float>> xDataResized = xData;
    QVector<float> yDataResized = yData;

    if (yData.size() > 2e6) {
        const int resizeCount = yDataResized.size() / 2;
        xDataResized.resize(resizeCount);
        yDataResized.resize(resizeCount);
    }
    
    return d_ptr->trainCandidateClassifier(
            xDataResized,
            yDataResized,
            epochsMax,
            batchSize,
            learningRate,
            focalLossGamma,
            seed,
            nodeFraction,
            verbosity
            );
}

bool CandidateClassifier::predict(
        const QVector<QVector<float>> &xData,
        QVector<float> *predictions
        ) const {
    return d_ptr->predict(xData, predictions);
}
