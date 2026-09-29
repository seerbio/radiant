#include "PeptideFamilyNeuralNet.h"

#include <QCoreApplication>
#include <QThreadPool>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>

namespace {
void require(bool condition, const char *message) {
    if (!condition) {
        std::cerr << message << '\n';
        std::exit(20);
    }
}
}

int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    require(PeptideFamilyNeuralNet::familyHash("_PEPI(UniMod:35)DE_")
            == PeptideFamilyNeuralNet::familyHash("PEPLDE"), "Family equivalence");
    require(PeptideFamilyNeuralNet::familyHash("PEPLDE")
            != PeptideFamilyNeuralNet::familyHash("PEPTDE"), "Distinct family hashes");
    QVector<QVector<float>> values;
    QVector<quint32> decoys, families;
    for (int row = 0; row < 600; ++row) {
        QVector<float> full(FeaturesSize, 1.0f);
        for (int column = 0; column < FeaturesSize; ++column)
            full[column] = 1.0f + .15f * std::sin(float(row) * .031f + float(column) * .17f);
        full[CosineSimSpectrum] = .9f - .65f * (row % 2) + .04f * std::sin(float(row));
        full[CosineSim100MS1] = .9f - .6f * (row % 2);
        values.push_back(full);
        decoys.push_back(row % 2);
        // Six target/decoy/view rows for each origin family.
        families.push_back(row / 6);
    }
    const auto original = values;
    const auto originalLabels = decoys;
    const auto features = PeptideFamilyNeuralNet::nonTimsFeatures();
    PeptideFamilyNeuralNet::Settings settings;
    settings.epochs = 2;
    settings.networks = 4;
    settings.threads = 1;
    PeptideFamilyNeuralNet::Predictions serial, parallel, changed;
    QThreadPool::globalInstance()->setMaxThreadCount(7);
    require(PeptideFamilyNeuralNet::score(values, decoys, families, features, settings, &serial) == 0,
            "Serial score");
    settings.threads = 12;
    require(PeptideFamilyNeuralNet::score(values, decoys, families, features, settings, &parallel) == 0,
            "Parallel score");
    require(serial.perNetwork == parallel.perNetwork, "Serial/parallel determinism");
    require(serial.meanDecoyProbability == parallel.meanDecoyProbability, "Serial/parallel means");
    require(serial.heldOutFold == parallel.heldOutFold, "Serial/parallel folds");
    for (int threads : {2, 3, 7}) {
        settings.threads = threads;
        PeptideFamilyNeuralNet::Predictions limited;
        require(PeptideFamilyNeuralNet::score(values, decoys, families, features, settings, &limited) == 0,
                "Limited worker score");
        require(serial.perNetwork == limited.perNetwork, "Limited worker determinism");
    }
    settings.threads = 12;
    require(QThreadPool::globalInstance()->maxThreadCount() == 7, "Global pool unchanged");
    for (int row = 0; row < values.size(); ++row)
        if (families[row] % 3 == 0) decoys[row] = 1 - decoys[row];
    require(PeptideFamilyNeuralNet::score(values, decoys, families, features, settings, &changed) == 0,
            "Perturbed-label score");
    bool otherFoldChanged = false;
    for (int row = 0; row < values.size(); ++row) {
        require(parallel.heldOutFold[row] == families[row] % 3, "Family fold assignment");
        double total = 0;
        for (int net = 0; net < settings.networks; ++net)
            total += parallel.perNetwork[row * settings.networks + net];
        const float mean = float(total / settings.networks);
        require(mean == parallel.meanDecoyProbability[row], "Ensemble mean");
        for (int net = 0; net < settings.networks; ++net) {
            const auto index = row * settings.networks + net;
            if (families[row] % 3 == 0)
                require(parallel.perNetwork[index] == changed.perNetwork[index], "Held-out label leakage");
            else
                otherFoldChanged |= parallel.perNetwork[index] != changed.perNetwork[index];
        }
    }
    require(otherFoldChanged, "Label perturbation must affect other folds");
    require(values == original, "Evidence inputs unchanged");
    decoys = originalLabels;
    settings.normalizationGroups = {600};
    PeptideFamilyNeuralNet::Predictions singleGroup;
    require(PeptideFamilyNeuralNet::score(values, decoys, families, features, settings, &singleGroup) == 0,
            "Single normalization group");
    require(singleGroup.perNetwork == parallel.perNetwork, "Global/single-group scale parity");
    settings.normalizationGroups.clear();

    PeptideFamilyNeuralNet::Predictions unchanged;
    unchanged.perNetwork = {123.0f};
    const auto rejects = [&](const auto &x, const auto &y, const auto &groups,
                             const auto &featureList, const auto &options) {
        require(PeptideFamilyNeuralNet::score(x, y, groups, featureList, options, &unchanged)
                == Error::eValueError, "Invalid input accepted");
        require(unchanged.perNetwork == QVector<float>{123.0f}, "Error modified output");
    };
    auto invalidValues = values;
    invalidValues[0].removeLast();
    rejects(invalidValues, decoys, families, features, settings);
    invalidValues = values;
    invalidValues[0][CosineSimSpectrum] = std::numeric_limits<float>::quiet_NaN();
    rejects(invalidValues, decoys, families, features, settings);
    auto invalidLabels = decoys;
    invalidLabels[0] = 2;
    rejects(values, invalidLabels, families, features, settings);
    rejects(values, decoys, QVector<quint32>(600, 0), features, settings);
    auto invalidFeatures = features;
    invalidFeatures.push_back(features.front());
    rejects(values, decoys, families, invalidFeatures, settings);
    auto invalidSettings = settings;
    invalidSettings.normalizationGroups = {300, 299};
    rejects(values, decoys, families, features, invalidSettings);
    invalidSettings = settings;
    invalidSettings.learningRate = std::numeric_limits<double>::quiet_NaN();
    rejects(values, decoys, families, features, invalidSettings);
    invalidSettings = settings;
    invalidSettings.folds = 1;
    rejects(values, decoys, families, features, invalidSettings);
    std::cout << "Family NN: grouping, held-label isolation, serial/parallel and scale parity, "
                 "input/output preservation, and eight invalid-input checks passed\n";
}
