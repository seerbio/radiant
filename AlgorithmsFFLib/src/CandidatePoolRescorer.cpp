#include "CandidatePoolRescorer.h"

#include <QSet>
#include <algorithm>
#include <cmath>
#include <limits>

namespace {
struct Pool {
    QVector<QVector<float>> features;
    QVector<quint32> labels, families;

    void append(const CandidatePoolRescorer::Run &run, int index) {
        features.push_back(run.rawFeatures[index]);
        labels.push_back(run.identities[index].isDecoy);
        families.push_back(PeptideFamilyNeuralNet::familyHash(run.identities[index].originPeptide));
    }
    Error::Err score(const PeptideFamilyNeuralNet::Settings &settings, QVector<float> *probability) const {
        PeptideFamilyNeuralNet::Predictions predictions;
        const auto error = PeptideFamilyNeuralNet::score(
            features, labels, families, PeptideFamilyNeuralNet::nonTimsFeatures(), settings, &predictions);
        if (error == Error::eNoError) *probability = std::move(predictions.meanDecoyProbability);
        return error;
    }
};
}

double CandidatePoolRescorer::equalLogit(double first, double second) {
    first = std::clamp(first, 1e-7, 1.0 - 1e-7);
    second = std::clamp(second, 1e-7, 1.0 - 1e-7);
    return 1.0 / (1.0 + std::exp(-.5 * (std::log(first) - std::log1p(-first)
                                      + std::log(second) - std::log1p(-second))));
}

Error::Err CandidatePoolRescorer::score(
    const QVector<Run> &runs, const Settings &settings, QVector<Result> *results) {
    if (results == nullptr || runs.isEmpty() || settings.originalCandidateLimit < 2
        || settings.expandedCandidateLimit < settings.originalCandidateLimit
        || settings.minimumSharedFragments < 2 || settings.minimumSharedFragments > 12
        || !settings.neuralNet.normalizationGroups.isEmpty()) return Error::eValueError;
    QSet<QString> names;
    qint64 total = 0;
    for (const auto &run : runs) {
        if (run.name.isEmpty() || names.contains(run.name) || run.identities.size() < 2
            || run.identities.size() != run.evidence.size()
            || run.identities.size() != run.rawFeatures.size()) return Error::eValueError;
        names.insert(run.name);
        total += std::min(run.identities.size(), settings.expandedCandidateLimit);
        if (total > std::numeric_limits<int>::max()) return Error::eValueError;
        for (int index = 0; index < run.identities.size(); ++index) {
            const auto &identity = run.identities[index];
            const auto &evidence = run.evidence[index];
            if (identity.reportedPeptide != evidence.reportedPeptide
                || identity.charge != evidence.charge || identity.apex != evidence.apex
                || identity.isDecoy != evidence.isDecoy
                || run.rawFeatures[index].size() != FeaturesSize) return Error::eValueError;
        }
    }
    QVector<Result> output;
    output.reserve(runs.size());
    Pool originalPool, expandedPool;
    // Selection is run-local; pooled normalization and family splitting apply
    // across every selected run, in caller-supplied run order.
    for (const auto &run : runs) {
        Result result;
        result.name = run.name;
        auto error = CandidatePoolSelection::rankedUnique(
            run.identities, settings.expandedCandidateLimit, &result.selectedInputIndices);
        if (error != Error::eNoError) return error;
        Pool originalRun;
        const int originalCount = std::min(settings.originalCandidateLimit, result.selectedInputIndices.size());
        for (int row = 0; row < result.selectedInputIndices.size(); ++row) {
            const int index = result.selectedInputIndices[row];
            expandedPool.append(run, index);
            if (row < originalCount) {
                originalRun.append(run, index);
                originalPool.append(run, index);
            }
        }
        error = originalRun.score(settings.neuralNet, &result.originalPerRunProbability);
        if (error != Error::eNoError) return error;
        output.push_back(std::move(result));
    }
    QVector<float> originalProbability, expandedProbability;
    auto error = originalPool.score(settings.neuralNet, &originalProbability);
    if (error != Error::eNoError) return error;
    error = expandedPool.score(settings.neuralNet, &expandedProbability);
    if (error != Error::eNoError) return error;
    int originalOffset = 0, expandedOffset = 0;
    for (int runIndex = 0; runIndex < runs.size(); ++runIndex) {
        const auto &run = runs[runIndex];
        auto &result = output[runIndex];
        const int originalCount = result.originalPerRunProbability.size();
        const int expandedCount = result.selectedInputIndices.size();
        result.originalPooledProbability = originalProbability.mid(originalOffset, originalCount);
        result.expandedPooledProbability = expandedProbability.mid(expandedOffset, expandedCount);
        QVector<FragmentCompetition::Evidence> selectedEvidence;
        QVector<QString> selectedOrigins;
        selectedEvidence.reserve(expandedCount);
        selectedOrigins.reserve(expandedCount);
        result.combinedProbability.reserve(expandedCount);
        for (int row = 0; row < expandedCount; ++row) {
            selectedEvidence.push_back(run.evidence[result.selectedInputIndices[row]]);
            selectedOrigins.push_back(run.identities[result.selectedInputIndices[row]].originPeptide);
            double probability = result.expandedPooledProbability[row];
            if (row < originalCount) {
                const double original = equalLogit(result.originalPooledProbability[row],
                                                   result.originalPerRunProbability[row]);
                probability = equalLogit(original, probability);
            }
            result.combinedProbability.push_back(probability);
        }
        error = CandidatePoolSelection::finalize(
            selectedEvidence, selectedOrigins, result.combinedProbability,
            settings.minimumSharedFragments, &result.confidence);
        if (error != Error::eNoError) return error;
        for (int &index : result.confidence.inputIndices) index = result.selectedInputIndices[index];
        originalOffset += originalCount;
        expandedOffset += expandedCount;
    }
    *results = std::move(output);
    return Error::eNoError;
}
