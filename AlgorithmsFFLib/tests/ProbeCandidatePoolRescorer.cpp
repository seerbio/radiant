#include "CandidatePoolRescorer.h"

#include <QCoreApplication>
#include <QDebug>
#include <cmath>
#include <limits>

#define CHECK(expression) do { if (!(expression)) { qCritical() << "Failed" << #expression << __LINE__; return 1; } } while (false)

int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    QVector<CandidatePoolRescorer::Run> runs(2);
    for (int sample = 0; sample < runs.size(); ++sample) {
        auto &run = runs[sample];
        run.name = QString::number(sample);
        // Deliberately reversed input view order; view zero wins duplicates.
        for (int view : {1, 0}) {
            for (int family = 0; family < 200; ++family) {
                for (int label = 0; label < 2; ++label) {
                    CandidatePoolSelection::Identity id;
                    id.originPeptide = "PEPTIDE" + QString::number(family);
                    id.reportedPeptide = (label ? "DECOY" : "") + id.originPeptide;
                    id.targetKey = QString::number(family);
                    id.charge = 2;
                    id.view = view;
                    id.viewRank = family * 2 + label;
                    id.isDecoy = label;
                    id.apex = 60 + label * 40;
                    run.identities.push_back(id);
                    FragmentCompetition::Evidence evidence;
                    evidence.reportedPeptide = id.reportedPeptide;
                    evidence.mass = 1000 + family * 10;
                    evidence.charge = id.charge;
                    evidence.apex = id.apex;
                    evidence.width = 10;
                    evidence.isDecoy = id.isDecoy;
                    evidence.searched = {300, 400, 500, 600};
                    evidence.intensity.fill(100);
                    evidence.cosine.fill(1);
                    run.evidence.push_back(evidence);
                    QVector<float> features(FeaturesSize);
                    for (int feature = 0; feature < FeaturesSize; ++feature)
                        features[feature] = float((family * 31 + feature * 13 + sample * 7 + label * 3) % 97) / 97;
                    features[Mass] = evidence.mass;
                    features[Charge] = evidence.charge;
                    run.rawFeatures.push_back(std::move(features));
                }
            }
        }
    }
    CandidatePoolRescorer::Settings settings;
    settings.originalCandidateLimit = 300;
    settings.expandedCandidateLimit = 400;
    settings.neuralNet.networks = 2;
    settings.neuralNet.epochs = 3;
    settings.neuralNet.threads = 2;
    QVector<CandidatePoolRescorer::Result> first, changed;
    CHECK(CandidatePoolRescorer::score(runs, settings, &first) == Error::eNoError);
    CHECK(first.size() == 2);
    for (const auto &result : first) {
        CHECK(result.selectedInputIndices.size() == 400);
        CHECK(result.originalPerRunProbability.size() == 300);
        CHECK(result.originalPooledProbability.size() == 300);
        CHECK(result.expandedPooledProbability.size() == 400);
        CHECK(!result.confidence.inputIndices.isEmpty());
        for (int index = 0; index < 400; ++index) {
            CHECK(result.selectedInputIndices[index] == 400 + index);
            CHECK(std::isfinite(result.combinedProbability[index]));
            if (index >= 300) CHECK(result.combinedProbability[index] == result.expandedPooledProbability[index]);
        }
        for (int index : result.confidence.inputIndices) CHECK(index >= 400 && index < 800);
    }
    auto perturbed = runs;
    for (auto &run : perturbed) {
        for (int index = 0; index < run.identities.size(); ++index) {
            if (PeptideFamilyNeuralNet::familyHash(run.identities[index].originPeptide) % 3 == 0) {
                run.identities[index].isDecoy = !run.identities[index].isDecoy;
                run.evidence[index].isDecoy = run.identities[index].isDecoy;
            }
        }
    }
    CHECK(CandidatePoolRescorer::score(perturbed, settings, &changed) == Error::eNoError);
    bool otherPredictionChanged = false;
    for (int sample = 0; sample < 2; ++sample) {
        CHECK(first[sample].selectedInputIndices == changed[sample].selectedInputIndices);
        for (int index = 0; index < 400; ++index) {
            const int original = first[sample].selectedInputIndices[index];
            const auto hash = PeptideFamilyNeuralNet::familyHash(runs[sample].identities[original].originPeptide);
            if (hash % 3 == 0) {
                CHECK(first[sample].combinedProbability[index] == changed[sample].combinedProbability[index]);
            } else {
                otherPredictionChanged |= first[sample].combinedProbability[index] != changed[sample].combinedProbability[index];
            }
        }
    }
    CHECK(otherPredictionChanged);
    const auto saved = first;
    perturbed = runs;
    perturbed.last().name = runs.first().name;
    CHECK(CandidatePoolRescorer::score(perturbed, settings, &first) != Error::eNoError);
    CHECK(first.first().combinedProbability == saved.first().combinedProbability);
    perturbed = runs;
    perturbed.last().evidence.last().reportedPeptide = "MISMATCH";
    CHECK(CandidatePoolRescorer::score(perturbed, settings, &first) != Error::eNoError);
    CHECK(first.first().combinedProbability == saved.first().combinedProbability);
    settings.neuralNet.normalizationGroups = {800};
    CHECK(CandidatePoolRescorer::score(runs, settings, &first) != Error::eNoError);
    CHECK(first.first().combinedProbability == saved.first().combinedProbability);
    CHECK(std::abs(CandidatePoolRescorer::equalLogit(.1, .2)
                   - CandidatePoolRescorer::equalLogit(.2, .1)) < 1e-15);
    qInfo() << "Complete native pool scoring: selection, output mapping, cross-run family isolation and invalid inputs passed";
    return 0;
}
