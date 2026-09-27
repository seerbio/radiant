#include "CandidatePoolSelection.h"
#include "FragmentCompetitionTestData.h"

#include <QCoreApplication>
#include <QDebug>
#include <limits>

#define CHECK(expression) do { if (!(expression)) { qCritical() << "Failed" << #expression << __LINE__; return 1; } } while (false)

FragmentCompetition::Evidence row(QString peptide, double mass = 1000, bool decoy = false) {
    FragmentCompetition::Evidence e;
    e.reportedPeptide = peptide;
    e.mass = mass;
    e.charge = 2;
    e.apex = 60;
    e.width = 20;
    e.priority = 1;
    e.isDecoy = decoy;
    e.searched = {300, 400, 500, 600, 700};
    e.intensity.fill(100);
    e.cosine.fill(1);
    return e;
}

int main(int argc, char **argv) {
    QCoreApplication app(argc, argv);
    using Identity = CandidatePoolSelection::Identity;
    QVector<Identity> identities = {
        {"A", "A", "1", 2, 1, 0, false, 10},
        {"B", "B", "2", 2, 0, 0, false, 20},
        {"A", "A", "1", 2, 0, 1, false, 10},
        {"C", "C", "3", 2, 1, 1, true, 30},
        {"A", "A", "1", 2, 0, 2, true, 10},
        {"A", "A", "1", 2, 1, 2, false, 11},
    };
    QVector<int> selected = {999};
    CHECK(CandidatePoolSelection::rankedUnique(identities, 5, &selected) == Error::eNoError);
    CHECK(selected == QVector<int>({1, 0, 3, 4, 5}));
    CHECK(CandidatePoolSelection::rankedUnique(identities, 2, &selected) == Error::eNoError);
    CHECK(selected == QVector<int>({1, 0}));
    auto invalid = identities;
    invalid.last().apex = std::numeric_limits<double>::quiet_NaN();
    CHECK(CandidatePoolSelection::rankedUnique(invalid, 2, &selected) != Error::eNoError);
    CHECK(selected == QVector<int>({1, 0}));
    CHECK(CandidatePoolSelection::rankedUnique({}, 1, &selected) == Error::eNoError && selected.isEmpty());

    // Both endpoints lie outside the library reader's 200..1500 filter.
    // Captured searched evidence must keep them without filtering or mutation.
    auto a = row("A"), b = row("B", 1000, true);
    a.searched = {100, 1600, 800};
    b.searched = {100, 1600, 900};
    a.intensity[2] = 200;
    CHECK(FragmentCompetition::retainEvidence({a, b}, 2, &selected) == Error::eNoError);
    CHECK(selected == QVector<int>({0}));
    a.intensity[2] = 100;
    CHECK(FragmentCompetition::retainEvidence({a, b}, 2, &selected) == Error::eNoError);
    CHECK(selected == QVector<int>({1})); // equal evidence/priority prefers decoy
    CHECK(FragmentCompetition::retainEvidence({b, a}, 2, &selected) == Error::eNoError);
    CHECK(selected == QVector<int>({0}));
    b.width = 0;
    CHECK(FragmentCompetition::retainEvidence({a, b}, 2, &selected) != Error::eNoError);
    CHECK(selected == QVector<int>({0}));

    // Independently constructed native and captured rows share one outcome.
    CompetitionCandidateFixture first("A", {300, 400, 800}), second("B", {300, 400, 900});
    first.scores.integrations[2] = 200;
    QVector<CandidateScores*> native = {&first.scores, &second.scores};
    CHECK(FragmentCompetition::removeCompetingCandidates(true, &native, 2) == Error::eNoError);
    a = row("A"); b = row("B");
    a.searched = {300, 400, 800}; b.searched = {300, 400, 900};
    a.intensity[2] = 200;
    CHECK(FragmentCompetition::retainEvidence({a, b}, 2, &selected) == Error::eNoError);
    CHECK(native == QVector<CandidateScores*>({&first.scores}) && selected == QVector<int>({0}));

    // Peak selection precedes physical competition: the duplicate weak peak
    // cannot make the best peak lose its unique fragment to itself.
    auto duplicate = a;
    duplicate.apex += 1;
    CandidatePoolSelection::Confidence confidence;
    CHECK(CandidatePoolSelection::finalize({a, duplicate, b}, {.1, .2, .3}, 2, &confidence) == Error::eNoError);
    CHECK(confidence.inputIndices == QVector<int>({0}));
    CHECK(confidence.qValues == QVector<double>({1}));
    // Exact ties form one confidence block, including the decoy first.
    QVector<FragmentCompetition::Evidence> independent = {
        row("_A_", 1000), row("A", 1000), row("B", 1100), row("D", 1200, true),
        row("C", 1300), row("E", 1400)};
    CHECK(CandidatePoolSelection::finalize(independent, {.1, .2, .1, .1, .3, .4}, 2, &confidence) == Error::eNoError);
    CHECK(confidence.inputIndices == QVector<int>({3, 0, 2, 4, 5}));
    CHECK(confidence.qValues == QVector<double>({.5, .5, .5, .5, .5}));
    const auto previous = confidence;
    CHECK(CandidatePoolSelection::finalize(independent, {.1}, 2, &confidence) != Error::eNoError);
    CHECK(confidence.inputIndices == previous.inputIndices && confidence.qValues == previous.qValues);
    CHECK(CandidatePoolSelection::finalize({}, {}, 2, &confidence) == Error::eNoError);
    CHECK(confidence.inputIndices.isEmpty() && confidence.qValues.isEmpty());

    // Synthetic decoy display sequences can collide for distinct origins,
    // including hypotheses of different masses. Preserve both null units,
    // while collapsing repeat peaks of the same origin before competition.
    QVector<FragmentCompetition::Evidence> collisions = {
        row("ALFDLFENLK", 1235.631592, true),
        row("ALFDLFENLK", 1207.625366, true),
        row("_ALFDLFENLK_", 1235.631592, true)};
    collisions[2].apex += 1;
    QVector<QString> origins = {"APFDLFENRK", "APFDLFENKK", "_APFDLFENRK_"};
    QVector<double> probabilities = {.1, .1, .2};
    for (int index = 0; index < 8; ++index) {
        const auto peptide = QString("TARGET%1").arg(index);
        collisions.push_back(row(peptide, 1500 + 100 * index));
        origins.push_back(peptide);
        probabilities.push_back(.1);
    }
    CHECK(CandidatePoolSelection::finalize(collisions, origins, probabilities, 2, &confidence) == Error::eNoError);
    CHECK(confidence.inputIndices == QVector<int>({0, 1, 3, 4, 5, 6, 7, 8, 9, 10}));
    CHECK(confidence.qValues == QVector<double>(10, 3.0 / 8.0));
    const auto corrected = confidence;
    origins[1] = "_";
    CHECK(CandidatePoolSelection::finalize(collisions, origins, probabilities, 2, &confidence) != Error::eNoError);
    CHECK(confidence.inputIndices == corrected.inputIndices && confidence.qValues == corrected.qValues);
    CHECK(CandidatePoolSelection::finalize(collisions, {}, probabilities, 2, &confidence) != Error::eNoError);
    CHECK(confidence.inputIndices == corrected.inputIndices && confidence.qValues == corrected.qValues);
    // Legacy unique-name callers retain their behavior and cannot silently
    // substitute the display name for a supplied origin.
    CHECK(CandidatePoolSelection::finalize(collisions, probabilities, 2, &confidence) == Error::eNoError);
    CHECK(confidence.inputIndices == QVector<int>({0, 3, 4, 5, 6, 7, 8, 9, 10}));
    CHECK(confidence.qValues == QVector<double>(9, 2.0 / 8.0));
    // A library decoy and a synthetic decoy can share an origin but represent
    // different reported sequences. The origin alone must not merge them.
    origins[1] = "APFDLFENKK";
    auto libraryDecoy = row("APFDLFENRK", 1235.631592, true);
    libraryDecoy.apex = 160;
    collisions.push_back(libraryDecoy);
    origins.push_back("APFDLFENRK");
    probabilities.push_back(.1);
    CHECK(CandidatePoolSelection::finalize(collisions, origins, probabilities, 2, &confidence) == Error::eNoError);
    CHECK(confidence.inputIndices == QVector<int>({0, 1, 11, 3, 4, 5, 6, 7, 8, 9, 10}));
    CHECK(confidence.qValues == QVector<double>(11, 4.0 / 8.0));
    CHECK(CandidatePoolSelection::finalize({}, {}, {}, 2, &confidence) == Error::eNoError);
    CHECK(confidence.inputIndices.isEmpty() && confidence.qValues.isEmpty());
    qInfo() << "Native selection, shared physical evidence, decoy symmetry, ties and confidence probes passed";
    return 0;
}
