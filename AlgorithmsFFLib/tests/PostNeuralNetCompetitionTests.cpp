#include "PostNeuralNetCompetition.h"
#include "FragmentCompetitionTestData.h"
#include "MathUtils.h"
#include "AminoAcids.h"

#include <QtTest/QtTest>
#include <limits>
#include <memory>

class PostNeuralNetCompetitionTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void disabledIsExactNoOp();
    void rejectsInvalidInputsWithoutMutation();
    void restoresLdaOnPhysicalCompetitionError();
    void physicalEvidencePrecedesNeuralScore();
    void neuralScoreBreaksPhysicalTiesAndRestoresLda();
    void exactProbabilityTiesShareConservativeConfidence();
    void removingAnApproximateTieStillProducesExactScoreOrder();
    void duplicatePeaksCannotInflateConfidence();
    void libraryAndSyntheticDecoysBothCount();
    void distinctDecoyOriginsWithIdenticalDisplaySequences();
    void distinctChargesAndModificationsRemainDistinct();
    void emptyOrUnresolvedInput();
};

void PostNeuralNetCompetitionTests::disabledIsExactNoOp() {
    CompetitionCandidateFixture fixture;
    fixture.scores.classifierScore = 0.1;
    fixture.scores.qValue = 0.123;
    fixture.scores.precursorQValue = 0.234;
    QVector<CandidateScores*> rows{nullptr, &fixture.scores, &fixture.scores};
    const auto original = rows;
    QCOMPARE(PostNeuralNetCompetition::apply(0, &rows), eNoError);
    QCOMPARE(rows, original);
    QCOMPARE(fixture.scores.qValue, 0.123);
    QCOMPARE(fixture.scores.precursorQValue, 0.234);
}

void PostNeuralNetCompetitionTests::rejectsInvalidInputsWithoutMutation() {
    CompetitionCandidateFixture fixture;
    fixture.scores.classifierScore = 0.1;
    fixture.scores.qValue = 0.123;
    QVector<CandidateScores*> rows{&fixture.scores};
    const auto original = rows;
    QCOMPARE(PostNeuralNetCompetition::apply(2, nullptr), eValueError);
    for (const int threshold : {-1, 1, 13}) {
        QCOMPARE(PostNeuralNetCompetition::apply(threshold, &rows), eValueError);
        QCOMPARE(rows, original);
    }
    for (const double score : {-0.1, 1.1, std::numeric_limits<double>::quiet_NaN()}) {
        fixture.scores.classifierScore = score;
        QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eValueError);
        QCOMPARE(rows, original);
        QCOMPARE(fixture.scores.discriminantScore, 1.0);
        QCOMPARE(fixture.scores.qValue, 0.123);
    }
    fixture.scores.classifierScore = 0.1;
    rows.push_back(&fixture.scores);
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eValueError);
    QCOMPARE(rows.size(), 2);
}

void PostNeuralNetCompetitionTests::restoresLdaOnPhysicalCompetitionError() {
    CompetitionCandidateFixture fixture;
    fixture.scores.classifierScore = 0.1;
    fixture.scores.discriminantScore = 42.0;
    fixture.scores.featuresArray[Mass] = std::numeric_limits<float>::quiet_NaN();
    QVector<CandidateScores*> rows{&fixture.scores};
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eValueError);
    QCOMPARE(rows, QVector<CandidateScores*>({&fixture.scores}));
    QCOMPARE(fixture.scores.discriminantScore, 42.0);
}

void PostNeuralNetCompetitionTests::physicalEvidencePrecedesNeuralScore() {
    CompetitionCandidateFixture weak;
    CompetitionCandidateFixture strong("PEPTIDER", {300, 400, 500, 600, 800});
    weak.scores.classifierScore = 0.001;
    strong.scores.classifierScore = 0.8;
    weak.scores.integrations[4] = 1.0f;
    strong.scores.integrations[4] = 100.0f;
    strong.scores.proteinGroup = "entrapment";
    const auto features = strong.scores.featuresArray;
    const auto intensities = strong.scores.integrations;
    QVector<CandidateScores*> rows{&weak.scores, &strong.scores};
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&strong.scores}));
    QCOMPARE(strong.scores.featuresArray, features);
    QCOMPARE(strong.scores.integrations, intensities);
    QCOMPARE(strong.scores.classifierScore, 0.8);
    QCOMPARE(strong.scores.qValue, 1.0); // The +1 cannot vanish with no decoys.
}

void PostNeuralNetCompetitionTests::neuralScoreBreaksPhysicalTiesAndRestoresLda() {
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture second("PEPTIDER", {300, 400, 500, 600, 800});
    first.scores.discriminantScore = 100.0;
    second.scores.discriminantScore = -100.0;
    first.scores.classifierScore = 0.2;
    second.scores.classifierScore = 0.1;
    QVector<CandidateScores*> rows{&first.scores, &second.scores};
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&second.scores}));
    QCOMPARE(first.scores.discriminantScore, 100.0);
    QCOMPARE(second.scores.discriminantScore, -100.0);
    QCOMPARE(second.scores.classifierScore, 0.1);
}

void PostNeuralNetCompetitionTests::exactProbabilityTiesShareConservativeConfidence() {
    CompetitionCandidateFixture first("PEPTIDEK");
    CompetitionCandidateFixture second("PEPTIDER");
    CompetitionCandidateFixture tiedDecoy("PEPTIDEA");
    CompetitionCandidateFixture last("PEPTIDEC");
    CompetitionCandidateFixture laterDecoy("PEPTIDEG");
    QVector<CandidateScores*> rows{&first.scores, &second.scores, &tiedDecoy.scores,
                                   &last.scores, &laterDecoy.scores};
    for (int index = 0; index < rows.size(); ++index) rows[index]->scanTime += 100 * index;
    first.scores.classifierScore = 0.0;
    second.scores.classifierScore = tiedDecoy.scores.classifierScore = 0.1;
    last.scores.classifierScore = 0.2;
    laterDecoy.scores.classifierScore = 1.0;
    tiedDecoy.scores.isDecoy = laterDecoy.scores.isDecoy = true;
    std::reverse(rows.begin(), rows.end());
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&first.scores, &tiedDecoy.scores,
        &second.scores, &last.scores, &laterDecoy.scores}));
    for (auto *candidate : {&first.scores, &second.scores, &tiedDecoy.scores, &last.scores}) {
        QCOMPARE(candidate->qValue, 2.0 / 3.0);
        QCOMPARE(candidate->precursorQValue, candidate->qValue);
        QCOMPARE(candidate->isBestPrecursorCandidate, 1);
    }
    QCOMPARE(laterDecoy.scores.qValue, 1.0);
}

void PostNeuralNetCompetitionTests::removingAnApproximateTieStillProducesExactScoreOrder() {
    CompetitionCandidateFixture first("PEPTIDEK");
    CompetitionCandidateFixture middle("PEPTIDER");
    CompetitionCandidateFixture last("PEPTIDEA");
    CompetitionCandidateFixture competitor("PEPTIDEC", {300, 400, 500, 600, 800});
    const double tolerance = S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL;
    first.scores.classifierScore = static_cast<float>(0.001 + 1.5 * tolerance);
    middle.scores.classifierScore = static_cast<float>(0.001 + 0.75 * tolerance);
    last.scores.classifierScore = static_cast<float>(0.001);
    competitor.scores.classifierScore = 0.01;
    first.scores.discriminantScore = 3.0;
    middle.scores.discriminantScore = 2.0;
    last.scores.discriminantScore = 1.0;
    competitor.scores.discriminantScore = 0.0;
    middle.library.mass = competitor.library.mass = 1100.0;
    middle.scores.featuresArray[Mass] = competitor.scores.featuresArray[Mass] = 1100.0f;
    last.library.mass = 1200.0;
    last.scores.featuresArray[Mass] = 1200.0f;
    middle.scores.integrations[4] = 1.0f;
    const auto approximate = [](const CandidateScores *left, const CandidateScores *right) {
        if (MathUtils::tSame(left->classifierScore, right->classifierScore,
                            S_GLOBAL_SETTINGS.ROUNDING_PRECISION_DECIMAL))
            return left->discriminantScore > right->discriminantScore;
        return left->classifierScore < right->classifierScore;
    };
    QVector<CandidateScores*> rows{
        &first.scores, &middle.scores, &last.scores, &competitor.scores};
    QVERIFY(std::is_sorted(rows.begin(), rows.end(), approximate));
    const QVector<CandidateScores*> subsequence{
        &first.scores, &last.scores, &competitor.scores};
    QVERIFY(!std::is_sorted(subsequence.begin(), subsequence.end(), approximate));
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({
        &last.scores, &first.scores, &competitor.scores}));
    QVERIFY(std::is_sorted(rows.begin(), rows.end(),
        [](const CandidateScores *left, const CandidateScores *right) {
            return left->classifierScore < right->classifierScore;
        }));
}

void PostNeuralNetCompetitionTests::duplicatePeaksCannotInflateConfidence() {
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture duplicate;
    CompetitionCandidateFixture second("PEPTIDER");
    duplicate.scores.scanTime = 160.0f;
    second.scores.scanTime = 260.0f;
    first.scores.classifierScore = 0.2;
    duplicate.scores.classifierScore = 0.1;
    second.scores.classifierScore = 0.3;
    QVector<CandidateScores*> rows{&first.scores, &second.scores, &duplicate.scores};
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&duplicate.scores, &second.scores}));
    QCOMPARE(second.scores.qValue, 0.5);
    QCOMPARE(duplicate.scores.qValue, 0.5);
}

void PostNeuralNetCompetitionTests::libraryAndSyntheticDecoysBothCount() {
    CompetitionCandidateFixture target;
    CompetitionCandidateFixture synthetic;
    CompetitionCandidateFixture library;
    synthetic.scores.scanTime = 160.0f;
    library.scores.scanTime = 260.0f;
    target.scores.classifierScore = 0.1;
    synthetic.scores.classifierScore = library.scores.classifierScore = 0.2;
    synthetic.scores.isDecoy = true;
    library.library.isDecoy = true;
    QVector<CandidateScores*> rows{&target.scores, &synthetic.scores, &library.scores};
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QCOMPARE(rows.size(), 3); // Library decoy must not collapse with the target.
    for (const auto *candidate : rows) QCOMPARE(candidate->qValue, 1.0);
}

void PostNeuralNetCompetitionTests::distinctDecoyOriginsWithIdenticalDisplaySequences() {
    CompetitionCandidateFixture first("APFDLFENRK"), second("APFDLFENKK"), duplicate("APFDLFENRK");
    QCOMPARE(AminoAcids::mutatePenultimatePeptideResidues(first.pair.peptideStringWithMods()),
             AminoAcids::mutatePenultimatePeptideResidues(second.pair.peptideStringWithMods()));
    QVector<CandidateScores*> rows{&first.scores, &second.scores, &duplicate.scores};
    for (auto *candidate : rows) {
        candidate->isDecoy = true;
        candidate->classifierScore = .1;
    }
    second.scores.scanTime = 160;
    duplicate.scores.scanTime = 260;
    duplicate.scores.classifierScore = .2;
    std::vector<std::unique_ptr<CompetitionCandidateFixture>> targets;
    for (const auto &sequence : {"PEPTIDEK", "PEPTIDER", "PEPTIDEA", "PEPTIDEC",
                                 "PEPTIDEG", "PEPTIDEL", "PEPTIDEV", "PEPTIDEW"}) {
        auto target = std::make_unique<CompetitionCandidateFixture>(sequence);
        target->scores.scanTime = 360 + 100 * targets.size();
        target->scores.classifierScore = .1;
        rows.push_back(&target->scores);
        targets.push_back(std::move(target));
    }
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QCOMPARE(rows.size(), 10);
    QVERIFY(rows.contains(&first.scores));
    QVERIFY(rows.contains(&second.scores));
    QVERIFY(!rows.contains(&duplicate.scores));
    for (const auto *candidate : rows) QCOMPARE(candidate->qValue, 3.0 / 8.0);
}

void PostNeuralNetCompetitionTests::distinctChargesAndModificationsRemainDistinct() {
    CompetitionCandidateFixture unmodified("PEPTIDEK");
    CompetitionCandidateFixture modified("PEPT(Oxidation)IDEK");
    CompetitionCandidateFixture charged("PEPTIDEK");
    charged.library.precursorCharge = 3;
    charged.scores.featuresArray[Charge] = 3;
    modified.scores.scanTime = 160.0f;
    charged.scores.scanTime = 260.0f;
    QVector<CandidateScores*> rows{&unmodified.scores, &modified.scores, &charged.scores};
    for (auto *candidate : rows) candidate->classifierScore = 0.1;
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QCOMPARE(rows.size(), 3);
    for (const auto *candidate : rows) QCOMPARE(candidate->qValue, 1.0 / 3.0);
}

void PostNeuralNetCompetitionTests::emptyOrUnresolvedInput() {
    QVector<CandidateScores*> empty;
    QCOMPARE(PostNeuralNetCompetition::apply(2, &empty), eNoError);
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture second;
    first.scores.classifierScore = second.scores.classifierScore = 0.1;
    QVector<CandidateScores*> rows{&first.scores, &second.scores};
    QCOMPARE(PostNeuralNetCompetition::apply(2, &rows), eNoError);
    QVERIFY(rows.isEmpty());
}

QTEST_MAIN(PostNeuralNetCompetitionTests)
#include "PostNeuralNetCompetitionTests.moc"
