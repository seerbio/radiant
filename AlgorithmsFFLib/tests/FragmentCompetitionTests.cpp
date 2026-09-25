#include "FragmentCompetition.h"
#include "FragmentCompetitionTestData.h"

#include <QtTest/QtTest>

#include <limits>
#include <memory>
#include <vector>

class FragmentCompetitionTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void unsharedEvidenceBeatsScoresAndLabels();
    void requiresMatchingPrecursorAndPeak_data();
    void requiresMatchingPrecursorAndPeak();
    void acceptsBoundaryMatchesAcrossScans();
    void needsFourDistinctSupportedFragments_data();
    void needsFourDistinctSupportedFragments();
    void rejectsUnresolvedGroups();
    void resolvesTransitiveComponents();
    void unsupportedCompetitorFragmentsPreventUniqueness();
    void cosineWeightedEvidenceAndDeduplication();
    void usesLdaThenDecoyForTies_data();
    void usesLdaThenDecoyForTies();
    void disabledIsExactNoOp();
    void skipsNoPeakPlaceholders();
    void rejectsMalformedObservedRowsWithoutMutation();
    void matchesPythonReference();
};

void FragmentCompetitionTests::unsharedEvidenceBeatsScoresAndLabels() {
    CompetitionCandidateFixture weak("PEPTIDEK", {300, 400, 500, 600, 700});
    CompetitionCandidateFixture strong("PEPTIDER", {300, 400, 500, 600, 800});
    weak.scores.integrations[4] = 0.0f;
    weak.scores.discriminantScore = 100.0;
    weak.scores.classifierScore = 0.001;
    weak.scores.qValue = 0.001;
    weak.scores.proteinGroup = "human";
    strong.scores.integrations[4] = 80.0f;
    strong.scores.featuresArray[IntensityFoundMax5] = 0.0f;
    strong.scores.classifierScore = 0.99;
    strong.scores.qValue = 0.99;
    strong.scores.proteinGroup = "entrapment";
    const auto originalFeatures = strong.scores.featuresArray;
    const auto originalIntegrations = strong.scores.integrations;
    QVector<CandidateScores*> rows = {&weak.scores, &strong.scores};
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&strong.scores}));
    QCOMPARE(strong.scores.featuresArray, originalFeatures);
    QCOMPARE(strong.scores.integrations, originalIntegrations);
    QCOMPARE(strong.scores.classifierScore, 0.99);
    QCOMPARE(strong.scores.qValue, 0.99);
    QCOMPARE(weak.scores.discriminantScore, 100.0);
    std::swap(weak.scores.proteinGroup, strong.scores.proteinGroup);
    rows = {&strong.scores, &weak.scores};
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&strong.scores}));
}

void FragmentCompetitionTests::requiresMatchingPrecursorAndPeak_data() {
    QTest::addColumn<int>("change");
    QTest::newRow("different-charge") << 0;
    QTest::newRow("outside-precursor-ppm") << 1;
    QTest::newRow("outside-apex-window") << 2;
    QTest::newRow("narrower-width-controls-window") << 3;
}

void FragmentCompetitionTests::requiresMatchingPrecursorAndPeak() {
    QFETCH(int, change);
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture second("PEPTIDER", {300, 400, 500, 600, 800});
    if (change == 0) second.scores.featuresArray[Charge] = 3;
    if (change == 1) second.scores.featuresArray[Mass] = 1000.01f;
    if (change == 2) second.scores.scanTime = 70.01f;
    if (change == 3) {
        second.scores.scanTime = 62.0f;
        second.scores.scanTimeStart = 61.0f;
        second.scores.scanTimeEnd = 63.0f;
    }
    QVector<CandidateScores*> rows = {&first.scores, &second.scores};
    const auto original = rows;
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, original);
}

void FragmentCompetitionTests::acceptsBoundaryMatchesAcrossScans() {
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture second(
        "PEPTIDER", {300.0059f, 400.0079f, 500.0099f, 600.0119f, 800});
    CompetitionCandidateFixture singleton("PEPTIDEA", {810, 820, 830, 840, 850});
    second.scores.featuresArray[Mass] = 1000.0048828125f;
    second.scores.scanNumber = 101;
    second.scores.scanTime = 70.0f; // Exactly half of the narrower width.
    second.scores.integrations[4] = 1;
    QVector<CandidateScores*> rows = {&singleton.scores, &second.scores, &first.scores};
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&singleton.scores, &first.scores}));
}

void FragmentCompetitionTests::needsFourDistinctSupportedFragments_data() {
    QTest::addColumn<int>("change");
    QTest::newRow("only-three-shared") << 0;
    QTest::newRow("duplicate-fragments-count-once") << 1;
    QTest::newRow("low-trace-cosine") << 2;
    QTest::newRow("zero-intensity") << 3;
    QTest::newRow("nonfinite-trace") << 4;
    QTest::newRow("outside-fragment-ppm") << 5;
}

void FragmentCompetitionTests::needsFourDistinctSupportedFragments() {
    QFETCH(int, change);
    const QVector<float> firstMzs = change == 1
        ? QVector<float>({300, 300.001f, 400, 500, 700})
        : QVector<float>({300, 400, 500, 600, 700});
    QVector<float> secondMzs = firstMzs;
    secondMzs[4] = 800;
    if (change == 0) secondMzs[3] = 650;
    if (change == 5) secondMzs[3] = 600.02f;
    CompetitionCandidateFixture first("PEPTIDEK", firstMzs);
    CompetitionCandidateFixture second("PEPTIDER", secondMzs);
    if (change == 2) second.scores.featuresArray[CosineSimToAnchor4] = 0.49f;
    if (change == 3) second.scores.integrations[3] = 0;
    if (change == 4) second.scores.featuresArray[CosineSimToAnchor4] = std::numeric_limits<float>::quiet_NaN();
    QVector<CandidateScores*> rows = {&first.scores, &second.scores};
    const auto original = rows;
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, original);
}

void FragmentCompetitionTests::rejectsUnresolvedGroups() {
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture second;
    CompetitionCandidateFixture singleton("PEPTIDER", {810, 820, 830, 840});
    QVector<CandidateScores*> rows = {&singleton.scores, &second.scores, &first.scores};
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&singleton.scores}));
}

void FragmentCompetitionTests::resolvesTransitiveComponents() {
    CompetitionCandidateFixture first("PEPTIDEK", {300, 400, 500, 600, 700});
    CompetitionCandidateFixture middle("PEPTIDER", {300, 400, 500, 600, 810, 820, 830, 840, 850});
    CompetitionCandidateFixture last("PEPTIDEA", {810, 820, 830, 840, 900});
    first.scores.integrations[4] = 20;
    middle.scores.integrations[8] = 30;
    last.scores.integrations[4] = 80;
    QVector<CandidateScores*> rows = {&first.scores, &last.scores, &middle.scores};
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&last.scores}));
}

void FragmentCompetitionTests::unsupportedCompetitorFragmentsPreventUniqueness() {
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture second("PEPTIDER", {300, 400, 500, 600, 700, 800});
    first.scores.integrations[4] = 1000;
    second.scores.integrations[4] = 0;
    second.scores.integrations[5] = 1;
    QVector<CandidateScores*> rows = {&first.scores, &second.scores};
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&second.scores}));
}

void FragmentCompetitionTests::cosineWeightedEvidenceAndDeduplication() {
    CompetitionCandidateFixture first("PEPTIDEK", {300, 400, 500, 600, 700, 700.001f});
    CompetitionCandidateFixture second("PEPTIDER", {300, 400, 500, 600, 800});
    first.scores.integrations[4] = 100;
    first.scores.integrations[5] = 200;
    first.scores.featuresArray[CosineSimToAnchor5] = 0.5f;
    first.scores.featuresArray[CosineSimToAnchor6] = 0.5f;
    second.scores.integrations[4] = 60;
    second.scores.featuresArray[CosineSimToAnchor5] = 2.0f; // Clip to one.
    QVector<CandidateScores*> rows = {&first.scores, &second.scores};
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({&second.scores}));
}

void FragmentCompetitionTests::usesLdaThenDecoyForTies_data() {
    QTest::addColumn<bool>("equalLda");
    QTest::addColumn<bool>("syntheticDecoy");
    QTest::newRow("lda-before-decoy") << false << true;
    QTest::newRow("synthetic-decoy-tie") << true << true;
    QTest::newRow("library-decoy-tie") << true << false;
}

void FragmentCompetitionTests::usesLdaThenDecoyForTies() {
    QFETCH(bool, equalLda);
    QFETCH(bool, syntheticDecoy);
    CompetitionCandidateFixture target;
    CompetitionCandidateFixture decoy("PEPTIDER", {300, 400, 500, 600, 800});
    decoy.scores.isDecoy = syntheticDecoy;
    decoy.library.isDecoy = !syntheticDecoy;
    target.scores.discriminantScore = equalLda ? 1.0 : 2.0;
    target.scores.classifierScore = 1.0;
    decoy.scores.classifierScore = 0.0;
    QVector<CandidateScores*> rows = {&target.scores, &decoy.scores};
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows, QVector<CandidateScores*>({equalLda ? &decoy.scores : &target.scores}));
}

void FragmentCompetitionTests::disabledIsExactNoOp() {
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture second;
    CandidateScores invalid;
    QVector<CandidateScores*> rows = {&second.scores, nullptr, &invalid, &first.scores};
    const auto original = rows;
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(false, &rows), eNoError);
    QCOMPARE(rows, original);
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, nullptr), eValueError);
}

void FragmentCompetitionTests::skipsNoPeakPlaceholders() {
    CandidateScores placeholder;
    placeholder.initFeaturesArray();
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture second("PEPTIDER", {300, 400, 500, 600, 800});
    first.scores.integrations[4] = 0;
    QVector<CandidateScores*> rows(200000, &placeholder);
    rows.push_back(nullptr);
    rows.push_back(&first.scores);
    rows.push_back(&second.scores);
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    QCOMPARE(rows.size(), 200001);
    QCOMPARE(rows.front(), &placeholder);
    QCOMPARE(rows.back(), &second.scores);
}

void FragmentCompetitionTests::rejectsMalformedObservedRowsWithoutMutation() {
    CompetitionCandidateFixture first;
    CompetitionCandidateFixture bad;
    bad.scores.scanTimeEnd = bad.scores.scanTimeStart;
    QVector<CandidateScores*> rows = {&first.scores, &bad.scores};
    const auto original = rows;
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eValueError);
    QCOMPARE(rows, original);
}

void FragmentCompetitionTests::matchesPythonReference() {
    const QString location = qEnvironmentVariable(
        "RADIANT_COMPETITION_REFERENCE",
        QDir(qApp->applicationDirPath()).filePath("fragment_competition_reference.json"));
    QFile file(location);
    QVERIFY(file.open(QIODevice::ReadOnly));
    QJsonParseError parseError;
    const QJsonDocument document = QJsonDocument::fromJson(file.readAll(), &parseError);
    QCOMPARE(parseError.error, QJsonParseError::NoError);
    QVERIFY(document.isObject());
    const QJsonObject reference = document.object();
    const QJsonArray input = reference["rows"].toArray();
    QVERIFY(!input.isEmpty());
    QVERIFY(reference["keep"].isArray());
    std::vector<std::unique_ptr<CompetitionCandidateFixture>> fixtures;
    QVector<CandidateScores*> rows;
    for (const QJsonValue &value : input) {
        const QJsonObject row = value.toObject();
        QVector<float> mzs;
        for (const QJsonValue &mz : row["mzs"].toArray()) mzs.push_back(mz.toDouble());
        auto fixture = std::make_unique<CompetitionCandidateFixture>(row["peptide"].toString(), mzs);
        fixture->library.isDecoy = row["decoy"].toBool();
        fixture->scores.featuresArray[Mass] = row["mass"].toDouble();
        fixture->scores.featuresArray[Charge] = row["charge"].toDouble();
        fixture->scores.scanTime = row["apex"].toDouble();
        fixture->scores.scanTimeStart = row["start"].toDouble();
        fixture->scores.scanTimeEnd = row["end"].toDouble();
        fixture->scores.discriminantScore = row["lda"].toDouble();
        const QJsonArray intensities = row["intensities"].toArray();
        const QJsonArray cosines = row["cosines"].toArray();
        for (int ion = 0; ion < 12; ++ion) {
            fixture->scores.integrations[ion] = intensities[ion].toDouble();
            fixture->scores.featuresArray[CosineSimToAnchor1 + ion] = cosines[ion].toDouble();
        }
        rows.push_back(&fixture->scores);
        fixtures.push_back(std::move(fixture));
    }
    const auto original = rows;
    QCOMPARE(FragmentCompetition::removeCompetingCandidates(true, &rows), eNoError);
    const QJsonArray keep = reference["keep"].toArray();
    QCOMPARE(rows.size(), keep.size());
    for (int index = 0; index < rows.size(); ++index) {
        QCOMPARE(rows.at(index), original.at(keep[index].toInt()));
    }
}

QTEST_MAIN(FragmentCompetitionTests)
#include "FragmentCompetitionTests.moc"
