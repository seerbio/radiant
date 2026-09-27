#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.h"
#include "PythiaDIAFFWorkflowAlgos/PythiaDIAFFWorkflowSharedMethods.h"
#include "../../AlgorithmsFFLib/tests/FragmentCompetitionTestData.h"

#include <QtTest/QtTest>

class MsCalibratomaticSettertronTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void competitionToggle_data();
    void competitionToggle();
    void emptyCalibrationPool();
    void ms1AnchorUsesExtractedOriginMass_data();
    void ms1AnchorUsesExtractedOriginMass();
};

void MsCalibratomaticSettertronTests::competitionToggle_data() {
    QTest::addColumn<bool>("enabled");
    QTest::newRow("enabled") << true;
    QTest::newRow("disabled") << false;
}

void MsCalibratomaticSettertronTests::competitionToggle() {
    QFETCH(bool, enabled);
    CompetitionCandidateFixture weak;
    CompetitionCandidateFixture strong("PEPTIDER", {300, 400, 500, 600, 800});
    weak.scores.integrations[4] = 0;
    weak.scores.discriminantScore = 100.0;
    PythiaParameters parameters;
    parameters.competitionEnabled = enabled;
    MsCalibratomaticSettertron calibration;
    calibration.m_pythiaParameters = &parameters;
    QVector<CandidateScores*> candidates = {&strong.scores, nullptr, &weak.scores};
    QVector<CandidateScores*> selected;
    QCOMPARE(calibration.selectCalibrationCandidates(&candidates, 1000, &selected), eNoError);
    QCOMPARE(selected.size(), enabled ? 1 : 2);
    QVERIFY(!selected.contains(nullptr));
    QVERIFY(selected.contains(&strong.scores));
    QCOMPARE(selected.contains(&weak.scores), !enabled);
    QCOMPARE(candidates.size(), 2);
    QCOMPARE(candidates.front(), &weak.scores);
}

void MsCalibratomaticSettertronTests::emptyCalibrationPool() {
    PythiaParameters parameters;
    MsCalibratomaticSettertron calibration;
    calibration.m_pythiaParameters = &parameters;
    for (bool enabled : {false, true}) {
        parameters.competitionEnabled = enabled;
        QVector<CandidateScores*> candidates = {nullptr};
        QVector<CandidateScores*> selected;
        QCOMPARE(calibration.selectCalibrationCandidates(&candidates, 50, &selected), eNoError);
        QVERIFY(selected.isEmpty());
    }
}

void MsCalibratomaticSettertronTests::ms1AnchorUsesExtractedOriginMass_data() {
    QTest::addColumn<bool>("decoy");
    QTest::addColumn<int>("charge");
    QTest::newRow("target-charge2") << false << 2;
    QTest::newRow("decoy-charge2") << true << 2;
    QTest::newRow("target-charge3") << false << 3;
    QTest::newRow("decoy-charge3") << true << 3;
}

void MsCalibratomaticSettertronTests::ms1AnchorUsesExtractedOriginMass() {
    QFETCH(bool, decoy);
    QFETCH(int, charge);
    CompetitionCandidateFixture fixture;
    fixture.library.precursorCharge = charge;
    // A nonzero decoy shift is essential: a zero-shift fixture hides the bug.
    TargetDecoyCandidatePair pair(PeptideStringWithMods("PEPTIDEK"), 20.0f);
    pair.setFragLibReaderRowPntr(&fixture.library);
    QVERIFY(pair.mz(false) != pair.mz(true));
    fixture.scores.targetDecoyCandidatePair = &pair;
    fixture.scores.isDecoy = decoy;
    fixture.scores.featuresArray[Ms1MzMeanFound100] = pair.mz(false) * (1.0f + 1.25e-6f);
    fixture.scores.featuresArray[Ms1MzStDevFound100] = .0001f;
    fixture.scores.featuresArray[Ms1IntensityFound100] = 1000.0f;

    QVector<MsCalibarationReaderRow> rows;
    QCOMPARE(PythiaDIAFFWorkflowSharedMethods::buildMsCalibrationReaderRows(
                 MSLevelEnum::MS1, {&fixture.scores}, 0, &rows), eNoError);
    QCOMPARE(rows.size(), 1);
    QCOMPARE(rows.front().mzSearchedVec.front(), pair.mz(false));
    QCOMPARE(rows.front().mzFoundMeanVec.front(), fixture.scores.featuresArray[Ms1MzMeanFound100]);
    QCOMPARE(rows.front().mzFoundStDevVec.front(), .0001f);
    QCOMPARE(rows.front().intensityFoundMaxVec.front(), 1000.0f);
    QCOMPARE(fixture.scores.isDecoy, decoy);
}

QTEST_MAIN(MsCalibratomaticSettertronTests)
#include "MsCalibratomaticSettertronTests.moc"
