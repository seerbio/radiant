#include "PythiaDIAFFWorkflowAlgos/MsCalibratomaticSettertron.h"
#include "../../AlgorithmsFFLib/tests/FragmentCompetitionTestData.h"

#include <QtTest/QtTest>

class MsCalibratomaticSettertronTests : public QObject {
    Q_OBJECT
private Q_SLOTS:
    void competitionToggle_data();
    void competitionToggle();
    void emptyCalibrationPool();
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

QTEST_MAIN(MsCalibratomaticSettertronTests)
#include "MsCalibratomaticSettertronTests.moc"
