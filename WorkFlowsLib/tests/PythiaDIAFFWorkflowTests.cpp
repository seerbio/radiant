#include "PythiaDIAFFWorkflow.h"

#include "Error.h"
#include "../../AlgorithmsFFLib/tests/FragmentCompetitionTestData.h"
#include "MsReaderParquet.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

#include <QtTest>

class PythiaDIAFFWorkflowTests : public QObject
{
    Q_OBJECT

public:

    PythiaDIAFFWorkflowTests() = default;
    ~PythiaDIAFFWorkflowTests() override = default;

private slots:

    void initTest();

    void buildUniqueInfoScanKeyVsTargetDecoyCandidatePointersTest();

    void competitionRejectedRowsAreNotRestored_data();
    void competitionRejectedRowsAreNotRestored();
    void competitionToggle_data();
    void competitionToggle();
    void candidatePairCompletionHandlesEmptySelection();


};

void PythiaDIAFFWorkflowTests::initTest() {
    ERR_INIT

    const QString &testFragLibFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibFF");

    const QString &testFastaFilePath
            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    PythiaDIAFFWorkflow pythiaDiaffWorkflow;
    e = pythiaDiaffWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFragLibFilePath,
            testFastaFilePath,
            ""
            );
    QCOMPARE(e, eNoError);

    e = pythiaDiaffWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            "kalliope.fragLibDF",
            testFastaFilePath,
            ""
    );
    QCOMPARE(e, eFileError);

    e = pythiaDiaffWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFragLibFilePath,
            "bellatrix.fasta",
            ""
    );
    QCOMPARE(e, eFileError);

}

void PythiaDIAFFWorkflowTests::buildUniqueInfoScanKeyVsTargetDecoyCandidatePointersTest() {

    ERR_INIT

    const QString &testFragLibFilePath
            = QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibFF");

    const QString &testFastaFilePath
            = QDir(qApp->applicationDirPath()).filePath("human_plasma_entrapment_super_trunc.fasta");

    const QString &testMsFilePath
            = QDir(qApp->applicationDirPath()).filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");

    MsReaderParquet msReaderParquet;
    e = msReaderParquet.openFile(testMsFilePath);
    QCOMPARE(e, eNoError);

    const QVector<MsScanInfo> uniqueMsScanInfos = msReaderParquet.getUniqueTandemMsScanInfos();

    PythiaDIAFFWorkflow pythiaDiaffWorkflow;
    e = pythiaDiaffWorkflow.init(
            PythiaParameterReader::genericPythiaParametersForTests(),
            testFragLibFilePath,
            testFastaFilePath,
            ""
            );
    QCOMPARE(e, eNoError);

}


void PythiaDIAFFWorkflowTests::competitionRejectedRowsAreNotRestored_data() {
    QTest::addColumn<bool>("rejectedIsDecoy");
    QTest::newRow("rejected-decoy-stays-excluded") << true;
    QTest::newRow("rejected-target-stays-excluded") << false;
}

void PythiaDIAFFWorkflowTests::competitionRejectedRowsAreNotRestored() {
    QFETCH(bool, rejectedIsDecoy);

    PythiaDIAFFWorkflow workflow;
    workflow.m_candidateScorePairs.resize(3);
    CompetitionCandidateFixture loser;
    CompetitionCandidateFixture winner("PEPTIDER", {300, 400, 500, 600, 800});
    CompetitionCandidateFixture unrelated("PEPTIDEA", {810, 820, 830, 840, 850});
    CompetitionCandidateFixture *fixtures[] = {&loser, &winner, &unrelated};
    QVector<CandidateScores*> availableRows;
    for (int pairIndex = 0; pairIndex < workflow.m_candidateScorePairs.size(); ++pairIndex) {
        auto &scorePair = workflow.m_candidateScorePairs[pairIndex];
        scorePair.first = fixtures[pairIndex]->scores;
        scorePair.second = fixtures[pairIndex]->scores;
        scorePair.first.scanTime = 200.0f + pairIndex * 100;
        scorePair.second.scanTime = 500.0f + pairIndex * 100;
        scorePair.first.discriminantScore = 10.0 - pairIndex * 2;
        scorePair.second.discriminantScore = 9.0 - pairIndex * 2;
        scorePair.second.isDecoy = true;
        availableRows.push_back(&scorePair.first);
        availableRows.push_back(&scorePair.second);
    }

    auto &competingPair = workflow.m_candidateScorePairs[0];
    CandidateScores *rejected = rejectedIsDecoy ? &competingPair.second : &competingPair.first;
    CandidateScores *retained = rejectedIsDecoy ? &competingPair.first : &competingPair.second;
    rejected->discriminantScore = 1.0;
    retained->discriminantScore = 10.0;
    rejected->scanTime = 60.0f;
    rejected->integrations[4] = 0;
    workflow.m_candidateScorePairs[1].first.scanTime = 60.0f;

    const Err error = workflow.applyFragmentCompetition(&availableRows);
    QCOMPARE(error, eNoError);
    QCOMPARE(availableRows.size(), 5);
    QVERIFY(!availableRows.contains(rejected));

    auto &selectedPair = workflow.m_candidateScorePairs[1];
    auto &unselectedPair = workflow.m_candidateScorePairs[2];
    QVector<CandidateScores*> trainingRows = {retained, &selectedPair.first};
    QVector<CandidateScores*> inferenceRows = trainingRows;

    workflow.completeCandidateRowsToTargetDecoyPairs(availableRows, &trainingRows);
    workflow.completeCandidateRowsToTargetDecoyPairs(availableRows, &inferenceRows);

    for (const QVector<CandidateScores*> &selectedRows : {trainingRows, inferenceRows}) {
        QCOMPARE(selectedRows.size(), 3);
        QVERIFY(selectedRows.contains(retained));
        QVERIFY(selectedRows.contains(&selectedPair.first));
        QVERIFY(selectedRows.contains(&selectedPair.second));
        QVERIFY(!selectedRows.contains(rejected));
        QVERIFY(!selectedRows.contains(&unselectedPair.first));
        QVERIFY(!selectedRows.contains(&unselectedPair.second));
        for (CandidateScores *selected : selectedRows) {
            QVERIFY(availableRows.contains(selected));
        }
    }

    workflow.completeCandidateRowsToTargetDecoyPairs(availableRows, &trainingRows);
    QCOMPARE(trainingRows.size(), 3);
    QCOMPARE(workflow.m_candidateScorePairs.size(), 3);
    QCOMPARE(rejected->discriminantScore, 1.0);
}

void PythiaDIAFFWorkflowTests::competitionToggle_data() {
    QTest::addColumn<bool>("enabled");
    QTest::newRow("enabled") << true;
    QTest::newRow("disabled") << false;
}

void PythiaDIAFFWorkflowTests::competitionToggle() {
    QFETCH(bool, enabled);
    CompetitionCandidateFixture weak;
    CompetitionCandidateFixture strong("PEPTIDER", {300, 400, 500, 600, 800});
    weak.scores.integrations[4] = 0;
    PythiaDIAFFWorkflow workflow;
    workflow.m_pythiaParameters.competitionEnabled = enabled;
    QVector<CandidateScores*> rows = {&weak.scores, &strong.scores};
    QCOMPARE(workflow.applyFragmentCompetition(&rows), eNoError);
    QCOMPARE(rows.size(), enabled ? 1 : 2);
    QCOMPARE(rows.contains(&weak.scores), !enabled);
    QVERIFY(rows.contains(&strong.scores));
}

void PythiaDIAFFWorkflowTests::candidatePairCompletionHandlesEmptySelection() {
    PythiaDIAFFWorkflow workflow;
    QVector<CandidateScores*> availableRows;
    QVector<CandidateScores*> selectedRows;

    workflow.completeCandidateRowsToTargetDecoyPairs(availableRows, nullptr);
    workflow.completeCandidateRowsToTargetDecoyPairs(availableRows, &selectedRows);

    QVERIFY(selectedRows.isEmpty());
}

QTEST_MAIN(PythiaDIAFFWorkflowTests)

#include "PythiaDIAFFWorkflowTests.moc"
