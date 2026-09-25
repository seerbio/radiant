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

    void competitionToggle_data();
    void competitionToggle();


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

QTEST_MAIN(PythiaDIAFFWorkflowTests)

#include "PythiaDIAFFWorkflowTests.moc"
