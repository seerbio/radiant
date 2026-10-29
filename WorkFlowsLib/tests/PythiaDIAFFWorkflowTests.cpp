#include "PythiaDIAFFWorkflow.h"

#include "Error.h"
#include "../../AlgorithmsFFLib/tests/FragmentCompetitionTestData.h"
#include "MsReaderParquet.h"
#include "PythiaParameterReader.h"
#include "TargetDecoyCandidatePairManager.h"

#include <QtTest>
#include <QTemporaryDir>

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

    void repeatInitializationReplacesLibrary();
    void failedFileResetsStateAndRetainsLibrary();
    void featureSettersUpdateFileBaseline();
    void processFileRequiresInitialization();

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

void PythiaDIAFFWorkflowTests::repeatInitializationReplacesLibrary() {
    const QDir testFiles(qApp->applicationDirPath());
    const QString library = testFiles.filePath("FragLibReaderTests.fragLibFF");
    const QString fasta = testFiles.filePath("human_plasma_entrapment_super_trunc.fasta");
    PythiaDIAFFWorkflow workflow;
    const PythiaParameters parameters = PythiaParameterReader::genericPythiaParametersForTests();

    QCOMPARE(workflow.init(parameters, library, fasta, ""), eNoError);
    const int candidateCount = workflow.m_targetDecoyPairPntrs.size();
    const int rowCount = workflow.m_fragLibReaderRows.size();
    QVERIFY(candidateCount > 0);

    QCOMPARE(workflow.init(parameters, library, "", ""), eNoError);
    QCOMPARE(workflow.m_targetDecoyPairPntrs.size(), candidateCount);
    QCOMPARE(workflow.m_targetDecoyCandidatePairManager.targetsCount(), candidateCount);
    QCOMPARE(workflow.m_fragLibReaderRows.size(), rowCount);
    QVERIFY(workflow.m_fastaUri.isEmpty());
    QVERIFY(workflow.m_isInitialized);
}

void PythiaDIAFFWorkflowTests::failedFileResetsStateAndRetainsLibrary() {
    QTemporaryDir temporaryFiles;
    QVERIFY(temporaryFiles.isValid());
    const QString library = temporaryFiles.filePath("resident.fragLibFF");
    QVERIFY(QFile::copy(
        QDir(qApp->applicationDirPath()).filePath("FragLibReaderTests.fragLibFF"),
        library));
    PythiaParameters parameters = PythiaParameterReader::genericPythiaParametersForTests();
    parameters.alignMs1ScanTimes = true;
    PythiaDIAFFWorkflow workflow;
    QCOMPARE(workflow.init(parameters, library, "", ""), eNoError);
    const QVector<TargetDecoyCandidatePair*> libraryPointers = workflow.m_targetDecoyPairPntrs;
    QVERIFY(!libraryPointers.isEmpty());
    const QVector<MS2Ion> ions = libraryPointers.first()->ms2IonsTarget();
    QVERIFY(!ions.isEmpty());
    QVERIFY(QFile::remove(library));

    // Mimic parameters, calibration and result buffers left by a completed file.
    workflow.m_pythiaParameters.ms1ExtractionWidthPPM = 3.0;
    workflow.m_pythiaParameters.ms2ExtractionWidthPPM = 4.0;
    workflow.m_pythiaParameters.alignMs1ScanTimes = false;
    workflow.m_msCalibratomatic.setScanTimeStDev(12.0);
    workflow.m_msCalibratomatic.setIonMobilityStDev(0.05);
    workflow.m_msCalibratomatic.setMzStDevMS2(2.0);
    workflow.m_weights = {1.0f, 2.0f};
    workflow.m_targetDecoyCandidatePairsTopScores = libraryPointers;
    workflow.m_entered.insert(libraryPointers.first(), true);
    workflow.m_candidateScorePairs.push_back({CandidateScores(), CandidateScores()});
    workflow.m_resultsSummary.psmCountNeuralNet = 42;
    workflow.m_resultsSummary.psmCountLDA = 43;

    QCOMPARE(workflow.processFile(temporaryFiles.filePath("missing.prqFF")), eFileError);
    QCOMPARE(workflow.m_pythiaParameters.ms1ExtractionWidthPPM, parameters.ms1ExtractionWidthPPM);
    QCOMPARE(workflow.m_pythiaParameters.ms2ExtractionWidthPPM, parameters.ms2ExtractionWidthPPM);
    QCOMPARE(workflow.m_pythiaParameters.alignMs1ScanTimes, true);
    QCOMPARE(workflow.m_msCalibratomatic.scanTimeStDev(), -1.0f);
    QCOMPARE(workflow.m_msCalibratomatic.ionMobilityStDev(), -1.0f);
    QCOMPARE(workflow.m_msCalibratomatic.mzStDevMS2(), -1.0f);
    QVERIFY(!workflow.m_msCalibratomatic.isInitRT());
    QVERIFY(!workflow.m_msCalibratomatic.isInitIM());
    QVERIFY(!workflow.m_msCalibratomatic.isInitCalMS1());
    QVERIFY(!workflow.m_msCalibratomatic.isInitCalMS2());
    QVERIFY(workflow.m_weights.isEmpty());
    QVERIFY(workflow.m_targetDecoyCandidatePairsTopScores.isEmpty());
    QVERIFY(workflow.m_entered.isEmpty());
    QVERIFY(workflow.m_candidateScorePairs.isEmpty());
    QCOMPARE(workflow.resultsSummary().psmCountNeuralNet, -1);
    QCOMPARE(workflow.resultsSummary().psmCountLDA, -1);
    QVERIFY(!workflow.m_targetDecoyCandidatePairScoretron.isInit());
    QCOMPARE(workflow.m_targetDecoyPairPntrs, libraryPointers);

    const QVector<MS2Ion> residentIons = libraryPointers.first()->ms2IonsTarget();
    QCOMPARE(residentIons.size(), ions.size());
    for (int index = 0; index < ions.size(); ++index) {
        QCOMPARE(residentIons.at(index).mz, ions.at(index).mz);
        QCOMPARE(residentIons.at(index).intensity, ions.at(index).intensity);
    }

    // Reader-open failures must also be retryable without re-reading the library.
    QFile unsupported(temporaryFiles.filePath("invalid.unsupported"));
    QVERIFY(unsupported.open(QIODevice::WriteOnly));
    unsupported.write("invalid MS data");
    unsupported.close();
    QVERIFY(workflow.processFile(unsupported.fileName()) != eNoError);
    QVERIFY(!workflow.m_targetDecoyCandidatePairScoretron.isInit());
    QCOMPARE(workflow.m_targetDecoyPairPntrs, libraryPointers);
    QCOMPARE(workflow.processFile(temporaryFiles.filePath("missing-again.prqFF")), eFileError);
    QCOMPARE(workflow.m_targetDecoyPairPntrs, libraryPointers);
}

void PythiaDIAFFWorkflowTests::featureSettersUpdateFileBaseline() {
    PythiaDIAFFWorkflow workflow;
    const QVector<Features> calibration = {CosineSimSum100, IonMobilityDeltaAbs};
    const QVector<Features> optimization = {CosineSim100MS1, IonMobilityPdAbs};
    const QVector<Features> neuralNet = {CosineSimSpectrumOverTime};
    QCOMPARE(workflow.setCalibratomaticFeatures(calibration), eNoError);
    QCOMPARE(workflow.setPPMOptimizationFeatures(optimization), eNoError);
    QCOMPARE(workflow.setNeuralNetFeatures(neuralNet), eNoError);

    workflow.m_calibratomaticFeatures.clear();
    workflow.m_ppmOptimizationFeatures.clear();
    workflow.m_neuralNetFeatures.clear();
    workflow.resetFileState();
    QCOMPARE(workflow.m_calibratomaticFeatures, calibration);
    QCOMPARE(workflow.m_ppmOptimizationFeatures, optimization);
    QCOMPARE(workflow.m_neuralNetFeatures, neuralNet);

    const QVector<Features> replacement = {CosineSimSum100Top12};
    QCOMPARE(workflow.setNeuralNetFeatures(replacement), eNoError);
    workflow.m_neuralNetFeatures.clear();
    workflow.resetFileState();
    QCOMPARE(workflow.m_neuralNetFeatures, replacement);
}

void PythiaDIAFFWorkflowTests::processFileRequiresInitialization() {
    PythiaDIAFFWorkflow workflow;
    const QString dataFile = QDir(qApp->applicationDirPath())
        .filePath("EXP22092_2022ms0742X32_A.raw.mzML.trunc.prqFF");
    QVERIFY(workflow.processFile(dataFile) != eNoError);
    QVERIFY(workflow.m_targetDecoyPairPntrs.isEmpty());
    QCOMPARE(workflow.resultsSummary().psmCountNeuralNet, -1);
}

QTEST_MAIN(PythiaDIAFFWorkflowTests)

#include "PythiaDIAFFWorkflowTests.moc"
