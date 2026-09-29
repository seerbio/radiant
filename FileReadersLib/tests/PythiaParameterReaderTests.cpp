//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "PythiaParameterReader.h"

#include <QTemporaryFile>
#include <QtTest/QtTest>
#include <QVariant>

class PythiaParameterReaderTests : public QObject
{
    Q_OBJECT

public:
    PythiaParameterReaderTests() = default;
    ~PythiaParameterReaderTests() override = default;

private Q_SLOTS:

    static void readFileTest();
    static void competitionEnabled_data();
    static void competitionEnabled();
    static void sensitivityOptions();
    static void rejectsInvalidSensitivityOptions();
    static void candidateBundleOption();
    static void candidateOnlyOption();

};

void PythiaParameterReaderTests::readFileTest() {

    ERR_INIT

    const QString filePath
            = QDir(qApp->applicationDirPath()).filePath("test_params_wide_window.pythiaConfig");

    PythiaParameterReader reader;

    PythiaParameters pythiaParameters;
    PythiaParameterReader::buildPythiaParameters(filePath, &pythiaParameters);
    QCOMPARE(e, eNoError);

    QCOMPARE(pythiaParameters.threadCount, 16);
    QCOMPARE(pythiaParameters.verbosity, 1);
    QCOMPARE(pythiaParameters.writeRadiantDIA, true);
    QCOMPARE(pythiaParameters.chargeStateMin, 1);
    QCOMPARE(pythiaParameters.chargeStateMax, 4);
    QCOMPARE(pythiaParameters.mzMinMS2, 200.0);
    QCOMPARE(pythiaParameters.mzMaxMS2, 1500.0);
    QCOMPARE(pythiaParameters.peptideLengthMin, 8);
    QCOMPARE(pythiaParameters.peptideLengthMax, 31);
    QCOMPARE(pythiaParameters.trancheSizeMax, 1e4);
    QCOMPARE(pythiaParameters.useAlternativeDecoys, false);
    QCOMPARE(pythiaParameters.precursorExtractionWindowThomsons, 0.5);
    QCOMPARE(pythiaParameters.ms1ExtractionWidthPPM, 20.0);
    QCOMPARE(pythiaParameters.filterLengthIntegration, 6);
    QCOMPARE(pythiaParameters.filterLengthMS2, 4);
    QCOMPARE(pythiaParameters.competitionEnabled, true);
    QCOMPARE(pythiaParameters.ionsSharedToReject, 2);
    QCOMPARE(pythiaParameters.ms2ExtractionWidthPPM, 21.0);
    QCOMPARE(pythiaParameters.minMs2FragCount, 3);
    QCOMPARE(pythiaParameters.scanTimeWindowStDevs, 4);
    QCOMPARE(pythiaParameters.subtractShadows, false);
    QCOMPARE(pythiaParameters.smoothCountMS2, 2);
    QCOMPARE(pythiaParameters.stopThresholdFractionMS2, 0.666f);
    QCOMPARE(pythiaParameters.timsMainCandidateBudgetPerTargetKey, 4000);
    QCOMPARE(pythiaParameters.timsStratifyCandidateBudget, true);
    QCOMPARE(pythiaParameters.timsHighEvidenceMinCosineSimSum100, 3.8f);
    QCOMPARE(pythiaParameters.timsHighEvidenceMinCosineSimSpectrumOverTimeCubed, 0.2f);
    QCOMPARE(pythiaParameters.timsHighEvidenceMaxScanTimeDeltaAbs, 90.0f);
    QCOMPARE(pythiaParameters.timsHighEvidenceFilterSweep, true);
    QCOMPARE(pythiaParameters.timsHighEvidenceFilterEnabled, false);
    QCOMPARE(pythiaParameters.timsSecondStageCandidateRowLimit, 24000);
    QCOMPARE(pythiaParameters.timsSecondStageUniquePrecursorLimit, 12000);
    QCOMPARE(pythiaParameters.timsLocalFdrRtBinSeconds, 150.0);
    QCOMPARE(pythiaParameters.percentFDR, 2.0);
    QCOMPARE(pythiaParameters.reportDecoys, true);
    QCOMPARE(pythiaParameters.filterLength, 4);
    QCOMPARE(pythiaParameters.sigma, 1.1);
    QCOMPARE(pythiaParameters.signalToNoiseRatio, 2.1);
    QCOMPARE(pythiaParameters.smoothCount, 3);
    QCOMPARE(pythiaParameters.minScanCount, 4);
    QCOMPARE(pythiaParameters.skipScanCount, 4);

    pythiaParameters.print();

}

void PythiaParameterReaderTests::competitionEnabled_data() {
    QTest::addColumn<QByteArray>("config");
    QTest::addColumn<bool>("expected");
    QTest::newRow("missing-section-defaults-on")
        << QByteArray("[General]\nthreadCount = 1\n") << true;
    QTest::newRow("missing-option-defaults-on")
        << QByteArray("[MS2Params]\nionsSharedToReject = 7\n") << true;
    QTest::newRow("explicitly-enabled")
        << QByteArray("[MS2Params]\ncompetitionEnabled = true\nionsSharedToReject = 7\n") << true;
    QTest::newRow("explicitly-disabled")
        << QByteArray("[MS2Params]\ncompetitionEnabled = false\nionsSharedToReject = 7\n") << false;
}

void PythiaParameterReaderTests::competitionEnabled() {
    QFETCH(QByteArray, config);
    QFETCH(bool, expected);

    PythiaParameters parameters;
    QVERIFY(parameters.competitionEnabled);
    // Loading a new config with an omitted option must restore the default.
    parameters.competitionEnabled = false;

    QTemporaryFile file;
    QVERIFY(file.open());
    QCOMPARE(file.write(config), static_cast<qint64>(config.size()));
    QVERIFY(file.flush());

    const Err error = PythiaParameterReader::buildPythiaParameters(
        file.fileName(), &parameters);
    QCOMPARE(error, eNoError);
    QCOMPARE(parameters.competitionEnabled, expected);
    if (config.contains("ionsSharedToReject")) {
        QCOMPARE(parameters.ionsSharedToReject, 7);
    }
}


void PythiaParameterReaderTests::sensitivityOptions() {
    const QByteArray config =
        "[MS1Params]\nalignMs1ScanTimes = true\n"
        "[MS2Params]\nmainMinSimultaneousFragments = 3\nionsSharedToReject = 3\n"
        "postNeuralNetSharedFragments = 2\n"
        "[NeuralNetParams]\nneuralNetEnsembleSize = 4\n"
        "neuralNetGroupPeptideFamilies = true\nneuralNetLogIntensities = true\n"
        "neuralNetShuffleEachEpoch = true\n";
    QTemporaryFile file;
    QVERIFY(file.open());
    QCOMPARE(file.write(config), static_cast<qint64>(config.size()));
    QVERIFY(file.flush());
    PythiaParameters parameters;
    QCOMPARE(PythiaParameterReader::buildPythiaParameters(file.fileName(), &parameters), eNoError);
    QVERIFY(parameters.alignMs1ScanTimes);
    QCOMPARE(parameters.mainMinSimultaneousFragments, 3);
    QCOMPARE(parameters.postNeuralNetSharedFragments, 2);
    QCOMPARE(parameters.neuralNetEnsembleSize, 4);
    QVERIFY(parameters.neuralNetGroupPeptideFamilies);
    QVERIFY(parameters.neuralNetLogIntensities);
    QVERIFY(parameters.neuralNetShuffleEachEpoch);
    QVERIFY(file.resize(0));
    QVERIFY(file.seek(0));
    QVERIFY(file.write("[General]\nthreadCount = 1\n") > 0);
    QVERIFY(file.flush());
    QCOMPARE(PythiaParameterReader::buildPythiaParameters(file.fileName(), &parameters), eNoError);
    QVERIFY(!parameters.alignMs1ScanTimes);
    QCOMPARE(parameters.mainMinSimultaneousFragments, 4);
    QCOMPARE(parameters.postNeuralNetSharedFragments, 0);
    QCOMPARE(parameters.neuralNetEnsembleSize, 1);
    QVERIFY(!parameters.neuralNetGroupPeptideFamilies);
    QVERIFY(!parameters.neuralNetLogIntensities);
    QVERIFY(!parameters.neuralNetShuffleEachEpoch);
}

void PythiaParameterReaderTests::rejectsInvalidSensitivityOptions() {
    auto parameters = PythiaParameterReader::genericPythiaParametersForTests();
    QVERIFY(parameters.isValid());
    parameters.neuralNetEnsembleSize = 0;
    QVERIFY(!parameters.isValid());
    parameters.neuralNetEnsembleSize = 4;
    parameters.mainMinSimultaneousFragments = 2;
    QVERIFY(!parameters.isValid());
    parameters.mainMinSimultaneousFragments = 3;
    parameters.ionsSharedToReject = 1;
    QVERIFY(!parameters.isValid());
    parameters.ionsSharedToReject = 3;
    QVERIFY(parameters.isValid());
    for (const int value : {-1, 1, 13}) {
        parameters.postNeuralNetSharedFragments = value;
        QVERIFY(!parameters.isValid());
    }
    for (const int value : {0, 2, 12}) {
        parameters.postNeuralNetSharedFragments = value;
        QVERIFY(parameters.isValid());
    }
}

void PythiaParameterReaderTests::candidateBundleOption() {
    QTemporaryFile file;
    QVERIFY(file.open());
    PythiaParameters parameters;
    for (const QByteArray value : {QByteArray("0"), QByteArray("2"), QByteArray("200000"), QByteArray("1000000")}) {
        QVERIFY(file.resize(0)); QVERIFY(file.seek(0));
        const QByteArray config = "[NeuralNetParams]\ncandidateBundleLimit = " + value + "\n";
        QCOMPARE(file.write(config), qint64(config.size())); QVERIFY(file.flush());
        QCOMPARE(PythiaParameterReader::buildPythiaParameters(file.fileName(), &parameters), eNoError);
        QCOMPARE(parameters.candidateBundleLimit, value.toInt());
    }
    for (const QByteArray value : {QByteArray("-1"), QByteArray("1"), QByteArray("1000001"), QByteArray("true"), QByteArray("2.5"), QByteArray("\"200000\"")}) {
        QVERIFY(file.resize(0)); QVERIFY(file.seek(0));
        const QByteArray config = "[NeuralNetParams]\ncandidateBundleLimit = " + value + "\n";
        QCOMPARE(file.write(config), qint64(config.size())); QVERIFY(file.flush());
        QVERIFY(PythiaParameterReader::buildPythiaParameters(file.fileName(), &parameters) != eNoError);
    }
    QVERIFY(file.resize(0)); QVERIFY(file.seek(0));
    QVERIFY(file.write("[General]\nthreadCount = 1\n") > 0); QVERIFY(file.flush());
    parameters.candidateBundleLimit = 200000;
    QCOMPARE(PythiaParameterReader::buildPythiaParameters(file.fileName(), &parameters), eNoError);
    QCOMPARE(parameters.candidateBundleLimit, 0);
}

void PythiaParameterReaderTests::candidateOnlyOption() {
    QTemporaryFile file;
    QVERIFY(file.open());
    PythiaParameters parameters;
    const auto read = [&](const QByteArray &configuration) {
        file.resize(0);
        file.seek(0);
        file.write(configuration);
        file.flush();
        return PythiaParameterReader::buildPythiaParameters(file.fileName(), &parameters);
    };
    QCOMPARE(read("[NeuralNetParams]\ncandidateBundleLimit = 200000\ncandidateBundleOnly = true\n"),
             eNoError);
    QVERIFY(parameters.candidateBundleOnly);
    QCOMPARE(read("[NeuralNetParams]\ncandidateBundleOnly = true\n"), eValueError);
    for (const QByteArray value : {QByteArray("1"), QByteArray("\"true\""), QByteArray("1.0")}) {
        QVERIFY(read("[NeuralNetParams]\ncandidateBundleLimit = 200000\ncandidateBundleOnly = "
                     + value + "\n") != eNoError);
    }
    QCOMPARE(read("[General]\nthreadCount = 1\n"), eNoError);
    QVERIFY(!parameters.candidateBundleOnly);
    parameters = PythiaParameterReader::genericPythiaParametersForTests();
    parameters.candidateBundleOnly = true;
    parameters.candidateBundleLimit = 0;
    QVERIFY(!parameters.isValid());
    parameters.candidateBundleLimit = 200000;
    QVERIFY(parameters.isValid());
}

QTEST_MAIN(PythiaParameterReaderTests)
#include "PythiaParameterReaderTests.moc"
