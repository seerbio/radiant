#include "CommandLineParser.h"

#include <QDir>
#include <QFile>
#include <QTemporaryDir>
#include <QtTest/QtTest>

namespace {
QStringList sharedArguments() {
    return {
        QStringLiteral("RadiantDIA"),
        QStringLiteral("library.tsv"),
        QStringLiteral("background.fasta"),
        QStringLiteral("parameters.radiantConfig")
    };
}
}

class RadiantCommandLineParserTests : public QObject
{
    Q_OBJECT

private Q_SLOTS:
    void singleFileArgumentsRemainCompatible();
    void acceptsAllDataPathsInInputOrder_data();
    void acceptsAllDataPathsInInputOrder();
    void rejectsInvalidSharedArguments_data();
    void rejectsInvalidSharedArguments();
    void perFileLogsFlagAppliesToBatch_data();
    void perFileLogsFlagAppliesToBatch();
    void perFileLogsOptionIsAdvertised();
    void defersIndividualDataPathValidation();
    void sharedOptionsApplyToEntireBatch();
    void endOfOptionsAllowsHyphenLeadingInput();
    void existingOutputFolderHasTrailingSlash_data();
    void existingOutputFolderHasTrailingSlash();
    void rejectsUncreatableOutputFolder();
    void rejectsEmptyOutputFolder();
};

void RadiantCommandLineParserTests::singleFileArgumentsRemainCompatible() {
    const QString input = QStringLiteral("sample.mzML");
    QStringList args = sharedArguments();
    args.append(input);

    CommandLineParser parser;
    QVERIFY(parser.validateArguments(args));
    const auto &parameters = parser.getCliParams();
    QCOMPARE(parameters.fragLibFilePath, QStringLiteral("library.tsv"));
    QCOMPARE(parameters.fastaFilePath, QStringLiteral("background.fasta"));
    QCOMPARE(parameters.pythiaParametersFilePath, QStringLiteral("parameters.radiantConfig"));
    QCOMPARE(parameters.msDataFiles, QStringList{input});
    QVERIFY(parameters.outputFolderPath.isEmpty());
    QVERIFY(parameters.thermoRawReaderPath.isEmpty());
    QVERIFY(!parameters.perFileLogs);
}

void RadiantCommandLineParserTests::acceptsAllDataPathsInInputOrder_data() {
    QTest::addColumn<int>("inputCount");
    QTest::newRow("two-files") << 2;
    QTest::newRow("seventeen-files") << 17;
    QTest::newRow("sixty-four-files") << 64;
}

void RadiantCommandLineParserTests::acceptsAllDataPathsInInputOrder() {
    QFETCH(int, inputCount);
    QStringList inputs;
    for (int index = 0; index < inputCount; ++index) {
        // Reverse names catch accidental sorting as well as truncation.
        inputs.append(QStringLiteral("sample folder/run_%1.mzML").arg(inputCount - index));
    }
    QStringList args = sharedArguments();
    args.append(inputs);

    CommandLineParser parser;
    QVERIFY(parser.validateArguments(args));
    QCOMPARE(parser.getCliParams().msDataFiles, inputs);
}

void RadiantCommandLineParserTests::rejectsInvalidSharedArguments_data() {
    QTest::addColumn<QStringList>("arguments");
    QTest::newRow("missing-data") << sharedArguments();

    QStringList missingConfig = sharedArguments();
    missingConfig.removeLast();
    missingConfig.append(QStringLiteral("sample.mzML"));
    QTest::newRow("missing-config") << missingConfig;

    QStringList arguments = sharedArguments();
    arguments.append(QStringLiteral("sample.mzML"));
    QStringList wrongLibrary = arguments;
    wrongLibrary[1] = QStringLiteral("library.txt");
    QTest::newRow("wrong-library-extension") << wrongLibrary;

    QStringList wrongFasta = arguments;
    wrongFasta[2] = QStringLiteral("background.txt");
    QTest::newRow("wrong-fasta-extension") << wrongFasta;

    QStringList wrongConfig = arguments;
    wrongConfig[3] = QStringLiteral("parameters.txt");
    QTest::newRow("wrong-config-extension") << wrongConfig;
}

void RadiantCommandLineParserTests::rejectsInvalidSharedArguments() {
    QFETCH(QStringList, arguments);
    CommandLineParser parser;
    QVERIFY(!parser.validateArguments(arguments));
}

void RadiantCommandLineParserTests::perFileLogsFlagAppliesToBatch_data() {
    QTest::addColumn<bool>("flagBeforeInputs");
    QTest::newRow("before-inputs") << true;
    QTest::newRow("after-inputs") << false;
}

void RadiantCommandLineParserTests::perFileLogsFlagAppliesToBatch() {
    QFETCH(bool, flagBeforeInputs);
    const QStringList inputs = {
        QStringLiteral("first.mzML"),
        QStringLiteral("second.raw")
    };
    QStringList args = sharedArguments();
    if (flagBeforeInputs) {
        args.append(QStringLiteral("--per-file-logs"));
    }
    args.append(inputs);
    if (!flagBeforeInputs) {
        args.append(QStringLiteral("--per-file-logs"));
    }

    CommandLineParser parser;
    QVERIFY(parser.validateArguments(args));
    QVERIFY(parser.getCliParams().perFileLogs);
    QCOMPARE(parser.getCliParams().msDataFiles, inputs);
}

void RadiantCommandLineParserTests::perFileLogsOptionIsAdvertised() {
    CommandLineParser parser;
    QVERIFY(parser.helpText().contains(QStringLiteral("--per-file-logs")));
}

void RadiantCommandLineParserTests::defersIndividualDataPathValidation() {
    const QStringList inputs = {
        QStringLiteral("unsupported.txt"),
        QStringLiteral("missing.mzML"),
        QStringLiteral("thermo.RAW"),
        QStringLiteral("bruker.d/"),
        QStringLiteral("missing.mzML")
    };
    QStringList args = sharedArguments();
    args.append(inputs);

    CommandLineParser parser;
    QVERIFY(parser.validateArguments(args));
    QCOMPARE(parser.getCliParams().msDataFiles, inputs);
}

void RadiantCommandLineParserTests::sharedOptionsApplyToEntireBatch() {
    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());
    const QString outputFolder = temporaryDir.filePath(QStringLiteral("output folder/nested"));
    const QString readerPath = QStringLiteral("/opt/reader folder/RadiantThermoRawReader");
    const QStringList inputs = {
        QStringLiteral("first.raw"),
        QStringLiteral("second.mzML")
    };
    QStringList args = sharedArguments();
    args.append(inputs[0]);
    args.append(QStringLiteral("--thermo-raw-reader"));
    args.append(readerPath);
    args.append(inputs[1]);
    args.append(QStringLiteral("--output-folder"));
    args.append(outputFolder);

    CommandLineParser parser;
    QVERIFY(parser.validateArguments(args));
    const auto &parameters = parser.getCliParams();
    QCOMPARE(parameters.msDataFiles, inputs);
    QCOMPARE(parameters.thermoRawReaderPath, readerPath);
    QCOMPARE(parameters.outputFolderPath, outputFolder + QStringLiteral("/"));
    QVERIFY(QDir(outputFolder).exists());
}

void RadiantCommandLineParserTests::endOfOptionsAllowsHyphenLeadingInput() {
    const QString input = QStringLiteral("-sample.mzML");
    QStringList args = sharedArguments();
    args.append(QStringLiteral("--"));
    args.append(input);

    CommandLineParser parser;
    QVERIFY(parser.validateArguments(args));
    QCOMPARE(parser.getCliParams().msDataFiles, QStringList{input});
}

void RadiantCommandLineParserTests::existingOutputFolderHasTrailingSlash_data() {
    QTest::addColumn<bool>("hasTrailingSlash");
    QTest::newRow("without-slash") << false;
    QTest::newRow("with-slash") << true;
}

void RadiantCommandLineParserTests::existingOutputFolderHasTrailingSlash() {
    QFETCH(bool, hasTrailingSlash);
    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());
    const QString outputFolder = temporaryDir.path();
    QStringList args = sharedArguments();
    args.append(QStringLiteral("sample.mzML"));
    args.append(QStringLiteral("--output-folder"));
    args.append(outputFolder + (hasTrailingSlash ? QStringLiteral("/") : QString{}));

    CommandLineParser parser;
    QVERIFY(parser.validateArguments(args));
    QCOMPARE(parser.getCliParams().outputFolderPath, outputFolder + QStringLiteral("/"));
}

void RadiantCommandLineParserTests::rejectsUncreatableOutputFolder() {
    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());
    const QString blockerPath = temporaryDir.filePath(QStringLiteral("ordinary-file"));
    QFile blocker(blockerPath);
    QVERIFY(blocker.open(QIODevice::WriteOnly));
    blocker.close();
    QStringList args = sharedArguments();
    args.append(QStringLiteral("sample.mzML"));
    args.append(QStringLiteral("--output-folder"));
    args.append(blockerPath + QStringLiteral("/results"));

    CommandLineParser parser;
    QVERIFY(!parser.validateArguments(args));
}

void RadiantCommandLineParserTests::rejectsEmptyOutputFolder() {
    QStringList args = sharedArguments();
    args.append(QStringLiteral("sample.mzML"));
    args.append(QStringLiteral("--output-folder"));
    args.append(QString{});

    CommandLineParser parser;
    QVERIFY(!parser.validateArguments(args));
}

QTEST_GUILESS_MAIN(RadiantCommandLineParserTests)
#include "CommandLineParserTests.moc"
