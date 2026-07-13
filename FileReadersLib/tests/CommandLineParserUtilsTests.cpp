//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "CommandLineParserUtils.h"

#include <QDir>
#include <QFile>
#include <QTemporaryDir>
#include <QtTest/QtTest>

class CommandLineParserUtilsTests : public QObject
{
    Q_OBJECT

public:
    CommandLineParserUtilsTests() = default;
    ~CommandLineParserUtilsTests() override = default;

private Q_SLOTS:
    void checkFileNameExtensionTest();
    void isMassSpectrometryDataPathTest();
    void getDataFilesFromDirectoryTest();

};

void CommandLineParserUtilsTests::checkFileNameExtensionTest() {

    const QString fileName = QStringLiteral("bellatrixAndKallipe.cat");

    const bool trueTest = CommandLineParserUtils::checkFileNameExtensions(
            fileName,
            {"cat"}
            );
    QCOMPARE(trueTest, true);

    const bool falseTest = CommandLineParserUtils::checkFileNameExtensions(
            fileName,
            {"catt"}
    );
    QCOMPARE(falseTest, false);

    const QString directoryPathTrailingSlash = QStringLiteral("/tmp/example_run.d/");
    const bool directorySuffixTrailingSlashTest = CommandLineParserUtils::checkFileNameExtensions(
            directoryPathTrailingSlash,
            {"d"}
            );
    QCOMPARE(directorySuffixTrailingSlashTest, true);
}

void CommandLineParserUtilsTests::isMassSpectrometryDataPathTest() {

    const QString sidecarPathBySuffix = QStringLiteral("/tmp/example_run.d.idx");
    QCOMPARE(CommandLineParserUtils::isMassSpectrometryDataPath(sidecarPathBySuffix), true);
    QCOMPARE(CommandLineParserUtils::isMassSpectrometryDataPath(sidecarPathBySuffix + QStringLiteral("/")), true);

    const QString brukerPathBySuffix = QStringLiteral("/tmp/example_run.d");
    QCOMPARE(CommandLineParserUtils::isMassSpectrometryDataPath(brukerPathBySuffix), true);
    QCOMPARE(CommandLineParserUtils::isMassSpectrometryDataPath(brukerPathBySuffix + QStringLiteral("/")), true);

    QTemporaryDir temporaryDir;
    QVERIFY(temporaryDir.isValid());

    const QString metadataSidecarRoot = QDir(temporaryDir.path()).filePath("custom_sidecar_root");
    QVERIFY(QDir().mkpath(metadataSidecarRoot));

    QFile metadataFile(QDir(metadataSidecarRoot).filePath("metadata.json"));
    QVERIFY(metadataFile.open(QIODevice::WriteOnly | QIODevice::Text));
    metadataFile.write("{}");
    metadataFile.close();

    QCOMPARE(CommandLineParserUtils::isMassSpectrometryDataPath(metadataSidecarRoot), true);

    const QString plainDirectoryPath = QDir(temporaryDir.path()).filePath("plain_directory");
    QVERIFY(QDir().mkpath(plainDirectoryPath));
    QCOMPARE(CommandLineParserUtils::isMassSpectrometryDataPath(plainDirectoryPath), false);
}

void CommandLineParserUtilsTests::getDataFilesFromDirectoryTest() {

    ERR_INIT

    const QString &fastaFilePath = qApp->applicationDirPath();

    QStringList dataFiles;
    e = CommandLineParserUtils::getDataFilesFromDirectory(fastaFilePath, &dataFiles);
    QCOMPARE(e, eNoError);

}

QTEST_MAIN(CommandLineParserUtilsTests)
#include "CommandLineParserUtilsTests.moc"
