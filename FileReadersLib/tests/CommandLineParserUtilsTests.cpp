//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "CommandLineParserUtils.h"

#include <QtTest/QtTest>

class CommandLineParserUtilsTests : public QObject
{
    Q_OBJECT

public:
    CommandLineParserUtilsTests() = default;
    ~CommandLineParserUtilsTests() override = default;

private Q_SLOTS:
    void checkFileNameExtensionTest();
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
