//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "JsonParametersReader.h"

#include <QtTest/QtTest>
#include <QVariant>

class JsonParameterReaderTests : public QObject
{
    Q_OBJECT

public:
    JsonParameterReaderTests() = default;
    ~JsonParameterReaderTests() = default;

private Q_SLOTS:

    void readFileTest();

};

void JsonParameterReaderTests::readFileTest() {

    ERR_INIT

    const QString &jsonFilePath
            = QDir(qApp->applicationDirPath()).filePath("DigestParamsTest.dgst");


    JsonParametersReader reader;
    e = reader.readFile(jsonFilePath);

    QCOMPARE(e, eNoError);

   const QMap<QString, QVariant> &data = reader.jsonContents();

    QVERIFY(data.contains(QStringLiteral("allowedMissedCleavages")));
    QCOMPARE(data.value("allowedMissedCleavages"), QVariant(1));

}


QTEST_MAIN(JsonParameterReaderTests)
#include "JsonParameterReaderTests.moc"
