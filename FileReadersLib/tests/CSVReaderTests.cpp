//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"
#include "CSVReader.h"

#include <QtTest/QtTest>

class CSVReaderTests : public QObject
{
    Q_OBJECT

public:
    CSVReaderTests() = default;
    ~CSVReaderTests() override = default;

private Q_SLOTS:
    void readWriteTest();

};


namespace TestCSVNamespace {

    const QString X = QStringLiteral("x");
    const QString Y = QStringLiteral("y");

    const QStringList keysToCheck = {
            X,
            Y
    };

}//namespace
void CSVReaderTests::readWriteTest() {

    ERR_INIT

    const QString &csvFilePath = QDir(qApp->applicationDirPath()).filePath("test.csv");

    struct TestCSV: public CSVReaderInputBase {

        int x = -1;
        int y = -1;

        QMap<QString, QVariant> map() override {

            using namespace TestCSVNamespace;

            return {
                    {X, x},
                    {Y, y}
            };
        }

        Err initFromRead(const CSVReaderInputBase &row) override {

            using namespace TestCSVNamespace;

            ERR_INIT

            const QMap<QString, QVariant> &dataMap = row.dataMap();
            const bool allKeysPresent = CSVReaderInputBase::checkIfExpectedKeysArePresent(
                    dataMap,
                    keysToCheck
            );

            if (!allKeysPresent) {
                qDebug() << dataMap;
                qDebug() << keysToCheck;
            }

            e = ErrorUtils::isTrue(allKeysPresent); ree;

            x = dataMap.value(X).toInt();
            y = dataMap.value(Y).toInt();

            ERR_RETURN
        }
    };

    TestCSV testCsv;
    testCsv.x = 666;
    testCsv.y = 777;

    QVector<TestCSV> testCSVs = {testCsv, testCsv};
    e = CSVReader::write(testCSVs, csvFilePath);
    QCOMPARE(e, eNoError);

    e = ErrorUtils::fileExists(csvFilePath);
    QCOMPARE(e, eNoError);

    testCSVs.clear();
    e = CSVReader::read(csvFilePath, &testCSVs);
    QCOMPARE(e, eNoError);

    QCOMPARE(testCSVs.size(), 2);
    QCOMPARE(testCSVs.front().x, 666);
    QCOMPARE(testCSVs.front().y, 777);

}

QTEST_MAIN(CSVReaderTests)
#include "CSVReaderTests.moc"
