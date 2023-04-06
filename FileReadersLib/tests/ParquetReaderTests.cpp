//
// Created by anichols on 11/07/2021.
//

#include "ErrorUtils.h"

#include "ParquetReader.h"

#include <QtTest/QtTest>

class ParquetReaderTests : public QObject
{
    Q_OBJECT

public:
    ParquetReaderTests() = default;
    ~ParquetReaderTests() override = default;

private Q_SLOTS:

    void qVectorToStdStringAndbytesStdStringToQVectorTest();
    void readWriteCombinedTest();

};

void ParquetReaderTests::qVectorToStdStringAndbytesStdStringToQVectorTest() {

    const QVector<double> vecs = {666.6, 333.3, 666.6};

    const std::string bytesStdString
            = ParquetReaderInputBase::qVectorToQByteArray(vecs).toStdString();

    const QVector<double> vecDecoded
            = ParquetReaderInputBase::bytesStdStringToQVector<double>(bytesStdString);

    QCOMPARE(vecDecoded.size(), vecs.size());
    QCOMPARE(vecDecoded, vecs);
}

void ParquetReaderTests::readWriteCombinedTest() {

    struct TestRow : public ParquetReaderInputBase {
        QString testString;
        QVector<double> testVecDouble;
        QVector<int> testVecInt;
        QVector<float> testVecFloat;
        QVector<bool> testVecBool;
        QStringList testStringList;
        double d;
        int i;
        float f;
        bool b;
        QString s;

        QMap<QString, QVariant> map() override {

            return {
                {"testString", QVariant(testString)},
                {"testVecDouble", QVariant(qVectorToQByteArray(testVecDouble))},
                {"testVecInt", QVariant(qVectorToQByteArray(testVecInt))},
                {"testVecFloat", QVariant(qVectorToQByteArray(testVecFloat))},
                {"testVecBool", QVariant(qVectorToQByteArray(testVecBool))},
                {"testStringList", QVariant(joinQStringList(testStringList))},
                {"d", QVariant(d)},
                {"i", QVariant(i)},
                {"f", QVariant(f)},
                {"b", QVariant(b)},
                {"s", QVariant(s)}
            };
        }

    };

    TestRow testRow;
    testRow.testString = QStringLiteral("Bellatrix and Kalliope are the best kittehs eva");
    testRow.testVecDouble = {666.6, 66.6, 6.6};
    testRow.testVecInt = {666, 666, 666};
    testRow.testVecFloat = {666.6, 66.6, 6.6};;
    testRow.testVecBool = {true, false};
    testRow.testStringList = QStringList({"Chauncy", "Flops"});
    testRow.d = 666.6;
    testRow.i = 666;
    testRow.f = 666.6;
    testRow.b = false;
    testRow.s = QStringLiteral("RIP Thor");

    QVector<TestRow> testRows(10, testRow);

    const QString &testFilePath
            = QDir(qApp->applicationDirPath()).filePath("test_write.parquet");

    ERR_INIT

    ParquetReader parquetReader;
    e = parquetReader.writeDataToParquet(testFilePath, testRows);
    QCOMPARE(e, eNoError);

}


QTEST_MAIN(ParquetReaderTests)
#include "ParquetReaderTests.moc"
