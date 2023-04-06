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

    const QByteArray bytesString
            = ParquetReaderInputBase::qVectorToQByteArray(vecs);

    const QVector<double> vecDecoded
            = ParquetReaderInputBase::bytesArrayToQVector<double>(bytesString);

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
            };
        }

        Err initFromRead(const ParquetReaderInputBase &row) {

            ERR_INIT

            const QMap<QString, QVariant> &dataMap = row.dataMap();

            testString = dataMap.value("testString").toString();
            testVecDouble << bytesArrayToQVector<double>(dataMap.value("testVecDouble").toByteArray());
            testVecInt = bytesArrayToQVector<int>(dataMap.value("testVecInt").toByteArray());
            testVecFloat = bytesArrayToQVector<float>(dataMap.value("testVecFloat").toByteArray());
            testVecBool = bytesArrayToQVector<bool>(dataMap.value("testVecBool").toByteArray());
            testStringList = dataMap.value("testStringList").toString().split(S_GLOBAL_SETTINGS.SEPARATOR);
            d = dataMap.value("d").toDouble();
            i = dataMap.value("i").toInt();
            f = dataMap.value("f").toFloat();
            b = dataMap.value("b").toBool();

            ERR_RETURN
        }

    };

    TestRow testRow;
    testRow.testString = QStringLiteral("Bellatrix and Kalliope are the best kittehs eva");
    testRow.testVecDouble = {666.6, 66.6, 6.6};
    testRow.testVecInt = {666, 666, 666, 6666666};
    testRow.testVecFloat = {666.6, 66.6, 6.6};;
    testRow.testVecBool = {true, false};
    testRow.testStringList = QStringList({"Chauncy", "Flops"});
    testRow.d = 666.6;
    testRow.i = 666;
    testRow.f = 666.6;
    testRow.b = false;

    QVector<TestRow> testRows(10, testRow);

    const QString &testFilePath
            = QDir(qApp->applicationDirPath()).filePath("test_write.parquet");

    ERR_INIT

    QVector<QSharedPointer<ParquetReaderInputBase>> ptrs;
    for (const TestRow &tr : testRows) {
        QSharedPointer<ParquetReaderInputBase> ptr(new TestRow(tr));
        ptrs.push_back(ptr);
    }

    ParquetReader parquetReader;
    e = parquetReader.writeDataToParquet(
            testFilePath,
            ptrs
            );
    QCOMPARE(e, eNoError);

    qDebug() << testFilePath;

    QVector<ParquetReaderInputBase> ptrsRead;
    e = parquetReader.readDataFromParquet(
            testFilePath,
            &ptrsRead
            );
    QCOMPARE(e, eNoError);

    const ParquetReaderInputBase &r = ptrsRead.first();
    TestRow readRow;
    e = readRow.initFromRead(r);
    QCOMPARE(e, eNoError);

    QCOMPARE(readRow.testString, testRow.testString);
    QCOMPARE(readRow.testStringList, testRow.testStringList);
    QCOMPARE(readRow.testVecDouble, testRow.testVecDouble);
    QCOMPARE(readRow.testVecInt, testRow.testVecInt);
    QCOMPARE(readRow.testVecBool, testRow.testVecBool);
    QCOMPARE(readRow.testVecFloat, testRow.testVecFloat);
    QCOMPARE(readRow.d, testRow.d);
    QCOMPARE(readRow.i, testRow.i);
    QCOMPARE(readRow.f, testRow.f);
    QCOMPARE(readRow.b, testRow.b);

}


QTEST_MAIN(ParquetReaderTests)
#include "ParquetReaderTests.moc"
