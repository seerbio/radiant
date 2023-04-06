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
        QStringList ionLabels;
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
                {"ionLabels", QVariant(joinQStringList(ionLabels))},
                {"d", QVariant(d)},
                {"i", QVariant(i)},
                {"f", QVariant(f)},
                {"b", QVariant(b)},
                {"s", QVariant(s)}
            };
        }



    };

}


QTEST_MAIN(ParquetReaderTests)
#include "ParquetReaderTests.moc"
