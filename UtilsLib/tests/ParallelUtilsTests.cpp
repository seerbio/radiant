//
// Created by anichols on 11/07/2021.
//

#include "ParallelUtils.h"

#include <QtTest/QtTest>

class ParallelUtilsTests : public QObject
{
    Q_OBJECT

public:
    ParallelUtilsTests();
    ~ParallelUtilsTests() override = default;

private Q_SLOTS:
    void trancheVectorForParallelizationTest();
    void trancheVectorForParallelizationInOrderTest();


private:

    const QVector<int> m_testData = {1,2,3,4,5,6,7,9,10,11};
    const QVector<int> m_testData2 = {1,2,3,4,5,6,7,8,9,10,11};

};


ParallelUtilsTests::ParallelUtilsTests()
{
}


void ParallelUtilsTests::trancheVectorForParallelizationTest()
{
    ERR_INIT

    const QVector<QVector<int>> expectedResult = {{1, 4, 7, 11}, {2, 5, 9}, {3, 6, 10}};
    QVector<QVector<int>> output;

    e = ParallelUtils::trancheVectorForParallelization(m_testData, 3, &output);
    QCOMPARE(e, eNoError);
    QCOMPARE(output.size(), expectedResult.size());

    for (int i = 0; i < output.size(); i++) {
        QCOMPARE(output.at(i), expectedResult.at(i));
    }

    const QVector<QVector<int>> expectedResult2 = {{1, 3, 5, 7, 10}, {2, 4, 6, 9, 11}};
    QVector<QVector<int>> output2;

    e = ParallelUtils::trancheVectorForParallelization(m_testData, 2, &output2);
    QCOMPARE(e, eNoError);
    QCOMPARE(output2.size(), expectedResult2.size());
    QCOMPARE(output2, expectedResult2);

}


void ParallelUtilsTests::trancheVectorForParallelizationInOrderTest() {

    ERR_INIT

    QVector<QVector<int>> output;
    e = ParallelUtils::trancheVectorForParallelizationInOrder(m_testData2, 12, 0, &output);
    QCOMPARE(e, eNoError);
    QCOMPARE(output.size(), m_testData2.size());

    e = ParallelUtils::trancheVectorForParallelizationInOrder(m_testData2, 3, 0, &output);
    QCOMPARE(e, eNoError);

    const QVector<QVector<int>> expectedResult = {{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11}};
    qDebug() << "DSD" << output;
    QCOMPARE(output, expectedResult);

    e = ParallelUtils::trancheVectorForParallelizationInOrder(m_testData2, 3, 2, &output);
    QCOMPARE(e, eNoError);

    const QVector<QVector<int>> expectedResultBuffer = {{1, 2, 3, 4, 5, 6}, {5, 6, 7, 8, 9, 10}, {9, 10, 11}};
    QCOMPARE(output, expectedResultBuffer);

    const QVector<int> truncData = m_testData2.mid(0, 5);
    e = ParallelUtils::trancheVectorForParallelizationInOrder(truncData, 8, 1, &output);
    QCOMPARE(e, eNoError);

    const QVector<QVector<int>> expectedResultBufferTrunc = {{1, 2}, {2, 3}, {3, 4}, {4, 5}, {5}};
    QCOMPARE(output, expectedResultBufferTrunc);

    QVector<int> bigger(35);
    std::iota(bigger.begin(), bigger.end(), 0);

//    e = ParallelUtils::trancheVectorForParallelizationInOrder(bigger, 8, 1, &output);
//    QCOMPARE(e, eNoError);
//    qDebug() << "SDFDS" << output;
//    qDebug() << "output size" << output.size();

}


QTEST_MAIN(ParallelUtilsTests)
#include "ParallelUtilsTests.moc"
