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

ParallelUtilsTests::ParallelUtilsTests(){}

void ParallelUtilsTests::trancheVectorForParallelizationTest() {
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

    const QVector<QVector<int>> expectedResultBigger = {
            QVector<int>({0, 1, 2, 3}),
            QVector<int>({4, 5, 6, 7}),
            QVector<int>({8, 9, 10, 11}),
            QVector<int>({12, 13, 14, 15}),
            QVector<int>({16, 17, 18, 19}),
            QVector<int>({20, 21, 22, 23}),
            QVector<int>({24, 25, 26, 27}),
            QVector<int>({28, 29, 30, 31}),
            QVector<int>({32, 33, 34})
    };
    e = ParallelUtils::trancheVectorForParallelizationInOrder(bigger, 8, 0, &output);
    QCOMPARE(e, eNoError);
    QCOMPARE(output.size(), 9);
    QCOMPARE(output, expectedResultBigger);

    output.clear();
    const QVector<QVector<int>> expectedResultBiggerBuffer = {
            QVector({0, 1, 2, 3, 4}),
            QVector({4, 5, 6, 7, 8}),
            QVector({8, 9, 10, 11, 12}),
            QVector({12, 13, 14, 15, 16}),
            QVector({16, 17, 18, 19, 20}),
            QVector({20, 21, 22, 23, 24}),
            QVector({24, 25, 26, 27, 28}),
            QVector({28, 29, 30, 31, 32}),
            QVector({32, 33, 34})
    };
    e = ParallelUtils::trancheVectorForParallelizationInOrder(bigger, 8, 1, &output);
    QCOMPARE(e, eNoError);
    QCOMPARE(output.size(), 9);
    QCOMPARE(output, expectedResultBiggerBuffer);

}


QTEST_MAIN(ParallelUtilsTests)
#include "ParallelUtilsTests.moc"
