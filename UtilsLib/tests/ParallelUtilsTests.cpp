//
// Created by anichols on 11/07/2021.
//

#include "ParallelUtils.h"

#include <QtTest/QtTest>

class ParallelUtilsTests : public QObject
{
    Q_OBJECT

public:
    ParallelUtilsTests() = default;
    ~ParallelUtilsTests() override = default;

private Q_SLOTS:
    void trancheVectorForParallelizationTest();
    void trancheVectorForParallelizationInOrderTest();
    static void trancheMapValueVectorsByKeyForParallelization();
    static void convertMapToVectorPairsTest();
    static void zipTest();
    static void unZipTest();
    static void convertVectorToMapTest();
    static void convertMapToPointsTest();

private:

    const QVector<int> m_testData = {1,2,3,4,5,6,7,9,10,11};
    const QVector<int> m_testData2 = {1,2,3,4,5,6,7,8,9,10,11};

};

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

void ParallelUtilsTests::trancheMapValueVectorsByKeyForParallelization() {

    ERR_INIT

    QMap<int, QVector<double>> testMap;
    testMap[1] = QVector<double>({1.0, 2.0, 3.0, 4.0});
    testMap[2] = QVector<double>({5.0, 6.0, 7.0, 8.0});
    testMap[3] = QVector<double>({9.0, 10.0, 11.0, 12.0});

    QVector<QMap<int, QVector<double>>> result;
    e = ParallelUtils::trancheMapValueVectorsByKeyForParallelization(
            testMap,
            2,
            &result
            );

    QCOMPARE(e, eNoError);
    QCOMPARE(result.size(), 2);
    QCOMPARE(result.at(0).value(1), QVector<double>({1.0, 3.0}));
    QCOMPARE(result.at(0).value(2), QVector<double>({5.0, 7.0}));
    QCOMPARE(result.at(0).value(3), QVector<double>({9.0, 11.0}));
    QCOMPARE(result.at(1).value(1), QVector<double>({2.0, 4.0}));
    QCOMPARE(result.at(1).value(2), QVector<double>({6.0, 8.0}));
    QCOMPARE(result.at(1).value(3), QVector<double>({10.0, 12.0}));

}

void ParallelUtilsTests::convertMapToVectorPairsTest() {

    QMap<int, int> testMap;
    testMap[1] = 10;
    testMap[2] = 20;
    testMap[3] = 30;

    QVector<QPair<int, int>> resultVector = ParallelUtils::convertMapToVectorPairs(testMap);

    QCOMPARE(resultVector.size(), testMap.size());

    for (const auto &pair : resultVector) {
        QVERIFY(testMap.contains(pair.first));
        QCOMPARE(pair.second, testMap.value(pair.first));
    }
}

void ParallelUtilsTests::zipTest() {

    QMap<int, int> testMap;
    testMap.insert(1, 100);
    testMap.insert(2, 200);
    testMap.insert(3, 300);

    QVector<QPair<int, int>> result = ParallelUtils::convertMapToVectorPairs<int, int>(testMap);
    QVERIFY(result.size() == testMap.size());

    for(int i = 0; i < result.size(); i++) {
        QCOMPARE(result[i].first, testMap.keys()[i]);
        QCOMPARE(result[i].second, testMap.value(result[i].first));
    }

    QVector<float> z1{1, 2, 3};
    QVector<float> z2{4, 5, 6};
    QVector<QPointF> zipResult;

    ERR_INIT
    e = ParallelUtils::zip(z1, z2, &zipResult);
    QCOMPARE(e, eNoError);

    for(int i = 0; i < zipResult.size(); ++i) {
        QCOMPARE(zipResult[i].x(), z1[i]);
        QCOMPARE(zipResult[i].y(), z2[i]);
    }

    QVector<QPair<float, float>> zipResultPair;
    e = ParallelUtils::zip(z1, z2, &zipResultPair);
    QCOMPARE(e, eNoError);


    for(int i = 0; i < zipResult.size(); ++i) {
        QCOMPARE(zipResultPair[i].first, z1[i]);
        QCOMPARE(zipResultPair[i].second, z2[i]);
    }

}

void ParallelUtilsTests::unZipTest() {

    QVector<QPair<int, int>> input{qMakePair(1, 2), qMakePair(3, 4), qMakePair(5, 6)};

    const QPair<QVector<int>, QVector<int>> output = ParallelUtils::unZip(input);
    QVector<int> expectedFirst{1, 3, 5};
    QVector<int> expectedSecond{2, 4, 6};

    QCOMPARE(output.first, expectedFirst);
    QCOMPARE(output.second, expectedSecond);

}

void ParallelUtilsTests::convertVectorToMapTest() {

    QVector<int> input{1, 2, 3, 4, 5};

    QMap<int, int> output = ParallelUtils::convertVectorToMap(input);

    QCOMPARE(output.size(), input.size());

    for(int i = 0; i < output.size(); ++i) {
        QCOMPARE(output[i], input[i]);
    }

}

void ParallelUtilsTests::convertMapToPointsTest() {

    QMap<float, float> testMap;
    testMap.insert(1, 10);
    testMap.insert(2, 20);
    testMap.insert(3, 30);

    const QVector<QPointF> resultPoints = ParallelUtils::convertMapToPoints(testMap);
    QCOMPARE(resultPoints.size(), testMap.size());


    for(int i = 0; i < resultPoints.size(); i++) {
        QPointF point = resultPoints[i];
        const float key = testMap.keys()[i];
        const float value = testMap[key];

        QCOMPARE(point.x(), key);
        QCOMPARE(point.y(), value);
    }
}


QTEST_MAIN(ParallelUtilsTests)
#include "ParallelUtilsTests.moc"
