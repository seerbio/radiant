//
// Created by anichols on 11/07/2021.
//

#include "MathUtils.h"

#include <QDebug>
#include <QList>
#include <QtTest/QtTest>
#include <QVector>

#include <vector>


class MathUtilsTests : public QObject
{
    Q_OBJECT

public:
    MathUtilsTests();
    ~MathUtilsTests() override = default;


private Q_SLOTS:

    void initTestCase();
    void medianTest();
    void meanTest();
    void pRoundTest();
    void hashDecimalTest();
    void unHashDecimalTest();
    void factorialTest();
    void closestTest();
    void cleanupTestCase();

};

MathUtilsTests::MathUtilsTests() : QObject()
{
}

void MathUtilsTests::initTestCase()
{
}

void MathUtilsTests::medianTest()
{
    const std::vector<int> vec = {2,3,2,7,4,7};
    const double vecMedian = MathUtils::median(vec);
    QCOMPARE(vecMedian, 3.5);

    const QList<int> list = {2,3,2,3,9,7,6};
    double listMedian = MathUtils::median(list);
    QCOMPARE(listMedian, 3.);

    const QVector<int> qvec = {3,1,2};
    double qvecMedian = MathUtils::median(qvec);
    QCOMPARE(qvecMedian, 2.);

    const std::vector<double> vec2Test = {4,5};
    const double vec2TestMedian = MathUtils::median(vec2Test);
    QCOMPARE(vec2TestMedian, 4.5);
}

void MathUtilsTests::meanTest()
{
    const QVector<int> v1 = {};
    const QVector<int> v2 = {1,2,3,4,5};
    const std::vector<int> v1Std = v1.toStdVector();

    const double meanV1 = MathUtils::mean(v1);
    const double meanV2 = MathUtils::mean(v2);
    const double meanV1Std = MathUtils::mean(v1Std);

    QCOMPARE(meanV1, 0.0);
    QCOMPARE(meanV2, 3.0);
    QCOMPARE(meanV1Std, 0.0);
}

void MathUtilsTests::pRoundTest() {

    const double roundVal1 = MathUtils::pRound(3.135, 2);
    QCOMPARE(roundVal1, 3.14);

    const double roundVal2 = MathUtils::pRound(3.134, 2);
    QCOMPARE(roundVal2, 3.13);
}

void MathUtilsTests::cleanupTestCase()
{
}

void MathUtilsTests::hashDecimalTest() {

    const double testVal = 666.66;
    const int hashedTestVal1 = MathUtils::hashDecimal(testVal, 4);
    QCOMPARE(hashedTestVal1, 6666600);

    const int hashedTestVal2 = MathUtils::hashDecimal(testVal, 1);
    QCOMPARE(hashedTestVal2, 6667);

}

void MathUtilsTests::unHashDecimalTest() {

    const int testVal = 66666;
    const auto unHashedTestVal1 = MathUtils::unHashDecimal<double>(testVal, 1);
    QCOMPARE(QString::number(unHashedTestVal1), QString::number(6666.6));

    const auto unHashedTestVal2 = MathUtils::unHashDecimal<double>(testVal, 0);
    QCOMPARE(QString::number(unHashedTestVal2), QString::number(66666));
}

void MathUtilsTests::factorialTest() {

    QCOMPARE(MathUtils::factorial(30), 9223372036854775807);
    QCOMPARE(MathUtils::factorial(0), 1);
    QCOMPARE(MathUtils::factorial(-1), 0);
    QCOMPARE(MathUtils::factorial(3), 6);

}

void MathUtilsTests::closestTest() {

    const QVector<double> vec = {1.2, 1.3, 1.4, 1.45, 1.5};

    int closestIndex = MathUtils::closest(vec, 1.349);
    QCOMPARE(closestIndex, 1);

    int closestIndex1 = MathUtils::closest(vec, 1.49);
    QCOMPARE(closestIndex1, 4);
}


QTEST_MAIN(MathUtilsTests)
#include "MathUtilsTests.moc"
