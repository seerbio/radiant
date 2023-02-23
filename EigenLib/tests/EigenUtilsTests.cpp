//
// Created by anichols on 11/07/2021.
//

#include "EigenUtils.h"

#include <QDebug>
#include <QMap>
#include <QtTest/QtTest>


class EigenUtilsTests : public QObject
{
    Q_OBJECT

public:
    EigenUtilsTests();
    ~EigenUtilsTests() override = default;


private:
    Eigen::VectorX<int> m_testVecInt;
    Eigen::VectorX<double> m_testVecDouble;
    Eigen::MatrixX<double> m_testVecDoubleRowMajor;
    Eigen::VectorX<int> m_testMatInt;
    Eigen::MatrixX<int> m_testMatIntRowMajor;
    Eigen::VectorX<int> m_testMatIntEmpty;


private Q_SLOTS:

    void initTestCase();
    void cleanupTestCase();

    void apexesTest();
    void troughtsTest();

};


EigenUtilsTests::EigenUtilsTests() : QObject()
{}

void EigenUtilsTests::initTestCase() {
    m_testVecDouble.resize(100);
    m_testVecDouble.coeffRef(2) = 10.0;
    m_testVecDouble.coeffRef(5) = -1.5;
    m_testVecDouble.coeffRef(22) = 666.6;
    m_testVecDouble.coeffRef(88) = -666.6;

    m_testVecInt = m_testVecDouble.cast<int>();
    m_testVecDoubleRowMajor = m_testVecDouble;

    m_testMatInt.resize(100, 100);
    m_testMatInt.coeffRef(2,2) = 10;
    m_testMatInt.coeffRef(5,5) = -1;
    m_testMatInt.coeffRef(22,22) = 666;
    m_testMatInt.coeffRef(88,88) = -666;

    m_testMatIntRowMajor = m_testMatInt;
}


void EigenUtilsTests::cleanupTestCase()
{}

void EigenUtilsTests::apexesTest() {

    Eigen::VectorX<double> vec(8);
    vec.coeffRef(0) = 0.999;
    vec.coeffRef(1) = 1.000;
    vec.coeffRef(2) = 1.000;
    vec.coeffRef(3) = 0.999;
    vec.coeffRef(4) = 0.999;
    vec.coeffRef(5) = 1.001;
    vec.coeffRef(6) = 1.000;
    vec.coeffRef(7) = 0.999;

    QMap<int, double> apexResult = EigenUtils::apexes(vec);

    QCOMPARE(apexResult.size(), 2);
    QCOMPARE(apexResult.firstKey(), 1);
    QCOMPARE(apexResult.lastKey(), 5);
    QCOMPARE(static_cast<int>(apexResult.first()), 1);
    QCOMPARE(static_cast<int>(apexResult.last()), 1);
}

void EigenUtilsTests::troughtsTest() {

    Eigen::VectorX<double> vec(8);
    vec.coeffRef(0) = 0.999;
    vec.coeffRef(1) = 1.000;
    vec.coeffRef(2) = 1.000;
    vec.coeffRef(3) = 0.999;
    vec.coeffRef(4) = 0.999;
    vec.coeffRef(5) = 1.001;
    vec.coeffRef(6) = 1.000;
    vec.coeffRef(7) = 0.999;

    QMap<int, double> troughsResult = EigenUtils::troughs(vec);

    QCOMPARE(troughsResult.size(), 3);
    QCOMPARE(troughsResult.firstKey(), 0);
    QCOMPARE(troughsResult.lastKey(), 7);
    QCOMPARE(static_cast<int>(troughsResult.first()), 0);
    QCOMPARE(static_cast<int>(troughsResult.last()), 0);
}


QTEST_MAIN(EigenUtilsTests)
#include "EigenUtilsTests.moc"
