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

    void rowWiseCosineSimilarOfMatricesTest();
    void rowWiseKLDivergeneceTest();

    void removeRowsBelowThresholdTest();

    void minMaxScaleVectorTest();
    void minMaxScaleMatrixTest();

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


    Eigen::VectorX<double> vec2(14);
    vec2.coeffRef(0) = 9465.27;
    vec2.coeffRef(1) = 0.0;
    vec2.coeffRef(2) = 54665.6;
    vec2.coeffRef(3) = 539886.0;
    vec2.coeffRef(4) = 47320.7;
    vec2.coeffRef(5) = 7205.26;
    vec2.coeffRef(6) = 19682.5 ;
    vec2.coeffRef(7) = 13817.7;
    vec2.coeffRef(8) = 23367.9;
    vec2.coeffRef(9) = 20809.4;
    vec2.coeffRef(10) = 0;
    vec2.coeffRef(11) = 0;
    vec2.coeffRef(12) = 7581.3;
    vec2.coeffRef(13) = 0;

    QMap<int, double> apexResult2 = EigenUtils::apexes(vec2);

    QCOMPARE(apexResult2.size(), 5);
    QCOMPARE(apexResult2.firstKey(), 0);
    QCOMPARE(apexResult2.lastKey(), 12);
    QCOMPARE(static_cast<int>(apexResult2.first()), 9465);
    QCOMPARE(static_cast<int>(apexResult2.last()), 7581);

    qDebug() << apexResult2;

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

void EigenUtilsTests::rowWiseCosineSimilarOfMatricesTest() {

    Eigen::MatrixX<double> mat1(2,2);
    Eigen::MatrixX<double> mat2(2,2);

    mat1 << 1,2,3,4;
    mat2 << 2,4,6,8;

    Eigen::VectorX<double> cosineSim = EigenUtils::rowWiseCosineSimilarOfMatrices(mat1, mat2);

    QCOMPARE(cosineSim.coeff(0), 1.0);
    QCOMPARE(cosineSim.coeff(1), 1.0);

}

void EigenUtilsTests::rowWiseKLDivergeneceTest() {

    Eigen::MatrixX<double> mat1(3,3);
    Eigen::MatrixX<double> mat2(3,3);

    mat1 << 1,2,3,4,5,6,7,8,9;
    mat2 << 1,2,3,4,2.5,3,0,0,0;

    Eigen::VectorX<double> res = EigenUtils::rowWiseKLDivergence(mat1, mat2);


    EigenUtils::klDivergence(mat1.row(1), mat2.row(1));

    //TODO make these resluts into a test.
    for (int i = 0; i < 3; i++) {
        std::cout << EigenUtils::klDivergence(mat1.row(i), mat2.row(i)) << std::endl;
    }

    std::cout << res << std::endl;

}

void EigenUtilsTests::removeRowsBelowThresholdTest() {

    Eigen::MatrixX<double> mat(4,2);

    mat << 1.0, 2.0, 3.0 , 4.0 , 5.0 , 6.0 , 7.0 , 8.0;
    EigenUtils::removeRowsAboveOrBelowThreshold<double>(
            3.5,
            1,
            EigenUtils::ThresholderDirection::Below,
            &mat
    );

    QCOMPARE(mat.rows(), 3);
    QCOMPARE(mat.coeffRef(0,0), 3.0);
    QCOMPARE(mat.coeffRef(2,0), 7.0);

    Eigen::MatrixX<int> matInt(4,2);

    matInt << 1, 2, 3 , 4 , 5 , 6 , 7 , 8;
    EigenUtils::removeRowsAboveOrBelowThreshold(
            4,
            0,
            EigenUtils::ThresholderDirection::Below,
            &matInt
    );

    QCOMPARE(matInt.rows(), 2);
    QCOMPARE(matInt.coeffRef(0,0), 5);
    QCOMPARE(matInt.coeffRef(1,0), 7);

    Eigen::MatrixX<int> matIntAbove(4,2);

    matIntAbove << 1, 2, 3 , 4 , 5 , 6 , 7 , 8;
    EigenUtils::removeRowsAboveOrBelowThreshold(
            4,
            0,
            EigenUtils::ThresholderDirection::Above,
            &matIntAbove
    );

    QCOMPARE(matIntAbove.rows(), 2);
    QCOMPARE(matIntAbove.coeffRef(0,0), 1);
    QCOMPARE(matIntAbove.coeffRef(1,0), 3);

}

void EigenUtilsTests::minMaxScaleVectorTest() {

    Eigen::VectorX<double> vec(4);
    vec << 1, 2, 3, 4;

    EigenUtils::minMaxScaleVector(&vec);

    QVERIFY(MathUtils::tSame(vec.coeff(2), 0.66666));
}

void EigenUtilsTests::minMaxScaleMatrixTest() {

    Eigen::MatrixX<double> mat(4, 2);
    mat << 1, 1, 2, 2, 3, 3, 4, 4;

    EigenUtils::minMaxScaleMatrix(&mat);

    QVERIFY(MathUtils::tSame(mat.coeff(2,0), 0.66666));
    QVERIFY(MathUtils::tSame(mat.coeff(2,1), 0.66666));
}


QTEST_MAIN(EigenUtilsTests)
#include "EigenUtilsTests.moc"
