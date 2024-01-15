//
// Created by anichols on 11/07/2021.
//

#include "EigenKernelUtils.h"
#include "EigenUtils.h"

#include <QDebug>
#include <QMap>
#include <QtTest/QtTest>


class EigenKernelUtilsTests : public QObject
{
    Q_OBJECT

public:
    EigenKernelUtilsTests();
    ~EigenKernelUtilsTests() override = default;


//TODO write tests for the rest of the methods.
private Q_SLOTS:

    static void buildSavitzkyGolayKernelTest();
    static void lineSpaceTest();
    static void buildGaussianFilter1DTest();
    static void buildMexicanHatFilter1DTest();
    static void calculateNumberOfStridesTest();
    static void addPaddingToMatrixRowWiseTest();

    static void applyKernelToEachColumnInMatrix();

};

EigenKernelUtilsTests::EigenKernelUtilsTests() : QObject() {}

void EigenKernelUtilsTests::buildSavitzkyGolayKernelTest() {

    ERR_INIT

    const int windowSize = 7;
    const int order = 2;
    const int derivative = 0;
    const int rate = 1;

    Eigen::MatrixX<double> savitskyGolayKernel;
    e = EigenKernelUtils::buildSavitzkyGolayKernel(
            windowSize,
            order,
            derivative,
            rate,
            &savitskyGolayKernel
    );
    QCOMPARE(e, eNoError);

    const double threshold = 0.00001;
    QVERIFY((savitskyGolayKernel(0) - -0.0952381 < threshold));
    QVERIFY((savitskyGolayKernel(1) - 0.142857 < threshold));
    QVERIFY((savitskyGolayKernel(2) - 0.285714 < threshold));
    QVERIFY((savitskyGolayKernel(3) - 0.333333 < threshold));
    QVERIFY((savitskyGolayKernel(4) - 0.285714 < threshold));
    QVERIFY((savitskyGolayKernel(5) - 0.142857 < threshold));
    QVERIFY((savitskyGolayKernel(6) - -0.0952381 < threshold));

}

void EigenKernelUtilsTests::lineSpaceTest() {

    const Eigen::VectorX<double> lineSpaceVec = EigenKernelUtils::lineSpace<double>(
            -1.0,
            1.0,
            11,
            true
            );

    const QVector<double> expected = {
            -1,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1
    };

    QCOMPARE(lineSpaceVec.size(), expected.size());
    for (int i = 0; i < expected.size(); i++) {
        QCOMPARE(lineSpaceVec.coeff(i), expected.at(i));
    }

}

void EigenKernelUtilsTests::buildGaussianFilter1DTest() {

    const int filterLen = 3;
    const double sigma = 2.0;
    Eigen::VectorX<double> gaussianFilter = EigenKernelUtils::buildGaussianFilter1D<double>(
            filterLen,
            sigma
    );

    QVector<double> expected = {
            0.02154,
            0.15915,
            0.02154
    };

    QCOMPARE(gaussianFilter.size(), expected.size());
    for (int i = 0; i < expected.size(); i++) {
        QCOMPARE(MathUtils::pRound(gaussianFilter.coeff(i), 5), expected.at(i));
    }

}

void EigenKernelUtilsTests::buildMexicanHatFilter1DTest() {

    const Eigen::VectorX<double> mexHatVec = EigenKernelUtils::buildMexicanHatFilter1D<double>(
            11,
            2.0
            );

    const QVector<double> expected = {
            -0.141,
            -0.249,
            -0.249,
            0.0,
            0.406,
            0.613,
            0.406,
            0.0,
            -0.249,
            -0.249,
            -0.141
    };

    QCOMPARE(mexHatVec.size(), expected.size());

    for (int i = 0; i < expected.size(); i++) {
        QCOMPARE(MathUtils::pRound(mexHatVec.coeff(i), 3), expected.at(i));
    }

}

void EigenKernelUtilsTests::calculateNumberOfStridesTest() {
    const int numberOfStrides = EigenKernelUtils::calculateNumberOfStrides(11, 3, 1);
    QCOMPARE(numberOfStrides, 9);
}

void EigenKernelUtilsTests::addPaddingToMatrixRowWiseTest() {

    Eigen::MatrixX<double> mat(7,3);
    mat.setOnes();

    const int filterLength = 3;

    const Eigen::MatrixX<double> paddedMatrix = EigenKernelUtils::addPaddingToMatrixRowWise(
            mat,
            filterLength
            );

    QCOMPARE(static_cast<int>(paddedMatrix.rows()), 9);
    QCOMPARE(static_cast<int>(paddedMatrix.cols()), 3);
    QCOMPARE(paddedMatrix.coeff(0,0), 0.0);
    QCOMPARE(paddedMatrix.coeff(1,0), 1.0);
    QCOMPARE(paddedMatrix.coeff(7,0), 1.0);
    QCOMPARE(paddedMatrix.coeff(8,0), 0.0);


}

void EigenKernelUtilsTests::applyKernelToEachColumnInMatrix() {

    Eigen::VectorX<double> kernel(3);
    kernel.setOnes();

    Eigen::MatrixX<double> mat(5,2);
    mat << 1,0,2,0,3,0,4,0,5,0;

    const Eigen::MatrixX<double> matSmoothed = EigenKernelUtils::applyKernelToEachColumnInMatrix(
            mat,
            kernel
            ) ;

    QCOMPARE(static_cast<int>(matSmoothed.rows()), 5);
    QCOMPARE(static_cast<int>(matSmoothed.cols()), 2);
    QCOMPARE(matSmoothed.coeff(0,0), 3.0);
    QCOMPARE(matSmoothed.coeff(1,0), 6.0);
    QCOMPARE(matSmoothed.coeff(2,0), 9.0);
    QCOMPARE(matSmoothed.coeff(3,0), 12.0);
    QCOMPARE(matSmoothed.coeff(4,0), 9.0);
    QCOMPARE(matSmoothed.coeff(0,1), 0.0);
    QCOMPARE(matSmoothed.coeff(1,1), 0.0);
    QCOMPARE(matSmoothed.coeff(4,1), 0.0);

}


QTEST_MAIN(EigenKernelUtilsTests)
#include "EigenKernelUtilsTests.moc"
