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

    void buildSavitzkyGolayKernelTest();

    void buildGaussianFilterTest();

    void addPaddingToMatrixRowWiseTest();

    void applyKernelRowumnWiseToMatrixTest();

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

void EigenKernelUtilsTests::buildGaussianFilterTest() {

    const int filterLen = 3;
    const double sigma = 2.0;
    Eigen::VectorX<double> gaussianFilter = EigenKernelUtils::buildGaussianFilter1D(
            filterLen,
            sigma
    );

    std::cout << gaussianFilter << std::endl;

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

void EigenKernelUtilsTests::applyKernelRowumnWiseToMatrixTest() {

    Eigen::VectorX<double> kernel(3);
    kernel.setOnes();

    Eigen::MatrixX<double> mat(5,2);
    mat << 1,1,2,2,3,3,4,4,5,5;

    const Eigen::MatrixX<double> matSmoothed = EigenKernelUtils::applyKernelRowumnWiseToMatrix(
            mat,
            kernel
            ) ;

    QCOMPARE(static_cast<int>(matSmoothed.rows()), 5);
    QCOMPARE(static_cast<int>(matSmoothed.cols()), 2);
    QCOMPARE(matSmoothed.coeff(0,0), 3.0);
    QCOMPARE(matSmoothed.coeff(1,0), 6.0);
    QCOMPARE(matSmoothed.coeff(4,0), 9.0);

}


QTEST_MAIN(EigenKernelUtilsTests)
#include "EigenKernelUtilsTests.moc"
