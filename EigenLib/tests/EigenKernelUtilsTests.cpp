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



private Q_SLOTS:

    void buildSavitzkyGolayKernelTest();

    void buildGaussianFilterTest();



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


QTEST_MAIN(EigenKernelUtilsTests)
#include "EigenKernelUtilsTests.moc"
