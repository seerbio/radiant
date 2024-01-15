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
    static void convolveVectorWithKernelTest();
    static void addPaddingToSparseVectorTest();
    static void addPaddingToSparseMatrixRowWiseTest();
    static void addPaddingToMatrixRowWiseTest();
    static void addPaddingToSparseMatrixColWiseTest();
    static void applyKernelToEachColumnInMatrixTest();
    static void applyKernelToEachRowInMatrixTest();
    static void savitskyGolaySmoothTest();

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

void EigenKernelUtilsTests::convolveVectorWithKernelTest() {

    Eigen::VectorX<double> kernel(3);
    kernel.setOnes();

    Eigen::VectorX<double> vec(10);
    vec.setZero();
    vec.coeffRef(4) = 666.6;

    const Eigen::VectorX<double> vecSmoothed = EigenKernelUtils::convolveVectorWithKernel(
            vec,
            kernel
            );

    QCOMPARE(vecSmoothed.coeff(3), 666.6);
    QCOMPARE(vecSmoothed.coeff(4), 666.6);
    QCOMPARE(vecSmoothed.coeff(5), 666.6);

    Eigen::SparseVector<double> vecSparse(10);
    vecSparse.setZero();
    vecSparse.coeffRef(4) = 666.6;

    const Eigen::SparseVector<double> vecSparseSmoothed = EigenKernelUtils::convolveVectorWithKernel(
            vecSparse,
            kernel
    );

    QCOMPARE(vecSparseSmoothed.coeff(3), 666.6);
    QCOMPARE(vecSparseSmoothed.coeff(4), 666.6);
    QCOMPARE(vecSparseSmoothed.coeff(5), 666.6);
}

void EigenKernelUtilsTests::addPaddingToSparseVectorTest() {

    Eigen::SparseVector<double> vecSparse(10);
    vecSparse.setZero();
    vecSparse.coeffRef(4) = 666.6;

    const Eigen::SparseVector<double> vecSparsePadded = EigenKernelUtils::addPaddingToSparseVector(
            vecSparse,
            3
            );

    QCOMPARE(vecSparse.coeff(4), 666.6);
    QCOMPARE(vecSparsePadded.coeff(5), 666.6);
    QCOMPARE(vecSparsePadded.size(), vecSparse.size() + 2);

}

void EigenKernelUtilsTests::addPaddingToSparseMatrixRowWiseTest() {

    Eigen::SparseMatrix<double, Eigen::RowMajor> mat(7,3);
    for (int r = 0; r < mat.rows(); r++) {
        for (int c = 0; c < mat.cols(); c++) {
            mat.insert(r,c) = 1.0;
        }
    }

    const int filterLength = 3;

    const Eigen::SparseMatrix<double, Eigen::RowMajor> paddedMatrix = EigenKernelUtils::addPaddingToSparseMatrixRowWise(
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

void EigenKernelUtilsTests::addPaddingToSparseMatrixColWiseTest() {
    Eigen::SparseMatrix<double> mat(5000,20);
    mat.coeffRef(0, 0) = 666.6;
    mat.coeffRef(0, 19) = 666.6;

    const Eigen::MatrixX<double> matSmoothed = EigenKernelUtils::addPaddingToSparseMatrixColWise<double>(
            mat,
            3
    ) ;

    QCOMPARE(static_cast<int>(matSmoothed.rows()), 5000);
    QCOMPARE(static_cast<int>(matSmoothed.cols()), 22);
    QCOMPARE(matSmoothed.coeff(0,0), 0.0);
    QCOMPARE(matSmoothed.coeff(0,1), 666.6);
    QCOMPARE(matSmoothed.coeff(0,20), 666.6);
    QCOMPARE(matSmoothed.coeff(0,21), 0.0);
}

void EigenKernelUtilsTests::applyKernelToEachRowInMatrixTest() {

    Eigen::SparseMatrix<double> mat(10, 10);
    mat.setZero();
    mat.coeffRef(4, 4) = 666.6;

    Eigen::VectorX<double> kernel(3);
    kernel.setOnes();

    const Eigen::SparseMatrix<double> matSmoothed = EigenKernelUtils::applyKernelToEachRowInMatrix(mat, kernel);

    QCOMPARE(matSmoothed.rows(), 10);
    QCOMPARE(matSmoothed.cols(), 10);
    QCOMPARE(matSmoothed.coeff(4, 3), 666.6);
    QCOMPARE(matSmoothed.coeff(4, 4), 666.6);
    QCOMPARE(matSmoothed.coeff(4, 5), 666.6);
    QCOMPARE(matSmoothed.coeff(3, 3), 0.0);
    QCOMPARE(matSmoothed.coeff(3, 4), 0.0);
    QCOMPARE(matSmoothed.coeff(3, 5), 0.0);
    QCOMPARE(matSmoothed.coeff(5, 3), 0.0);
    QCOMPARE(matSmoothed.coeff(5, 4), 0.0);
    QCOMPARE(matSmoothed.coeff(5, 5), 0.0);


}

void EigenKernelUtilsTests::applyKernelToEachColumnInMatrixTest() {

    Eigen::MatrixX<double> mat(10, 10);
    mat.setZero();
    mat.coeffRef(4, 4) = 666.6;

    Eigen::VectorX<double> kernel(3);
    kernel.setOnes();

    const Eigen::MatrixX<double> matSmoothed = EigenKernelUtils::applyKernelToEachColumnInMatrix(mat, kernel);

    QCOMPARE(matSmoothed.rows(), 10);
    QCOMPARE(matSmoothed.cols(), 10);
    QCOMPARE(matSmoothed.coeff(3, 4), 666.6);
    QCOMPARE(matSmoothed.coeff(4, 4), 666.6);
    QCOMPARE(matSmoothed.coeff(5, 4), 666.6);
    QCOMPARE(matSmoothed.coeff(3, 3), 0.0);
    QCOMPARE(matSmoothed.coeff(4, 3), 0.0);
    QCOMPARE(matSmoothed.coeff(5, 3), 0.0);
    QCOMPARE(matSmoothed.coeff(3, 5), 0.0);
    QCOMPARE(matSmoothed.coeff(4, 5), 0.0);
    QCOMPARE(matSmoothed.coeff(5, 5), 0.0);


}

void EigenKernelUtilsTests::savitskyGolaySmoothTest() {

    ERR_INIT

    Eigen::VectorX<double> vecOG(10);
    vecOG.setZero();
    vecOG.coeffRef(4) = 666.6;

    Eigen::VectorX<double> vec1 = vecOG;
    e = EigenKernelUtils::savitskyGolaySmooth(
            3,
            1,
            0,
            1,
            &vec1
            );
    QCOMPARE(e, eNoError);
    QCOMPARE(vec1.size(), vecOG.size());
    QCOMPARE(vec1.coeff(3), 222.2);
    QCOMPARE(vec1.coeff(4), 222.2);
    QCOMPARE(vec1.coeff(5), 222.2);

    Eigen::VectorX<double> vec2 = vecOG;
    e = EigenKernelUtils::savitskyGolaySmooth(
            5,
            2,
            0,
            1,
            &vec2
    );
    QCOMPARE(e, eNoError);
    QCOMPARE(vec2.size(), vecOG.size());
    QCOMPARE(MathUtils::pRound(vec2.coeff(2), 3) , -57.137);
    QCOMPARE(MathUtils::pRound(vec2.coeff(3), 3) , 228.549);
    QCOMPARE(MathUtils::pRound(vec2.coeff(4), 3) , 323.777);
    QCOMPARE(MathUtils::pRound(vec2.coeff(5), 3) , 228.549);
    QCOMPARE(MathUtils::pRound(vec2.coeff(6), 3) , -57.137);

    Eigen::VectorX<double> vec3 = vecOG;
    e = EigenKernelUtils::savitskyGolaySmooth(
            5,
            2,
            1,
            1,
            &vec3
    );

    QCOMPARE(e, eNoError);
    QCOMPARE(MathUtils::pRound(vec3.coeff(2), 3) , 133.32);
    QCOMPARE(MathUtils::pRound(vec3.coeff(3), 3) , 66.66);
    QCOMPARE(MathUtils::pRound(vec3.coeff(4), 3) , MathUtils::pRound(-2.27346e-14, 3));
    QCOMPARE(MathUtils::pRound(vec3.coeff(5), 3) , -66.66);
    QCOMPARE(MathUtils::pRound(vec3.coeff(6), 3) , -133.32);

}


QTEST_MAIN(EigenKernelUtilsTests)
#include "EigenKernelUtilsTests.moc"
