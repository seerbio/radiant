//
// Created by anichols on 11/07/2021.
//

#include "EigenSparseUtils.h"

#include <QDebug>
#include <QMap>
#include <QtTest/QtTest>

#include <Eigen/IterativeLinearSolvers>


class EigenSparseUtilsTests : public QObject
{
    Q_OBJECT

public:
    EigenSparseUtilsTests();
    ~EigenSparseUtilsTests() override = default;


private:
    Eigen::SparseVector<int> m_testVecInt;
    Eigen::SparseVector<double> m_testVecDouble;
    Eigen::SparseVector<double, Eigen::RowMajor> m_testVecDoubleRowMajor;
    Eigen::SparseMatrix<int> m_testMatInt;
    Eigen::SparseMatrix<int, Eigen::RowMajor> m_testMatIntRowMajor;
    Eigen::SparseMatrix<int> m_testMatIntEmpty;


private Q_SLOTS:

    void initTestCase();
    void maxTest();
    void minTest();
    void meanTest();
    void stDevTest();
    void isValidTest();
    void medianTest();
    void removeElementBelowThresholdTest();
    void rollTest();
    void apexesTest();
    void buildCombFilterTest();
    void markVectorApexesTest();
    void markMatrixApexesTest();
    void solverTest();
    void cleanupTestCase();

};


EigenSparseUtilsTests::EigenSparseUtilsTests() : QObject()
{}

void EigenSparseUtilsTests::initTestCase() {
    m_testVecDouble.resize(100);
    m_testVecDouble.insert(2) = 10.0;
    m_testVecDouble.insert(5) = -1.5;
    m_testVecDouble.insert(22) = 666.6;
    m_testVecDouble.insert(88) = -666.6;

    m_testVecInt = m_testVecDouble.cast<int>();
    m_testVecDoubleRowMajor = m_testVecDouble;

    m_testMatInt.resize(100, 100);
    m_testMatInt.insert(2,2) = 10;
    m_testMatInt.insert(5,5) = -1;
    m_testMatInt.insert(22,22) = 666;
    m_testMatInt.insert(88,88) = -666;

    m_testMatIntRowMajor = m_testMatInt;
}


void EigenSparseUtilsTests::maxTest() {
    QCOMPARE(EigenSparseUtils::max(m_testVecInt), 666);
    QCOMPARE(EigenSparseUtils::max(m_testVecDouble), 666.6);
    QCOMPARE(EigenSparseUtils::max(m_testVecDoubleRowMajor), 666.6);
    QCOMPARE(EigenSparseUtils::max(m_testMatInt), 666);
    QCOMPARE(EigenSparseUtils::max(m_testMatIntRowMajor), 666);
    QCOMPARE(EigenSparseUtils::max(m_testMatIntEmpty), 0);
}


void EigenSparseUtilsTests::minTest() {
    QCOMPARE(EigenSparseUtils::min(m_testVecInt), -666);
    QCOMPARE(EigenSparseUtils::min(m_testVecDouble), -666.6);
    QCOMPARE(EigenSparseUtils::min(m_testVecDoubleRowMajor), -666.6);
    QCOMPARE(EigenSparseUtils::min(m_testMatInt), -666);
    QCOMPARE(EigenSparseUtils::min(m_testMatIntRowMajor), -666);
    QCOMPARE(EigenSparseUtils::min(m_testMatIntEmpty), 0);
}


void EigenSparseUtilsTests::meanTest() {
    QCOMPARE(EigenSparseUtils::mean(m_testVecInt), 2.25);
    QCOMPARE(EigenSparseUtils::mean(m_testVecDouble), 2.125);
    QCOMPARE(EigenSparseUtils::mean(m_testVecDoubleRowMajor), 2.125);
    QCOMPARE(EigenSparseUtils::mean(m_testMatInt), 2.25);
    QCOMPARE(EigenSparseUtils::mean(m_testMatIntRowMajor), 2.25);
    QCOMPARE(EigenSparseUtils::mean(m_testMatIntEmpty), 0.0);
}


void EigenSparseUtilsTests::stDevTest() {
    QCOMPARE(std::floor(EigenSparseUtils::stDev(m_testVecInt)), std::floor(470.955));
    QCOMPARE(std::floor(EigenSparseUtils::stDev(m_testVecDouble)), std::floor(471.38));
    QCOMPARE(std::floor(EigenSparseUtils::stDev(m_testVecDoubleRowMajor)), std::floor(471.38));
    QCOMPARE(std::floor(EigenSparseUtils::stDev(m_testMatInt)), std::floor(470.5));
    QCOMPARE(std::floor(EigenSparseUtils::stDev(m_testMatIntRowMajor)), std::floor(470.5));
    QCOMPARE(std::floor(EigenSparseUtils::stDev(m_testMatIntEmpty)), 0.0);
}


void EigenSparseUtilsTests::isValidTest() {
    QCOMPARE(EigenSparseUtils::validIndex(m_testVecInt, 99),  true);
    QCOMPARE(EigenSparseUtils::validIndex(m_testVecInt, 100),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testVecInt, 0),  true);
    QCOMPARE(EigenSparseUtils::validIndex(m_testVecInt, -1),  false);

    QCOMPARE(EigenSparseUtils::validIndex(m_testVecDoubleRowMajor, 99),  true);
    QCOMPARE(EigenSparseUtils::validIndex(m_testVecDoubleRowMajor, 100),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testVecDoubleRowMajor, 0),  true);
    QCOMPARE(EigenSparseUtils::validIndex(m_testVecDoubleRowMajor, -1),  false);

    QCOMPARE(EigenSparseUtils::validIndex(m_testMatInt, 99, 99),  true);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatInt, 99, 100),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatInt, 100, 99),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatInt, 100, 100),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatInt, 0, 0),  true);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatInt, 0, -1),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatInt, -1, 0),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatInt, -1, -1), false);

    QCOMPARE(EigenSparseUtils::validIndex(m_testMatIntRowMajor, 99, 99),  true);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatIntRowMajor, 99, 100),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatIntRowMajor, 100, 99),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatIntRowMajor, 100, 100),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatIntRowMajor, 0, 0),  true);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatIntRowMajor, 0, -1),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatIntRowMajor, -1, 0),  false);
    QCOMPARE(EigenSparseUtils::validIndex(m_testMatIntRowMajor, -1, -1), false);
}


void EigenSparseUtilsTests::medianTest() {
    QCOMPARE(EigenSparseUtils::median(m_testVecInt), 4.5);
    QCOMPARE(EigenSparseUtils::median(m_testVecDouble), 4.25);
    QCOMPARE(EigenSparseUtils::median(m_testVecDoubleRowMajor), 4.25);
    QCOMPARE(EigenSparseUtils::median(m_testMatInt), 4.5);
    QCOMPARE(EigenSparseUtils::median(m_testMatIntRowMajor), 4.5);
    QCOMPARE(EigenSparseUtils::median(m_testMatIntEmpty), 0.0);
}


void EigenSparseUtilsTests::removeElementBelowThresholdTest() {

    Eigen::SparseVector<int> testVecInt = m_testVecInt;
    Eigen::SparseVector<double> testVecDouble = m_testVecDouble;
    Eigen::SparseVector<double, Eigen::RowMajor> testVecDoubleRowMajor = m_testVecDoubleRowMajor;
    Eigen::SparseMatrix<int> testMatInt = m_testMatInt;
    Eigen::SparseMatrix<int, Eigen::RowMajor> testMatIntRowMajor = m_testMatIntRowMajor;

    const int threshold1 = -1;
    EigenSparseUtils::removeElementsBelowThreshold(threshold1, &testVecInt);
    const int minAfterThreshold1 = EigenSparseUtils::min(testVecInt);
    QCOMPARE(minAfterThreshold1, threshold1);

    const double threshold2 = -1.5;
    EigenSparseUtils::removeElementsBelowThreshold(threshold2, &testVecDouble);
    const double minAfterThreshold2 = EigenSparseUtils::min(testVecDouble);
    QCOMPARE(minAfterThreshold2, threshold2);

    const double threshold3 = -1.5;
    EigenSparseUtils::removeElementsBelowThreshold(threshold3, &testVecDoubleRowMajor);
    const double minAfterThreshold3 = EigenSparseUtils::min(testVecDoubleRowMajor);
    QCOMPARE(minAfterThreshold3, threshold3);

    const int threshold4 = -1;
    EigenSparseUtils::removeElementsBelowThreshold(threshold4, &testMatInt);
    const int minAfterThreshold4 = EigenSparseUtils::min(testMatInt);
    QCOMPARE(minAfterThreshold4, threshold4);

    const int threshold5 = -1;
    EigenSparseUtils::removeElementsBelowThreshold(threshold5, &testMatIntRowMajor);
    const int minAfterThreshold5 = EigenSparseUtils::min(testMatIntRowMajor);
    QCOMPARE(minAfterThreshold5, threshold5);
}


void EigenSparseUtilsTests::rollTest() {

    const int rollDistance = 20;

    const Eigen::SparseVector<double> rolledVec = EigenSparseUtils::roll(m_testVecDouble,
                                                                         rollDistance);
    QCOMPARE(rolledVec.coeffs().size(), 3);
    QCOMPARE(rolledVec.coeff(22), m_testVecDouble.coeff(2));
    QCOMPARE(rolledVec.coeff(25), m_testVecDouble.coeff(5));
    QCOMPARE(rolledVec.coeff(42), m_testVecDouble.coeff(22));


    const Eigen::SparseVector<double, Eigen::RowMajor> rolledVecRowMajor = EigenSparseUtils::roll(m_testVecDoubleRowMajor,
                                                                                                  rollDistance);
    QCOMPARE(rolledVecRowMajor.coeffs().size(), 3);
    QCOMPARE(rolledVecRowMajor.coeff(22), m_testVecDoubleRowMajor.coeff(2));
    QCOMPARE(rolledVecRowMajor.coeff(25), m_testVecDoubleRowMajor.coeff(5));
    QCOMPARE(rolledVecRowMajor.coeff(42), m_testVecDoubleRowMajor.coeff(22));
}


void EigenSparseUtilsTests::apexesTest() {

    Eigen::SparseVector<double, Eigen::RowMajor> vec(10000);
    vec.coeffRef(999) = 0.999;
    vec.coeffRef(1000) = 1.000;
    vec.coeffRef(1001) = 1.000;
    vec.coeffRef(1002) = 0.999;
    vec.coeffRef(1003) = 0.999;
    vec.coeffRef(1004) = 1.001;
    vec.coeffRef(1005) = 1.000;
    vec.coeffRef(1006) = 0.999;

    QMap<int, double> apexResult = EigenSparseUtils::apexes(vec);

    QCOMPARE(apexResult.size(), 2);
    QCOMPARE(apexResult.firstKey(), 1000);
    QCOMPARE(apexResult.lastKey(), 1004);
    QCOMPARE(static_cast<int>(apexResult.first()), 1);
    QCOMPARE(static_cast<int>(apexResult.last()), 1);

    Eigen::SparseMatrix<double, Eigen::RowMajor> mat(10000, 10000);

    mat.coeffRef(999, 999) = 1.999;
    mat.coeffRef(1000, 999) = 1.000;
    mat.coeffRef(1001, 999) = 1.000;
    mat.coeffRef(999, 1000) = 0.999;
    mat.coeffRef(1000, 1000) = 2.999;
    mat.coeffRef(1001, 1000) = 1.001;
    mat.coeffRef(999, 1001) = 0.999;
    mat.coeffRef(1000, 1001) = 0.999;
    mat.coeffRef(1001, 1001) = 1.001;

    mat.coeffRef(999, 1999) = 0.999;
    mat.coeffRef(1000, 1999) = 1.000;
    mat.coeffRef(1001, 1999) = 1.000;
    mat.coeffRef(999, 2000) = 0.999;
    mat.coeffRef(1000, 2000) = 3.999;
    mat.coeffRef(1001, 2000) = 1.001;
    mat.coeffRef(999, 2001) = 0.999;
    mat.coeffRef(1000, 2001) = 0.999;
    mat.coeffRef(1001, 2001) = 1.001;

    QVector<EigenSparseUtils::SparseMatrixPoint> apexes =EigenSparseUtils::apexes(mat);
    QCOMPARE(apexes.size(), 2);
    QCOMPARE(apexes.front().row,1000);
    QCOMPARE(apexes.front().col,2000);
    QCOMPARE(apexes.back().row,1000);
    QCOMPARE(apexes.back().col,1000);

    Eigen::SparseMatrix<double> matTest(4,10);

    QVector<double> vec1 = {0, 0, 0.000493879, 0.0140574, 0.163132, 1.00056, 3.54482, 7.52937, 9.68593, 7.52937};
    QVector<double> vec2 = {0.040279, 0.491388, 3.25722, 13.1678, 35.7941, 74.3548, 132.886, 199.942, 221.167, 161.062};
    QVector<double> vec3 = {0.809024, 9.8697, 65.2242, 258.826, 653.294, 1090.8, 1242.56, 985.914, 544.366, 204.988};
    QVector<double> vec4 = {0.040279, 0.491388, 3.2473, 12.8855, 32.5175, 54.258, 61.6867, 48.7109, 26.6201, 9.83087};

    int row = 0;
    for (const QVector<double> &v : {vec1, vec2, vec3, vec4}) {
        for (int c = 0; c < v.size(); c++) {
            matTest.coeffRef(row, c) = v.at(c);
        }
        row++;
    }

    const QVector<EigenSparseUtils::SparseMatrixPoint> apexBlock = EigenSparseUtils::apexes(matTest);
    QCOMPARE(apexBlock.size(), 1);
    QCOMPARE(apexBlock.front().row, 2);
    QCOMPARE(apexBlock.front().col, 6);
}


void EigenSparseUtilsTests::buildCombFilterTest() {

    const QVector<double> vals = {10.0, 20.0, 30.0, 40.0, 50.0};

    const int precision = 1;
    const double tolerance = 1.0;
    const double buffer = 1.0;
    const double maxLength = 50.0 + buffer;
    const int hashedMaxLength = MathUtils::hashDecimal(maxLength, precision);

    Eigen::SparseVector<double> vec(hashedMaxLength);
    for (double v : vals) {

        const int vHash = MathUtils::hashDecimal(v, precision);
        vec.coeffRef(vHash) = v;
    }

    const QVector<double> extractVals = {9.0, 28.9, 40.0, 49.1};
    Eigen::SparseMatrix<double> combFilter
                        = EigenSparseUtils::buildCombFilter(extractVals, tolerance, maxLength, precision);

    const Eigen::VectorXd result = combFilter * vec;
    const QVector<int> expectedResult = {10, 0, 40, 50};

    QCOMPARE(static_cast<int>(result.size()) ,expectedResult.size());
    for (int i = 0; i < expectedResult.size(); i++) {
        QCOMPARE( static_cast<int>(result.coeff(i)), expectedResult.at(i));
    }

//TODO figure out why you commented this out and left it in.

//    std::uniform_real_distribution<double> unif(100.0, 1800.0);
//    std::mt19937 re(std::random_device{}());
//    auto generator = bind(unif, std::ref(re));
//
////    std::srand(i);
//    const int n = 50; //rand() % 100 + 10;
//    QVector<double> v(n);
//    std::generate(v.begin(), v.end(), std::ref(generator) );
//
//
//    for (int i = 0; i < 10 * 100; i++) {
//
//        if (i % 10 == 0) {
//            qDebug() << i / 10;
//        }
//
//        EigenSparseUtils::buildCombFilter<double>(v, 0.01, 2000.0, 3);
////        break;
//    }

}


void EigenSparseUtilsTests::markVectorApexesTest() {

    Eigen::SparseVector<int> vecNormal(1000);
    vecNormal.insert(500) = 1;
    vecNormal.insert(501) = 2;
    vecNormal.insert(502) = 3;
    vecNormal.insert(503) = 1;
    vecNormal.insert(504) = 2;
    vecNormal.insert(505) = 1;

    const Eigen::SparseVector<int> resultNormal
            = EigenSparseUtils::markVectorApexes(vecNormal);

    QCOMPARE(resultNormal.coeff(502), 1);
    QCOMPARE(resultNormal.coeff(504), 1);


    Eigen::SparseVector<int> vecTriple(1000);
    vecTriple.insert(500) = 1;
    vecTriple.insert(501) = 3;
    vecTriple.insert(502) = 3;
    vecTriple.insert(503) = 3;
    vecTriple.insert(504) = 2;
    vecTriple.insert(505) = 1;

    const Eigen::SparseVector<int> resultTriple
            = EigenSparseUtils::markVectorApexes(vecTriple);

    QCOMPARE(resultTriple.coeff(501), 1);

    Eigen::SparseVector<int> vecLeftEdge(1000);
    vecLeftEdge.insert(0) = 11;
    vecLeftEdge.insert(1) = 3;

    const Eigen::SparseVector<int> resultLeftEdge
            = EigenSparseUtils::markVectorApexes(vecLeftEdge);

    QCOMPARE(resultLeftEdge.coeff(0), 1);


    Eigen::SparseVector<int> vecRightEdge(1000);
    vecRightEdge.insert(999) = 11;
    vecRightEdge.insert(998) = 3;

    const Eigen::SparseVector<int> resultRightEdge
            = EigenSparseUtils::markVectorApexes(vecRightEdge);

    QCOMPARE(resultRightEdge.coeff(999), 1);

}


void EigenSparseUtilsTests::markMatrixApexesTest() {

    Eigen::SparseMatrix<int> mat(1000, 1000);

    mat.insert(500, 499) = 9;
    mat.insert(500, 500) = 10;
    mat.insert(500, 501) = 9;

    mat.insert(499, 499) = 8;
    mat.insert(499, 500) = 9;
    mat.insert(499, 501) = 8;

    mat.insert(501, 499) = 8;
    mat.insert(501, 500) = 9;
    mat.insert(501, 501) = 8;

    Eigen::SparseMatrix<int> resultMat = EigenSparseUtils::markMatrixApexes(mat);
//    std::cout << resultMat.block(495, 495, 11, 11) << std::endl;
    QCOMPARE(resultMat.coeff(500, 500), 2);

    QVector<EigenSparseUtils::SparseMatrixPoint> smp = EigenSparseUtils::findMatrixApexes(mat);

    const EigenSparseUtils::SparseMatrixPoint p = smp.front();
    QCOMPARE(p.row, 500);
    QCOMPARE(p.col, 500);
    QCOMPARE(p.value, 10);

}


void EigenSparseUtilsTests::cleanupTestCase()
{}

void EigenSparseUtilsTests::solverTest() {

    int m=10000, n = 3;
    Eigen::SparseVector<double> x(n), b(m);
    Eigen::SparseMatrix<double> A(m,n);

    b.coeffRef(100) = 3.0;
    b.coeffRef(200) = 2.0;
    b.coeffRef(300) = 1.0;
    b.coeffRef(400) = 2.0;
    b.coeffRef(500) = 1.0;

    A.coeffRef(100, 0) = 1.0;
    A.coeffRef(200, 0) = 1.0;
    A.coeffRef(300, 0) = 1.0;
    A.coeffRef(100, 1) = 1.0;
    A.coeffRef(400, 1) = 1.0;
    A.coeffRef(400, 2) = 1.0;
    A.coeffRef(500, 2) = 1.0;

    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor> > lscg;
    lscg.setMaxIterations(10);
    lscg.setTolerance(1e-10);

    lscg.compute(A);
    x = lscg.solve(b);
    std::cout << "#iterations:     " << lscg.iterations() << std::endl;
    std::cout << "estimated error: " << lscg.error()      << std::endl;
    std::cout << x << std::endl;
// update b, and solve again
    x = lscg.solve(b);
    std::cout << x << std::endl;

}


QTEST_MAIN(EigenSparseUtilsTests)
#include "EigenSparseUtilsTests.moc"
