//
// Created by anichols on 11/07/2021.
//

#include "TandemSpectraDeconvolvotron.h"
#include "FragLibReader.h"

#include "ParallelUtils.h"

#include <QElapsedTimer>
#include <QDebug>
#include <QtTest/QtTest>

#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>

class TandemSpectraDeconvolvotronTests : public QObject
{

Q_OBJECT

public:
    TandemSpectraDeconvolvotronTests();
    ~TandemSpectraDeconvolvotronTests() override = default;

private Q_SLOTS:

    void deconvolveTandemSpectraTest();
    void deconvolveSuperToy();

};

TandemSpectraDeconvolvotronTests::TandemSpectraDeconvolvotronTests() : QObject(){}

void TandemSpectraDeconvolvotronTests::deconvolveTandemSpectraTest() {

    ERR_INIT

    const MS2Ion i1(100.0, 1, "");
    const MS2Ion i2(200.0, 1, "");
    const MS2Ion i3(300.0, 1, "");

    const MS2Ion j1(100.0, 1, "");
    const MS2Ion j2(200.0, 1, "");
    const MS2Ion j3(600.0, 1, "");

    const QVector<QPointF> points = {
            {100, 1.5},
            {200, 1.5},
            {300, 1},
            {600, 0.5},
    };

    const QVector<MS2Ion> cand1 = {
            i1, i2, i3
    };

    const QVector<MS2Ion> cand2 = {
            j1, j2, j3
    };

    QMap<PeptideSequenceChargeKey, QVector<MS2Ion>> input = {
            {"C1", cand1},
            {"C2", cand2}
    };

    QMap<PeptideSequenceChargeKey, TandemDeconvolverResult> result;

    TandemSpectraDeconvolvotron deconvolvotron;
    e = deconvolvotron.init(3, 1000, 20, 0.000000000000000001, 0.05);
    QCOMPARE(e, eNoError);

    double fStat;
    double pValFTest;
    e = deconvolvotron.deconvolveTandemSpectra(
            points,
            input,
            &result
            );
    QCOMPARE(e, eNoError);

    const QStringList expectedKeys = {"C1", "C2"};

    QCOMPARE(result.keys(), expectedKeys);
    QCOMPARE(QString::number(result.value("C1").discScore), QString::number(0.666666666667));
    QCOMPARE(QString::number(result.value("C2").discScore), QString::number(0.333333333333));

}

void TandemSpectraDeconvolvotronTests::deconvolveSuperToy() {

    QSKIP("Leave this for troubleshooting");

    Eigen::MatrixX<double> Xd(5, 3);
    Xd << 1,1,1,1,2,3,1,3,6,1,5,14,1,17,2;

    Eigen::SparseMatrix<double> A = Xd.sparseView();

    Eigen::VectorX<double> b(5);
    b << 2,3,4,6,18;

    Eigen::VectorX<double> x(3);
    x << 1,1,0;

    const double residualSumOfSquares = (b - A * x).squaredNorm();

    int n = static_cast<int>(b.size());
    int p = static_cast<int>(A.cols());
    int df1 = p - 1;
    int df2 = n - p;

    //NOTE: not using numeric_limits because that produces inf below.
    const double nearZero = 1e-33;

    double mse = residualSumOfSquares / df2;
    mse = MathUtils::tZero(mse) ? nearZero : mse;

    const double fStat = (((x.transpose() * A.transpose() * A * x) / p) / mse)(0);
    const double p_value_f = 1 - boost::math::cdf(boost::math::fisher_f(df1, df2), fStat);

    Eigen::SparseMatrix<double> AtA = A.transpose() * A;
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(AtA);
    Eigen::MatrixX<double> cov = mse * solver.solve(Eigen::SparseMatrix<double>(AtA)).rhs();
    Eigen::VectorX<double> se = cov.diagonal().array().sqrt();
    Eigen::VectorX<double> t = x.array() / se.array();

    std::vector<double> p_values_t(p);
    for (int i = 0; i < p; i++) {
        double abs_t = std::abs(t(i));
        p_values_t[i] = 2 * (1 - boost::math::cdf(boost::math::students_t(df2), abs_t));
    }

    // Print results
    std::cout << "F-statistic: " << fStat << std::endl;
    std::cout << "p-value (F-test): " << p_value_f << std::endl;
    for (int i = 0; i < p; i++) {
        std::cout << "t-test for coefficient " << i << ": " << t(i) << std::endl;
        std::cout << "p-value for coefficient " << i << ": " << p_values_t[i] << std::endl;
    }

}

QTEST_MAIN(TandemSpectraDeconvolvotronTests)
#include "TandemSpectraDeconvolvotronTests.moc"
