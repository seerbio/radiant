//
// Created by anichols on 4/21/23.
//
#include "TandemSpectraDeconvolvotron.h"

#include "EigenUtils.h"
#include "ErrorUtils.h"
#include "MathUtils.h"
#include "FragLibReader.h"

#include <Eigen/Sparse>

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>


TandemSpectraDeconvolvotron::TandemSpectraDeconvolvotron()
: m_precision(-1)
, m_mzMax(-1.0)
, m_iterationsMax(-1)
, m_stopTolerance(-1.0)
, m_isInit(false)
, m_pValThreshold(-1.0)
{}

Err TandemSpectraDeconvolvotron::init(
        int precision,
        double mzMax,
        int iterationsMax,
        double stopTolerance,
        double pValThreshold
        ) {

    ERR_INIT

    const int precisionMin = 1;
    e = ErrorUtils::isAboveThreshold(
            precision,
            precisionMin,
            ErrorUtilsParam::IncludeThreshold
    ); ree;
    m_precision = precision;

    const double mzMaxMin = 100.0;
    e = ErrorUtils::isAboveThreshold(
            mzMax,
            mzMaxMin,
            ErrorUtilsParam::IncludeThreshold
    ); ree;
    m_mzMax = mzMax;

    const int iterationsMaxMin = 1;
    e = ErrorUtils::isAboveThreshold(
            iterationsMax,
            iterationsMaxMin,
            ErrorUtilsParam::IncludeThreshold
    ); ree;
    m_iterationsMax = iterationsMax;

    const double stopTolMin = 1e-100;
    m_stopTolerance = std::max(stopTolerance, stopTolMin);

    m_pValThreshold = pValThreshold;

    m_isInit = true;

    ERR_RETURN
}

namespace {

    Eigen::SparseMatrix<double> buildDeconvolveMatrix(
            const QList<QVector<MS2Ion>> &targetsMs2Ions,
            int precision,
            double mzMax
            ) {

        const int rows = MathUtils::hashDecimal(mzMax, precision);
        const int cols = targetsMs2Ions.size();

        Eigen::SparseMatrix<double> mat(rows, cols);

        int currentTargetIndex = 0;
        for (const QVector<MS2Ion> &targetIons : targetsMs2Ions) {

            for (const MS2Ion &ion : targetIons) {

                const int mzHashed = MathUtils::hashDecimal(ion.mz, precision);

                if (mzHashed >= rows || currentTargetIndex >= cols){
                    continue;
                }

                mat.insert(mzHashed, currentTargetIndex) = ion.intensity;

            }

            currentTargetIndex++;
        }

        return mat;
    }

    Err deconvolveStats(
            const Eigen::SparseMatrix<double> &A,
            const Eigen::VectorX<double> &x,
            const Eigen::VectorX<double> &b,
            double *fStat,
            double *pValFTest,
            QVector<double> *coeffsTTests,
            QVector<double> *coeffsPVals
    ) {

        ERR_INIT

        e = ErrorUtils::isAboveThreshold(
                static_cast<int>(A.nonZeros()),
                0,
                ErrorUtilsParam::ExcludeThreshold
                ); ree;

        const double residualSumOfSquares = (b - A * x).squaredNorm();

        const int n = static_cast<int>(b.size());
        const int p = static_cast<int>(A.cols());
        const int df1 = p - 1;
        const int df2 = n - p;

        //NOTE: not using numeric_limits::min() because that produces inf below.
        const double nearZero = 1e-33;

        double mse = residualSumOfSquares / df2;
        mse = MathUtils::tZero(mse) ? nearZero : mse;

        *fStat = (((x.transpose() * A.transpose() * A * x) / p) / mse)(0);
        *pValFTest = 1 - boost::math::cdf(boost::math::fisher_f(df1, df2), *fStat);

        const Eigen::SparseMatrix<double> AtA = A.transpose() * A;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(AtA);

        const Eigen::MatrixX<double> cov = mse * solver.solve(Eigen::SparseMatrix<double>(AtA)).rhs();
        const Eigen::VectorX<double> se = cov.diagonal().array().sqrt();
        const Eigen::VectorX<double> t = x.array() / se.array();

        std::vector<double> p_values_t(p);
        for (int i = 0; i < p; i++) {
            double abs_t = std::abs(t(i));
            p_values_t[i] = 2 * (1 - boost::math::cdf(boost::math::students_t(df2), abs_t));
        }

        for (int i = 0; i < p; i++) {
            coeffsTTests->push_back(t(i));
            coeffsPVals->push_back(p_values_t[i]);
        }

        ERR_RETURN
    }

}//namespace
Err TandemSpectraDeconvolvotron::deconvolveTandemSpectra(
        const ScanPoints &scanPoints,
        const QMap<PeptideStringWithMods, QVector<MS2Ion>> &tandemPredictions,
        QMap<PeptideStringWithMods, TandemDeconvolverResult> *pepSeqVsWeight
        ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    e = ErrorUtils::isNotEmpty(tandemPredictions); ree;

    pepSeqVsWeight->clear();

    e = ErrorUtils::isTrue(m_isInit); ree;
    e = ErrorUtils::isNotEmpty(scanPoints); ree;
    e = ErrorUtils::isNotEmpty(tandemPredictions); ree;

    const Eigen::SparseMatrix<double> mat = buildDeconvolveMatrix(
            tandemPredictions.values(),
            m_precision,
            m_mzMax
            );

    Eigen::VectorX<double> vecScanPoints = EigenUtils::convertQPointFVecToEigen(
            scanPoints,
            m_precision,
            m_mzMax
            );

    vecScanPoints /= vecScanPoints.maxCoeff();

    Eigen::VectorX<double> x(mat.cols());
    x.setZero();

    Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> lscg;
    lscg.setMaxIterations(m_iterationsMax);
    lscg.setTolerance(m_stopTolerance);

    lscg.compute(mat);
    x = lscg.solve(vecScanPoints);
    const double error = lscg.error();

    const int minScanSize = 2;
    if (tandemPredictions.size() < minScanSize) {

        TandemDeconvolverResult tdr;
        tdr.discScore = x.coeff(0);
        tdr.tTestVal = 1000;
        tdr.pVal = 0;
        tdr.frameFStat = -1.0;
        tdr.pValFrameFtest = -1.0;
        tdr.frameError = error;

        const PeptideStringWithMods &peptideSequenceChargeKey = tandemPredictions.firstKey();
        pepSeqVsWeight->insert(peptideSequenceChargeKey, tdr);

        ERR_RETURN
    }

    QVector<double> coeffsTTests;
    QVector<double> coeffsPVal;
    double fStat;
    double pValFTest;
    e = deconvolveStats(
            mat,
            x,
            vecScanPoints,
            &fStat,
            &pValFTest,
            &coeffsTTests,
            &coeffsPVal
            ); ree;

    const QList<PeptideStringWithMods> &keys = tandemPredictions.keys();
    e = ErrorUtils::isEqual(coeffsTTests.size(), coeffsPVal.size()); ree;
    e = ErrorUtils::isEqual(coeffsTTests.size(), keys.size()); ree;

    for (int i = 0; i < keys.size(); i++) {

        const PeptideStringWithMods &peptideSequenceChargeKey = keys.at(i);

        TandemDeconvolverResult tdr;
        tdr.discScore = x.coeff(i);
        tdr.tTestVal = coeffsTTests.at(i);
        tdr.pVal = coeffsPVal.at(i);
        tdr.frameFStat = fStat;
        tdr.pValFrameFtest = pValFTest;
        tdr.frameError = error;

        pepSeqVsWeight->insert(peptideSequenceChargeKey, tdr);
    }

//#define DEBUG_MAT_SOLVE
#ifdef DEBUG_MAT_SOLVE
    qDebug() << "#iterations:" << lscg.iterations()
             << "estimated error: " << lscg.error();
#endif

    ERR_RETURN
}
