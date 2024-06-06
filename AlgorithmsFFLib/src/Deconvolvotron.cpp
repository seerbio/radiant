//
// Created by anichols on 6/5/24.
//

#include "Deconvolvotron.h"

#include <unsupported/Eigen/NNLS>

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>

#include "EigenUtils.h"

Deconvolvotron::Deconvolvotron() : m_precision(-1) {}

Err Deconvolvotron::init(int precision) {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(precision, 0, ErrorUtilsParam::ExcludeThreshold); ree;
    m_precision = precision;

    ERR_RETURN
}

namespace {

    Err buildAMatrixAndBVec(
        const QMap<IdStr, QVector<QPointF>>& aMatrixPoints,
        const QVector<QPointF>& bVecPoints,
        const QVector<int> &hashedXVals,
        int precision,
        Eigen::MatrixX<double> *aMat,
        Eigen::VectorX<double> *bVec
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
        e = ErrorUtils::isNotEmpty(hashedXVals); ree;

        const QVector<IdStr> keys = aMatrixPoints.keys().toVector();

        QVector<int> hashedXValsSorted = hashedXVals;
        std::sort(hashedXValsSorted.begin(), hashedXValsSorted.end());

        int columnCounter = 0;
        QHash<int, int> hashedXValVsMatrixColumn;
        for (int hashedXVal : hashedXValsSorted) {
            hashedXValVsMatrixColumn.insert(hashedXVal, columnCounter++);
        }

        aMat->resize(hashedXValsSorted.size(), aMatrixPoints.size());
        aMat->setZero();
        for (int col = 0; col < aMatrixPoints.size(); col++) {
            const QVector<QPointF> &pnts = aMatrixPoints.value(keys.at(col));
            for (const QPointF &pnt : pnts) {
                const int hashedXVal = MathUtils::hashDecimal(pnt.x(), precision);

                e = ErrorUtils::contains(hashedXVal, hashedXValVsMatrixColumn); ree;
                const int row = hashedXValVsMatrixColumn.value(hashedXVal);

                aMat->coeffRef(row, col) = pnt.y();
            }
        }

        bVec->resize(hashedXValsSorted.size());
        bVec->setZero();
        for (const QPointF &pnt : bVecPoints) {
            const int hashedXVal = MathUtils::hashDecimal(pnt.x(), precision);

            e = ErrorUtils::contains(hashedXVal, hashedXValVsMatrixColumn); ree;
            const int row = hashedXValVsMatrixColumn.value(hashedXVal);
            bVec->coeffRef(row) = pnt.y();
        }

        ERR_RETURN
    }

    void deconvolveStats(
            const Eigen::MatrixX<double> &A,
            const Eigen::VectorX<double> &x,
            const Eigen::VectorX<double> &b,
            double *pValFTest,
            QVector<double> *coeffsPVals
            ) {

        const int n = static_cast<int>(b.size());
        const int p = static_cast<int>(A.cols());
        const int df1 = p - 1;
        const int df2 = n - p;

        const Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const Eigen::VectorXd beta = svd.solve(b);

        const Eigen::VectorXd residuals = b - A * beta;

        const double variance = residuals.squaredNorm() / (static_cast<double>(A.rows()) - static_cast<double>(A.cols()));
        const Eigen::VectorXd se = variance * (A.transpose() * A).inverse().diagonal().cwiseSqrt();

        // Calculate mean square of residuals (MSE):
        const double mse = residuals.array().square().sum() / (static_cast<double>(A.rows()) - static_cast<double>(A.cols()));
        const double tss = (b.array() - b.array().mean()).square().sum();
        const double rss = tss - mse * (A.rows() - 1.0);

        const double fStatistic = (rss / A.cols()) / mse;
        const boost::math::fisher_f distF(df1, df2);
        *pValFTest = 1.0 - boost::math::cdf(distF, fStatistic);

        const Eigen::VectorXd t = beta.cwiseQuotient(se);

        const boost::math::students_t dist(static_cast<double>(A.rows()) - static_cast<double>(A.cols()));
        Eigen::VectorXd pVal(A.cols());
        for (int i = 0; i < A.cols(); ++i){
            pVal[i] = static_cast<double>(2.0 * static_cast<double>(1.0 - boost::math::cdf(dist, fabs(t[i]))));
        }

        *coeffsPVals = EigenUtils::convertEigenVectorToQVector(pVal);
    }

}//namespace
Err Deconvolvotron::deconvolve(
    const QMap<IdStr, QVector<QPointF>>& aMatrixPoints,
    const QVector<QPointF>& bVecPoints,
    QVector<QPair<IdStr, DeconvolvotronResult>> *idStrVsScore
    ) const  {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(m_precision, 0, ErrorUtilsParam::ExcludeThreshold); ree;
    e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
    e = ErrorUtils::isNotEmpty(bVecPoints); ree;

    idStrVsScore->clear();

    const QVector<IdStr> idStrs = aMatrixPoints.keys().toVector();

    QVector<int> hashedXVals;
    e = buildXValsSet(
        aMatrixPoints,
        bVecPoints,
        &hashedXVals
        ); ree;

    Eigen::MatrixX<double> aMat;
    Eigen::VectorX<double> bVec;
    e = buildAMatrixAndBVec(
        aMatrixPoints,
        bVecPoints,
        hashedXVals,
        m_precision,
        &aMat,
        &bVec
        ); ree;

    Eigen::NNLS<Eigen::MatrixXd> nnls(aMat);
    Eigen::VectorXd x = nnls.solve(bVec);
    x /= x.sum();

    double pValFTest;
    QVector<double> coeffsPVals;
    deconvolveStats(
        aMat,
        x,
        bVec,
        &pValFTest,
        &coeffsPVals
        );

    e = ErrorUtils::isEqual(aMatrixPoints.size(), static_cast<int>(x.size())); ree;
    e = ErrorUtils::isEqual(aMatrixPoints.size(), coeffsPVals.size()); ree;
    idStrVsScore->resize(aMatrixPoints.size());

    for (int i = 0; i < idStrVsScore->size(); i++) {
        DeconvolvotronResult deconvolvotronResult;
        deconvolvotronResult.discScore = x.coeff(i);
        deconvolvotronResult.pVal = coeffsPVals.at(i);
        deconvolvotronResult.pValFrameFtest = pValFTest;
        (*idStrVsScore)[i] = {idStrs.at(i), deconvolvotronResult};
    }

    ERR_RETURN
}

Err Deconvolvotron::buildXValsSet(
        const QMap<IdStr, QVector<QPointF>>& aMatrixPoints,
        const QVector<QPointF>& bVecPoints,
        QVector<int> *hashedXValsReturn
        ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
    e = ErrorUtils::isNotEmpty(bVecPoints); ree;

    QVector<int> hashedXVals;
    for (const QVector<QPointF> &aMatPoints : aMatrixPoints) {
        for (const QPointF &pnt : aMatPoints) {
            const int hashedXVal = MathUtils::hashDecimal(pnt.x(), m_precision);
            hashedXVals.push_back(hashedXVal);
        }
    }

    for (const QPointF &pnt : bVecPoints) {
        const int hashedXVal = MathUtils::hashDecimal(pnt.x(), m_precision);
        hashedXVals.push_back(hashedXVal);
    }

    QSet<int> hashedXValsSet = {hashedXVals.begin(), hashedXVals.end()};

    *hashedXValsReturn = {hashedXValsSet.begin(), hashedXValsSet.end()};

    ERR_RETURN
}