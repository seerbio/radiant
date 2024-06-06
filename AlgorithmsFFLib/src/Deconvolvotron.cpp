//
// Created by anichols on 6/5/24.
//

#include "Deconvolvotron.h"

#include "MsUtils.h"

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
        const QVector<QPair<IdStr, QVector<QPointF>>> &aMatrixPoints,
        const QVector<QPointF>& bVecPoints,
        const QVector<double> &xVals,
        double ppmExtractTol,
        Eigen::MatrixX<double> *aMat,
        Eigen::VectorX<double> *bVec
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
        e = ErrorUtils::isNotEmpty(xVals); ree;

        QVector<QVector<QPointF>> points;
        std::transform(
            aMatrixPoints.begin(),
            aMatrixPoints.end(),
            std::back_inserter(points),
            [](const QPair<IdStr, QVector<QPointF>> &pr){return pr.second;}
            );

        QVector<double> xValsSorted = xVals;
        std::sort(xValsSorted.begin(), xValsSorted.end());

        QVector<QPointF> bVecExtractPoints = MsUtils::extractPointsFromPoints(
            bVecPoints,
            xValsSorted,
            ppmExtractTol
            );

        aMat->resize(bVecExtractPoints.size(), aMatrixPoints.size());
        aMat->setZero();
        for (int col = 0; col < aMatrixPoints.size(); col++) {
            const QVector<QPointF> &pnts = points.at(col);
            for (const QPointF &pnt : pnts) {
                const int row = MathUtils::closest(xValsSorted, pnt.x());

                if (row < 0 || row >= bVecExtractPoints.size()) {
                    continue;
                }

                aMat->coeffRef(row, col) = pnt.y();
            }
        }

        bVec->resize(bVecExtractPoints.size());
        bVec->setZero();
        for (int i = 0; i < bVecExtractPoints.size(); i++) {
            const QPointF &pnt = bVecExtractPoints.at(i);
            bVec->coeffRef(i) = pnt.y();
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

        const double mse = residuals.array().square().sum() / (static_cast<double>(A.rows()) - static_cast<double>(A.cols()));
        const double tss = (b.array() - b.array().mean()).square().sum();
        const double rss = tss - mse * (A.rows() - 1.0);

        const double fStatistic = (rss / A.cols()) / mse;
        const boost::math::fisher_f distF(df1, df2);
        *pValFTest = 1.0 - boost::math::cdf(distF, std::max(fStatistic, 0.00));

        const Eigen::VectorXd t = beta.cwiseQuotient(se);

        const boost::math::students_t dist(static_cast<double>(A.rows()) - static_cast<double>(A.cols()));
        Eigen::VectorXd pVal(A.cols());
        for (int i = 0; i < A.cols(); ++i){
            const double fbs = std::isnan(fabs(t[i])) ? 0.0 : fabs(t[i]);
            pVal[i] = static_cast<double>(2.0 * static_cast<double>(1.0 - boost::math::cdf(dist, fbs)));
        }

        *coeffsPVals = EigenUtils::convertEigenVectorToQVector(pVal);
    }

}//namespace
Err Deconvolvotron::deconvolve(
    const QVector<QPair<IdStr, QVector<QPointF>>> &aMatrixPoints,
    const QVector<QPointF>& bVecPoints,
    QVector<QPair<IdStr, DeconvolvotronResult>> *idStrVsScore
    ) const  {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(m_precision, 0, ErrorUtilsParam::ExcludeThreshold); ree;
    e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
    e = ErrorUtils::isNotEmpty(bVecPoints); ree;

    idStrVsScore->clear();

    QVector<IdStr> idStrs;
    std::transform(
        aMatrixPoints.begin(),
        aMatrixPoints.end(),
        std::back_inserter(idStrs),
        [](const QPair<IdStr, QVector<QPointF>> &pr){return pr.first;}
        );

    QVector<double> xVals;
    e = buildXValsSet(
        aMatrixPoints,
        bVecPoints,
        &xVals
        ); ree;

    Eigen::MatrixX<double> aMat;
    Eigen::VectorX<double> bVec;
    e = buildAMatrixAndBVec(
        aMatrixPoints,
        bVecPoints,
        xVals,
        m_precision,
        &aMat,
        &bVec
        ); ree;

    if (MathUtils::tZero(bVec.maxCoeff())) {
        ERR_RETURN
    }

    Eigen::NNLS<Eigen::MatrixXd> nnls(aMat);
    Eigen::VectorXd x = nnls.solve(bVec);
    x /= x.sum();

    // double pValFTest;
    // QVector<double> coeffsPVals;
    // deconvolveStats(
    //     aMat,
    //     x,
    //     bVec,
    //     &pValFTest,
    //     &coeffsPVals
    //     );

    // e = ErrorUtils::isEqual(aMatrixPoints.size(), static_cast<int>(x.size())); ree;
    // e = ErrorUtils::isEqual(aMatrixPoints.size(), coeffsPVals.size()); ree;

    for (int i = 0; i < x.size(); i++) {

        if (MathUtils::tZero(x.coeff(i))) {
            continue;
        }

        DeconvolvotronResult deconvolvotronResult;
        deconvolvotronResult.discScore = x.coeff(i);
        // deconvolvotronResult.pVal = coeffsPVals.at(i);
        // deconvolvotronResult.pValFrameFtest = pValFTest;

        idStrVsScore->push_back({idStrs.at(i), deconvolvotronResult});
    }

    ERR_RETURN
}

Err Deconvolvotron::buildXValsSet(
        const QVector<QPair<IdStr, QVector<QPointF>>> &aMatrixPoints,
        const QVector<QPointF>& bVecPoints,
        QVector<double> *xValsReturn
        ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
    e = ErrorUtils::isNotEmpty(bVecPoints); ree;

    QVector<int> hashedXVals;
    for (const QPair<IdStr, QVector<QPointF>> &aMatPointsPr : aMatrixPoints) {
        for (const QPointF &pnt : aMatPointsPr.second) {
            const int hashedXVal = MathUtils::hashDecimal(pnt.x(), m_precision);
            hashedXVals.push_back(hashedXVal);
        }
    }

    QSet<int> hashedXValsSet = {hashedXVals.begin(), hashedXVals.end()};

    QVector<double> xValsSet;
    std::transform(
        hashedXValsSet.begin(),
        hashedXValsSet.end(),
        std::back_inserter(xValsSet),
        [&](int i){return MathUtils::unHashDecimal<double>(i, m_precision);}
        );

    *xValsReturn = xValsSet;

    ERR_RETURN
}
