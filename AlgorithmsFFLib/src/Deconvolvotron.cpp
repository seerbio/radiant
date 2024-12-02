//
// Created by anichols on 6/5/24.
//

#include "Deconvolvotron.h"

#include "CandidateScores.h"
#include "MsUtils.h"

#include <unsupported/Eigen/NNLS>

#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>

#include "EigenUtils.h"

Deconvolvotron::Deconvolvotron()
: m_mzGroupingPrecision(-1)
, m_ppmExtractionTolerance(-1.0)
{}

Err Deconvolvotron::init(
    int mzGroupingPrecision,
    double ppmExtractionTolerance
    ) {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(mzGroupingPrecision, 0, ErrorUtilsParam::ExcludeThreshold); ree;
    m_mzGroupingPrecision = mzGroupingPrecision;

    ERR_RETURN
}

namespace {

    Err buildAMatrixAndBVec(
        const QVector<QPair<CandidateScores*, QVector<QPointF>>> &aMatrixPoints,
        const QVector<QPointF>& bVecPoints,
        const QMap<int, QVector<double>> &hashedXValsVsMzValsGrouped,
        double ppmExtractTol,
        Eigen::MatrixX<double> *aMat,
        Eigen::VectorX<double> *bVec
        ) {

        ERR_INIT

        e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
        e = ErrorUtils::isNotEmpty(hashedXValsVsMzValsGrouped); ree;

        QVector<QVector<QPointF>> points;
        std::transform(
            aMatrixPoints.begin(),
            aMatrixPoints.end(),
            std::back_inserter(points),
            [](const QPair<CandidateScores*, QVector<QPointF>> &pr){return pr.second;}
            );

        QVector<double> xValMeans;
        std::transform(
            hashedXValsVsMzValsGrouped.begin(),
            hashedXValsVsMzValsGrouped.end(),
            std::back_inserter(xValMeans),
            [](const QVector<double> &xVals){return MathUtils::mean(xVals);}
            );
        std::sort(xValMeans.begin(), xValMeans.end());

        QVector<QPointF> bVecExtractPoints = MsUtils::extractPointsFromPoints(
            bVecPoints,
            xValMeans,
            ppmExtractTol,
            false
            );

        aMat->resize(bVecExtractPoints.size(), aMatrixPoints.size());
        aMat->setZero();

        for (int col = 0; col < aMatrixPoints.size(); col++) {
            const QVector<QPointF> &pnts = points.at(col);
            for (const QPointF &pnt : pnts) {
                const int row = MathUtils::closest(xValMeans, pnt.x());

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
        const int df1 = std::max(p - 1, 1);
        const int df2 = std::max(n - p, 1);

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
    const QVector<QPair<CandidateScores*, QVector<QPointF>>> &aMatrixPoints,
    const QVector<QPointF>& bVecPoints,
    QVector<QPair<CandidateScores*, DeconvolvotronResult>> *candScoresVsScore
    ) const  {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(m_mzGroupingPrecision, 0, ErrorUtilsParam::ExcludeThreshold); ree;
    e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
    e = ErrorUtils::isNotEmpty(bVecPoints); ree;

    candScoresVsScore->clear();

    QVector<CandidateScores*> candidateScoreses;
    std::transform(
        aMatrixPoints.begin(),
        aMatrixPoints.end(),
        std::back_inserter(candidateScoreses),
        [](const QPair<CandidateScores*, QVector<QPointF>> &pr){return pr.first;}
        );

    QMap<int, QVector<double>> hashedXValsVsMzValsGrouped;
    e = buildXValsSet(
        aMatrixPoints,
        &hashedXValsVsMzValsGrouped
        ); ree;

    Eigen::MatrixX<double> aMat;
    Eigen::VectorX<double> bVec;
    e = buildAMatrixAndBVec(
        aMatrixPoints,
        bVecPoints,
        hashedXValsVsMzValsGrouped,
        20,
        &aMat,
        &bVec
        ); ree;

    if (MathUtils::tZero(bVec.maxCoeff())) {
        ERR_RETURN
    }

    Eigen::NNLS<Eigen::MatrixXd> nnls(aMat);
    Eigen::VectorXd x = nnls.solve(bVec);
    x /= x.sum();

// #define PRINT_MATS
#ifdef PRINT_MATS
std::cout << "aMat = [" << std::endl;
for (int row = 0; row < aMat.rows(); row++) {
    for (int col = 0; col < aMat.cols(); col++) {
        std::cout << aMat.coeff(row, col) << ",";
    }
    std::cout << std::endl;
}
std::cout << "]" << std::endl;

std::cout << "bVec = [" << std::endl;
for (int i = 0; i < bVec.size(); i++) {
    std::cout << bVec.coeff(i) << ",";
}
std::cout << std::endl;
std::cout << "]" << std::endl;
#endif

    double pValFTest;
    QVector<double> coeffsPVals;
    // deconvolveStats(
    //     aMat,
    //     x,
    //     bVec,
    //     &pValFTest,
    //     &coeffsPVals
    //     );

    e = ErrorUtils::isEqual(aMatrixPoints.size(), static_cast<int>(x.size())); ree;
    // e = ErrorUtils::isEqual(aMatrixPoints.size(), coeffsPVals.size()); ree;

    for (int i = 0; i < x.size(); i++) {

        if (MathUtils::tZero(x.coeff(i))) {
            continue;
        }

        DeconvolvotronResult deconvolvotronResult;
        deconvolvotronResult.discScore = x.coeff(i);
        // deconvolvotronResult.pVal = coeffsPVals.at(i);
        // deconvolvotronResult.pValFrameFtest = pValFTest;

        candScoresVsScore->push_back({candidateScoreses.at(i), deconvolvotronResult});
    }

    ERR_RETURN
}

Err Deconvolvotron::buildXValsSet(
        const QVector<QPair<CandidateScores*, QVector<QPointF>>> &aMatrixPoints,
        QMap<int, QVector<double>> *hashedXValsVsMzValsGrouped
        ) const {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;

    for (const QPair<CandidateScores*, QVector<QPointF>> &aMatPointsPr : aMatrixPoints) {
        for (const QPointF &pnt : aMatPointsPr.second) {
            const int hashedXVal = MathUtils::hashDecimal(pnt.x(), m_mzGroupingPrecision);
            (*hashedXValsVsMzValsGrouped)[hashedXVal].push_back(pnt.x());
        }
    }

    ERR_RETURN
}

Err Deconvolvotron::buildAMatrixAndBVecTestAccess(
    const QVector<QPair<CandidateScores*, QVector<QPointF>>>& aMatrixPoints,
    const QVector<QPointF>& bVecPoints,
    const QMap<int, QVector<double>>& hashedXValsVsMzValsGrouped,
    double ppmExtractTol,
    QVector<QVector<double>> *aMat,
    QVector<double> *bVec
    ) {

    ERR_INIT

    e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
    e = ErrorUtils::isNotEmpty(bVecPoints); ree;
    e = ErrorUtils::isNotEmpty(hashedXValsVsMzValsGrouped); ree;
    e = ErrorUtils::isTrue(ppmExtractTol > 0);

    Eigen::MatrixX<double> aMatEigen;
    Eigen::VectorX<double> bVecEigen;
    e = buildAMatrixAndBVec(
        aMatrixPoints,
        bVecPoints,
        hashedXValsVsMzValsGrouped,
        ppmExtractTol,
        &aMatEigen,
        &bVecEigen
        ); ree;

    std::cout << aMatEigen << std::endl;
    std::cout << "****" << std::endl;
    std::cout << bVecEigen << std::endl;

    *aMat = EigenUtils::convertEigenMatrixToQVectors(aMatEigen);
    *bVec = EigenUtils::convertEigenVectorToQVector(bVecEigen);

    ERR_RETURN
}
