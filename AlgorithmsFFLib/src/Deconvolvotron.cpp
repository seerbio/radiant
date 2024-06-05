//
// Created by anichols on 6/5/24.
//

#include "Deconvolvotron.h"

#include "Eigen/Dense"
#include <unsupported/Eigen/NNLS>

#include <iostream>

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



}//namespace
Err Deconvolvotron::deconvolve(
    const QMap<IdStr, QVector<QPointF>>& aMatrixPoints,
    const QVector<QPointF>& bVecPoints,
    QVector<QPair<IdStr, Score>>* idStrVsScore
    ) {

    ERR_INIT

    e = ErrorUtils::isAboveThreshold(m_precision, 0, ErrorUtilsParam::ExcludeThreshold); ree;
    e = ErrorUtils::isNotEmpty(aMatrixPoints); ree;
    e = ErrorUtils::isNotEmpty(bVecPoints); ree;

    idStrVsScore->clear();

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

    std::cout << aMat << std::endl;
    std::cout << "---" << std::endl;
    std::cout << bVec << std::endl;
    std::cout << "***" << std::endl;
    std::cout << x << std::endl;


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